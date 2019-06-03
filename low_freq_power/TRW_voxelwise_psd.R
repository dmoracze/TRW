########################### LIBRARIES ###########################
library(bspec)
library(abind)
library(ggplot2)
library(reshape2)

########################### FUNCTIONS ###########################
# function to find AFNI's R functions
first.in.path <- function(file) {
   ff <- paste(strsplit(Sys.getenv('PATH'),':')[[1]],'/', file, sep='')
   ff<-ff[lapply(ff,file.exists)==TRUE];
   return(gsub('//','/',ff[1], fixed=TRUE)) 
}

# function to read in the all participants time series
# also,
# 1. cut the timeseries to remove the ramp up and end baseline
# 2. scales by unit variance
# reads in int1 and int2 for both hemis
read.ts <- function(dir,subj) {
	source(first.in.path('AFNIio.R')) # source AFNI's R functions
	fin <- list() # create empty list container for output
	# loop through conditions
	for (cc in levels(as.factor(c('int1','int2')))) {
		cat('--------\n',cc,'\n')
		hemi <- list() # create empty list container for hemi
		# loop though hemispheres
		for (hh in levels(as.factor(c('lh','rh')))) {
			cat(' ',hh,'\n')
			# 271=timepoints, 36002=nodes
			ts <- array(0,dim=c(271,36002,length(subj))) # create output container for timeseries brick
			# loop through subjects
			for (ss in 1:length(subj)) {
				cat('   ',subj[ss],'\n')
				# read the lh data into R
				temp <- read.AFNI(paste0(dir,'/',subj[ss],'/bold/',subj[ss],'_',cc,'.pb06.',hh,'.smooth.niml.dset'),forcedset=TRUE)$brk
				temp <- t(temp[,,,])[14:284,] # reduce dimensions, transpose, cut timeseries
				ts[,,ss] <- scale(temp) # scale by unit variance
			}
			hemi[[hh]] <- ts
		}
		fin[[cc]] <- hemi
	}
	return(fin)
}

# function to remove the subject's global mean from the nodewise timeseries
rm.global <- function(data) {
	fin <- list() # create empty container for results
	# loop through conditions, want to keep int1 and int2 separate still
	for (cc in levels(as.factor(c('int1','int2')))) {
		cat(' ',cc,'\n')
		# concatenate the left and right hemis
		temp <- abind(data[[cc]]$lh,data[[cc]]$rh,along=2)
		# find the global mean of all nodes' timeseries
		glb <- apply(temp,c(1,3),mean)
		# subtract the global mean from timeseries
		glb.rm <- sweep(temp,c(1,3),glb,'-')
		# split back into lh and rh
		fin[[cc]] <- list(lh=unname(glb.rm[,1:36002,]),rh=unname(glb.rm[,36003:72004,]))
	}
	return(fin)
}

# function to find the power spectrum at the node level
voxel.indPSD <- function(data,win,olap) {
	fin <- list() # create empty container for results
	# find the number of power spectra in a sample to create container later
	nP <- length(welchPSD(ts(data[[1]][[1]][,1,1]),seglength=win,windowfun=tukeywindow,two.sided=FALSE,r=olap)$power)
	# loop through conditions, want to keep int1 and int2 separate still
	for (cc in levels(as.factor(c('int1','int2')))) {
		cat(cc,'\n')
		hemi <- list() # create empty container for hemi results
		for (hh in levels(as.factor(c('lh','rh')))) {
			cat(' ',hh,'\n')
			vox <- dim(data[[cc]][[hh]])[2] # number of voxels (well, nodes)
			subs <- dim(data[[cc]][[hh]])[3] # number of subjects
			out <- array(NA, dim=c(nP,vox,subs)) # output container
			# for each voxel (SOO SLOW, you can do better...)
			for (vv in 1:vox) {
				# for each subject
				for (ss in 1:subs) {
					temp <- ts(data[[cc]][[hh]][,vv,ss]) # turn voxel's data into timeseries object
					# find the power spectrum
					out[,vv,ss] <- welchPSD(temp,seglength=win,windowfun=tukeywindow,two.sided=FALSE,r=olap)$power
				}
			}
			hemi[[hh]] <- out
		}
		fin[[cc]] <- hemi
	}
	return(fin)
}

# function to calculate the mean PSD for each voxel (and normalize)
voxel.meanPSD <- function(data) {
	fin <- list() # create empty container for results
	# loop through conditions, want to keep int1 and int2 separate still
	for (cc in levels(as.factor(c('int1','int2')))) {
		cat(cc,'\n')
		hemi <- list() # create empty container for hemi results
		for (hh in levels(as.factor(c('lh','rh')))) {
			temp <- apply(data[[cc]][[hh]],c(1,2),mean)
			hemi[[hh]] <- apply(temp,2,function(x) x/sum(x))
		}
		fin[[cc]] <- hemi
	}
	return(fin)
}

# function to read ROIs (labels must be in order for the ROI index numbers in the filename)
read.ROIS <- function(dir,labels) {
	source(first.in.path('AFNIio.R')) # source AFNI's R functions
	fin <- list() # create empty container for results
	for (rr in 1:length(labels)) {
		cat(labels[rr],'\n')
		hemi <- list() # create empty container for hemi results
		for (hh in levels(as.factor(c('lh','rh')))) {
			filename <- paste0(dir,'/schurz_ROIs.',hh,'.lim5.',rr,'.1D') # roi filename
			# check if the file exists
			if (file.exists(filename)) {
				# read in the ROI file
				hemi[[hh]] <- read.AFNI(filename,forcedset=TRUE)$brk[,,,][,2]
			} else {
				# if it doesn't exists, then a vector of 0
				hemi[[hh]] <- vector('numeric',length=36002)
			}
		}
		fin[[labels[rr]]] <- hemi
	}
	return(fin)
}

# function to get the average PSD at the subject level across different ROIs
roi.indPSD <- function(data,rois,df) {
	labels <- names(rois) # names of ROIs
	nR <- length(labels) # number of ROIs
	nf <- dim(data[[1]][[1]])[1] # number of frequency bins
	nS <- dim(data[[1]][[1]])[3] # number of subjects
	fin <- array(NA, dim=c(nf,nR,nS)) # output container
	for (ss in 1:nS) {
		for (rr in 1:nR) {
			int1.lh <- data[['int1']]$lh[,which(rois[[rr]]$lh!=0),ss]
			int1.rh <- data[['int1']]$rh[,which(rois[[rr]]$rh!=0),ss]
			int2.lh <- data[['int2']]$lh[,which(rois[[rr]]$lh!=0),ss]
			int2.rh <- data[['int2']]$rh[,which(rois[[rr]]$rh!=0),ss]
			temp <- cbind(int1.lh,int1.rh,int2.lh,int2.rh)
			temp <- apply(temp,2,function(x) x/sum(x))
			temp <- rowMeans(temp)
			fin[,rr,ss] <- temp
		}
	}
	return(fin)
}

# function to take the mean of mean voxel PSDs across an ROI
# UGH this is messy
roi.meanPSD <- function(data,rois,df) {
	fin <- list() # create empty container for results
	labels <- names(rois)
	for (rr in 1:length(rois)) {
		cat(labels[rr],'\n')
		int1.lh <- data[['int1']]$lh[,which(rois[[rr]]$lh!=0)]
		int1.rh <- data[['int1']]$rh[,which(rois[[rr]]$rh!=0)]
		int2.lh <- data[['int2']]$lh[,which(rois[[rr]]$lh!=0)]
		int2.rh <- data[['int2']]$rh[,which(rois[[rr]]$rh!=0)]
		temp <- cbind(int1.lh,int1.rh,int2.lh,int2.rh)
		fin[[rr]] <- rowMeans(temp)
		fin[[rr]] <- fin[[rr]]/df
	}
	fin <- data.frame(fin[[1]],fin[[2]],fin[[3]],fin[[4]],fin[[5]])
	names(fin) <- names(rois)
	return(fin)
}

# function to compute a bootstrapped mean on the ROI PSD data
bootPSD <- function(data,rois,nBoot,df) {
	fin <- list() # create empty container for results
	dims <- dim(data[[1]][[1]])
	nSubj <- dims[3] # grab number of subjects from the data
	labels <- names(rois) # grab ROI names
	# loop through the ROIs
	for (rr in 1:length(rois)) {
		cat(labels[rr],'\n')
		ROIidx <- list() # container for ROI indices
		# loop through hemis to grab ROI values
		for (hh in levels(as.factor(c('lh','rh')))) {
			ROIidx[[hh]] <- which(rois[[rr]][[hh]]!=0) # find ROI nodes
		}
		# run the bootstrap trials
		bootDat <- matrix(NA,nrow=dims[1], ncol=nBoot)
		for (bb in 1:nBoot) {
			SUBidx <- sample(1:nSubj, replace=TRUE) # sample the subjects for this trial
			int1.lh <- apply(data[['int1']]$lh[,ROIidx[['lh']],SUBidx],c(1,2),mean)
			int1.rh <- apply(data[['int1']]$rh[,ROIidx[['rh']],SUBidx],c(1,2),mean)
			int2.lh <- apply(data[['int2']]$lh[,ROIidx[['lh']],SUBidx],c(1,2),mean)
			int2.rh <- apply(data[['int2']]$rh[,ROIidx[['rh']],SUBidx],c(1,2),mean)
			temp <- cbind(int1.lh,int1.rh,int2.lh,int2.rh)
			bootDat[,bb] <- rowMeans(temp)
			bootDat[,bb] <- bootDat[,bb]/sum(bootDat[,bb])
		}
		fin[[rr]] <- apply(bootDat,1,sd)
		fin[[rr]] <- fin[[rr]]/df
	}
	fin <- data.frame(fin)
	names(fin) <- names(rois)
	return(fin)
}

# # function to calculate the proportion of low-frequency power in different ROIS, given a frequency band
# lowfreq.prop <- function(data,rois,band) {
# 	fin <- list() # create empty container for results
# 	for (cc in levels(as.factor(c('int1','int2')))) {
# 		cat(cc,'\n')
# 		hemi <- list() # create empty container for hemi results
# 		for (hh in levels(as.factor(c('lh','rh')))) {
# 			cat(' ',hh,'\n')
# 			vox <- dim(data[[cc]][[hh]])[2] # number of voxels (well, nodes)
# 			subs <- dim(data[[cc]][[hh]])[3] # number of subjects
# 			out <- array(NA, dim=c(nP,vox,subs)) # output container
# 			# for each voxel (SOO SLOW, you can do better...)
# 			for (vv in 1:vox) {
# 				# for each subject
# 				for (ss in 1:subs) {
# 					temp <- ts(data[[cc]][[hh]][,vv,ss]) # turn voxel's data into timeseries object
# 					# find the power spectrum
# 					out[,vv,ss] <- welchPSD(temp,seglength=win,windowfun=tukeywindow,two.sided=FALSE,r=olap)$power
# 				}
# 			}
# 			hemi[[hh]] <- out
# 		}
# 		fin[[cc]] <- hemi
# 	}
# 	return(fin)


# }


########################### ANALYSIS ###########################

data.dir <- '/export/data/neuron/TRW'
out <- '/export/data/neuron/TRW/reprocess/psd.Rdata'

roi.dir <- '/export/data/neuron/TRW/reprocess/roi_txt/schurzROIs'
rois <- c('vis','aud','tpj','dmPFC','precun') # must be in order according to filename labels
r <- read.ROIS(roi.dir,rois)


adult <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041','RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')


setwd(out)

a.ts <- read.ts(data.dir,adult)
save(a.ts,file=paste0(out,'/a.ts.RData'))

a.glb <- rm.global(a.ts)
save(a.glb,file=paste0(out,'/a.glb.RData'))

a.voxPSD <- voxel.indPSD(a.glb,100,0.5)
save(a.voxPSD,file=paste0(out,'/a.voxPSD.RData'))

a.mPSD <- voxel.meanPSD(a.voxPSD)
save(a.mPSD,file=paste0(out,'/a.mPSD.RData'))

a.roiPSD <- roi.meanPSD(a.mPSD,r,0.01)
save(a.roiPSD,file=paste0(out,'/a.roiPSD.RData'))

a.bootPSD <- bootPSD(a.voxPSD,r,10,0.01)
save(a.bootPSD,file=paste0(out,'/a.bootPSD.RData'))


child <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173', 'RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')

c.ts <- read.ts(data.dir,child)
save(c.ts,file=paste0(out,'/c.ts.RData'))

c.glb <- rm.global(c.ts)
save(c.glb,file=paste0(out,'/c.glb.RData'))

c.voxPSD <- voxel.indPSD(c.glb,100,0.5)
save(c.voxPSD,file=paste0(out,'/c.voxPSD.RData'))

c.mPSD <- voxel.meanPSD(c.voxPSD)
save(c.mPSD,file=paste0(out,'/c.mPSD.RData'))

c.roiPSD <- roi.meanPSD(c.mPSD,r,0.01)
save(c.roiPSD,file=paste0(out,'/c.roiPSD.RData'))

c.bootPSD <- bootPSD(c.voxPSD,r,10,0.01)
save(c.bootPSD,file=paste0(out,'/c.bootPSD.RData'))

######################################

load('a.ts.RData')
load('a.glb.RData')
load('a.voxPSD.RData')
load('a.mPSD.RData')
load('a.roiPSD.RData')
load('a.bootPSD.RData')

load('c.ts.RData')
load('c.glb.RData')
load('c.voxPSD.RData')
load('c.mPSD.RData')
load('c.roiPSD.RData')
load('c.bootPSD.RData')
 

 # p <- region/sum(region)

## plot ##
a.sub <- data.frame(a.roiPSD[2:34,])
a.sub$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
a.sub$grp <- 'adult'

a.se <- data.frame(a.bootPSD[2:34,])
a.se$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
a.se$grp <- 'adult'
ma.se <- melt(a.se,id.var=c('freq','grp'))

ma <- melt(a.sub,id.var=c('freq','grp'))
ma$se <- ma.se$value

ggplot(ma, aes(freq,value,group=variable)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se), fill='grey70', color='grey70') +
	geom_line() + 
	geom_vline(xintercept = 0.04) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	ggtitle('Adult')


c.sub <- c.roiPSD[2:34,]
c.sub$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
c.sub$grp <- 'child'

c.se <- data.frame(c.bootPSD[2:34,])
c.se$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
c.se$grp <- 'child'
mc.se <- melt(c.se,id.var=c('freq','grp'))

mc <- melt(c.sub,id.var=c('freq','grp'))
mc$se <- mc.se$value


ggplot(mc, aes(freq,value,color=variable)) + 
	geom_line() + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se)) +
	geom_vline(xintercept = 0.04) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	ggtitle('Child')

d <- rbind(ma,mc)
write.csv(d, 'all.boot.se.csv',row.names=FALSE,quote=FALSE)

ggplot(d, aes(freq,value,color=grp)) + 
	geom_line() + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se)) +
	geom_vline(xintercept=0.04) + 
	scale_x_continuous(trans='log10') + 
	facet_wrap(~variable)



# test <- ts(g[['int1']]$lh[,1,1])
# test <- welchPSD(test,seglength=100,windowfun=tukeywindow,two.sided=FALSE,r=0.5)
# s$freq <- test$frequency

# ms <- melt(s,id.var='freq')



# ggplot(ms, aes(freq,value,color=variable)) + geom_line()


###############

boot <- read.csv('iCloud/TRW/second_submission/psd_analysis/all.boot.se.csv')
boot$variable <- factor(boot$variable,levels=c('dmPFC','tpj','precun','aud','vis'))

ggplot(boot, aes(freq,value,fill=grp)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),alpha=0.75) +
	geom_line() + 
	geom_vline(xintercept=0.04) + 
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) + 
	facet_wrap(~variable) +
	scale_fill_manual(values=c("darkorchid4","darkslategray4"))



########### proportion analysis

f <- seq(0,0.5,0.01)

