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
			fin[,rr,ss] <- temp/df
		}
	}
	return(fin)
}

roi.indPSD.sep <- function(data,rois,df) {
	labels <- names(rois) # names of ROIs
	nR <- 2 # number of ROIs (only one ROI at a time, left and right)
	nf <- dim(data[[1]][[1]])[1] # number of frequency bins
	nS <- dim(data[[1]][[1]])[3] # number of subjects
	fin <- array(NA, dim=c(nf,nR,nS)) # output container
	for (ss in 1:nS) {
		int1.lh <- data[['int1']]$lh[,which(rois[[1]]$lh!=0),ss]
		int2.lh <- data[['int2']]$lh[,which(rois[[1]]$lh!=0),ss]
		temp <- cbind(int1.lh,int2.lh)
		temp <- apply(temp,2,function(x) x/sum(x))
		temp <- rowMeans(temp)
		fin[,1,ss] <- temp/df

		int1.rh <- data[['int1']]$rh[,which(rois[[1]]$rh!=0),ss]
		int2.rh <- data[['int2']]$rh[,which(rois[[1]]$rh!=0),ss]
		temp <- cbind(int1.rh,int2.rh)
		temp <- apply(temp,2,function(x) x/sum(x))
		temp <- rowMeans(temp)
		fin[,2,ss] <- temp/df
	}
	return(fin)
}                  

# function calculate proportion of PSD in a frequency band
roi.band.indPSD <- function(data,rois,df,band) {
	labels <- names(rois)
	# find the indices of the frequency band
	f <- seq(0,0.5,df)
	f.min <- which(f==band[1])
	f.max <- which(f==band[2])
	temp <- data[f.min:f.max,,]
	temp <- apply(temp,c(2,3),sum)
	fin <- data.frame(t(temp))
	names(fin) <- labels
	return(fin)
}

# function calculate proportion of PSD in a frequency band
roi.band.indPSD.sep <- function(data,rois,df,band) {
	labels <- c('lh','rh')
	# find the indices of the frequency band
	f <- seq(0,0.5,df)
	f.min <- which(f==band[1])
	f.max <- which(f==band[2])
	temp <- data[f.min:f.max,,]
	temp <- apply(temp,c(2,3),sum)
	fin <- data.frame(t(temp))
	names(fin) <- labels
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

# function to take the mean of mean voxel PSDs across an ROI
# UGH this is messy
roi.meanPSD.sep <- function(data,rois,df) {
	temp.lh <- list()
	temp.rh <- list()
	labels <- names(rois)
	for (rr in 1:length(rois)) {
		cat(labels[rr],'\n')
		int1.lh <- data[['int1']]$lh[,which(rois[[rr]]$lh!=0)]
		int1.rh <- data[['int1']]$rh[,which(rois[[rr]]$rh!=0)]
		int2.lh <- data[['int2']]$lh[,which(rois[[rr]]$lh!=0)]
		int2.rh <- data[['int2']]$rh[,which(rois[[rr]]$rh!=0)]

		temp1 <- cbind(int1.lh,int2.lh)
		temp.lh[[rr]] <- rowMeans(temp1)
		temp.lh[[rr]] <- temp.lh[[rr]]/df

		temp2 <- cbind(int1.rh,int2.rh)
		temp.rh[[rr]] <- rowMeans(temp2)
		temp.rh[[rr]] <- temp.rh[[rr]]/df
	}
	names(temp.lh) <- names(rois)
	names(temp.rh) <- names(rois)
	fin <- list('lh'=temp.lh, 'rh'=temp.rh)
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
			#temp <- cbind(int1.lh,int1.rh,int2.lh,int2.rh)
			temp <- cbind(int1.lh,int1.rh,int2.lh)
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

# function to compute a bootstrapped mean on the ROI PSD data
bootPSD.sep <- function(data,rois,nBoot,df) {
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
		bootDat.lh <- matrix(NA,nrow=dims[1], ncol=nBoot)
		bootDat.rh <- matrix(NA,nrow=dims[1], ncol=nBoot)
		for (bb in 1:nBoot) {
			SUBidx <- sample(1:nSubj, replace=TRUE) # sample the subjects for this trial

			int1.lh <- apply(data[['int1']]$lh[,ROIidx[['lh']],SUBidx],c(1,2),mean)
			int2.lh <- apply(data[['int2']]$lh[,ROIidx[['lh']],SUBidx],c(1,2),mean)
			temp.lh <- cbind(int1.lh,int2.lh)
			bootDat.lh[,bb] <- rowMeans(temp.lh)
			bootDat.lh[,bb] <- bootDat.lh[,bb]/sum(bootDat.lh[,bb])

			int1.rh <- apply(data[['int1']]$rh[,ROIidx[['rh']],SUBidx],c(1,2),mean)
			int2.rh <- apply(data[['int2']]$rh[,ROIidx[['rh']],SUBidx],c(1,2),mean)
			temp.rh <- cbind(int1.rh,int2.rh)
			bootDat.rh[,bb] <- rowMeans(temp.rh)
			bootDat.rh[,bb] <- bootDat.rh[,bb]/sum(bootDat.rh[,bb])
		}
		fin.lh <- apply(bootDat.lh,1,sd)
		fin.lh <- fin.lh/df

		fin.rh <- apply(bootDat.rh,1,sd)
		fin.rh <- fin.rh/df
	}

	fin <- list('lh'=fin.lh, 'rh'=fin.rh)
	return(fin)
}


########################### ANALYSIS ###########################

data.dir <- '/export/data/neuron/TRW'
out <- '/export/data/neuron/TRW/reprocess/psd.Rdata'
age <- read.table('/export/data/neuron/TRW/scripts/trw_age_final_sample.txt', header=TRUE)
mot <- read.csv('/export/data/neuron/TRW/scripts/TRW_final_sample_motion.csv')
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

a.roi.indPSD <- roi.indPSD(a.voxPSD,r,0.01)
save(a.roi.indPSD,file=paste0(out,'/a.roi.indPSD.RData'))

a.roi.bandPSD <- roi.band.indPSD(a.roi.indPSD,r,0.01,c(0.01,0.04))
save(a.roi.bandPSD,file=paste0(out,'/a.roi.bandPSD.RData'))

a.roiPSD <- roi.meanPSD(a.mPSD,r,0.01)
save(a.roiPSD,file=paste0(out,'/a.roiPSD.RData'))

a.bootPSD <- bootPSD(a.voxPSD,r,10,0.01)
save(a.bootPSD,file=paste0(out,'/a.bootPSD.RData'))


#child <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173', 'RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
child <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_117', 'RED_TRW_118', 'RED_TRW_119', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_140', 'RED_TRW_141', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_155', 'RED_TRW_156', 'RED_TRW_173', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')

c.ts <- read.ts(data.dir,child)
save(c.ts,file=paste0(out,'/c.ts.RData'))

c.glb <- rm.global(c.ts)
save(c.glb,file=paste0(out,'/c.glb.RData'))

c.voxPSD <- voxel.indPSD(c.glb,100,0.5)
save(c.voxPSD,file=paste0(out,'/c.voxPSD.RData'))

c.mPSD <- voxel.meanPSD(c.voxPSD)
save(c.mPSD,file=paste0(out,'/c.mPSD.RData'))

c.roi.indPSD <- roi.indPSD(c.voxPSD,r,0.01)
save(c.roi.indPSD,file=paste0(out,'/c.roi.indPSD.RData'))

c.roi.bandPSD <- roi.band.indPSD(c.roi.indPSD,r,0.01,c(0.01,0.04))
save(c.roi.bandPSD,file=paste0(out,'/c.roi.bandPSD.RData'))

c.roiPSD <- roi.meanPSD(c.mPSD,r,0.01)
save(c.roiPSD,file=paste0(out,'/c.roiPSD.RData'))

c.bootPSD <- bootPSD(c.voxPSD,r,10,0.01)
save(c.bootPSD,file=paste0(out,'/c.bootPSD.RData'))

######################################

load('a.ts.RData')
load('a.glb.RData')
load('a.voxPSD.RData')
load('a.mPSD.RData')
load('a.roi.indPSD.RData')
load('a.roi.bandPSD.RData')
load('a.roiPSD.RData')
load('a.bootPSD.RData')

load('c.ts.RData')
load('c.glb.RData')
load('c.voxPSD.RData')
load('c.mPSD.RData')
load('c.roi.indPSD.RData')
load('c.roi.bandPSD.RData')
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
library(lme4)
a.band <- a.roi.bandPSD
a.band$group <- 'adult'
a.band$subject <- adult

c.band <- c.roi.bandPSD
c.band$group <- 'child'
c.band$subject <- child

both <- rbind(a.band,c.band)
int.mot <- subset(mot,mot$con=='intFD')
new.m <- merge(both,int.mot,by=c('subject','group'))


m <- melt(new.m,id.vars=c('subject','group','meanFD','con','age'))
names(m) <- c('subject','group','meanFD','con','age','roi','alpha')
m <- subset(m,m$roi!='vis')



write.csv(m,'all.band.prop.csv',row.names=FALSE,quote=FALSE)

ggplot(m,aes(group,alpha,color=roi)) + geom_point() + facet_wrap(~roi)


m.adult <- subset(m,m$group=='adult')
pairwise.t.test(m.adult$alpha,m.adult$roi,p.adjust='bonferroni')

m.child <- subset(m,m$group=='child')
pairwise.t.test(m.child$alpha,m.child$roi,p.adjust='bonferroni')

# group comparisons
tpj <- subset(m,m$roi=='tpj')
t.test(alpha~group,data=tpj)
mod <- lm(alpha ~ group + meanFD, data = tpj)
summary(mod)

dmpfc <- subset(m,m$roi=='dmPFC')
t.test(alpha~group,data=dmpfc)
mod <- lm(alpha ~ group + meanFD, data = dmpfc)
summary(mod)

precun <- subset(m,m$roi=='precun')
t.test(alpha~group,data=precun)
mod <- lm(alpha ~ group + meanFD, data = precun)
summary(mod)

aud <- subset(m,m$roi=='aud')
t.test(alpha~group,data=aud)
mod <- lm(alpha ~ group + meanFD, data = aud)
summary(mod)



p <- c(0.0072,0.073,0.087,0.2112)
p.adjust(p,method='bonferroni')


#### individual differences in age
mot.dat <- subset(mot, mot$group=='child' & mot$con=='intFD')
ind.dat <- data.frame(c.roi.bandPSD,mot.dat)
m.dat <- melt(ind.dat, id.vars=c('subject','age','group','con','meanFD'))
m.dat <- subset(m.dat,m.dat$variable!='vis')

ggplot(m.dat, aes(age,value, color=variable)) + geom_point() + geom_smooth(method='lm') + facet_wrap(~variable)


aud <- lm(aud ~ age+meanFD,data=ind.dat)
summary(aud)

tpj <- lm(tpj ~ age+meanFD,data=ind.dat)
summary(tpj)

dmpfc <- lm(dmPFC ~ age+meanFD,data=ind.dat)
summary(dmpfc)

precun <- lm(precun ~ age+meanFD,data=ind.dat)
summary(precun)

p <- c(0.04,0.73,0.99,0.35)
p.adjust(p,method='bonferroni')

write.csv(m.dat, 'age.lowfreq.power.csv',row.names=FALSE,quote=FALSE)



################ Separate left and right TPJ

#### Left and right TPJ by group
# adult
a.sep.mPSD <- roi.meanPSD.sep(a.mPSD,r,0.01)
a.tpj.mPSD <- data.frame(a.sep.mPSD[['lh']]['tpj'],a.sep.mPSD[['rh']]['tpj'])
names(a.tpj.mPSD) <- c('lh','rh')

# child
c.sep.mPSD <- roi.meanPSD.sep(c.mPSD,r,0.01)
c.tpj.mPSD <- data.frame(c.sep.mPSD[['lh']]['tpj'],c.sep.mPSD[['rh']]['tpj'])
names(c.tpj.mPSD) <- c('lh','rh')

# grab only the tpj ROIS
tpjs <- list(tpj=r[[c('tpj')]])

# do the bootsrap
a.tpj.bootPSD <- bootPSD.sep(a.voxPSD,tpjs,10,0.01)
a.tpj.bootPSD <- data.frame(a.tpj.bootPSD$lh, a.tpj.bootPSD$rh)
names(a.tpj.bootPSD) <- c('lh','rh')

c.tpj.bootPSD <- bootPSD.sep(c.voxPSD,tpjs,10,0.01)
c.tpj.bootPSD <- data.frame(c.tpj.bootPSD$lh, c.tpj.bootPSD$rh)
names(c.tpj.bootPSD) <- c('lh','rh')

## put the tpj data together for plotting
# adult
## plot ##
a.sub <- data.frame(a.tpj.mPSD[2:34,])
a.sub$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
a.sub$grp <- 'adult'

a.se <- data.frame(a.tpj.bootPSD[2:34,])
a.se$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
a.se$grp <- 'adult'
ma.se <- melt(a.se,id.var=c('freq','grp'))

ma <- melt(a.sub,id.var=c('freq','grp'))
ma$se <- ma.se$value

# child
## plot ##
c.sub <- data.frame(c.tpj.mPSD[2:34,])
c.sub$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
c.sub$grp <- 'child'

c.se <- data.frame(c.tpj.bootPSD[2:34,])
c.se$freq <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33)
c.se$grp <- 'child'
mc.se <- melt(c.se,id.var=c('freq','grp'))

mc <- melt(c.sub,id.var=c('freq','grp'))
mc$se <- ma.se$value

d <- rbind(ma,mc)
write.csv(d, 'tpj.boot.se.csv',row.names=FALSE,quote=FALSE)

ggplot(d, aes(freq,value,fill=grp)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04,size=1) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	facet_wrap(~variable) +
	theme_classic() +
	scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10)) 

	+
	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14))


### grab the band from the separate TPJs
# adult
a.tpj.indPSD <- roi.indPSD.sep(a.voxPSD,tpjs,0.01)
a.tpj.bandPSD <- roi.band.indPSD.sep(a.tpj.indPSD,tpjs,0.01,c(0.01,0.04))

# child
c.tpj.indPSD <- roi.indPSD.sep(c.voxPSD,tpjs,0.01)
c.tpj.bandPSD <- roi.band.indPSD.sep(c.tpj.indPSD,tpjs,0.01,c(0.01,0.04))


# put together with age
ind.dat <- data.frame(c.tpj.bandPSD,mot.dat)
m.dat <- melt(ind.dat, id.vars=c('subject','age','group','con','meanFD'))





ggplot(m.dat, aes(age,value, color=variable)) + geom_point() + geom_smooth(method='lm') + facet_wrap(~variable)

lh <- subset(m.dat, m.dat$variable=='lh')
lh <- lm(value ~ age+meanFD,data=lh)
summary(lh)

rh <- subset(m.dat, m.dat$variable=='rh')
rh <- lm(value ~ age+meanFD,data=rh)
summary(rh)

write.csv(m.dat, 'tpj.age.lowfreq.power.csv',row.names=FALSE,quote=FALSE)





library(lme4)
a.band <- a.tpj.bandPSD
a.band$group <- 'adult'
a.band$subject <- adult

c.band <- c.tpj.bandPSD
c.band$group <- 'child'
c.band$subject <- child

both <- rbind(a.band,c.band)
int.mot <- subset(mot,mot$con=='intFD')
new.m <- merge(both,int.mot,by=c('subject','group'))


m <- melt(new.m,id.vars=c('subject','group','meanFD','con','age'))
names(m) <- c('subject','group','meanFD','con','age','roi','alpha')
=


write.csv(m,'tpj.band.prop.csv',row.names=FALSE,quote=FALSE)

ggplot(m,aes(group,alpha,color=roi)) + geom_point() + facet_wrap(~roi)


m.adult <- subset(m,m$group=='adult')
pairwise.t.test(m.adult$alpha,m.adult$roi,p.adjust='bonferroni')

m.child <- subset(m,m$group=='child')
pairwise.t.test(m.child$alpha,m.child$roi,p.adjust='bonferroni')

# group comparisons
lh <- subset(m,m$roi=='lh')
t.test(alpha~group,data=lh)
mod <- lm(alpha ~ group + meanFD, data = lh)
summary(mod)

rh <- subset(m,m$roi=='rh')
t.test(alpha~group,data=rh)
mod <- lm(alpha ~ group + meanFD, data = rh)
summary(mod)



