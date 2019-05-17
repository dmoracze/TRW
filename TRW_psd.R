library(ggplot2)
library(bspec)
dir <- '~/iCloud/TRW/second_submission/roi_timeseries/'

# read in data
dat <- read.table(paste0(dir,'RED_TRW_001.int1.ROI5.lh.1D'),stringsAsFactors=FALSE)[-1,-1]
dat <- unlist(lapply(dat,as.numeric))
dat <- data.frame(ts=dat)

# plot time series
ggplot(dat, aes(1:271,ts)) + geom_line()

dat.ts <- ts(dat$ts)
w <- welchPSD(dat.ts, seglength=100, windowfun=tukeywindow, r=0.5)

fin <- data.frame(freq=w$frequency[2:34],power=w$power[2:34])

ggplot(fin, aes(freq,power)) + geom_line() + scale_x_continuous(trans='log10')



read.ts <- function(dir,subj,run,roi,hemi) {
	dat <- read.table(paste0(dir,subj,'.',run,'.',roi,'.',hemi,'.1D'),stringsAsFactors=FALSE)[-1,-1]
	dat <- unlist(lapply(dat,as.numeric))
	return(dat)
}

ts.mat <- function(dir,subjs,roi,hemi) {
	n <- length(read.ts(dir,subjs[1],'int1',roi,hemi))
	mat1 <- matrix(NA,nrow=n, ncol=length(subjs))
	mat2 <- matrix(NA,nrow=n, ncol=length(subjs))
	for (ss in 1:length(subjs)) {
		mat1[,ss] <- scale(read.ts(dir,subjs[ss],'int1',roi,hemi))
		mat2[,ss] <- scale(read.ts(dir,subjs[ss],'int2',roi,hemi))
		dat <- list(int1=mat1,int2=mat2)
	}
	fin <- list(data=dat,roi=roi,hemi=hemi)
	return(fin)
}

psd <- function(data) {
	# mean, sd, power matrix, freq
	nT <- length(data$data$int1[,1])
	nS <- length(data$data$int1[1,])
	mat <- matrix(NA,nrow=33,ncol=nS)
	for (ss in 1:nS) {
		w1 <- welchPSD(ts(data$data$int1[,ss]),seglength=100,windowfun=tukeywindow,r=0.5)
		w2 <- welchPSD(ts(data$data$int2[,ss]),seglength=100,windowfun=tukeywindow,r=0.5)
		mat[,ss] <- rowMeans(data.frame(w1=w1$power,w2=w2$power))[2:34]
	}
	m <- rowMeans(mat)
	se <- apply(mat,1,function(x) sd(x)/length(x))
	freq <- w1$frequency[2:34]
	fin <- list(summary=data.frame(mean=m,se=se,freq=freq,roi=data$roi,hemi=data$hemi),mat=mat)







	return(fin)
}



new <- rbind(v1.lh.a$summary,v1.rh.a$summary)

ggplot(new, aes(freq,mean,color=hemi)) + geom_line() + scale_x_continuous(trans='log10')



adults <- as.factor(c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041', 'RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040'))
kids <- as.factor(c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173', 'RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180'))

v1.lh.a <- psd(ts.mat(dir,adults,'ROI1','lh'))
v1.lh.c <- psd(ts.mat(dir,kids,'ROI1','lh'))

v1.rh.a <- psd(ts.mat(dir,adults,'ROI1','rh'))
v1.rh.c <- psd(ts.mat(dir,kids,'ROI1','rh'))

a1.lh.a <- psd(ts.mat(dir,adults,'ROI2','lh'))
a1.lh.c <- psd(ts.mat(dir,kids,'ROI2','lh'))

a1.rh.a <- psd(ts.mat(dir,adults,'ROI2','rh'))
a1.rh.c <- psd(ts.mat(dir,kids,'ROI2','rh'))

tpj.lh.a <- psd(ts.mat(dir,adults,'ROI3','lh'))
tpj.lh.c <- psd(ts.mat(dir,kids,'ROI3','lh'))

tpj.rh.a <- psd(ts.mat(dir,adults,'ROI3','rh'))
tpj.rh.c <- psd(ts.mat(dir,kids,'ROI3','rh'))

dmpfc.lh.a <- psd(ts.mat(dir,adults,'ROI4','lh'))
dmpfc.lh.c <- psd(ts.mat(dir,kids,'ROI4','lh'))

pre.lh.a <- psd(ts.mat(dir,adults,'ROI5','lh'))
pre.lh.c <- psd(ts.mat(dir,kids,'ROI5','lh'))















