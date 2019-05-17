##### Libraries
library(ggplot2)
library(lme4)
library(dplyr)
library(lmerTest)

##### Function
r2z <- function(r) 0.5*(log(1+r)-log(1-r))
z2r <- function(z) (exp(2*z)-1)/(exp(2*z)+1)
dataTable <- function(dir,adult,child,condition,episode) {
	len <- length(adult)*length(child)*2
	isc.vec <- vector('numeric',len) # ISC
	sub1.vec <- vector('character',len) # child subject
	sub2.vec <- vector('character',len) # adult subject
	con.vec <- vector('character',len) # condition
	ep.vec <- vector('character',len) # episode
	m.vec <- vector('numeric',len) # child motion
	Mcomp.vec <- vector('numeric',len) # child mental comprehension
	NMcomp.vec <- vector('numeric',len) # child nonmental comprehension
	age.vec <- vector('numeric',len) # child age
	count <- 0
	for (aa in 1:length(adult)) {
		for (cc in 1:length(child)) {
			count <- count+1
			file <- paste0(dir,'/dmPFC_pairwise_stats/',adult[aa],'.',child[cc],'.',condition,'.dmPFC.rh.1D')
			isc.vec[count] <- r2z(as.numeric(read.table(file, stringsAsFactors=FALSE)[-1,-1]))
			sub1.vec[count] <- child[cc]
			sub2.vec[count] <- adult[aa]
			con.vec[count] <- condition
			ep.vec[count] <- episode
			m.vec[count] <- mot[[paste0(condition,'FD')]][grep(child[cc],mot$subj)]
			temp <- comp %>% filter(subj==child[cc] & con==condition)
			Mcomp.vec[count] <- temp$percM
			NMcomp.vec[count] <- temp$percNM
			age.vec[count] <- samp$age[grep(child[cc],samp$subj)]
		}
	}
	for (cc in 1:length(child)) {
		for (aa in 1:length(adult)) {
			count <- count+1
			file <- paste0(dir,'/dmPFC_pairwise_stats/',adult[aa],'.',child[cc],'.',condition,'.dmPFC.rh.1D')
			isc.vec[count] <- r2z(as.numeric(read.table(file, stringsAsFactors=FALSE)[-1,-1]))
			sub2.vec[count] <- child[cc]
			sub1.vec[count] <- adult[aa]
			con.vec[count] <- condition
			ep.vec[count] <- episode
			m.vec[count] <- mot[[paste0(condition,'FD')]][grep(child[cc],mot$subj)]
			temp <- comp %>% filter(subj==child[cc] & con==condition)
			Mcomp.vec[count] <- temp$percM
			NMcomp.vec[count] <- temp$percNM
			age.vec[count] <- samp$age[grep(child[cc],samp$subj)]
		}
	}
	fin <- data.frame(subj1=sub1.vec,subj2=sub2.vec,isc=isc.vec,con=con.vec,ep=ep.vec,fd=m.vec,M=Mcomp.vec,NM=NMcomp.vec,age=age.vec)
	return(fin)
}

avg.kid <- function(dat,subs) {
	len <- length(subs)*2
	res <- vector('numeric',len)
	con <- vector('character',len)
	mOVER <- vector('numeric',len)
	child <- vector('character',len)
	count <- 0
	for (cond in levels(dat$Condition)) {
		new <- dat %>% filter(Condition==cond)
		for (cc in 1:length(subs)) {
			count <- count+1
			temp <- new %>% filter(subj1==subs[cc] | subj2==subs[cc])
			res[count] <- mean(temp$resd)
			con[count] <- cond
			mOVER[count] <- temp$mOVER[1]
			child[count] <- subs[cc]
		}
	}
	fin <- data.frame(Residuals=res,Condition=con,mOVER=mOVER,Child=child)
	return(fin)
}


##### Initial objects
adultA <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')
childA <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')
adultB <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
childB <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
con <- c('int','scr')
dir <- '~/iCloud/TRW/second_submission'
# for comprehension
comp <- read.table(paste0(dir,'/comp_dat_8-22-17.txt'),header=TRUE)
# for age
samp <- read.csv(paste0(dir,'/TRW_sample_summary.csv'),header=TRUE)
# for motion
mot <- read.csv('~/iCloud/TRW/first_submission/TRW_motion.csv')
mot$intFD <- (mot$int1.mean+mot$int2.mean)/2
mot$scrFD <- (mot$scr1.mean+mot$scr2.mean)/2


##### Create dataTable
A.int <- dataTable(dir,adultA,childA,'int','BB')
A.scr <- dataTable(dir,adultA,childA,'scr','MB')
B.int <- dataTable(dir,adultB,childB,'int','MB')
B.scr <- dataTable(dir,adultB,childB,'scr','BB')

all <- rbind(A.int,A.scr,B.int,B.scr)
all$con <- factor(all$con,levels=c('scr','int'))
all$ep <- factor(all$ep,levels=c('MB','BB'))
all$subj1 <- factor(all$subj1)
all$subj2 <- factor(all$subj2)
all$mOVER <- (all$M/(all$M+all$NM))*100

##### Run model and get ready to plot
mod <- lmer(isc~con*mOVER+age+ep+fd+(1|subj1)+(1|subj2),data=all)
summary(mod)

resd <- residuals(lmer(isc~age+ep+fd+(1|subj1)+(1|subj2),data=all))
resd.dat <- data.frame(resd=z2r(resd),Condition=all$con,mOVER=all$mOVER,subj1=all$subj1,subj2=all$subj2)
plt <- avg.kid(resd.dat,c(childA,childB))

ggplot(plt, aes(mOVER,Residuals,color=Condition)) + geom_point() + stat_smooth(method='lm')




