library(reshape2)
library(lme4)
library(ggplot2)
library(plyr)
library(lmerTest)
library(arm)




A <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041','RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')
B <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040','RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
con <- c('int','scr')
dir <- '~/iCloud/TRW/second_submission/roi_txt/roi_pairwise_stats/'
rois <- c('vis','aud','tpj','dmPFC','precun')

r2z <- function(r) 0.5*(log(1+r)-log(1-r))
z2r <- function(z) (exp(2*z)-1)/(exp(2*z)+1)
iskid <- function(s) grepl('RED_TRW_1',s)
isc.val <- function(dir,sub1,sub2,con,roi) {
	n1 <- as.numeric(strsplit(sub1,'_')[[1]][3])
	n2 <- as.numeric(strsplit(sub2,'_')[[1]][3])
	if (n1>n2) {
		s1 <- sub2
		s2 <- sub1
	} else if(n2>n1) {
		s1 <- sub1
		s2 <- sub2
	}
	if (roi==1 || roi==2 || roi==3) {
		lh.file <- paste0(dir,s1,'.',s2,'.',con,'.lh.roi',roi,'.1D')
		lh <- r2z(as.numeric(read.table(lh.file, stringsAsFactors=FALSE)[-1,-1]))
		rh.file <- paste0(dir,s1,'.',s2,'.',con,'.rh.roi',roi,'.1D')
		rh <- r2z(as.numeric(read.table(rh.file, stringsAsFactors=FALSE)[-1,-1]))
		fin <- (lh+rh)/2
		
	} else {
		lh.file <- paste0(dir,s1,'.',s2,'.',con,'.lh.roi',roi,'.1D')
		fin <- r2z(as.numeric(read.table(lh.file, stringsAsFactors=FALSE)[-1,-1]))
	}
	return(fin)
}
dataTable <- function(dir,subjects,condition,episode,rois) {
	n <- length(subjects)
	N <- ((n*n)-n)* length(rois)
	isc.vec <- vector('numeric',N) # ISC
	sub1.vec <- vector('character',N) # sub 1
	sub2.vec <- vector('character',N) # sub 2
	grp.vec <- vector('character',N) # group
	con.vec <- vector('character',N) # condition
	ep.vec <- vector('character',N) # episode
	roi.vec <- vector('character',N) # rois
	count <- 0
	for (ii in 1:n) {
		for (jj in 1:n) {
			if (ii!=jj) {
				for (rr in 1:length(rois)) {
					count <- count+1 # dataTable row number
					sub1.vec[count] <- subjects[ii]
					sub2.vec[count] <- subjects[jj]
					isc.vec[count] <- isc.val(dir,subjects[ii],subjects[jj],condition,rr)
					con.vec[count] <- condition
					ep.vec[count] <- episode
					roi.vec[count] <- rois[rr]

					# figure out which group
					check <- c(iskid(subjects[ii]),iskid(subjects[jj]))
					if (sum(check)==0) {
						grp.vec[count] <- 'adult'
					} else if (sum(check)==2) {
						grp.vec[count] <- 'child'
					} else if (sum(check)==1) {
						grp.vec[count] <- 'neumat'
						# figure out which kid (for neumat values)
						if (which(check)==1) {
							k <- subjects[ii]
						} else if (which(check)==2) {
							k <- subjects[jj]
						}
					}
				}
			}
		}
	}
	fin <- data.frame(sub1=sub1.vec,sub2=sub2.vec,group=grp.vec,ep=ep.vec,con=con.vec,roi=roi.vec,isc=isc.vec)
	return(fin)
}

avg.sub <- function(data,roi) {
	s <- unique(data$sub1)
	n <- length(s)
	c <- c('int','scr')
	isc <- vector('numeric',n*length(c))
	sub <- vector('character',n*length(c))
	con <- vector('character',n*length(c))
	r <- vector('character',n*length(c))
	count <- 0
	for (ss in 1:n) {
		temp <- subset(data,data$sub1==s[ss] | data$sub2==s[ss])
		for (cc in 1:length(c)) {
			count <- count+1
			new <- subset(temp,temp$con==c[cc])
			isc[count] <- mean(new$isc)
			sub[count] <- s[ss]
			con[count] <- c[cc]
			r[count] <- roi	
		}
	}
	fin <- data.frame(sub=sub,con=con,roi=r,isc=isc)
	return(fin)
}


A.int <- dataTable(dir,A,'int','BB',rois)
A.scr <-  dataTable(dir,A,'scr','MB',rois)
B.int <- dataTable(dir,B,'int','MB',rois)
B.scr <- dataTable(dir,B,'scr','BB',rois)
fin <- rbind(A.int,A.scr,B.int,B.scr)

##############
#### Strict motion (MAD(meanFD) < 2) ####
#### Subject IDS - Adult ####
#### Monkey Bars
# subject IDs for scrambled (A)
scrA.subj <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')
# subject IDs for intact (B)
intB.subj <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
#### Body Bus
# subject IDs for scrambled (B)
scrB.subj <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
# subject IDs for intact (A)
intA.subj <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')

### grab relevant numbers for each sample ###
# Monkey Bars
# Scrambled
scrA.n <- length(scrA.subj)
scrA.N <- scrA.n*(scrA.n-1)/2
# Intact 
intB.n <- length(intB.subj)
intB.N <- intB.n*(intB.n-1)/2

# Body Bus
# Scrambled
scrB.n <- length(scrB.subj)
scrB.N <- scrB.n*(scrB.n-1)/2
# Intact 
intA.n <- length(intA.subj)
intA.N <- intA.n*(intA.n-1)/2

# Whole sample
allSubj <- c(scrA.subj,intB.subj,scrB.subj,intA.subj)
nSubj <- length(allSubj)
nFile <- (scrA.N+intB.N+scrB.N+intA.N)*2 # total number of ISC files

# ADULT
df <- nSubj-2

a <- subset(fin,fin$group=='adult')
a$con <- factor(a$con,levels=c('scr','int'))

## AUD
a.aud <- subset(a,a$roi=='aud')
a.aud.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=a.aud)))
a.aud.fin <- data.frame(isc=a.aud.mod,sub1=a.aud$sub1,sub2=a.aud$sub2,con=a.aud$con,roi='aud')
a.aud.m <- avg.sub(a.aud.fin,'aud')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=a.aud)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.aud <- 2*pt(newt,df,lower=TRUE)

## TPJ
a.tpj <- subset(a,a$roi=='tpj')
a.tpj.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=a.tpj)))
a.tpj.fin <- data.frame(isc=a.tpj.mod,sub1=a.tpj$sub1,sub2=a.tpj$sub2,con=a.tpj$con,roi='tpj')
a.tpj.m <- avg.sub(a.tpj.fin,'tpj')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=a.tpj)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.tpj <- 2*pt(newt,df,lower=FALSE)

## DMPFC
a.dmpfc <- subset(a,a$roi=='dmPFC')
a.dmpfc.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=a.dmpfc)))
a.dmpfc.fin <- data.frame(isc=a.dmpfc.mod,sub1=a.dmpfc$sub1,sub2=a.dmpfc$sub2,con=a.dmpfc$con,roi='dmpfc')
a.dmpfc.m <- avg.sub(a.dmpfc.fin,'dmpfc')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=a.dmpfc)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.dmpfc <- 2*pt(newt,df,lower=FALSE)

# PRECUN
a.precun <- subset(a,a$roi=='precun')
a.precun.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=a.precun)))
a.precun.fin <- data.frame(isc=a.precun.mod,sub1=a.precun$sub1,sub2=a.precun$sub2,con=a.precun$con,roi='precun')
a.precun.m <- avg.sub(a.precun.fin,'precun')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=a.precun)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.precun <- 2*pt(newt,df,lower=FALSE)

ps <- c(p.aud,p.tpj,p.dmpfc,p.precun)
ps.adj <- p.adjust(ps, method='bonferroni')





all <- rbind(a.aud.m,a.tpj.m,a.dmpfc.m,a.precun.m)
all$con <- factor(all$con,levels=c('scr','int'))
all$roi <- factor(all$roi,levels=c('aud','precun','dmpfc','tpj'))

all$con <- revalue(all$con, c(int='Intact',scr='Scrambled'))
all$roi <- revalue(all$roi, c(aud='Auditory',precun='Precuneus',dmpfc='dmPFC',tpj='TPJ'))

aud <- subset(all,all$roi=='Auditory')
pre <- subset(all,all$roi=='Precuneus')
dmpfc <- subset(all,all$roi=='dmPFC')
tpj <- subset(all,all$roi=='TPJ')

ggplot(tpj, aes(con,isc,fill=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_boxplot(outlier.shape=NA,width=0.5) + 
	geom_jitter(width=0.15,alpha=0.25,size=0.75) +
	scale_fill_manual(values=c('goldenrod1','tomato')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	facet_wrap(~roi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))

scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +



########
#### Strict motion (MAD(meanFD) < 2) ####
#### Subject IDS - Child ####
#### Monkey Bars
# subject IDs for scrambled (A)
scrA.subj <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')
# subject IDs for intact (B)
intB.subj <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
#### Body Bus
# subject IDs for scrambled (B)
scrB.subj <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
# subject IDs for intact (A)
intA.subj <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')

### grab relevant numbers for each sample ###
# Monkey Bars
# Scrambled
scrA.n <- length(scrA.subj)
scrA.N <- scrA.n*(scrA.n-1)/2
# Intact 
intB.n <- length(intB.subj)
intB.N <- intB.n*(intB.n-1)/2

# Body Bus
# Scrambled
scrB.n <- length(scrB.subj)
scrB.N <- scrB.n*(scrB.n-1)/2
# Intact 
intA.n <- length(intA.subj)
intA.N <- intA.n*(intA.n-1)/2

# Whole sample
allSubj <- c(scrA.subj,intB.subj,scrB.subj,intA.subj)
nSubj <- length(allSubj)
nFile <- (scrA.N+intB.N+scrB.N+intA.N)*2 # total number of ISC files

# child
df <- nSubj-2

c <- subset(fin,fin$group=='child')
c$con <- factor(c$con,levels=c('scr','int'))

# Auditory
c.aud <- subset(c,c$roi=='aud')
c.aud.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=c.aud)))
c.aud.fin <- data.frame(isc=c.aud.mod,sub1=c.aud$sub1,sub2=c.aud$sub2,con=c.aud$con,roi='aud')
c.aud.m <- avg.sub(c.aud.fin,'aud')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=c.aud)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.aud <- 2*pt(newt,df,lower=TRUE)

# TPJ
c.tpj <- subset(c,c$roi=='tpj')
c.tpj.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=c.tpj)))
c.tpj.fin <- data.frame(isc=c.tpj.mod,sub1=c.tpj$sub1,sub2=c.tpj$sub2,con=c.tpj$con,roi='tpj')
c.tpj.m <- avg.sub(c.tpj.fin,'tpj')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=c.tpj)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.tpj <- 2*pt(newt,df,lower=FALSE)

# DMPFC
c.dmpfc <- subset(c,c$roi=='dmPFC')
c.dmpfc.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=c.dmpfc)))
c.dmpfc.fin <- data.frame(isc=c.dmpfc.mod,sub1=c.dmpfc$sub1,sub2=c.dmpfc$sub2,con=c.dmpfc$con,roi='dmpfc')
c.dmpfc.m <- avg.sub(c.dmpfc.fin,'dmpfc')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=c.dmpfc)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.dmpfc <- 2*pt(newt,df,lower=FALSE)

# PRECUN
c.precun <- subset(c,c$roi=='precun')
c.precun.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=c.precun)))
c.precun.fin <- data.frame(isc=c.precun.mod,sub1=c.precun$sub1,sub2=c.precun$sub2,con=c.precun$con,roi='precun')
c.precun.m <- avg.sub(c.precun.fin,'precun')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=c.precun)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 2
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.precun <- 2*pt(newt,df,lower=FALSE)


ps <- c(p.aud,p.tpj,p.dmpfc,p.precun)
p.adjust(ps, method='bonferroni')


c.all <- rbind(c.aud.m,c.tpj.m,c.dmpfc.m,c.precun.m)
c.all$con <- factor(c.all$con,levels=c('scr','int'))
c.all$roi <- factor(c.all$roi,levels=c('aud','precun','dmpfc','tpj'))

c.all$con <- revalue(c.all$con, c(int='Intact',scr='Scrambled'))
c.all$roi <- revalue(c.all$roi, c(aud='Auditory',precun='Precuneus',dmpfc='dmPFC',tpj='TPJ'))

aud <- subset(c.all,c.all$roi=='Auditory')
pre <- subset(c.all,c.all$roi=='Precuneus')
dmpfc <- subset(c.all,c.all$roi=='dmPFC')
tpj <- subset(c.all,c.all$roi=='TPJ')

ggplot(tpj, aes(con,isc,fill=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_boxplot(outlier.shape=NA,width=0.5) + 
	geom_jitter(width=0.15,alpha=0.25,size=0.75) +
	scale_fill_manual(values=c('goldenrod1','tomato')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	facet_wrap(~roi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))

ggplot(c.all, aes(con,isc,fill=roi)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_boxplot(outlier.shape=NA,width=0.5) + 
	geom_jitter(width=0.15,alpha=0.25,size=0.75) +
	facet_wrap(~roi,nrow=1) +
	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=18),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=16),
		axis.text.y=element_text(size=12))


avg.sub.g <- function(data,roi) {
	s <- unique(data$sub1)
	n <- length(s)
	c <- c('int','scr')
	isc <- vector('numeric',n*length(c))
	sub <- vector('character',n*length(c))
	con <- vector('character',n*length(c))
	r <- vector('character',n*length(c))
	grp <- vector('character',n*length(c))
	count <- 0
	for (ss in 1:n) {
		temp <- subset(data,data$sub1==s[ss] | data$sub2==s[ss])
		for (cc in 1:length(c)) {
			count <- count+1
			new <- subset(temp,temp$con==c[cc])
			isc[count] <- mean(new$isc)
			sub[count] <- s[ss]
			con[count] <- c[cc]
			r[count] <- roi
			grp[count] <- new$grp[1]
		}
	}
	fin <- data.frame(sub=sub,con=con,roi=r,isc=isc,grp=grp)
	return(fin)
}


#########
#### Strict motion (MAD(meanFD) < 2) ####
#### Subject IDS - Adult ####
#### Monkey Bars
# subject IDs for scrambled (A)
scrA.Asubj <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')
# subject IDs for intact (B)
intB.Asubj <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
#### Body Bus
# subject IDs for scrambled (B)
scrB.Asubj <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
# subject IDs for intact (A)
intA.Asubj <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')

#### Subject IDS - Child ####
#### Monkey Bars
# subject IDs for scrambled (A)
scrA.Csubj <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')
# subject IDs for intact (B)
intB.Csubj <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
#### Body Bus
# subject IDs for scrambled (B)
scrB.Csubj <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
# subject IDs for intact (A)
intA.Csubj <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')

#### grab relevant numbers for each sample - Adult ####
#### Monkey Bars
# Scrambled
scrA.An <- length(scrA.Asubj)
scrA.AN <- scrA.An*(scrA.An-1)/2
# Intact 
intB.An <- length(intB.Asubj)
intB.AN <- intB.An*(intB.An-1)/2
#### Body Bus
# Scrambled
scrB.An <- length(scrB.Asubj)
scrB.AN <- scrB.An*(scrB.An-1)/2
# Intact 
intA.An <- length(intA.Asubj)
intA.AN <- intA.An*(intA.An-1)/2

#### grab relevant numbers for each sample - Child ####
#### Monkey Bars
# Scrambled
scrA.Cn <- length(scrA.Csubj)
scrA.CN <- scrA.Cn*(scrA.Cn-1)/2
# Intact 
intB.Cn <- length(intB.Csubj)
intB.CN <- intB.Cn*(intB.Cn-1)/2
#### Body Bus
# Scrambled
scrB.Cn <- length(scrB.Csubj)
scrB.CN <- scrB.Cn*(scrB.Cn-1)/2
# Intact 
intA.Cn <- length(intA.Csubj)
intA.CN <- intA.Cn*(intA.Cn-1)/2

# Whole sample
allSubj <- c(scrA.Asubj,intB.Asubj,scrB.Asubj,intA.Asubj,scrA.Csubj,intB.Csubj,scrB.Csubj,intA.Csubj)
nSubj <- length(allSubj)
nFile <- scrA.AN+intB.AN+scrB.AN+intA.AN+scrA.CN+intB.CN+scrB.CN+intA.CN  # total number of ISC files

# adult - child
df <- nSubj-2

## Auditory
aud <- subset(fin,fin$roi=='aud' & fin$group!='neumat')
aud.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=aud)))
aud.fin <- data.frame(isc=aud.mod,sub1=aud$sub1,sub2=aud$sub2,con=aud$con,roi='aud',grp=aud$group)
aud.m <- avg.sub.g(aud.fin,'aud')

fit <- lmer(isc~con*group+ep+(1|sub1)+(1|sub2),data=aud)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 5
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.aud <- 2*pt(newt,df,lower=TRUE)

# TPJ
tpj <- subset(fin,fin$roi=='tpj' & fin$group!='neumat')
tpj.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj)))
tpj.fin <- data.frame(isc=tpj.mod,sub1=tpj$sub1,sub2=tpj$sub2,con=tpj$con,roi='tpj',grp=tpj$group)
tpj.m <- avg.sub.g(tpj.fin,'tpj')

fit <- lmer(isc~con*group+ep+(1|sub1)+(1|sub2),data=tpj)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 5
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.tpj <- 2*pt(newt,df,lower=FALSE)

# DMPFC
dmpfc <- subset(fin,fin$roi=='dmPFC' & fin$group!='neumat')
dmpfc.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=dmpfc)))
dmpfc.fin <- data.frame(isc=dmpfc.mod,sub1=dmpfc$sub1,sub2=dmpfc$sub2,con=dmpfc$con,roi='dmpfc',grp=dmpfc$group)
dmpfc.m <- avg.sub.g(dmpfc.fin,'dmpfc')

fit <- lmer(isc~con*group+ep+(1|sub1)+(1|sub2),data=dmpfc)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 5
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.dmpfc <- 2*pt(newt,df,lower=FALSE)

# Precuneus
precun <- subset(fin,fin$roi=='precun' & fin$group!='neumat')
precun.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=precun)))
precun.fin <- data.frame(isc=precun.mod,sub1=precun$sub1,sub2=precun$sub2,con=precun$con,roi='precun',grp=precun$group)
precun.m <- avg.sub.g(precun.fin,'precun')

fit <- lmer(isc~con*group+ep+(1|sub1)+(1|sub2),data=precun)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 5
quantile(sdat[,beta], c(0.025,0.975))
b <- summary(fit)$coefficients[beta,1]
se <- summary(fit)$coefficients[beta,2]
newt <- b*sqrt((nFile)-1)/(se*sqrt((2*nFile)-1))
newt
p.precun <- 2*pt(newt,df,lower=FALSE)


ps <- c(p.aud,p.tpj,p.dmpfc,p.precun)
ps.adj <- p.adjust(ps, method='bonferroni')


all <- rbind(aud.m,tpj.m,dmpfc.m,precun.m)
all$con <- factor(all$con,levels=c('scr','int'))
all$roi <- factor(all$roi,levels=c('aud','precun','dmpfc','tpj'))

all$con <- revalue(all$con, c(int='Intact',scr='Scrambled'))
all$roi <- revalue(all$roi, c(aud='Auditory',precun='Precuneus',dmpfc='dmPFC',tpj='TPJ'))
all$grp <- revalue(all$grp, c('1'='Adult','2'='Child'))

aud <- subset(all,all$roi=='Auditory')
pre <- subset(all,all$roi=='Precuneus')
dmpfc <- subset(all,all$roi=='dmPFC')
tpj <- subset(all,all$roi=='TPJ')

ggplot(dmpfc, aes(con,isc,fill=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_boxplot(outlier.shape=NA,width=0.5) + 
	geom_jitter(width=0.15,alpha=0.25,size=0.75) +
	scale_fill_manual(values=c('goldenrod1','tomato')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	facet_wrap(~grp) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))


ggplot(all, aes(grp,isc,fill=con)) + 
	geom_boxplot(outlier.shape=NA,color='black') +
	geom_point(alpha=0.6,size=0.5,position=position_jitterdodge()) +
	facet_wrap(~roi,nrow=1) +
	theme_classic() +
	scale_fill_manual(values=c("darkorchid4","darkslategray3")) +
	labs(x='',y='') +
	theme(
		strip.background=element_blank(),
		strip.text=element_text(size=18),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=16),
		axis.text.y=element_text(size=12))


# for comprehension
comp <- read.table(paste0(dir,'/comp_dat_8-22-17.txt'),header=TRUE)
# for age
samp <- read.csv(paste0(dir,'/TRW_sample_summary.csv'),header=TRUE)
# for motion
mot <- read.csv('~/iCloud/TRW/first_submission/TRW_motion.csv')
mot$intFD <- (mot$int1.mean+mot$int2.mean)/2
mot$scrFD <- (mot$scr1.mean+mot$scr2.mean)/2


p <- c(0.00000022,0.0001,0.0000000000055,0.0000000023,0.000000023,.412,.299,.000000023,.93,.02,.00000064,0.0000000000049)


##################
#### Subject IDS - Adult ####
#### Monkey Bars
# subject IDs for scrambled (A)
scrA.Asubj <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')
# subject IDs for intact (B)
intB.Asubj <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
#### Body Bus
# subject IDs for scrambled (B)
scrB.Asubj <- c('RED_TRW_002', 'RED_TRW_005', 'RED_TRW_006', 'RED_TRW_008', 'RED_TRW_014', 'RED_TRW_024', 'RED_TRW_025', 'RED_TRW_027', 'RED_TRW_033', 'RED_TRW_034', 'RED_TRW_035', 'RED_TRW_036', 'RED_TRW_040')
# subject IDs for intact (A)
intA.Asubj <- c('RED_TRW_001', 'RED_TRW_013', 'RED_TRW_015', 'RED_TRW_017', 'RED_TRW_022', 'RED_TRW_023', 'RED_TRW_026', 'RED_TRW_029', 'RED_TRW_037', 'RED_TRW_039', 'RED_TRW_041')

#### Subject IDS - Child ####
#### Monkey Bars
# subject IDs for scrambled (A)
scrA.Csubj <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')
# subject IDs for intact (B)
intB.Csubj <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
#### Body Bus
# subject IDs for scrambled (B)
scrB.Csubj <- c('RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')
# subject IDs for intact (A)
intA.Csubj <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173')

#### grab relevant numbers for each sample ####
#### Monkey Bars
# Scrambled
scrA.An <- length(scrA.Asubj)
scrA.AN <- scrA.An*(scrA.An-1)/2
# Intact 
intB.An <- length(intB.Asubj)
intB.AN <- intB.An*(intB.An-1)/2
#### Body Bus
# Scrambled
scrB.An <- length(scrB.Asubj)
scrB.AN <- scrB.An*(scrB.An-1)/2
# Intact 
intA.An <- length(intA.Asubj)
intA.AN <- intA.An*(intA.An-1)/2

#### grab relevant numbers for each sample - Child ####
#### Monkey Bars
# Scrambled
scrA.Cn <- length(scrA.Csubj)
scrA.CN <- scrA.Cn*(scrA.Cn-1)/2
# Intact 
intB.Cn <- length(intB.Csubj)
intB.CN <- intB.Cn*(intB.Cn-1)/2
#### Body Bus
# Scrambled
scrB.Cn <- length(scrB.Csubj)
scrB.CN <- scrB.Cn*(scrB.Cn-1)/2
# Intact 
intA.Cn <- length(intA.Csubj)
intA.CN <- intA.Cn*(intA.Cn-1)/2

#### grab relevant numbers for both Child + Adult ####
#### Monkey Bars
# Scrambled
scrA.subj <- c(scrA.Csubj,scrA.Asubj)
scrA.n <- length(c(scrA.Csubj,scrA.Asubj))
scrA.N <- scrA.n*(scrA.n-1)/2
# Intact 
intB.subj <- c(intB.Csubj,intB.Asubj)
intB.n <- length(c(intB.Csubj,intB.Asubj))
intB.N <- intB.n*(intB.n-1)/2
#### Body Bus
# Scrambled
scrB.subj <- c(scrB.Csubj,scrB.Asubj)
scrB.n <- length(c(scrB.Csubj,scrB.Asubj))
scrB.N <- scrB.n*(scrB.n-1)/2
# Intact 
intA.subj <- c(intA.Csubj,intA.Asubj)
intA.n <- length(c(intA.Csubj,intA.Asubj))
intA.N <- intA.n*(intA.n-1)/2

# Whole sample
nSubj <- scrA.n+intB.n+scrB.n+intA.n # total number of subjects
nFile <- scrA.N+intB.N+scrB.N+intA.N  # total number of ISC files
NN12 <- (scrA.Cn*scrA.An)+(intB.Cn*intB.An)+(scrB.Cn*scrB.An)+(intA.Cn*intA.An)

# neumat
neu <- subset(fin,fin$group=='neumat')
child <- c('RED_TRW_101', 'RED_TRW_103', 'RED_TRW_112', 'RED_TRW_113', 'RED_TRW_114', 'RED_TRW_115', 'RED_TRW_116', 'RED_TRW_119', 'RED_TRW_137', 'RED_TRW_139', 'RED_TRW_141', 'RED_TRW_150', 'RED_TRW_152', 'RED_TRW_156', 'RED_TRW_173', 'RED_TRW_117', 'RED_TRW_118', 'RED_TRW_120', 'RED_TRW_126', 'RED_TRW_127', 'RED_TRW_128', 'RED_TRW_130', 'RED_TRW_135', 'RED_TRW_140', 'RED_TRW_143', 'RED_TRW_144', 'RED_TRW_148', 'RED_TRW_155', 'RED_TRW_175', 'RED_TRW_177', 'RED_TRW_180')


# read sample info
samp <- read.csv('/Volumes/research$/redcay/DSCN lab/Experiments/TRW/data/TRW_sample_summary.csv',header=TRUE)
samp <- subset(samp,samp$allTRW==1 & samp$include_other==1 & samp$out==1)
samp <- subset(samp,samp$group!='adult')

comp <- read.table('/Volumes/research$/redcay/DSCN lab/Experiments/TRW/data/comp_dat_8-22-17.txt',header=TRUE)
comp<- subset(comp,comp$include==1)

mot <- read.csv('~/iCloud/TRW/first_submission/TRW_motion.csv')
mot$int <- (mot$int1.mean+mot$int2.mean)/2
mot$scr <- (mot$scr1.mean+mot$scr2.mean)/2
mot <- data.frame(subj=mot$subj,int=mot$int,scr=mot$scr)
mot <- melt(mot,id.var='subj')
names(mot) <- c('subj','con','fd')

# # merge
# new <- merge(neu,samp,by=c('subj'))
# # grab important columns
# good <- new[c('subj','group','con','percM','percNM','age','sex','order','condition','iq','allTRW','episode')]

age <- vector('numeric',nrow(neu))
m <- vector('numeric',nrow(neu))
nm <- vector('numeric',nrow(neu))
fd <- vector('numeric',nrow(neu))
for (ii in 1:nrow(neu)) {
	check <- c(iskid(neu$sub1[ii]),iskid(neu$sub2[ii]))
	ind <- which(check==1)
	kid.id <- neu[ii,ind]
	age[ii] <- samp$age[grep(kid.id,samp$subj)]
	c.ind <- grep(kid.id,comp$subj)
	temp <- comp[c.ind,]
	m[ii] <- temp$percM[grep(neu$con[ii],temp$con)]
	nm[ii] <- temp$percNM[grep(neu$con[ii],temp$con)]
	c.m.ind <- grep(kid.id,mot$subj)
	temp <- mot[c.m.ind,]
	fd[ii] <- temp$fd[grep(neu$con[ii],temp$con)]

}
all.neu <- data.frame(neu,age=age,ment=m,nonment=nm,fd=fd)

all.neu <- subset(all.neu,all.neu$roi!='vis')
all.neu$roi <- factor(all.neu$roi, levels=c('aud','precun','dmPFC','tpj'))
all.neu$roi <- revalue(all.neu$roi, c(aud='Auditory',precun='Precuneus',dmPFC='dmPFC',tpj='TPJ'))
all.neu$mOVER <- all.neu$ment/(all.neu$ment+all.neu$nonment)


# create avg neumat for plots
avg.sub.neu <- function(data,s,roi) {
	n <- length(s)
	c <- c('int','scr')
	isc <- vector('numeric',n*length(c))
	sub <- vector('character',n*length(c))
	con <- vector('character',n*length(c))
	r <- vector('character',n*length(c))
	age <- vector('numeric',n*length(c))
	m <- vector('numeric',n*length(c))
	nm <- vector('numeric',n*length(c))
	mOVER <- vector('numeric',n*length(c))
	count <- 0
	for (ss in 1:n) {
		temp <- subset(data,data$sub1==s[ss] | data$sub2==s[ss])
		for (cc in 1:length(c)) {
			count <- count+1
			new <- subset(temp,temp$con==c[cc])
			isc[count] <- mean(new$isc)
			sub[count] <- s[ss]
			con[count] <- c[cc]
			r[count] <- roi
			age[count] <- new$age[1]
			m[count] <- new$m[1]
			nm[count] <- new$nm[1]
			mOVER[count] <- new$mOVER[1]
		}
	}
	fin <- data.frame(sub=sub,con=con,roi=r,isc=isc,age=age,ment=m,nonment=nm,mOVER=mOVER)
	return(fin)
}
all.neu$mOVER <- all.neu$mOVER*100
all.neu$con <- factor(all.neu$con,levels=c('scr','int'))

# effect of age
df <- nSubj-4

pre <- subset(all.neu,all.neu$roi=='Precuneus')
fit <- lmer(isc~con*age+ep+fd+(1|sub1)+(1|sub2),data=pre)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 6
quantile(sdat[,beta], c(0.025,0.975))
beta <- summary(fit)$coefficients[6,1]
se <- summary(fit)$coefficients[6,2]
newt <- beta*sqrt((2*NN12)-1)/(se*sqrt((4*NN12)-1))
2*pt(newt,df,lower=FALSE)

pre.mod <- z2r(residuals(lmer(isc~ep+fd+(1|sub1)+(1|sub2),data=pre)))
pre.fin <- data.frame(isc=pre.mod,sub1=pre$sub1,sub2=pre$sub2,con=pre$con,age=pre$age,m=pre$ment,nm=pre$nonment,mOVER=pre$mOVER)
pre.m <- avg.sub.neu(pre.fin,child,'Precuneus')

ggplot(pre.m, aes(age,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(7,8,9,10,11,12,13)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))




aud <- subset(all.neu,all.neu$roi=='Auditory')
fit <- lmer(isc~con*age+ep+fd+(1|sub1)+(1|sub2),data=aud)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 6
quantile(sdat[,beta], c(0.025,0.975))
beta <- summary(fit)$coefficients[6,1]
se <- summary(fit)$coefficients[6,2]
newt <- beta*sqrt((2*NN12)-1)/(se*sqrt((4*NN12)-1))
2*pt(newt,df,lower=FALSE)

aud.mod <- z2r(residuals(lmer(isc~ep+fd+(1|sub1)+(1|sub2),data=aud)))
aud.fin <- data.frame(isc=aud.mod,sub1=aud$sub1,sub2=aud$sub2,con=aud$con,age=aud$age,m=aud$ment,nm=aud$nonment,mOVER=aud$mOVER)
aud.m <- avg.sub.neu(aud.fin,child,'Auditory')

ggplot(aud.m, aes(age,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(7,8,9,10,11,12,13)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))

dmpfc <- subset(all.neu,all.neu$roi=='dmPFC')
fit <- lmer(isc~con*age+ep+fd+(1|sub1)+(1|sub2),data=dmpfc)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 6
quantile(sdat[,beta], c(0.025,0.975))
beta <- summary(fit)$coefficients[6,1]
se <- summary(fit)$coefficients[6,2]
newt <- beta*sqrt((2*NN12)-1)/(se*sqrt((4*NN12)-1))
2*pt(newt,df,lower=FALSE)


dmpfc.mod <- z2r(residuals(lmer(isc~ep+fd+(1|sub1)+(1|sub2),data=dmpfc)))
dmpfc.fin <- data.frame(isc=dmpfc.mod,sub1=dmpfc$sub1,sub2=dmpfc$sub2,con=dmpfc$con,age=dmpfc$age,m=dmpfc$ment,nm=dmpfc$nonment,mOVER=dmpfc$mOVER)
dmpfc.m <- avg.sub.neu(dmpfc.fin,child,'dmpfc')

ggplot(dmpfc.m, aes(age,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(7,8,9,10,11,12,13)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))


tpj <- subset(all.neu,all.neu$roi=='TPJ')
fit <- lmer(isc~con*age+ep+fd+(1|sub1)+(1|sub2),data=tpj)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 6
quantile(sdat[,beta], c(0.025,0.975))

tpj.mod <- z2r(residuals(lmer(isc~ep+fd+(1|sub1)+(1|sub2),data=tpj)))
tpj.fin <- data.frame(isc=tpj.mod,sub1=tpj$sub1,sub2=tpj$sub2,con=tpj$con,age=tpj$age,m=tpj$ment,nm=tpj$nonment,mOVER=tpj$mOVER)
tpj.m <- avg.sub.neu(tpj.fin,child,'tpj')

ggplot(tpj.m, aes(age,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(7,8,9,10,11,12,13)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))



# effect of mOVER
pre <- subset(all.neu,all.neu$roi=='Precuneus')
fit <- lmer(isc~con*mOVER+age+ep+fd+(1|sub1)+(1|sub2),data=pre)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 7
quantile(sdat[,beta], c(0.025,0.975))

pre.mod <- z2r(residuals(lmer(isc~age+ep+fd+(1|sub1)+(1|sub2),data=pre)))
pre.fin <- data.frame(isc=pre.mod,sub1=pre$sub1,sub2=pre$sub2,con=pre$con,age=pre$age,m=pre$ment,nm=pre$nonment,mOVER=pre$mOVER)
pre.m <- avg.sub.neu(pre.fin,child,'Precuneus')

ggplot(pre.m, aes(mOVER,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(10,20,30,40,50,60,70)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))




aud <- subset(all.neu,all.neu$roi=='Auditory')
fit <- lmer(isc~con*mOVER+age+ep+fd+(1|sub1)+(1|sub2),data=aud)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 7
quantile(sdat[,beta], c(0.025,0.975))

aud.mod <- z2r(residuals(lmer(isc~age+ep+fd+(1|sub1)+(1|sub2),data=aud)))
aud.fin <- data.frame(isc=aud.mod,sub1=aud$sub1,sub2=aud$sub2,con=aud$con,age=aud$age,m=aud$ment,nm=aud$nonment,mOVER=aud$mOVER)
aud.m <- avg.sub.neu(aud.fin,child,'Auditory')

ggplot(aud.m, aes(mOVER,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(10,20,30,40,50,60,70)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))

dmpfc <- subset(all.neu,all.neu$roi=='dmPFC')
fit <- lmer(isc~con*mOVER+age+ep+fd+(1|sub1)+(1|sub2),data=dmpfc)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 7
quantile(sdat[,beta], c(0.025,0.975))

dmpfc.mod <- z2r(residuals(lmer(isc~age+ep+fd+(1|sub1)+(1|sub2),data=dmpfc)))
dmpfc.fin <- data.frame(isc=dmpfc.mod,sub1=dmpfc$sub1,sub2=dmpfc$sub2,con=dmpfc$con,age=dmpfc$age,m=dmpfc$ment,nm=dmpfc$nonment,mOVER=dmpfc$mOVER)
dmpfc.m <- avg.sub.neu(dmpfc.fin,child,'dmpfc')

ggplot(dmpfc.m, aes(mOVER,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(10,20,30,40,50,60,70)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))


tpj <- subset(all.neu,all.neu$roi=='TPJ')
fit <- lmer(isc~con*mOVER+age+ep+fd+(1|sub1)+(1|sub2),data=tpj)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 7
quantile(sdat[,beta], c(0.025,0.975))

tpj.mod <- z2r(residuals(lmer(isc~age+ep+fd+(1|sub1)+(1|sub2),data=tpj)))
tpj.fin <- data.frame(isc=tpj.mod,sub1=tpj$sub1,sub2=tpj$sub2,con=tpj$con,age=tpj$age,m=tpj$ment,nm=tpj$nonment,mOVER=tpj$mOVER)
tpj.m <- avg.sub.neu(tpj.fin,child,'tpj')

ggplot(tpj.m, aes(mOVER,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(10,20,30,40,50,60,70)) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))






