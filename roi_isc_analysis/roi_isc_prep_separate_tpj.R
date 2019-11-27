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
rois <- c('tpj')

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
	lh.file <- paste0(dir,s1,'.',s2,'.',con,'.lh.roi',roi,'.1D')
	lh <- r2z(as.numeric(read.table(lh.file, stringsAsFactors=FALSE)[-1,-1]))
	rh.file <- paste0(dir,s1,'.',s2,'.',con,'.rh.roi',roi,'.1D')
	rh <- r2z(as.numeric(read.table(rh.file, stringsAsFactors=FALSE)[-1,-1]))
	fin <- data.frame(lh=lh, rh=rh)
	return(fin)
}

test <- isc.val(dir,'RED_TRW_001','RED_TRW_015','int',3)

dataTable <- function(dir,subjects,condition,episode,rois) {
	n <- length(subjects)
	N <- ((n*n)-n)* length(rois)
	isc.vec <- vector('numeric',N) # ISC
	hemi.vec <- vector('character',N) # hemisphere
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
				# get isc values for lh and rh
				isc <- isc.val(dir,subjects[ii],subjects[jj],condition,3)

				# figure out which group
				check <- c(iskid(subjects[ii]),iskid(subjects[jj]))
				if (sum(check)==0) {
					grp <- 'adult'
				} else if (sum(check)==2) {
					grp <- 'child'
				} else if (sum(check)==1) {
					grp <- 'neumat'
					# figure out which kid (for neumat values)
					if (which(check)==1) {
						k <- subjects[ii]
					} else if (which(check)==2) {
						k <- subjects[jj]
					}
				}

				# lh
				count <- count+1 # dataTable row number
				sub1.vec[count] <- subjects[ii]
				sub2.vec[count] <- subjects[jj]
				isc.vec[count] <- isc$lh
				hemi.vec[count] <- 'lh'
				con.vec[count] <- condition
				ep.vec[count] <- episode
				roi.vec[count] <- rois[1]
				grp.vec[count] <- grp

				# rh
				count <- count+1 # dataTable row number
				sub1.vec[count] <- subjects[ii]
				sub2.vec[count] <- subjects[jj]
				isc.vec[count] <- isc$rh
				hemi.vec[count] <- 'rh'
				con.vec[count] <- condition
				ep.vec[count] <- episode
				roi.vec[count] <- rois[1]
				grp.vec[count] <- grp
			}
		}
	}
	fin <- data.frame(sub1=sub1.vec,sub2=sub2.vec,group=grp.vec,ep=ep.vec,con=con.vec,roi=roi.vec,hemi=hemi.vec,isc=isc.vec)
	return(fin)
}

avg.sub <- function(data,roi,hemi) {
	s <- unique(data$sub1)
	n <- length(s)
	c <- c('int','scr')
	isc <- vector('numeric',n*length(c))
	sub <- vector('character',n*length(c))
	con <- vector('character',n*length(c))
	r <- vector('character',n*length(c))
	h <- vector('character',n*length(c))
	count <- 0
	for (ss in 1:n) {
		temp <- subset(data,data$sub1==s[ss] | data$sub2==s[ss])
		for (cc in 1:length(c)) {
			count <- count+1
			new <- subset(temp,temp$con==c[cc])
			isc[count] <- mean(new$isc)
			sub[count] <- as.character(s[ss])
			con[count] <- c[cc]
			r[count] <- roi	
			h[count] <- hemi
		}
	}
	fin <- data.frame(sub=sub,con=con,roi=r,hemi=h,isc=isc)
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


## left TPJ
tpj.lh <- subset(a,a$hemi=='lh')
tpj.lh.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj.lh)))
tpj.lh.fin <- data.frame(isc=tpj.lh.mod,sub1=tpj.lh$sub1,sub2=tpj.lh$sub2,con=tpj.lh$con,roi='tpj')
tpj.lh.m <- avg.sub(tpj.lh.fin,'tpj','lh')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=tpj.lh)
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
p.lh <- 2*pt(newt,df,lower=FALSE)


## right TPJ
tpj.rh <- subset(a,a$hemi=='rh')
tpj.rh.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj.rh)))
tpj.rh.fin <- data.frame(isc=tpj.rh.mod,sub1=tpj.rh$sub1,sub2=tpj.rh$sub2,con=tpj.rh$con,roi='tpj',hemi='lh')
tpj.rh.m <- avg.sub(tpj.rh.fin,'tpj','rh')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=tpj.rh)
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
p.rh <- 2*pt(newt,df,lower=FALSE)

ps <- c(p.lh,p.rh)
p.adjust(ps, method='bonferroni')


all <- rbind(tpj.lh.m,tpj.rh.m)
all$con <- factor(all$con,levels=c('scr','int'))

all$con <- revalue(all$con, c(int='Intact',scr='Scrambled'))
all$roi <- revalue(all$roi, c(tpj='TPJ'))
all$hemi <- revalue(all$hemi, c(lh='Left', rh='Right'))

left <- subset(all,all$hemi=='Left')
right <- subset(all,all$hemi=='Right')

ggplot(right, aes(con,isc,fill=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_boxplot(outlier.shape=NA,width=0.5) + 
	geom_jitter(width=0.15,alpha=0.25,size=0.75) +
	scale_fill_manual(values=c('goldenrod1','tomato')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	facet_wrap(~hemi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))




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

## left TPJ
tpj.lh <- subset(c,c$hemi=='lh')
tpj.lh.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj.lh)))
tpj.lh.fin <- data.frame(isc=tpj.lh.mod,sub1=tpj.lh$sub1,sub2=tpj.lh$sub2,con=tpj.lh$con,roi='tpj')
tpj.lh.m <- avg.sub(tpj.lh.fin,'tpj','lh')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=tpj.lh)
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
p.lh <- 2*pt(newt,df,lower=FALSE)


## right TPJ
tpj.rh <- subset(c,c$hemi=='rh')
tpj.rh.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj.rh)))
tpj.rh.fin <- data.frame(isc=tpj.rh.mod,sub1=tpj.rh$sub1,sub2=tpj.rh$sub2,con=tpj.rh$con,roi='tpj',hemi='lh')
tpj.rh.m <- avg.sub(tpj.rh.fin,'tpj','rh')

fit <- lmer(isc~con+ep+(1|sub1)+(1|sub2),data=tpj.rh)
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
p.rh <- 2*pt(newt,df,lower=FALSE)

ps <- c(p.lh,p.rh)
p.adjust(ps, method='bonferroni')


all <- rbind(tpj.lh.m,tpj.rh.m)
all$con <- factor(all$con,levels=c('scr','int'))

all$con <- revalue(all$con, c(int='Intact',scr='Scrambled'))
all$roi <- revalue(all$roi, c(tpj='TPJ'))
all$hemi <- revalue(all$hemi, c(lh='Left', rh='Right'))

left <- subset(all,all$hemi=='Left')
right <- subset(all,all$hemi=='Right')

ggplot(left, aes(con,isc,fill=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_boxplot(outlier.shape=NA,width=0.5) + 
	geom_jitter(width=0.15,alpha=0.25,size=0.75) +
	scale_fill_manual(values=c('goldenrod1','tomato')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	facet_wrap(~hemi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.ticks=element_blank(),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))




avg.sub.g <- function(data,roi,hemi) {
	s <- unique(data$sub1)
	n <- length(s)
	c <- c('int','scr')
	isc <- vector('numeric',n*length(c))
	sub <- vector('character',n*length(c))
	con <- vector('character',n*length(c))
	r <- vector('character',n*length(c))
	h <- vector('character',n*length(c))
	grp <- vector('character',n*length(c))
	count <- 0
	for (ss in 1:n) {
		temp <- subset(data,data$sub1==s[ss] | data$sub2==s[ss])
		for (cc in 1:length(c)) {
			count <- count+1
			new <- subset(temp,temp$con==c[cc])
			isc[count] <- mean(new$isc)
			sub[count] <- as.character(s[ss])
			con[count] <- c[cc]
			r[count] <- roi
			h[count] <- hemi
			grp[count] <- as.character(new$grp[1])
		}
	}
	fin <- data.frame(sub=sub,con=con,roi=r,grp=grp,hemi=h,isc=isc)
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

# left TPJ
tpj.lh <- subset(fin,fin$hemi=='lh' & fin$group!='neumat')
tpj.lh.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj.lh)))
tpj.lh.fin <- data.frame(isc=tpj.lh.mod,sub1=tpj.lh$sub1,sub2=tpj.lh$sub2,con=tpj.lh$con,roi='tpj',grp=tpj.lh$group)
tpj.lh.m <- avg.sub.g(tpj.lh.fin,'tpj','lh')

fit <- lmer(isc~con*group+ep+(1|sub1)+(1|sub2),data=tpj.lh)
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
p.lh <- 2*pt(newt,df,lower=FALSE)


# right TPJ
tpj.rh <- subset(fin,fin$hemi=='rh' & fin$group!='neumat')
tpj.rh.mod <- z2r(residuals(lmer(isc~ep+(1|sub1)+(1|sub2),data=tpj.rh)))
tpj.rh.fin <- data.frame(isc=tpj.rh.mod,sub1=tpj.rh$sub1,sub2=tpj.rh$sub2,con=tpj.rh$con,roi='tpj',grp=tpj.rh$group)
tpj.rh.m <- avg.sub.g(tpj.rh.fin,'tpj','rh')

fit <- lmer(isc~con*group+ep+(1|sub1)+(1|sub2),data=tpj.rh)
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
p.rh <- 2*pt(newt,df,lower=FALSE)

ps <- c(p.lh,p.rh)
p.adjust(ps, method='bonferroni')


all <- rbind(tpj.lh.m,tpj.rh.m)
all$con <- factor(all$con,levels=c('scr','int'))

all$con <- revalue(all$con, c(int='Intact',scr='Scrambled'))
all$hemi <- revalue(all$hemi, c(lh='Left', rh='Right'))
all$grp <- revalue(all$grp, c(adult='Adult', child='Child'))

left <- subset(all,all$hemi=='Left')
right <- subset(all,all$hemi=='Right')

ggplot(left, aes(con,isc,fill=con)) + 
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



## Age related differences in age and comprehension 
# for comprehension
comp <- read.table('~/iCloud/TRW/second_submission/comp_dat_8-22-17.txt',header=TRUE)
# for age
samp <- read.csv('~/iCloud/TRW/second_submission/TRW_sample_summary.csv',header=TRUE)
# for motion
mot <- read.csv('~/iCloud/TRW/first_submission/TRW_motion.csv')
mot$int <- (mot$int1.mean+mot$int2.mean)/2
mot$scr <- (mot$scr1.mean+mot$scr2.mean)/2


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
#samp <- read.csv('/Volumes/research$/redcay/DSCN lab/Experiments/TRW/data/TRW_sample_summary.csv',header=TRUE)
samp <- subset(samp,samp$allTRW==1 & samp$include_other==1 & samp$out==1)
samp <- subset(samp,samp$group!='adult')

#comp <- read.table('/Volumes/research$/redcay/DSCN lab/Experiments/TRW/data/comp_dat_8-22-17.txt',header=TRUE)
comp<- subset(comp,comp$include==1)

#mot <- read.csv('~/iCloud/TRW/first_submission/TRW_motion.csv')
#mot$int <- (mot$int1.mean+mot$int2.mean)/2
#mot$scr <- (mot$scr1.mean+mot$scr2.mean)/2
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

all.neu$hemi <- revalue(all.neu$hemi, c(lh='Left', rh='Right'))
all.neu$mOVER <- all.neu$ment/(all.neu$ment+all.neu$nonment)
all.neu$mOVER <- all.neu$mOVER*100
all.neu$con <- factor(all.neu$con,levels=c('scr','int'))

# create avg neumat for plots
avg.sub.neu <- function(data,s,roi,hemi) {
	n <- length(s)
	c <- c('int','scr')
	isc <- vector('numeric',n*length(c))
	sub <- vector('character',n*length(c))
	con <- vector('character',n*length(c))
	r <- vector('character',n*length(c))
	h <- vector('character',n*length(c))
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
			h[count] <- hemi
			age[count] <- new$age[1]
			m[count] <- new$m[1]
			nm[count] <- new$nm[1]
			mOVER[count] <- new$mOVER[1]
		}
	}
	fin <- data.frame(sub=sub,con=con,roi=r,hemi=h,isc=isc,age=age,ment=m,nonment=nm,mOVER=mOVER)
	return(fin)
}


# effect of age
df <- nSubj-4
 

# Left TPJ
tpj.lh <- subset(all.neu,all.neu$hemi=='Left')
fit <- lmer(isc~con*age+ep+fd+(1|sub1)+(1|sub2),data=tpj.lh)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 6
quantile(sdat[,beta], c(0.025,0.975))

tpj.lh.mod <- z2r(residuals(lmer(isc~ep+fd+(1|sub1)+(1|sub2),data=tpj.lh)))
tpj.lh.fin <- data.frame(isc=tpj.lh.mod,sub1=tpj.lh$sub1,sub2=tpj.lh$sub2,con=tpj.lh$con,age=tpj.lh$age,m=tpj.lh$ment,nm=tpj.lh$nonment,mOVER=tpj.lh$mOVER)
tpj.lh.m <- avg.sub.neu(tpj.lh.fin,child,'tpj','Left')

ggplot(tpj.lh.m, aes(age,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(7,8,9,10,11,12,13)) +
	facet_wrap(~hemi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))

# Right TPJ
tpj.rh <- subset(all.neu,all.neu$hemi=='Right')
fit <- lmer(isc~con*age+ep+fd+(1|sub1)+(1|sub2),data=tpj.rh)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 6
quantile(sdat[,beta], c(0.025,0.975))

tpj.rh.mod <- z2r(residuals(lmer(isc~ep+fd+(1|sub1)+(1|sub2),data=tpj.rh)))
tpj.rh.fin <- data.frame(isc=tpj.rh.mod,sub1=tpj.rh$sub1,sub2=tpj.rh$sub2,con=tpj.rh$con,age=tpj.rh$age,m=tpj.rh$ment,nm=tpj.rh$nonment,mOVER=tpj.rh$mOVER)
tpj.rh.m <- avg.sub.neu(tpj.rh.fin,child,'tpj','Right')

ggplot(tpj.rh.m, aes(age,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(7,8,9,10,11,12,13)) +
	facet_wrap(~hemi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))



# effect of mOVER
# Left TPJ
tpj.lh <- subset(all.neu,all.neu$hemi=='Left')
fit <- lmer(isc~con*mOVER+age+ep+fd+(1|sub1)+(1|sub2),data=tpj.lh)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 7
quantile(sdat[,beta], c(0.025,0.975))

tpj.lh.mod <- z2r(residuals(lmer(isc~age+ep+fd+(1|sub1)+(1|sub2),data=tpj.lh)))
tpj.lh.fin <- data.frame(isc=tpj.lh.mod,sub1=tpj.lh$sub1,sub2=tpj.lh$sub2,con=tpj.lh$con,age=tpj.lh$age,m=tpj.lh$ment,nm=tpj.lh$nonment,mOVER=tpj.lh$mOVER)
tpj.lh.m <- avg.sub.neu(tpj.lh.fin,child,'tpj','Left')

ggplot(tpj.lh.m, aes(mOVER,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(10,20,30,40,50,60,70)) +
	facet_wrap(~hemi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))


# Right TPJ
tpj.rh <- subset(all.neu,all.neu$hemi=='Right')
fit <- lmer(isc~con*mOVER+age+ep+fd+(1|sub1)+(1|sub2),data=tpj.rh)
summary(fit)
s <- sim(fit,10000)
sdat <- data.frame(s@fixef)
# find confidence intervals
beta <- 7
quantile(sdat[,beta], c(0.025,0.975))

tpj.rh.mod <- z2r(residuals(lmer(isc~age+ep+fd+(1|sub1)+(1|sub2),data=tpj.rh)))
tpj.rh.fin <- data.frame(isc=tpj.rh.mod,sub1=tpj.rh$sub1,sub2=tpj.rh$sub2,con=tpj.rh$con,age=tpj.rh$age,m=tpj.rh$ment,nm=tpj.rh$nonment,mOVER=tpj.rh$mOVER)
tpj.rh.m <- avg.sub.neu(tpj.rh.fin,child,'tpj','Right')

ggplot(tpj.rh.m, aes(mOVER,isc,color=con)) + 
	geom_hline(yintercept=0,linetype='dashed',alpha=0.3) +
	geom_point(color='black',alpha=0.5) +
	stat_smooth(method='lm') +
	scale_color_manual(values=c('tomato','goldenrod1')) +
	scale_y_continuous(breaks=c(-0.04,-0.02,0,0.02,0.04),limits=c(-0.05,0.05)) +
	scale_x_continuous(breaks=c(10,20,30,40,50,60,70)) +
	facet_wrap(~hemi) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		strip.text=element_text(size=20),
		axis.text.x=element_text(size=18),
		axis.text.y=element_text(size=14))




