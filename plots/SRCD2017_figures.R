library(ggplot2)
library(lme4)
library(arm)
library(MASS)
library(lmerTest)
library(reshape2)
library(plyr)
library(grid)

dat <- read.table('~/Dropbox/DSCN/Experiments/TRW/data/comp_dat_strict_motion_Mar8.txt',header=T)


###############################
##### participant ageXsex #####
###############################

# ggplot(dat, aes(age, fill=sex)) + 
# 	geom_histogram(breaks=seq(6,14, by=0.99999),col="black",alpha=0.7) +
# 	scale_x_continuous(breaks=6:14) +
# 	labs(x="Age",y="Number of participants") +
# 	theme_bw() +
# 	scale_fill_manual(values=c("forestgreen","steelblue")) +
# 	theme(panel.grid.minor.x=element_blank(),
# 		panel.grid.minor.y=element_blank(),
# 		plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))

#####################################
##### mean+se for 2x2 factorial #####
#####################################

# set up average data
# create empty containers
condition <- matrix(NA,nrow=4)
type <- matrix(NA,nrow=4)
mean <- matrix(NA,nrow=4)
se <- matrix(NA,nrow=4)
count <- 1
# find the mean and se for each type/condition
for (c in levels(dat$con)) {
	print(c)
	for (t in levels(dat$type)) {
		print(t)
		temp <- subset(dat, dat$con == c & dat$type == t)
		condition[count] <- c
		type[count] <- t
		mean[count] <- mean(temp$perc)
		se[count] <- sd(temp$perc)/sqrt(length(temp$subj))
		count <- count + 1
	}
}
# concatenate average data
mdat <- data.frame(condition,type,mean,se)
mdat$mean <- mdat$mean*100
mdat$se <- mdat$se*100
mdat$condition <- revalue(mdat$condition, c(int='Intact', scr='Scrambled'))
mdat$type <- revalue(mdat$type, c(M='Mental', NM='Non-Mental'))
# plot the average data
dodge <- position_dodge(width=0.9)
ggplot(mdat, aes(condition,mean,fill=type)) +
	geom_bar(stat='identity',position=dodge) +
	geom_bar(stat='identity',position=dodge,color='black') +
	geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=dodge, width = .25) +
	labs(y='Percent',x="") +
	ylim(0,100) +
	scale_fill_manual(values=c("darkorchid4","darkslategray4")) +
	theme_bw() +
	theme(axis.title=element_blank(),
		legend.position='none',
		axis.text.y=element_text(size=20),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.background=element_blank(),
		strip.text=element_text(size=24),
		plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))

#########################
##### ageXcomp plot #####
#########################

dat$con <- revalue(dat$con, c(int='Intact', scr='Scrambled'))
dat$type <- revalue(dat$type, c(M='Mental', NM='Non-Mental'))

ggplot(dat, aes(age,perc, color=type)) + 
	geom_point(size=2.5,alpha=0.6) + 
	geom_point(size=2.5,alpha=0.95,shape=1) + 
	geom_smooth(method='lm',se=FALSE,size=2) + 
	facet_wrap(~con) +
	labs(x='Age',y='% correct') +
	scale_color_manual(values=c("darkorchid4","darkslategray4")) +
	scale_x_continuous(breaks=c(5,6,7,8,9,10,11,12,13,14)) +
	theme_bw() +
	theme(axis.title=element_blank(),
		legend.position='none',
		axis.text.y=element_text(size=20),
		axis.text.x=element_text(size=20),
		strip.background=element_blank(),
		strip.text=element_blank(),
		plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))


####################################
##### lmer model distributions #####
####################################

# define model

mod <- lmer(perc ~ con*type + age + (1|subj), data=dat)
display(mod)
summary(mod)


s <- sim(mod,100000)
sdat <- data.frame(s@fixef)
head(sdat)
names(sdat) <- c('Intercept','Scrambled','Non Mental','Age','Scrambled * Non Mental')
msdat <- melt(sdat[,2:ncol(sdat)])
msdat$variable <- factor(msdat$variable, levels=c('Age','Scrambled','Non Mental','Scrambled * Non Mental'))

ggplot(msdat, aes(value, fill=variable)) + 
	geom_vline(xintercept = 0,size=.25) +
	geom_density(alpha=0.85) + 
	scale_x_continuous(breaks=c(-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3)) +
	labs(x='',y='') +
	scale_fill_manual(values=c('peachpuff3','steelblue','forestgreen','firebrick')) +
	theme_bw() +
	theme(legend.position='none',
		axis.text.x=element_text(size=16),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank())
	

beta <- 2
quantile(sdat[,beta], c(0.025,0.975))



