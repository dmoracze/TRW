# libraries
library(ggplot2)
library(lme4)
library(nlme)
library(arm)
library(MASS)
library(lmerTest)
library(reshape2)
library(plyr)
library(grid)

# read data
dat <- read.table('/Volumes/research$/redcay/DSCN lab/Experiments/TRW/data/comp_dat_8-22-17.txt',header=TRUE)
dat <- subset(dat,dat$include==1)
# dat$m.nm <- dat$percM-dat$percNM
# dat$mOVER <- dat$percM/(dat$percM+dat$percNM)

# read sample info
samp <- read.csv('/Volumes/research$/redcay/DSCN lab/Experiments/TRW/data/TRW_sample_summary.csv',header=TRUE)
samp <- read.table('~/iCloud/TRW/first_submission/TRW_sample_summary.txt',header=TRUE)
samp$order <- as.factor(samp$order)
samp <- subset(samp,samp$allTRW==1 & samp$include_other==1 & samp$out==1)

# merge
new <- merge(dat,samp,by=c('subj','group'))
# grab important columns
good <- new[c('subj','group','con','percM','percNM','age','sex','order','condition','iq','allTRW','episode')]

# split open and melt
g <- melt(good, id.vars=c('subj','group','con','age','sex','order','condition','iq','allTRW','episode'))
#g <- subset(g, g$variable!='m.nm' & g$variable!='mOVER')

# factor
g$group <- factor(g$group,levels=c('child','adult'))
g$con <- factor(g$con,levels=c('int','scr'))
g$variable <- factor(g$variable,levels=c('percM','percNM'))

# subset groups
a <- subset(g,g$group=='adult')
k <- subset(g,g$group=='child')

# define ALL model
All.mod <- lmer(value ~ con*variable*group+episode+(1|subj), data=g)
summary(All.mod)
anova(All.mod)
all.s <- sim(All.mod,10000)
all.sdat <- data.frame(all.s@fixef)

beta <- 4
quantile(all.sdat[,beta], c(0.025,0.975))

# define adult model
amod <- lmer(value ~ con*variable+episode+(1|subj), data=a)
summary(amod)
anova(amod)
as <- sim(amod,10000)
asdat <- data.frame(as@fixef)

beta <- 5
quantile(asdat[,beta], c(0.025,0.975))

# define kid model
kmod <- lmer(value ~ con*variable+age+episode+(1|subj), data=k)
summary(kmod)
anova(kmod)
ks <- sim(kmod,10000)
ksdat <- data.frame(ks@fixef)

# find confidence intervals
beta <- 4
quantile(ksdat[,beta], c(0.025,0.975))


# plot kid distributions
ksdat <- data.frame(ksdat[,2:4],ksdat[,6])
names(ksdat) <- c('Scrambled','Non Mental','Age','Scrambled * Non Mental')
msdat <- melt(ksdat)
msdat$variable <- factor(msdat$variable, levels=c('Age','Scrambled','Non Mental','Scrambled * Non Mental'))
msdat$value <- msdat$value*100
ggplot(msdat, aes(value, fill=variable)) + 
	geom_vline(xintercept = 0,size=.25) +
	geom_density(alpha=0.85) + 
	scale_x_continuous(breaks=c(-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30)) +
	labs(x='',y='') +
	scale_fill_manual(values=c("forestgreen", "firebrick","steelblue", "palevioletred3")) +
	theme_bw() +
	theme(legend.position='none',
		axis.text.x=element_text(size=20),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank())

# plot adult distributions
asdat <- data.frame(asdat[,2:3],asdat[,5])
names(asdat) <- c('Scrambled','Non Mental','Scrambled * Non Mental')
msdat <- melt(asdat)
msdat$variable <- factor(msdat$variable, levels=c('Scrambled','Non Mental','Scrambled * Non Mental'))
msdat$value <- msdat$value*100
ggplot(msdat, aes(value, fill=variable)) + 
	geom_vline(xintercept = 0,size=.25) +
	geom_density(alpha=0.85) + 
	scale_x_continuous(breaks=c(-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30)) +
	labs(x='',y='') +
	scale_fill_manual(values=c("firebrick","steelblue", "palevioletred3")) +
	theme_bw() +
	theme(legend.position='none',
		axis.text.x=element_text(size=20),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank())


# set up average data - kid
# create empty containers
condition <- matrix(NA,nrow=4)
type <- matrix(NA,nrow=4)
mean <- matrix(NA,nrow=4)
se <- matrix(NA,nrow=4)
count <- 1
# find the mean and se for each type/condition
for (c in levels(k$con)) {
	print(c)
	for (t in levels(k$variable)) {
		print(t)
		temp <- subset(k, k$con == c & k$variable == t)
		condition[count] <- c
		type[count] <- t
		mean[count] <- mean(temp$value)*100
		se[count] <- sd(temp$value)/sqrt(length(temp$subj))*100
		count <- count + 1
	}
}
# concatenate average data
mdat.k <- data.frame(condition,type,mean,se)
mdat.k <- data.frame(condition,type,mean,se)
mdat.k$condition <- revalue(mdat.k$condition, c(int='Intact', scr='Scrambled'))
mdat.k$type <- revalue(mdat.k$type, c(percM='Mental', percNM='Non-Mental'))


# kid mean + se
dodge <- position_dodge(width=0.9)
ggplot(mdat.k, aes(type,mean,fill=condition)) + # set initial layout of graph
	geom_bar(stat='identity',position=dodge,color='black') +
	geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=dodge, width = .25) +
	labs(y='Percent',x="") +
	ylim(0,100) +
	scale_fill_manual(values=c("darkorchid4","darkslategray4")) +
	theme_bw() +
	theme(axis.title=element_blank(),
		axis.text.y=element_text(size=36),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.background=element_blank(),
		strip.text=element_text(size=24),
		legend.position='none',
		plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))

# set up average data - adult
# create empty containers
condition <- matrix(NA,nrow=4)
type <- matrix(NA,nrow=4)
mean <- matrix(NA,nrow=4)
se <- matrix(NA,nrow=4)
count <- 1
# find the mean and se for each type/condition
for (c in levels(a$con)) {
	print(c)
	for (t in levels(a$variable)) {
		print(t)
		temp <- subset(a, a$con == c & a$variable == t)
		condition[count] <- c
		type[count] <- t
		mean[count] <- mean(temp$value)*100
		se[count] <- sd(temp$value)/sqrt(length(temp$subj))*100
		count <- count + 1
	}
}
# concatenate average data
mdat.a <- data.frame(condition,type,mean,se)
mdat.a <- data.frame(condition,type,mean,se)
mdat.a$condition <- revalue(mdat.a$condition, c(int='Intact', scr='Scrambled'))
mdat.a$type <- revalue(mdat.a$type, c(percM='Mental', percNM='Non-Mental'))


# adult mean + se
dodge <- position_dodge(width=0.9)
ggplot(mdat.a, aes(type,mean,fill=condition)) + # set initial layout of graph
	geom_bar(stat='identity',position=dodge,color='black') +
	geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=dodge, width = .25) +
	labs(y='Percent',x="") +
	ylim(0,100) +
	scale_fill_manual(values=c("darkorchid4","darkslategray4")) +
	theme_bw() +
	theme(axis.title=element_blank(),
		axis.text.y=element_text(size=36),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.background=element_blank(),
		strip.text=element_text(size=24),
		legend.position='none',
		plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))


k$value <- k$value*100
# age
ggplot(k, aes(age,value, color=con)) + 
geom_point(size=2.5,alpha=0.6) + 
geom_point(size=2.5,alpha=0.95,shape=1) + 
geom_smooth(method='lm',se=FALSE,size=2) + 
facet_wrap(~variable) +
labs(x='Age',y='% correct') +
scale_color_manual(values=c("darkorchid4","darkslategray4")) +
scale_x_continuous(breaks=c(5,6,7,8,9,10,11,12,13,14)) +
theme_bw() +
theme(axis.title=element_blank(),
	legend.position='none',
	axis.text.y=element_text(size=20),
	axis.text.x=element_text(size=20),
	strip.background=element_blank(),
	strip.text = element_blank(),
	plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))


# put kids and adults bar graphs together
mdat.a$group <- 'Adult'
mdat.k$group <- 'Child'
all.mdat <- rbind(mdat.a,mdat.k)
# plot both kid and adult next to each other
dodge <- position_dodge(width=0.9)
ggplot(all.mdat, aes(type,mean,fill=group)) + # set initial layout of graph
	geom_bar(stat='identity',position=dodge,color='black') +
	geom_errorbar(aes(ymin=mean-se,ymax=mean+se), position=dodge, width = .25) +
	facet_wrap(~condition) +
	labs(y='Percent',x="") +
	ylim(0,100) +
	scale_fill_manual(values=c("darkorchid4","darkslategray4")) +
	theme_bw() +
	theme(axis.title=element_blank(),
		axis.text.y=element_text(size=36),
		strip.background=element_blank(),
		strip.text=element_text(size=24),
		plot.margin=unit(c(1.2,1.2,1.2,1.2), "cm"))

ggplot(all.mdat, aes(type,mean,fill=group)) + # set initial layout of graph
	geom_boxplot()




#### post hoc test in the kids
ks.M <- subset(k, k$con=='scr' & k$variable=='percM')
ks.NM <- subset(k, k$con=='scr' & k$variable=='percNM')
ki.M  <- subset(k, k$con=='int' & k$variable=='percM')
ki.NM <- subset(k, k$con=='int' & k$variable=='percNM')

# scrambled mental & scrambled non-mental
ksm.ksnm.t <- t.test(ks.M$value,ks.NM$value, paired=TRUE)
# intact mental & scrambled mental
kim.ksm.t <- t.test(ki.M$value,ks.M$value, paired=TRUE)
# intact non-mental & scrambled non-mental
kinm.ksnm.t <- t.test(ki.NM$value,ks.NM$value, paired=TRUE)
# intact mental & intact non-mental
kim.kinm.t <- t.test(ki.M$value,ki.NM$value, paired=TRUE)

kp <- c(ksm.ksnm.t$p.value, kim.ksm.t$p.value, kinm.ksnm.t$p.value ,kim.kinm.t$p.value)
p.adjust(kp,method='bonferroni')

#### post hoc test in the adults
as.M <- subset(a, a$con=='scr' & a$variable=='percM')
as.NM <- subset(a, a$con=='scr' & a$variable=='percNM')
ai.M  <- subset(a, a$con=='int' & a$variable=='percM')
ai.NM <- subset(a, a$con=='int' & a$variable=='percNM')

# scrambled mental & scrambled non-mental
asm.asnm.t <- t.test(as.M$value,as.NM$value, paired=TRUE)
# intact mental & scrambled mental
aim.asm.t <- t.test(ai.M$value,as.M$value, paired=TRUE)
# intact non-mental & scrambled non-mental
ainm.asnm.t <- t.test(ai.NM$value,as.NM$value, paired=TRUE)
# intact mental & intact non-mental
aim.ainm.t <- t.test(ai.M$value,ai.NM$value, paired=TRUE)

ap <- c(asm.asnm.t$p.value, aim.asm.t$p.value, ainm.asnm.t$p.value, aim.ainm.t$p.value)
p.adjust(ap,method='bonferroni')


# compare child and adult comprehension
# scrambled mental
asm.ksm.t <- t.test(as.M$value,ks.M$value, paired=FALSE)
# scrambled non-mental
asnm.ksnm.t <- t.test(as.NM$value,ks.NM$value, paired=FALSE)
# intact mental
aim.kim.t <- t.test(ai.M$value,ki.M$value, paired=FALSE)
# intact scrambled
ainm.kinm.t <- t.test(ai.NM$value,ki.NM$value, paired=FALSE)

akp <- c(asm.ksm.t$p.value, asnm.ksnm.t$p.value, aim.kim.t$p.value, ainm.kinm.t$p.value)
p.adjust(akp, method='bonferroni')
 



ggplot(g, aes(variable,value, fill=group)) + geom_boxplot() + facet_wrap(~con)





