library(ggplot2)
library(plyr)
library(lme4)
library(lmerTest)
library(arm)

dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/all.boot.se.csv')
#dat$grp <- factor(dat$grp, levels=c('child','adult'))
dat <- subset(dat,dat$variable!='vis')
dat$variable <- factor(dat$variable, levels=c('aud','precun','dmPFC','tpj'))
dat$variable <- revalue(dat$variable, c(aud='Auditory',precun='Precuneus',dmPFC='dmPFC',tpj='TPJ'))

a <- subset(dat,dat$grp=='adult')

c <- subset(dat,dat$grp=='child')


## adult group
ggplot(a, aes(freq,value,fill=variable,group=variable)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04,size=1) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	theme_classic() +
	scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10)) +
	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14))


	geom_line(color='black',alpha=0.9,linetype='longdash') +

aud <- subset(dat,dat$variable=='Auditory')
pre <- subset(dat,dat$variable=='Precuneus')
dmpfc <- subset(dat,dat$variable=='dmPFC')
tpj <- subset(dat,dat$variable=='TPJ')

# both groups 
ggplot(aud, aes(freq,value,fill=grp)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10)) +
	facet_wrap(~variable) +
	theme_classic() +
	scale_fill_manual(values=c("#9a989a","#e6e5e6")) +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))


# aud: "#9a989a","#e6e5e6"
# TPJ: "#1e651b","#3dc936"
# dmPFC: "#2c5987","#8bb3da"
# precun: "#742a89","#c98adb"




dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/all.band.prop.csv')
dat$alpha <- dat$alpha*.01
dat$roi <- factor(dat$roi,levels=c('aud','precun','dmPFC','tpj'))
dat$roi <- revalue(dat$roi, c(aud='Auditory',precun='Precuneus',dmPFC='dmPFC',tpj='TPJ'))
dat$group <- revalue(dat$group,c(adult='Adult',child='Child'))
dat$group <- factor(dat$group,levels=c('Child','Adult'))
# a <- subset(dat,dat$group=='adult')
# c <- subset(dat,dat$group=='child')

# ave <- vector('numeric',8)
# se <- vector('numeric',8)
# group <- vector('character',8)
# roi <- vector('character',8)

# count <- 0
# for (gg in levels(dat$group)) {
# 	for (rr in levels(dat$roi)) {
# 		count <- count+1
# 		temp <- subset(dat,dat$group==gg & dat$roi==rr)
# 		group[count] <- gg
# 		roi[count] <- rr
# 		ave[count] <- mean(temp$alpha)
# 		se[count] <- sd(temp$alpha)/sqrt(length(temp$alpha))
# 	}
# }
# fin <- data.frame(group,roi,ave,se)

# fin$roi <- factor(fin$roi,levels=c('aud','precun','dmPFC','tpj'))

# a <- subset(fin,fin$group=='adult')
# c <- subset(fin,fin$group=='child')

# # ggplot(c, aes(roi,ave,fill=roi)) + 
# # 	geom_errorbar(aes(ymin=ave-se,ymax=ave+se),width=0.1) +
# # 	geom_bar(stat='identity',color='black',width=.95) +
# # 	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
# # 	scale_y_continuous(breaks=c(0.0,0.1,0.2,0.3),limits=c(0,0.35)) +
# # 	theme_classic() +
# # 	labs(x='',y='') +
# # 	theme(legend.position='none',
# # 		axis.text=element_text(size=14),
# # 		axis.text.x=element_blank(),
# # 		axis.ticks.x=element_blank())

# ggplot(c, aes(roi,alpha,fill=roi)) + 
# 	geom_boxplot(outlier.shape=NA,color='black') +
# 	geom_jitter(width=0.15,alpha=0.6,size=0.5) +
# 	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
# 	scale_y_continuous(breaks=c(0.0,0.1,0.2,0.3),limits=c(0,0.35)) +
# 	theme_classic() +
# 	labs(x='',y='') +
# 	theme(legend.position='none',
# 		axis.text=element_text(size=14),
# 		axis.text.x=element_blank(),
# 		axis.ticks.x=element_blank())


# First test adults: is alpha different by ROIS?
a <- subset(dat, dat$group=='Adult')

mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=a)
anova(mod)

# Second test adults: Do the long frequency ROIs show greater alpha compared to Auditory?
tpj.aud <- subset(a, a$roi=='TPJ' | a$roi=='Auditory')
mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=tpj.aud)
summary(mod)
confint(mod)

# dmPFC and aud
dmpfc.aud <- subset(a, a$roi=='dmPFC' | a$roi=='Auditory')
mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=dmpfc.aud)
summary(mod)
confint(mod)

# precuneus and aud
precun.aud <- subset(a, a$roi=='Precuneus' | a$roi=='Auditory')
mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=precun.aud)
summary(mod)
confint(mod)


# First test child: is alpha different by ROIs?
c <- subset(dat, dat$group=='Child') 

mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=c)
anova(mod)

# Second test child: Do the long frequency ROIs show greater alpha compared to Auditory?
tpj.aud <- subset(c, c$roi=='TPJ' | c$roi=='Auditory')
mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=tpj.aud)
summary(mod)
confint(mod)

# dmPFC and aud
dmpfc.aud <- subset(c, c$roi=='dmPFC' | c$roi=='Auditory')
mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=dmpfc.aud)
summary(mod)
confint(mod)

# precuneus and aud
precun.aud <- subset(c, c$roi=='Precuneus' | c$roi=='Auditory')
mod <- lmer(alpha ~ roi + meanFD + (1|subject), data=precun.aud)
summary(mod)
confint(mod)


# Third test: Are there differences in alpha betwen groups?
# group comparisons
tpj <- subset(dat,dat$roi=='TPJ')
t.test(alpha~group,data=tpj)
mod <- lm(alpha ~ group + meanFD, data = tpj)
summary(mod)
confint(mod)

dmpfc <- subset(dat,dat$roi=='dmPFC')
t.test(alpha~group,data=dmpfc)
mod <- lm(alpha ~ group + meanFD, data = dmpfc)
summary(mod)
confint(mod)

precun <- subset(dat,dat$roi=='Precuneus')
t.test(alpha~group,data=precun)
mod <- lm(alpha ~ group + meanFD, data = precun)
summary(mod)
confint(mod)

aud <- subset(dat,dat$roi=='Auditory')
t.test(alpha~group,data=aud)
mod <- lm(alpha ~ group + meanFD, data = aud)
summary(mod)
confint(mod)

p <- c(0.0063,0.0647,0.098,0.223)
p.adjust(p,method='bonferroni')


# Fourth test: Are there age-related differences in alpha?
mod <- lmer(alpha ~ age*roi + meanFD + (1|subject), data=c)
summary(mod)
anova(mod)

tpj <- subset(c,c$roi=='TPJ')
mod <- lm(alpha ~ age + meanFD, data = tpj)
summary(mod)
confint(mod)

dmpfc <- subset(c,c$roi=='dmPFC')
mod <- lm(alpha ~ age + meanFD, data = dmpfc)
summary(mod)
confint(mod)

precun <- subset(c,c$roi=='Precuneus')
mod <- lm(alpha ~ age + meanFD, data = precun)
summary(mod)
confint(mod)

aud <- subset(c,c$roi=='Auditory')
mod <- lm(alpha ~ age + meanFD, data = aud)
summary(mod)
confint(mod)

p <- c(0.99, 0.043, 0.36, 0.72)
p.adjust(p,method='bonferroni')



fit <- lmer(alpha~group*roi+(1|subj),data=dat)
summary(fit)
anova(fit)

ggplot(tpj, aes(group,alpha,fill=group)) + 
	geom_boxplot(outlier.shape=NA,color='black') +
	geom_jitter(width=0.15,alpha=0.6,size=0.5) +
	facet_wrap(~roi) +
	scale_fill_manual(values=c("#1e651b","#3dc936")) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14),
		strip.text=element_blank(),
		panel.spacing = unit(2, "lines"),
		axis.ticks.x=element_blank(),
		axis.text.x=element_text(size=24))



a <- subset(dat,dat$group=='Adult')
pairwise.t.test(a$alpha,a$roi,p.adjust='bonferroni')

aud <- subset(a,a$roi=='Auditory')
pre <- subset(a,a$roi=='Precuneus')
dmpfc <- subset(a,a$roi=='dmPFC')
tpj <- subset(a,a$roi=='TPJ')

test <- rbind(aud,pre)
fit <- lmer(alpha~roi+(1|subj),data=test)
summary(fit)
fit.s <- sim(fit,10000)
sdat <- data.frame(fit.s@fixef)
quantile(sdat[,2], c(0.025,0.975))

test <- rbind(aud,dmpfc)
fit <- lmer(alpha~roi+(1|subj),data=test)
summary(fit)
fit.s <- sim(fit,10000)
sdat <- data.frame(fit.s@fixef)
quantile(sdat[,2], c(0.025,0.975))

test <- rbind(aud,tpj)
fit <- lmer(alpha~roi+(1|subj),data=test)
summary(fit)
fit.s <- sim(fit,10000)
sdat <- data.frame(fit.s@fixef)
quantile(sdat[,2], c(0.025,0.975))


fit <- lmer(alpha~roi+(1|subj),data=a)
summary(fit)
anova(fit)



c <- subset(dat,dat$group=='Child')
pairwise.t.test(c$alpha,c$roi,p.adjust='bonferroni')

aud <- subset(c,c$roi=='Auditory')
pre <- subset(c,c$roi=='Precuneus')
dmpfc <- subset(c,c$roi=='dmPFC')
tpj <- subset(c,c$roi=='TPJ')

test <- rbind(aud,pre)
fit <- lmer(alpha~roi+(1|subj),data=test)
summary(fit)
fit.s <- sim(fit,10000)
sdat <- data.frame(fit.s@fixef)
quantile(sdat[,2], c(0.025,0.975))

test <- rbind(aud,dmpfc)
fit <- lmer(alpha~roi+(1|subj),data=test)
summary(fit)
fit.s <- sim(fit,10000)
sdat <- data.frame(fit.s@fixef)
quantile(sdat[,2], c(0.025,0.975))

test <- rbind(aud,tpj)
fit <- lmer(alpha~roi+(1|subj),data=test)
summary(fit)
fit.s <- sim(fit,10000)
sdat <- data.frame(fit.s@fixef)
quantile(sdat[,2], c(0.025,0.975))


fit <- lmer(alpha~roi+(1|subj),data=c)
summary(fit)
anova(fit)





#### separated TPJ
dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/tpj.boot.se.csv')

dat$variable <- factor(dat$variable, levels=c('lh','rh'))
dat$variable <- revalue(dat$variable, c(lh='Left TPJ',rh='Right TPJ'))

lh <- subset(dat, dat$variable=='Left TPJ')
rh <- subset(dat, dat$variable=='Right TPJ')

## left tpj
ggplot(lh, aes(freq,value,fill=grp)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10)) +
	facet_wrap(~variable) +
	theme_classic() +
	scale_fill_manual(values=c("#1e651b","#3dc936")) +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))

ggplot(rh, aes(freq,value,fill=grp)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10)) +
	facet_wrap(~variable) +
	theme_classic() +
	scale_fill_manual(values=c("#1e651b","#3dc936")) +
	labs(x='',y='') +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))



dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/tpj.band.prop.csv')
dat$alpha <- dat$alpha*.01
dat$roi <- revalue(dat$roi, c(lh='Left TPJ',rh='Right TPJ'))
dat$group <- revalue(dat$group,c(adult='Adult',child='Child'))
dat$group <- factor(dat$group,levels=c('Child','Adult'))

lh <- subset(dat,dat$roi=='Left TPJ')
fit <- lm(alpha~group + meanFD,data=lh)
summary(fit)
confint(fit)

rh <- subset(dat,dat$roi=='Right TPJ')
fit <- lm(alpha~group + meanFD,data=rh)
summary(fit)
confint(fit)
 
p <- c(.002,.023)
p.adjust(p,method='bonferroni')


lh$group <- factor(lh$group,levels=c('Adult','Child'))
rh$group <- factor(rh$group,levels=c('Adult','Child'))

ggplot(rh, aes(group,alpha,fill=group)) + 
	geom_boxplot(outlier.shape=NA,color='black') +
	geom_jitter(width=0.15,alpha=0.6,size=0.5) +
	facet_wrap(~roi) +
	scale_fill_manual(values=c("#1e651b","#3dc936")) +
	theme_classic() +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14),
		strip.text=element_blank(),
		panel.spacing = unit(2, "lines"),
		axis.ticks.x=element_blank(),
		axis.text.x=element_text(size=24))





###### age and low frequency power
dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/age.lowfreq.power.csv')
dat$variable <- factor(dat$variable, levels=c('aud','precun','dmPFC','tpj'))
dat$variable <- revalue(dat$variable, c(aud='Auditory',precun='Precuneus',dmPFC='dmPFC',tpj='TPJ'))


aud <- subset(dat, dat$variable=='Auditory')
precun <- subset(dat, dat$variable=='Precuneus')
dmpfc <- subset(dat, dat$variable=='dmPFC')
tpj <- subset(dat, dat$variable=='TPJ')

ggplot(aud, aes(age,value, fill=variable)) + 
	geom_smooth(method='lm', color='grey30') + 
	geom_point(color='black',pch=21) + 
	scale_fill_manual(values=c('gainsboro')) +
	facet_wrap(~variable) +
	ylim(10,35) +
	labs(x='',y='') +
	theme_classic() +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))

ggplot(precun, aes(age,value, fill=variable)) + 
	geom_smooth(method='lm', color='grey30') + 
	geom_point(color='black',pch=21) + 
	scale_fill_manual(values=c('mediumorchid3')) +
	facet_wrap(~variable) +
	ylim(10,35) +
	labs(x='',y='') +
	theme_classic() +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))

ggplot(dmpfc, aes(age,value, fill=variable)) + 
	geom_smooth(method='lm', color='grey30') + 
	geom_point(color='black',pch=21) + 
	scale_fill_manual(values=c('steelblue3')) +
	facet_wrap(~variable) +
	ylim(10,35) +
	labs(x='',y='') +
	theme_classic() +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))


ggplot(tpj, aes(age,value, fill=variable)) + 
	geom_smooth(method='lm', color='grey30') + 
	geom_point(color='black',pch=21) + 
	scale_fill_manual(values=c('forestgreen')) +
	facet_wrap(~variable) +
	ylim(10,35) +
	labs(x='',y='') +
	theme_classic() +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))

ggplot(dat, aes(age,value, color=variable)) + geom_point() + geom_smooth(method='lm') + facet_wrap(~variable)



# scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
ggplot(a, aes(freq,value,fill=variable,group=variable)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04,size=1) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	theme_classic() +
	scale_y_continuous(breaks=c(0,2,4,6,8,10),limits=c(0,10)) +
	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14))


### with TPJ separate
dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/tpj.age.lowfreq.power.csv')
dat$value <- dat$value*0.01
dat$variable <- revalue(dat$variable, c(lh='Left TPJ',rh='Right TPJ'))

lh <- subset(dat, dat$variable=='Left TPJ')
fit <- lm(value ~ age + meanFD,data=lh)
summary(fit)
confint(fit)



rh <- subset(dat, dat$variable=='Right TPJ')
fit <- lm(value ~ age + meanFD,data=rh)
summary(fit)
confint(fit)



p <- c(.255,.64)
p.adjust(p,method='bonferroni')



ggplot(lh, aes(age,value, fill=variable)) + 
	geom_smooth(method='lm', color='grey30') + 
	geom_point(color='black',pch=21) + 
	scale_fill_manual(values=c('forestgreen')) +
	facet_wrap(~variable) +
	labs(x='',y='') +
	theme_classic() +
	theme(legend.position='none',
		strip.background=element_blank(),
		axis.text=element_text(size=14),
		strip.text=element_text(size=20),
		panel.spacing = unit(2, "lines"))