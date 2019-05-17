library(ggplot2)
library(wesanderson)

dat <- read.csv('~/iCloud/TRW/second_submission/psd_analysis/all.boot.se.csv')
#dat$grp <- factor(dat$grp, levels=c('child','adult'))
dat <- subset(dat,dat$variable!='vis')
dat$variable <- factor(dat$variable, levels=c('aud','precun','dmPFC','tpj'))

a <- subset(dat,dat$grp=='adult')
c <- subset(dat,dat$grp=='child')


## adult group
ggplot(c, aes(freq,value,fill=variable,group=variable)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04,size=1) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	theme_classic() +
	scale_fill_manual(values=c('gainsboro','mediumorchid3','steelblue3','forestgreen')) +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14))


	geom_line(color='black',alpha=0.9,linetype='longdash') +

dat$variable <- factor(dat$variable, levels=c('tpj','precun','dmPFC','aud'))
# both groups 
ggplot(dat, aes(freq,value,fill=grp)) + 
	geom_ribbon(aes(ymin=value-se,ymax=value+se),color='black',alpha=0.96) +
	geom_vline(xintercept = 0.04) +
	scale_x_continuous(trans='log10',breaks=c(0.01,0.04,0.1,0.33)) +
	facet_wrap(~variable) +
	theme_classic() +
	scale_fill_manual(values=c("darkorchid4","darkslategray3")) +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text=element_text(size=14),
		strip.text=element_blank())

