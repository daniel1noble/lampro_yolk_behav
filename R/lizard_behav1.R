#lizard performance and behavious as effect of incubation temperature and maternal effects

dat<-read.csv( "/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/dataset_reclean.csv")
str(dat)
dat$time_tohide<-as.numeric(dat$Time_hide_sec)  
dat$time_hiding<-as.numeric(dat$Time_snout_sec)
dat$time_active<-as.numeric(dat$Time_emerge_sec)
dat$temp<-as.factor(dat$temp)  
dat$maternal<-as.factor(dat$egg_treat)
dat$sp<-as.factor(dat$sp)
str(dat)

hist(dat$Distance.moved)#one outlier
shapiro.test(dat$Distance.moved)
hist(dat$time_tohide)#one clear outlier
hist(dat$time_hiding)
hist(dat$time_active)

unique(length (dat$id))

##checking repeatibility of behaviour responses before removing weird numbers 
#(weird numbers in terms of antipredation are those lonerg than 90 min or 5400sec)

#first take out animals without intact tail

dim(dat)
dat1<-dat[!(dat$Tail_intact=="No"),]
dim(dat1)#removed 15 lines, so 5 animals

#delicata

dim(dat1)
dat2<-dat1[!(dat1$sp=="Guich"),]
dim(dat2)
#dattrial<-dat2[!(dat2$speed_1m_s>30),]
#dim(dattrial)

#guich
dim(dat1)
dat3<-dat1[!(dat1$sp=="Deli"),]
dim(dat3)

library(rptR)
repbehavs <- rpt(time_active ~ (1 | id), grname = "id", data = dat1, 
              datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs)

repbehavs1 <- rpt(time_hiding ~ (1 | id), grname = "id", data = dat1, 
                 datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs1)

repbehavs2 <- rpt(speed_1m_s ~ (1 | id), grname = "id", data = dat1, 
                 datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs2)

repbehavs3 <- rpt(Distance.moved ~ (1 | id), grname = "id", data = dat1, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs3)

#one by one
#deli
repbehavs <- rpt(time_active ~ (1 | id), grname = "id", data = dat2, 
                 datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs)

repsprint <- rpt(burst_25cm ~ (1 | id), grname = "id", data = dat2, 
                 datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs)

repbehavs1 <- rpt(time_hiding ~ (1 | id), grname = "id", data = dat2, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs1)

repbehavs2 <- rpt(speed_1m_s ~ (1 | id), grname = "id", data = dat2, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs2)

repbehavs3 <- rpt(Distance.moved ~ (1 | id), grname = "id", data = dat2, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs3)
#gouch
repbehavs <- rpt(time_active ~ (1 | id), grname = "id", data = dat3, 
                 datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs)

repsprint <- rpt(burst_25cm ~ (1 | id), grname = "id", data = dat3, 
                 datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs)

repbehavs1 <- rpt(time_hiding ~ (1 | id), grname = "id", data = dat3, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs1)

repbehavs2 <- rpt(speed_1m_s ~ (1 | id), grname = "id", data = dat3, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs2)

repbehavs3 <- rpt(Distance.moved ~ (1 | id), grname = "id", data = dat3, 
                  datatype = "Gaussian", nboot = 1000, npermut = 0)
print(repbehavs3)

#Note: all values quite repeatable so can do the analyses with means or with animal ID as random
#I might check both for a couple of variables to see whether the results are similar

#also, check whether there is acclimation on the responses associated with the trial with the pass of time

m1m<-lmer(log(speed_1m_s)~day+(1|id), data=dat1)
hist(residuals(m1m))
summary(m1m)#no acclimation in time to run a m

mmove<-lmer(Distance.moved~day+(1|id), data=dat2)
hist(residuals(mmove))
shapiro.test(residuals(mmove))
summary(mmove)
boxplot(Distance.moved~day, data=dat1)

mmoveg<-lmer(Distance.moved~day+(1|id), data=dat3)
hist(residuals(mmoveg))
summary(mmoveg)
boxplot(Distance.moved~day, data=dat3)#guich habituates

mhide<-lmer(time_hiding~day+(1|id), data=dat2)
hist(residuals(mhide))
summary(mhide)

mhideg<-lmer(time_hiding~day+(1|id), data=dat3)
hist(residuals(mhideg))
summary(mhideg)
boxplot(time_hiding~day, data=dat3)#guich habituates

mtoactive<-lmer(time_active~day+(1|id), data=dat2)
hist(residuals(mtoactive))
summary(mtoactive)

mtoactiveg<-lmer(time_active~day+(1|id), data=dat3)
hist(residuals(mtoactiveg))
summary(mtoactiveg)#guich habituates
boxplot(time_active~day, data=dat3)

##Start analysis of running performance

(hist(dat1$speed_1m_s)) #there are 3 clear outliers. Check model and residuals
#log makes it much better so use that
hist(log(dat1$burst_25cm))

hist(dat1$SVL)
hist(dat1$Total)
plot(dat1$SVL,dat1$Total)

#delicata

dim(dat1)
dat2<-dat1[!(dat1$sp=="Guich"),]
dim(dat2)
#dattrial<-dat2[!(dat2$speed_1m_s>30),]
#dim(dattrial)

#guich
dim(dat1)
dat3<-dat1[!(dat1$sp=="Deli"),]
dim(dat3)

#with outliers
cor.test(dat2$SVL, dat2$Tail)

#SVL differs between temps so scale to add as covariate
dat2$svldeli <- scale(dat2$SVL)
dat3$svlguich <- scale(dat3$SVL)
mean2$svldelimean <- scale(mean2$SVL)
mean3$svlguichmean <- scale(mean3$SVL)


hist(dat2$speed_1m_s)
mveldel<-lmer(log(speed_1m_s)~temp*maternal+svldeli+(1|id), data=dat2)
hist(residuals(mveldel))
summary(mveldel)
plot(mveldel)

mveldel1<-lmer(log(speed_1m_s)~temp+maternal+svldeli+(1|id), data=dat2)
hist(residuals(mveldel1))
summary(mveldel1)
#plot(mveldel)

mveldel2<-lmer(log(burst_25cm)~temp*maternal+svldeli+(1|id), data=dat2)
hist(residuals(mveldel2))
summary(mveldel2)

mveldel3<-lmer(log(burst_25cm)~temp+maternal+svldeli+(1|id), data=dat2)
hist(residuals(mveldel3))
summary(mveldel3)

#without outliers
#hist(dattrial$speed_1m_s)
#mveldel1<-lmer(log(speed_1m_s)~temp*maternal+SVL+(1|id), data=dattrial)
#hist(residuals(mveldel1))
#summary(mveldel1)
#plot(mveldel1)#much better without outliers.  


#trying with average numbers since might be less extreme# USE THIS FOR MAIN TEXT

mean<-read.csv("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/average_all.csv")
str(mean)

mean$mean_dist<-as.numeric(mean$mean_dist)
mean$temp<-as.factor(mean$temp)  
mean$maternal<-as.factor(mean$egg_treat)
mean$sp<-as.factor(mean$sp)
mean$mean_tohide<-as.numeric(mean$mean_tohide)

dim(mean)
mean1<-mean[!(mean$Tail_intact=="No"),]
dim(mean1)#removed 5 animals

table(paste(mean1$maternal, mean1$temp, mean1$sp))

hist(mean1$mean_dist)
hist(mean1$mean_snout)
hist(mean1$mean_emerge)
hist(log(mean1$mean_1m))
hist(log(mean1$mean_25cm))

#deli
mean2<-mean1[!(mean1$sp=="Guich"),]


mmeanveldel<-lm(log(mean_1m)~temp*maternal+svldelimean, data=mean2)
hist(residuals(mmeanveldel))
shapiro.test(residuals(mmeanveldel))#not super good so try boxcox
summary(mmeanveldel)


library(car)
powerTransform(mean2$mean_1m)
mean2$bxcx1mdeli<-bcPower(mean2$mean_1m, -0.9150967)
hist(mean2$bxcx1mdeli)

mean2$maternal <- relevel(mean2$maternal, ref = "control")

mmeanveldel1<-lm(bxcx1mdeli~temp*maternal+svldelimean, data=mean2)
hist(residuals(mmeanveldel1))
shapiro.test(residuals(mmeanveldel1))
#plot(mmeanveldel1)
summary(mmeanveldel1)#no interaction. Simplify model. For ms
Anova(mmeanveldel1, 3)

mmeanveldel2<-lm(bxcx1mdeli~temp+maternal+svldelimean, data=mean2)
hist(residuals(mmeanveldel2))
shapiro.test(residuals(mmeanveldel2))
summary(mmeanveldel2)
boxplot(bxcx1mdeli~maternal, data=mean2)#ablated take longer to run a meter. For ms

mmeanveldel3<-lm(log(mean_25cm)~temp*maternal+svldelimean, data=mean2)
hist(residuals(mmeanveldel3))
shapiro.test(residuals(mmeanveldel3))
summary(mmeanveldel3)#not normal, use bxcxc

powerTransform(mean2$mean_25cm)
mean2$bxcx25deli<-bcPower(mean2$mean_25cm, -0.8680516)

mmeanveldel4<-lm(bxcx25deli~temp*maternal+svldelimean, data=mean2)
hist(residuals(mmeanveldel4))
shapiro.test(residuals(mmeanveldel4))
summary(mmeanveldel4)#this for ms

mmeanveldel5<-lm(bxcx25deli~temp+maternal+svldelimean, data=mean2)
hist(residuals(mmeanveldel5))
shapiro.test(residuals(mmeanveldel5))
summary(mmeanveldel5)#this for ms


#guich

mean3<-mean1[!(mean1$sp=="Deli"),]


mmean1mg<-lm(log(mean_1m)~temp*maternal+svlguichmean, data=mean3)
hist(residuals(mmean1mg))
shapiro.test(residuals(mmean1mg))#not super good so try boxcox
summary(mmean1mg)

powerTransform(mean3$mean_1m)
mean3$bxcx1guich<-bcPower(mean3$mean_1m, -0.9131911)

mmean1mg1<-lm(bxcx1guich~temp*maternal+svlguichmean, data=mean3)
hist(residuals(mmean1mg1))
shapiro.test(residuals(mmean1mg1))
summary(mmean1mg1)#no 2-way

mmean1mg1red<-lm(bxcx1guich~temp+maternal+svlguichmean, data=mean3)
hist(residuals(mmean1mg1red))
shapiro.test(residuals(mmean1mg1red))
summary(mmean1mg1red)

mmean25mg<-lm(log(mean_25cm)~temp*maternal+svlguichmean, data=mean3)
hist(residuals(mmean25mg))
shapiro.test(residuals(mmean25mg))#not super good so try boxcox
summary(mmean25mg)

powerTransform(mean3$mean_25cm)
mean3$bxcx25guich<-bcPower(mean3$mean_25cm, -0.8536025)

mmean25mg1<-lm(bxcx25guich~temp*maternal+svlguichmean, data=mean3)
hist(residuals(mmean25mg1))
shapiro.test(residuals(mmean25mg1))
summary(mmean25mg1)#for ms

mmean25mg2<-lm(bxcx25guich~temp+maternal+svlguichmean, data=mean3)
hist(residuals(mmean25mg2))
shapiro.test(residuals(mmean25mg2))
summary(mmean25mg2)#for ms

#Quick checK:does treatment affect morphological traits? in delicata

#scaling age to add as covariate
mean2$agesacledeli <- scale(mean2$age)
mean3$agesaclguic <- scale(mean3$age)

mtail<-lm(Tail~temp*maternal+agesacledeli, data=mean2)
hist(residuals(mtail))
summary(mtail)#intercation significant
lsmeans(mtail,pairwise~temp*maternal, adjust = "tukey")
boxplot(Tail~temp*maternal, data=mean2)
mtail1<-lm(Tail~temp+maternal+agesacledeli, data=mean2)
hist(residuals(mtail1))
summary(mtail1)

mtot<-lm(Total~temp*maternal+agesacledeli, data=mean2)
hist(residuals(mtot))
summary(mtot)
mtot1<-lm(Total~temp+maternal+agesacledeli, data=mean2)
hist(residuals(mtot1))
summary(mtot1)
boxplot(Total~temp, data=mean2)

msvl<-lm(SVL~temp*maternal+agesacledeli, data=mean2)
hist(residuals(msvl))
summary(msvl)
msvl1<-lm(SVL~temp+maternal+agesacledeli, data=mean2)
hist(residuals(msvl1))
summary(msvl1)
boxplot(SVL~temp*maternal, data=mean2)

mmass<-lm(Weigth~temp*maternal+agesacledeli, data=mean2)
hist(residuals(mmass))
summary(mmass)

#in guich
mtail<-lm(Tail~temp*maternal+agesaclguic, data=mean3)
hist(residuals(mtail))
summary(mtail)
mtail1<-lm(Tail~temp+maternal+agesaclguic, data=mean3)
hist(residuals(mtail1))
summary(mtail1)

mtot<-lm(Total~temp*maternal+agesaclguic, data=mean3)
hist(residuals(mtot))
summary(mtot)
mtot1<-lm(Total~temp+maternal+agesaclguic, data=mean3)
hist(residuals(mtot1))
summary(mtot1)
boxplot(Total~temp*maternal, data=mean3)

msvl<-lm(SVL~temp*maternal+agesaclguic, data=mean3)
hist(residuals(msvl))
summary(msvl)
msv1<-lm(SVL~temp+maternal+agesaclguic, data=mean3)
hist(residuals(msv1))
summary(msv1)
boxplot(SVL~temp*maternal, data=mean3)

mmass<-lm(Weigth~temp*maternal+agesaclguic, data=mean3)
hist(residuals(mmass))
summary(mmass)
mmass1<-lm(Weigth~temp+maternal+agesaclguic, data=mean3)
hist(residuals(mmass1))
summary(mmass1)


### antipredation#########################
#########################################

hist(mean2$mean_dist)
hist(dat2$Movement)
hist(mean2$mean_snout)
hist(dat2$Time_snout_sec)
hist(mean2$mean_emerge)
hist(dat2$Time_emerge_sec)
hist(log(mean2$mean_tohide))
hist(log(dat2$Time_hide_sec))

hist(mean3$mean_dist)
hist(dat3$Movement)
hist(mean3$mean_snout)
hist(dat3$Time_snout_sec)
hist(mean3$mean_emerge)
hist(dat3$Time_emerge_sec)
hist(mean3$mean_tohide)
hist(dat3$Time_hide_sec)



#mean
#deli
mactivity<-lm(mean_dist~temp*maternal+svldelimean, data=mean2)
hist(residuals(mactivity))
shapiro.test(residuals(mactivity))
#plot(mactivity)
summary (mactivity)#no interaction. Simplify model

mactivity1<-lm(mean_dist~temp+maternal+svldelimean, data=mean2)
hist(residuals(mactivity1))
shapiro.test(residuals(mactivity1))
summary (mactivity1)

msnout<-lm(mean_snout~temp*maternal+svldelimean+mean_dist, data=mean2)
hist(residuals(msnout))
shapiro.test(residuals(msnout))
#plot(msnout)
summary (msnout)#no interaction. Simplify model

msnout1<-lm(mean_snout~temp+maternal+svldelimean+mean_dist, data=mean2)
hist(residuals(msnout1))
shapiro.test(residuals(msnout1))
summary (msnout1)

memerge<-lm(mean_emerge~temp*maternal+svldelimean+mean_dist, data=mean2)
hist(residuals(memerge))
shapiro.test(residuals(memerge))
#plot(memerge)
summary (memerge)#no interaction. Simplify model

memerge1<-lm(mean_emerge~temp+maternal+svldelimean+mean_dist, data=mean2)
hist(residuals(memerge1))
shapiro.test(residuals(memerge1))
summary (memerge1)
plot(mean2$mean_emerge, mean2$svldelimean)


#ind as random

mactivity<-lmer(Distance.moved~temp*maternal+svldeli+(1|id), data=dat2)
hist(residuals(mactivity))
shapiro.test(residuals(mactivity))
plot(mactivity)
summary(mactivity)


mactivity.2<-lmer(Distance.moved~temp+maternal+svldeli+(1|id), data=dat2)
hist(residuals(mactivity.2))
shapiro.test(residuals(mactivity.2))
plot(mactivity.2)
summary(mactivity.2)

sim2 <- simulateResiduals(glmmTMB(Movement ~ temp*maternal+svldeli+(1|id), data=dat2, family="gaussian"))
testResiduals(sim2)#not bad

mactivity.3<-lmer(Distance.moved~temp+maternal+svldeli+(1|id), data=dat2)
hist(residuals(mactivity.3))
shapiro.test(residuals(mactivity.3))
plot(mactivity.3)
summary(mactivity.3)     

msnout<-lmer(Time_snout_sec~temp*maternal+svldeli+(1|id), data=dat2)
hist(residuals(msnout))
shapiro.test(residuals(msnout))
#plot(msnout)
summary (msnout)#no interaction. Simplify model

msnout1<-lmer(Time_snout_sec~temp+maternal+svldeli+(1|id), data=dat2)
hist(residuals(msnout1))
shapiro.test(residuals(msnout1))
summary (msnout1)

memerge<-lmer(Time_emerge_sec~temp*maternal+svldeli+(1|id), data=dat2)
hist(residuals(memerge))
shapiro.test(residuals(memerge))
#plot(memerge)
summary (memerge)#no interaction. Simplify model

memerge1<-lmer(Time_emerge_sec~temp+maternal+svldeli+(1|id), data=dat2)
hist(residuals(memerge1))
shapiro.test(residuals(memerge1))
summary (memerge1)


#####guich

mactivityg<-lm(log(mean_dist)~temp*maternal+svlguichmean, data=mean3)
hist(residuals(mactivityg))
shapiro.test(residuals(mactivityg))
#plot(mactivity)
summary (mactivityg)#no interaction. Simplify model

mactivity1g<-lm(log(mean_dist)~temp+maternal+svlguichmean, data=mean3)
hist(residuals(mactivity1g))
shapiro.test(residuals(mactivity1g))
summary (mactivity1g)

msnoutg<-lm(log(mean_snout)~temp*maternal+svlguichmean, data=mean3)
hist(residuals(msnoutg))
shapiro.test(residuals(msnoutg))
#plot(msnout)
summary (msnoutg)#no interaction. Simplify model

msnout1g<-lm(log(mean_snout)~temp+maternal+svlguichmean+mean_dist, data=mean3)
hist(residuals(msnout1g))
shapiro.test(residuals(msnout1g))
summary (msnout1g)

memergeg<-lm(log(mean_emerge)~temp*maternal+svlguichmean+mean_dist, data=mean3)
hist(residuals(memergeg))
shapiro.test(residuals(memergeg))
#plot(memerge)
summary (memergeg)#no interaction. Simplify model

memerge1g<-lm(log(mean_emerge)~temp+maternal+svlguichmean+mean_dist, data=mean3)
hist(residuals(memerge1g))
shapiro.test(residuals(memerge1g))
summary (memerge1g)


#ind as random

mactivityg<-lmer(Distance.moved~temp*maternal+svlguich+(1|id), data=dat3)
hist(residuals(mactivityg))
shapiro.test(residuals(mactivityg))
plot(mactivityg)
summary(mactivityg)


mactivity1g<-lmer(Distance.moved~temp+maternal+svlguich+(1|id), data=dat3)
hist(residuals(mactivity1g))
shapiro.test(residuals(mactivity1g))
plot(mactivity1g)
summary(mactivity1g)


msnoutg<-lmer(log(Time_snout_sec)~temp*maternal+svlguich+(1|id), data=dat3)
hist(residuals(msnoutg))
shapiro.test(residuals(msnoutg))
#plot(msnout)
summary (msnoutg)#no interaction. Simplify model

msnout1g<-lmer(log(Time_snout_sec)~temp+maternal+svlguich+(1|id), data=dat3)
hist(residuals(msnout1g))
shapiro.test(residuals(msnout1g))
summary (msnout1g)

memergeg<-lmer(Time_emerge_sec~temp*maternal+svlguich+(1|id), data=dat3)
hist(residuals(memergeg))
shapiro.test(residuals(memergeg))
#plot(memerge)
summary (memergeg)#no interaction. Simplify model

memerge1g<-lmer(Time_emerge_sec~temp+maternal+svlguich+(1|id), data=dat3)
hist(residuals(memerge1g))
shapiro.test(residuals(memerge1g))
summary (memerge1g)



####################################
###########figures################

#morphol deli
pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/taildeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=Tail, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Tail length (mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/svldeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=SVL, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Snout vent length (SVL, mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/massdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=Weigth, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Mass (gr)') + xlab('Temperature')+
  theme_bw()
dev.off()


#morphol gucih
pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/tailguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=Tail, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962'))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Tail length, mm') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/svlguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=SVL, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962'))+
  labs(fill = "ps") +
  ylab('Snout vent length (SVL, mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/massguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=Weigth, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962'))+
  labs(fill = "ps") +
  ylab('Mass (gr)') + xlab('Temperature')+
  theme_bw()
dev.off()


##performance
#deli

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/distancedeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=mean_dist, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Distance travelled (mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/1mdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=bxcx1mdeli, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Running velocity 1m (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/25cmdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=bxcx25deli, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Running velocity 25cm (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()

##guich

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/distanceguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=log(mean_dist), fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
  labs(fill = "ps") +
  ylab('Distance travelled (log, mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/1mguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=bxcx1guich, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
  labs(fill = "ps") +
  ylab('Running velocity 1m (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/25cmguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=bxcx25guich, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
  labs(fill = "ps") +
  ylab('Running velocity 25cm (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()


#behaviour

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/hiddingdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=mean_snout, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #labs(fill = "ps") +
  ylab('Time hiding (s)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/toactivedeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=mean_emerge, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
 # labs(fill = "ps") +
  ylab('Time to active (s)') + xlab('Temperature')+
  theme_bw()
dev.off()

##guich

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/hiddingguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=mean_snout, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
 # labs(fill = "ps") +
  ylab('Time hiding (log, s)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("/Users/maitxu/Documents/2023 EBD/papers/lizard_behav/figs/toactiveguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=mean_emerge, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
  labs(fill = "ps") +
  ylab('Time to active (log, s)') + xlab('Temperature')+
  theme_bw()
dev.off()
