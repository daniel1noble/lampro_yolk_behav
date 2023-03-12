#########################################################
# lizard performance and behavious as effect of incubation 
# temperature and maternal effects
#########################################################

## Load Packages
  pacman::p_load(rptR, brms, bayestestR, tidyverse)  

## Load Data
  dat<-read.csv( "./data/dataset_reclean.csv")

## Data cleaning and organisation
  str(dat)
  dat$time_tohide<-as.numeric(dat$Time_hide_sec)  
  dat$time_hiding<-as.numeric(dat$Time_snout_sec)
  dat$time_active<-as.numeric(dat$Time_emerge_sec)
  dat$temp<-as.factor(dat$temp)  
  dat$maternal<-as.factor(dat$egg_treat)
  dat$sp<-as.factor(dat$sp)
  str(dat)
#visualizing data
  hist(dat$Distance.moved)#one outlier
  shapiro.test(dat$Distance.moved)
  hist(dat$time_tohide)#one clear outlier
  hist(dat$time_hiding)
  hist(dat$time_active)
  hist(log(dat1$speed_1m_s)) #log makes it much better so use that
  hist(log(dat1$burst_25cm))
  hist(dat1$SVL)
  hist(dat1$Total)
  plot(dat1$SVL,dat1$Total)

  unique(length (dat$id))

  ##checking repeatibility of behavioural responses 

  #first take out animals without intact tail

  dim(dat)
  dat1<-dat[!(dat$Tail_intact=="No"),]
  dim(dat1)#removed 15 lines, so 5 animals

  #dataset for delicata

  dim(dat1)
  dat2<-dat1[!(dat1$sp=="Guich"),]
  dim(dat2)

  #dataset for guich
  dim(dat1)
  dat3<-dat1[!(dat1$sp=="Deli"),]
  dim(dat3)

## Repeatability
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

  #guich
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


#Note: all values quite repeatable

#checking whether there is acclimation on the responses associated with the trial with the pass of time

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

##Analysis of running performance
#####################################

#SVL differs between temps so scale to add as covariate
dat2$svldeli <- scale(dat2$SVL)
dat3$svlguich <- scale(dat3$SVL)
mean2$svldelimean <- scale(mean2$SVL)
mean3$svlguichmean <- scale(mean3$SVL)

hist (log(dat2$speed_1m_s))
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

### Antipredation#########################
#########################################

#deli

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
# Bayesian Multivariate models
####################################

# Repeatability
      # Transformations
        dat2 <- dat2 %>% mutate(logTimeSnout = log(Time_snout_sec),
                                logspeed_1m = log(speed_1m_s),
                                logspeed_burst = na_if(log(burst_25cm), -Inf),
                                logTime_emerge_sec = log(Time_emerge_sec)) 

        dat3 <- dat3 %>% mutate(logTimeSnout = log(Time_snout_sec),
                                logspeed_1m = log(speed_1m_s),
                                logspeed_burst = na_if(log(burst_25cm), -Inf),
                                logTime_emerge_sec = log(Time_emerge_sec)) 
        
        write.csv(dat2, file = "./output/dat2.csv")
        write.csv(dat3, file = "./output/dat3.csv")
      
        dat2  <- read.csv("./output/dat2.csv")
        dat3  <- read.csv("./output/dat3.csv")
    # Delicata
        tim_emerge_ap  <- bf(logTime_emerge_sec | mi()~ 1 + (1|id) + (1|clutch)) + gaussian()
         tim_snout_ap  <- bf(logTimeSnout | mi() ~ 1 + (1|id) + (1|clutch)) + gaussian()
         dist_move_ap  <- bf(Distance.moved | mi() ~ 1 + (1|id) + (1|clutch)) + gaussian()
             speed_per <- bf(logspeed_1m | mi() ~ 1 + (1|id) + (1|clutch)) + gaussian()
       speed_burst_per <- bf(logspeed_burst | mi() ~ 1 + (1|id) + (1|clutch)) + gaussian()

        deli_mv <- brms::brm(tim_emerge_ap + tim_snout_ap + dist_move_ap + speed_per + speed_burst_per + set_rescor(TRUE), iter = 2000, warmup = 1000, chains = 4, cores = 4, file = "./output/models/deli_mv", file_refit = "always", data = dat2)
    
    # Guichenoti

        guich_mv <- brms::brm(tim_emerge_ap + tim_snout_ap + dist_move_ap + speed_per + speed_burst_per + set_rescor(TRUE), iter = 2000, warmup = 1000, chains = 4, cores = 4, save_pars = save_pars(), file = "./output/models/deli_mv", file_refit = "on_change", control = list(adapt_delta = 0.98), data = dat3)
    
####################################
########### Figures ################
####################################
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

