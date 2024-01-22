#########################################################
# lizard performance and behavious as effect of incubation 
# temperature and maternal effects
#########################################################

#########################################################
## Load Packages and Visualising Data
#########################################################
    pacman::p_load(rptR, brms, bayestestR, tidyverse, patchwork)  

  ## Load Data
    dat<-read.csv( "./data/dataset_reclean_git.csv")
  
  ## Functions
  source("./R/func.R")

  ## Data cleaning and organisation
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
    hist(log(dat$X25fast))
    plot(dat1$SVL,dat1$Total)

    unique(length (dat$id)) 

    #first take out animals without intact tail

    dim(dat)
    dat1<-dat[!(dat$Tail_intact=="No"),]
    dim(dat1)#removed 15 lines, so 5 animals

    #dataset for delicata

    dim(dat1)
    dat2<-dat1[!(dat1$sp=="Guich"),]
    dim(dat2)
  hist(dat2$Time_snout_sec)
  shapiro.test(dat2$Time_snout_sec)
  hist(dat2$Time_emerge_sec)
  hist(dat2$speed_1m_s)
  dat22<-dat2[!(dat2$X25fast>8),]
  hist(dat2$Distance.moved)

  hist(dat22$X25fast)
    
    #dataset for guich
    dim(dat1)
    dat3<-dat1[!(dat1$sp=="Deli"),]
    dim(dat3)

    hist(log(dat3$Time_snout_sec))
    shapiro.test(dat3$Time_snout_sec)
    hist(log(dat3$Time_emerge_sec))
    hist(log(dat3$speed_1m_s))
    hist(dat3$Distance.moved)
    hist(log(dat3$X25fast))
    
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

############################################
# Morphology Data
############################################

# Use different dataset since we don't have repeated measures for this

          morph <- read.csv( "./data/morphol.csv")

     morph$temp <- as.factor(morph$temp)  
morph$egg_treat <- as.factor(morph$egg_treat)
       morph$sp <- as.factor(morph$sp)
 morph$scaleage <- scale(morph$age)

# Extract each species separately
  morph2<-morph[!(morph$sp=="Guich"),]
  morph3<-morph[!(morph$sp=="Deli"),]

############################################
# Morphology Models - delicata
############################################
  refit = FALSE
  ### INTERACTION MODEL: do temp and maternal interact? We can get everything we need (main effects etc) form this model.
  if(refit){
      svl   <- bf(SVL    ~ 1 + temp*egg_treat + (1|clutch)) + gaussian()
      mass  <- bf(Weigth ~ 1 + temp*egg_treat + (1|clutch)) + gaussian()
      tail  <- bf(Tail   ~ 1 + temp*egg_treat + (1|clutch)) + gaussian()

    deli_morph_int <- brms::brm(svl + mass + tail  + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "output/models/deli_morph_int", file_refit = "on_change", data = morph2, control = list(adapt_delta = 0.98))

  } else {
    
    deli_morph_int <- readRDS("./output/models/deli_morph_int.rds")
  }

  # Checking model fit
    brms_model_check(deli_morph_int)
    
  # With age

  if(refit){
      svl   <- bf(SVL    ~ 1 + temp*egg_treat  + scaleage + (1|clutch)) + gaussian()
      mass  <- bf(Weigth ~ 1 + temp*egg_treat  + scaleage + (1|clutch)) + gaussian()
      tail  <- bf(Tail   ~ 1 + temp*egg_treat  + scaleage + (1|clutch)) + gaussian()

    deli_morph_int_age <- brms::brm(svl + mass + tail  + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "output/models/deli_morph_int_age", file_refit = "on_change", data = morph2, control = list(adapt_delta = 0.98))

  } else {
    deli_morph_int_age <- readRDS("./output/models/deli_morph_int_age.rds")
  }

  # Checking model fit
    brms_model_check(deli_morph_int_age)

  # Extract the posteriors for each trait from the model and calculate the mean for each of the treatment groups and get the contrasts that are relevant. Note that if you want to use the age corrected models then you should be aware that the means are for an averaged aged animal.
         deli_svl <- extract_post(deli_morph_int, "SVL")
         contrast_post(deli_svl)
         contrast_post_main(deli_svl)
      deli_weight <- extract_post(deli_morph_int, "Weigth")
         contrast_post(deli_weight)
         contrast_post_main(deli_weight)
        deli_tail <- extract_post(deli_morph_int, "Tail")
         contrast_post(deli_tail)
         contrast_post_main(deli_tail)
  

############################################
# Morphology Models - guichenoti
############################################

### INTERACTION MODEL: do temp and maternal interact? We can also extract the main effects from this model. See below.
      
  if(refit){
    svl   <- bf(SVL     ~ 1 + temp*egg_treat + (1|clutch)) + gaussian()
    mass  <- bf(Weigth  ~ 1 + temp*egg_treat + (1|clutch)) + gaussian()
    tail  <- bf(Tail    ~ 1 + temp*egg_treat + (1|clutch)) + gaussian()

    guich_morph_int <- brms::brm(svl + mass + tail  + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "output/models/guich_morph_int", file_refit = "on_change", data = morph3, control = list(adapt_delta = 0.98))
    
  } else {
         guich_morph_int <- readRDS("./output/models/guich_morph_int.rds")
  }

    # Checking model fit
    brms_model_check(guich_morph_int)

  # Age model
  if(refit){
    svl   <- bf(SVL     ~ 1 + temp*egg_treat  + scaleage + (1|clutch)) + gaussian()
    mass  <- bf(Weigth  ~ 1 + temp*egg_treat  + scaleage + (1|clutch)) + gaussian()
    tail  <- bf(Tail    ~ 1 + temp*egg_treat  + scaleage + (1|clutch)) + gaussian()

    guich_morph_int_age <- brms::brm(svl + mass + tail  + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "output/models/guich_morph_int_age", file_refit = "on_change", data = morph3, control = list(adapt_delta = 0.98))
    
  } else {
         guich_morph_int_age <- readRDS("./output/models/guich_morph_int_age.rds")
  }

  # Checking model fit
    brms_model_check(guich_morph_int_age)

 # Extract the posteriors for each trait from the model and calculate the mean for each of the treatment groups and get the contrasts that are relevant. Note that if you want to use the age corrected models then you should be aware that the means are for an averaged aged animal.
         guich_svl <- extract_post(guich_morph_int, "SVL")
         contrast_post(guich_svl)
         contrast_post_main(guich_svl)
      guich_weight <- extract_post(guich_morph_int, "Weigth")
         contrast_post(guich_weight)
          contrast_post_main(guich_weight)
        guich_tail <- extract_post(guich_morph_int, "Tail")
         contrast_post(guich_tail)
         contrast_post_main(guich_tail)
         
  

#############################################
# Bayesian Multivariate models - Part I
############################################

# Transformations
        rerun = FALSE
        if(rerun){dat2 <- dat2 %>% mutate(logTimeSnout = log(Time_snout_sec),
                                logspeed_1m = log(speed_1m_s),
                                logspeed_burst = na_if(log(X25fast), -Inf),
                                logTime_emerge_sec = log(Time_emerge_sec),
                                z_svl = scale(SVL)) 

        dat3 <- dat3 %>% mutate(logTimeSnout = log(Time_snout_sec),
                                logspeed_1m = log(speed_1m_s),
                                logspeed_burst = na_if(log(X25fast), -Inf),
                                logTime_emerge_sec = log(Time_emerge_sec),
                                z_svl = scale(SVL)) 
        
        write.csv(dat2, file = "output/data/dat2.csv")
        write.csv(dat3, file = "output/data/dat3.csv")
        } else{
          dat2  <- read.csv("output/data/dat2.csv")
          dat3  <- read.csv("output/data/dat3.csv")}
colnames(dat2)[!colnames(dat2) %in% colnames(dat3)]
# Summarise missing data
  dat2_deli <- dat2  %>% select(logTime_emerge_sec, logTimeSnout, Distance.moved, logspeed_1m, logspeed_burst) %>% summarise_all(funs((sum(is.na(.))/length(.))*100))
 dat2_guich <- dat3  %>%  select(logTime_emerge_sec, logTimeSnout, Distance.moved, logspeed_1m, logspeed_burst)  %>% summarise_all(funs((sum(is.na(.))/length(.))*100))

# Missing data Table
  miss_data <- rbind(dat2_deli, dat2_guich)
  miss_data <- cbind(species=c("Deli", "Guich"), miss_data)  %>% mutate(across(where(is.numeric), round, 2))
  write.csv(miss_data, file = "output/tables/miss_data.csv")

#Transformations Maider (deli does no need to log-transform antipredatory behaviours)
    # The model. Intercept only controlling for ID and clutch. Most variables are approximately normal. Missing data will be dealt with during model fitting using data augmentation.

        tim_emerge_ap  <- bf(logTime_emerge_sec | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
         tim_snout_ap  <- bf(logTimeSnout       | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
         dist_move_ap  <- bf(Distance.moved  | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
             speed_per <- bf(logspeed_1m     | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
       speed_burst_per <- bf(logspeed_burst  | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()

    # Delicata
      if(rerun){
        deli_mv <- brms::brm(tim_emerge_ap + tim_snout_ap + dist_move_ap + speed_per + speed_burst_per + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "./output/models/deli_mv", file_refit = "on_change", data = dat2, control = list(adapt_delta = 0.98))
        saveRDS(deli_mv, file = "./output/models/deli_mv.rds")
        } else {
           deli_mv <- readRDS("./output/models/deli_mv.rds")
        }
      
        # Guichenoti
   
   tim_emerge_ap  <- bf(logTime_emerge_sec | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
    tim_snout_ap  <- bf(logTimeSnout       | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
    dist_move_ap  <- bf(Distance.moved     | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
        speed_per <- bf(logspeed_1m        | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
  speed_burst_per <- bf(logspeed_burst     | mi() ~ 1 + (1|q|id) + (1|clutch)) + gaussian()
   
        if(rerun){
        guich_mv <- brms::brm(tim_emerge_ap + tim_snout_ap + dist_move_ap + speed_per + speed_burst_per + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "output/models/guich_mv", file_refit = "on_change", data = dat3, control = list(adapt_delta = 0.98))
        saveRDS(guich_mv, file = "./output/models/guich_mv.rds")
        } else {
           guich_mv <- readRDS("./output/models/guich_mv.rds")
        }
        

############################################
# Repeatability
############################################

# Delicata
  # Extract posteriors
    posterior_deli <- posterior_samples(deli_mv, pars = c("^sd", "sigma"))
       clutch_deli <- posterior_deli[,grep("clutch", colnames(posterior_deli))]
           id_deli <- posterior_deli[,grep("id", colnames(posterior_deli))]
        sigma_deli <- posterior_deli[,grep("sigma", colnames(posterior_deli))]

    # Calculate the repeatability of each trait. Note the trait name needs to match how it's stored in brms
    R_deli_speed_burst <- repeatability(id_sd = id_deli, clutch_sd = clutch_deli, sigma_sd = sigma_deli, trait = "logspeedburst")
     R_deli_timeemerge <- repeatability(id_deli, clutch_deli, sigma_deli, trait = "logTimeemergesec")
        R_deli_speed1m <- repeatability(id_deli, clutch_deli, sigma_deli, trait = "logspeed1m")
      R_deli_timesnout <- repeatability(id_deli, clutch_deli, sigma_deli, trait = "logTimeSnout")
       R_deli_distmove <- repeatability(id_deli, clutch_deli, sigma_deli, trait = "Distancemoved")

    # Create an organised table to summarise all the repeatability. 
    table1 <- rbind(R_deli_speed_burst, R_deli_timeemerge, R_deli_speed1m, R_deli_timesnout, R_deli_timeemerge)
    traits <- c("Burst Speed (m/s - log)", "Time to Emerge (s - log)", "Sprint Speed - 1m (m/s - log)", "Time Snout Out (s - log)", "Distanced Moved (cm)")
    table1  <-  cbind(traits, table1)

# Guichenoti
    # Extract posteriors
    posterior_guich <- posterior_samples(guich_mv, pars = c("^sd", "sigma"))
       clutch_guich <- posterior_guich[,grep("clutch", colnames(posterior_guich))]
           id_guich <- posterior_guich[,grep("id", colnames(posterior_guich))]
        sigma_guich <- posterior_guich[,grep("sigma", colnames(posterior_guich))]

    # Calculate the repeatability of each trait. Note the trait name needs to match how it's stored in brms
    R_guich_speed_burst <- repeatability(id_guich, clutch_guich, sigma_guich, trait = "logspeedburst")
     R_guich_timeemerge <- repeatability(id_guich, clutch_guich, sigma_guich, trait = "logTimeemergesec")
        R_guich_speed1m <- repeatability(id_guich, clutch_guich, sigma_guich, trait = "logspeed1m")
      R_guich_timesnout <- repeatability(id_guich, clutch_guich, sigma_guich, trait = "logTimeSnout")
       R_guich_distmove <- repeatability(id_guich, clutch_guich, sigma_guich, trait = "Distancemoved")

    # Create an organised table to summarise all the repeatability. 
    table2 <- rbind(R_guich_speed_burst, R_guich_timeemerge, R_guich_speed1m, R_guich_timesnout, R_guich_timeemerge)
    traits <- c("Burst Speed (m/s - log)", "Time to Emerge (s - log)", "Sprint Speed - 1m (m/s - log)", "Time Snout Out (s - log)", "Distanced Moved (cm)")
    table2  <-  cbind(traits, table2)

############################################
# Bayesian Multivariate models - Part II
# Models controling for z_svl
############################################

# The interaction models first. 
        tim_emerge_ap_int  <- bf(Time_emerge_sec    | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
         tim_snout_ap_int  <- bf(Time_snout_sec     | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
         dist_move_ap_int  <- bf(Distance.moved     | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
             speed_per_int <- bf(logspeed_1m        | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
       speed_burst_per_int <- bf(logspeed_burst     | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()

    # Delicata
      if(rerun){
        
        deli_behav_int <- brms::brm(tim_emerge_ap_int + tim_snout_ap_int + dist_move_ap_int + speed_per_int + speed_burst_per_int + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, data = dat2, control = list(adapt_delta = 0.98))

        saveRDS(deli_behav_int, "./output/models/deli_behav_int.rds")
        
        } else{
          deli_behav_int <- readRDS("./output/models/deli_behav_int.rds")
        }
      
        deli_svl <- build_behav_table(deli_behav_int)
    
      # Model checking
        pp_check(deli_behav_int, resp = "Timeemergesec")
        brms_model_check_res(deli_behav_int)


    # Guichenoti

        #all but distance moved loged. Differnet from deli
         tim_emerge_ap_int  <- bf(Time_emerge_sec    | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
          tim_snout_ap_int  <- bf(Time_snout_sec     | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
          dist_move_ap_int  <- bf(Distance.moved     | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
              speed_per_int <- bf(logspeed_1m        | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
        speed_burst_per_int <- bf(logspeed_burst     | mi() ~ 1 + temp*egg_treat + z_svl + (1|q|id) + (1|clutch)) + gaussian()
        
        if(rerun){
          guich_mv_int <- brms::brm(tim_emerge_ap_int + tim_snout_ap_int + dist_move_ap_int + speed_per_int + speed_burst_per_int + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, save_pars = save_pars(), file = "output/models/guich_mv_int", file_refit = "on_change", control = list(adapt_delta = 0.98), data = dat3)

          saveRDS(guich_mv_int, "./output/models/guich_mv_int.rds")
          
        } else{
          guich_mv_int <- readRDS("./output/models/guich_mv_int.rds")
        }
      
       guich_svl <- build_behav_table(guich_mv_int)

       # Model checking
        pp_check(guich_mv_int, resp = "Timeemergesec")
        brms_model_check_res(guich_mv_int)


############################################
# Bayesian Multivariate models - Part III
# Models without controlling for z_svl
############################################

    # Deli
         tim_emerge_ap_int1  <- bf(Time_emerge_sec    | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
          tim_snout_ap_int1  <- bf(Time_snout_sec      | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
          dist_move_ap_int1  <- bf(Distance.moved     | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
              speed_per_int1 <- bf(logspeed_1m        | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
        speed_burst_per_int1 <- bf(logspeed_burst     | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
        
        if(rerun){
          deli_behav_int_nonSVL <- brms::brm(tim_emerge_ap_int1 + tim_snout_ap_int1 + dist_move_ap_int1 + speed_per_int1 + speed_burst_per_int1 + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, file = "output/models/deli_behav_int_nonSVL", file_refit = "on_change", data = dat2, control = list(adapt_delta = 0.98))
          saveRDS(deli_behav_int_nonSVL, "./output/models/deli_behav_int_nonSVL.rds")
        } else{
          deli_behav_int_nonSVL <- readRDS("./output/models/deli_behav_int_nonSVL.rds")
        }
        
        deli_behav <- build_behav_table(deli_behav_int_nonSVL)

        # Model checking
        pp_check(deli_behav_int_nonSVL, resp = "Timeemergesec")
        brms_model_check_res(deli_behav_int_nonSVL)

    #Guich
        #2-way
         tim_emerge_ap_int2  <- bf(Time_emerge_sec       | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
          tim_snout_ap_int2  <- bf(Time_snout_sec        | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
          dist_move_ap_int2  <- bf(Distance.moved        | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
              speed_per_int2 <- bf(logspeed_1m           | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
        speed_burst_per_int2 <- bf(logspeed_burst        | mi() ~ 1 + temp*egg_treat + (1|q|id) + (1|clutch)) + gaussian()
        
        if(rerun){
          guich_mv_int_nonSVL <- brms::brm(tim_emerge_ap_int2 + tim_snout_ap_int2 + dist_move_ap_int2 + speed_per_int2 + speed_burst_per_int2 + set_rescor(TRUE), iter = 4000, warmup = 1000, chains = 4, cores = 4, save_pars = save_pars(), file = "output/models/guich_mv_int_nonSVL", file_refit = "on_change", control = list(adapt_delta = 0.98), data = dat3)
          saveRDS(guich_mv_int_nonSVL, "./output/models/guich_mv_int_nonSVL.rds")
        } else{
          guich_mv_int_nonSVL <- readRDS("./output/models/guich_mv_int_nonSVL.rds")
        }
        
        guich_behav <- build_behav_table(guich_mv_int_nonSVL)  

        # Model checking
        pp_check(guich_mv_int_nonSVL, resp = "Timeemergesec")
        brms_model_check_res(guich_mv_int_nonSVL)                                    
        
############################################
########### Figures ########################
############################################
#morphol deli
pdf("output/figs/taildeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=Tail, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Tail length (mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/svldeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=SVL, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Snout vent length (SVL, mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/massdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=Weigth, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Mass (gr)') + xlab('Temperature')+
  theme_bw()
dev.off()


#morphol gucih
pdf("output/figs/tailguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=Tail, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962'))+
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Tail length, mm') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/svlguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=SVL, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962'))+
  labs(fill = "ps") +
  ylab('Snout vent length (SVL, mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/massguich.pdf", width = 4, height = 4.5, useDingbats = F)
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

pdf("output/figs/distancedeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=mean_dist, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Distance travelled (mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/1mdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=bxcx1mdeli, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Running velocity 1m (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/25cmdeli.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean2, aes(x=temp, y=bxcx25deli, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  labs(fill = "ps") +
  ylab('Running velocity 25cm (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()

##guich

pdf("output/figs/distanceguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=log(mean_dist), fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
  labs(fill = "ps") +
  ylab('Distance travelled (log, mm)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/1mguich.pdf", width = 4, height = 4.5, useDingbats = F)
ggplot(mean3, aes(x=temp, y=bxcx1guich, fill = maternal)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.05),alpha=0.8, size = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha=0.5) +
  scale_fill_manual(values=c('#287c8e', '#5ec962')) +
  labs(fill = "ps") +
  ylab('Running velocity 1m (boxcox)') + xlab('Temperature')+
  theme_bw()
dev.off()

pdf("output/figs/25cmguich.pdf", width = 4, height = 4.5, useDingbats = F)
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

