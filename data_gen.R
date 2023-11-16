packages <- c("survival", "nnet", "tidyverse", "data.table", "flexsurv", "parallel", 
              "doParallel", "geepack","here")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}


expit <- function(x) {1/(1+exp(-x))}

## This code generates data from a structural nested model 
## compatible with a marginal structural model.
## It is based on Jessica Young's algorithm, published here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3635680/
#
## This code, which extend's Young's algorithm to multiple outcomes
## was written by Erica Moodie, published here:
# https://www.ncbi.nlm.nih.gov/pubmed/24272681
##
### Data generation
##
n <- 1000 # Number of subjects
N <- 10 #number of intervals per subject
K <- 1 # Number of causes of death

## This is the matrix of parameters of interest, possibly different
## at each interval
psi.mat <- matrix(0, nrow=K, ncol=N+1)

##Here are the effect sizes for the K=1 cause
psi.mat[1, ] <- log(1)

##Here the (untreated) all-cause rate is set to lambda=0.01, with
##lambda/K per cause; muK=lambda is used in the algorithm.
lambda <- 0.2
gamma.vec <- rep(log(lambda/K))
muK <- sum(exp(gamma.vec)) #So this bit is necessary to deal with that inequality in algorithm?
A<-J<-M<-L<-ID<-Y<-Z<-Tv<-Int<-ALast<-LLast<-LFirst<-numeric()
T0.vec<-T.vec<-Y.vec<-Z.vec <- rep(0, n)

##Here are the coefficients determining the
##mediation and treatment assignment mechanisms.
bevec <- c(-2,rep(log(1), 3)) #Used to generate time-varying confounder L
alvec <- c(-1,rep(log(2), 3)) #Used to generate exposure (Intercept, L, LLast, ALast)

##cval is used as in Young's algorithm to introduce the confounding
cval <- 30

##Begin the data-generation loop
simulation <- function (exposure) {
    for (i in 1:n) {
      ##Generate the counterfactual (untreated) survival time
      T0 <- rexp(1, lambda) #Generate T0 from an exponential dist with constant rate=lamba
      Ival <- as.numeric(T0 < cval)
      ##Begin the interval-by-interval simulation
      m <- 0
      mu.tot <- 0
      A.vec<-L.vec<-ALast.vec<-LLast.vec<-LFirst.vec<-rep(0, N+1)
      ##Implement Young's algorithm with multiple causes
      ##Generate the survival time, then the cause
      while (muK*T0 > mu.tot & m <= N) {
        if (m == 0) {
          ##First interval
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*0 + bevec[4]*0
          pval <- 1 / (1 + exp(-eta))
          L.vec[m+1] <- rbinom(1, 1, pval)
          
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*0 + alvec[4]*0 #design matrix
          pval <- 1 / (1 + exp(-eta))
          if (is.null(exposure)) {A.vec[m + 1] <- rbinom(1, 1, pval)} #generating observed A based on binomial distr
          else {A.vec[m+1] <- exposure}
          
          ALast.vec[m+1] <- 0; LLast.vec[m+1] <- 0;
          LFirst.vec <- rep(L.vec[m+1], N + 1)
          
        } else {
          ##Subsequent intervals
          eta <- bevec[1] + bevec[2]*Ival + bevec[3]*A.vec[m] + bevec[4]*L.vec[m]
          pval <- 1 / (1 + exp(-eta))
          L.vec[m+1] <- rbinom(1, 1, pval) 
          
          eta <- alvec[1] + alvec[2]*L.vec[m+1] + alvec[3]*L.vec[m] + alvec[4]*A.vec[m]
          pval <- 1 / (1 + exp(-eta)) #A affected by L at this time point, last point and A at last time point
          if (is.null(exposure)&A.vec[m]==0) {A.vec[m+1] <- rbinom(1, 1, pval)} #changed exposure to absorbing state
          else if (is.null(exposure)&A.vec[m]==1) {A.vec[m+1] <- 1}
          else if (!is.null(exposure)) {A.vec[m+1] <- exposure}
          ALast.vec[m+1] <- A.vec[m]; LLast.vec[m+1] <- L.vec[m];
        }
        
        muval <- sum(exp(gamma.vec + A.vec[m+1]*psi.mat[ , m+1]))
        
        ##Tval is computed for each interval, but is overwritten
        ##until the final interval
        Tval <- m + (muK * T0 - mu.tot) / muval
        mu.tot <- mu.tot + muval
        m <- m + 1
      }
      
      ##After exiting the loop, the survival time has been generated as Tval
      ##Now need to generate the failure type.
      if (m > N) {
        ##In the case of censoring at tenth interval, no failure.
        Tval <- m - 1
        Z.vec[i] <- 0
      } else {
        ##In the case of failure, use the ratio hazards to define the
        ##relevant multinomial distribution on the K causes.
        Z.vec[i] <- sample(c(1:K), 1, prob = exp(gamma.vec + A.vec[m]*psi.mat[ ,m])) # I don't really get this step
      }
      
      ##Store the outcomes
      T0.vec[i] <- T0
      T.vec[i] <- Tval
      Y.vec[i] <- m - 1
      ID <- c(ID, rep(i,m)) #Individual
      Int <- c(Int, c(1:m)) #Time point
      A <- c(A, A.vec[1:m]) #Time-updated treatment
      L <- c(L, L.vec[1:m]) #Time-updated covariate L
      ALast <- c(ALast, ALast.vec[1:m]) #Treatment at last t point
      LLast <- c(LLast, LLast.vec[1:m]) #Covariate L at last t point
      LFirst <- c(LFirst, LFirst.vec[1:m]) #Baseline covariate L value
      Z <- c(Z, rep(0,m - 1), Z.vec[i]) #Outcome: Z>0 indicates outcome of some type, determined by value)
      tv <- c(1:m); tv[m] <- Tval
      Tv <- c(Tv, tv) #If event occurs, exact time at which it occurred; o.w. equal to Int)
      
    }
    
    DeathsK.df <- data.frame(ID, Int, Tv, A, ALast, L, LLast, LFirst, Z)
    
    ##Trim off the intervals beyond the Nth (loop goes one too far)
    DeathsK.df <- DeathsK.df[DeathsK.df$Int <= N, ]
    DeathsK.df$Int0 <- DeathsK.df$Int - 1
    
  
  return(DeathsK.df)
}


plot_dat <- simulation(exposure=NULL) %>% 
  mutate(first_exposure = as.numeric(A+ALast == 1),
         last_id = as.numeric(Int==10|Z==1))

p0 <- plot_dat %>% 
  filter(A==0) %>% 
  group_by(ID) %>% 
  mutate(max_time = max(Tv),
         min_time = min(Int0)) %>% 
  group_by(ID) %>% 
  mutate(last_id = !duplicated(ID, fromLast = T)) %>% 
  filter(last_id == 1) %>% 
  select(ID, A, Z, min_time, max_time)

p1 <- plot_dat %>% 
  filter(A==1) %>% 
  group_by(ID) %>% 
  mutate(max_time = max(Tv),
         min_time = min(Int0)) %>% 
  filter(last_id == 1) %>% 
  select(ID, A, Z, min_time, max_time)

p <- rbind(p0,p1) %>% 
  arrange(ID)  %>% 
  group_by(ID) %>% 
  mutate(last_id = !duplicated(ID, fromLast = T))

ggplot() + 
  geom_segment(data = p[p$ID<21,], 
               aes(x=min_time, 
                   xend=max_time, 
                   y = ID, 
                   yend=ID, color=factor(A))) +
  geom_point(data = p[p$ID<21&p$last_id==1,], 
             aes(x = max_time, 
                 y = ID, 
                 shape = factor(Z))) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed") +
  theme_classic() + 
  # scale_x_continuous(expand=c(0,0)) +
  # scale_y_continuous(expand=c(0,0)) + 
  xlab("Time on Study")

# create simulation to be used in analysis steps
s <- 50
sim <- s
all_res <- data.frame(Oracle = numeric(sim),
                      IPTW.PL = numeric(sim),
                      IPTW.COX = numeric(sim),
                      CCW1 = numeric(sim),
                      CCW2 = numeric(sim),
                      IPTW.TF1 = numeric(sim),
                      IPTW.TF2 = numeric(sim))
start <- Sys.time()

for (i in 1:sim){
  print(i)
  set.seed(1 + i)
  ###################################################
  # Oracle Effect
  ###################################################
  d1 <- simulation(exposure=1)
  d0 <- simulation(exposure=0)
  oracle_df <- bind_rows(d1, d0)
  oracle <- coxph(Surv(Int0, Tv, Z)  ~ A + ALast, data=oracle_df, ties="efron") 
  all_res[i, "Oracle"] <- coef(oracle)[1]
  
  ###################################################
  ## IPTW without bootstrap
  ###################################################
  #Denominator of weights
  d <- simulation(exposure=NULL)
  denominator <- rep(NA, nrow(d))
  logit <- predict(glm(A ~ ALast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=d))
  denominator[d$A == 1] <- expit(logit[d$A == 1])
  denominator[d$A == 0] <- 1 - expit(logit[d$A == 0])
  
  #Numerator of weights
  numerator <- rep(NA, nrow(d))
  logit <- predict(glm(A ~ ALast + as.factor(Int), family=binomial, data=d))
  numerator[d$A == 1] <- expit(logit[d$A == 1])
  numerator[d$A == 0] <- 1 - expit(logit[d$A == 0])
  wt <- unlist(tapply(numerator / denominator, d$ID, cumprod))
  
  #IP-weighted Pooled Logistic Model
  iptw.pl <- suppressWarnings(glm(Z ~ A + ALast + as.factor(Int), data=d, family=binomial(link="logit"), weights=wt))
  all_res[i, "IPTW.PL"] <- coef(iptw.pl)[1]
  
  #IP-weighted Cox Model
  iptw.cox <- coxph(Surv(Int0, Tv, Z)  ~ A + ALast, data=d, weights=wt, ties="efron", timefix=FALSE)
  all_res[i, "IPTW.COX"] <- coef(iptw.cox)[1]
  
  ###################################################
  ## Clone-Censor-Weighting Estimator: Approach 1
  ###################################################
  #Clone just those with a change in treatment from unexposed to exposed (0 to 1), 
      #create a clone that is set to 1 at all times prior to the interval where treatment(0-->1) (IDA)
  d_0 <- d %>% filter(Int == 1 & A == 0)
  
  d_clone01 <- d %>% filter(ID %in% d_0$ID) %>%
    group_by(ID) %>%
    mutate(Asum = cumsum(A),
           T01 = if_else(Asum == 1, Int, NA)) %>%   #T01 is the interval where treatment flips from 0 to 1
    filter(!is.na(T01)) %>% 
    select(ID, T01) %>%
    ungroup() %>%
    left_join(d, by = "ID") %>% #merge back to original data
    mutate(AClone = if_else(Int < T01, 1, A),
           IDClone = paste0(ID, "A"),
           ZClone = Z)
  
  #Generate clone & censored data (duplicate above, but when interval T01 is reached, censor)
  d_clone00 <- d_clone01 %>%
    group_by(ID) %>%
    filter(Int <= T01) %>%
    mutate(AClone = 0,
           IDClone = paste0(ID, "B"),
           ZClone = if_else(row_number()==n(), 2, ZClone)) %>%
    ungroup()
  
  #Bind data together
  d_cloned <- bind_rows(d_clone01, d_clone00)
  
  d_not_cloned <- d %>% filter(!(ID %in% d_cloned$ID)) %>%
    mutate(ZClone = Z)
  
  d_cloned_full <- bind_rows(d_cloned, d_not_cloned)
  
  d_cloned <- d_cloned_full %>% mutate(AClone = if_else(is.na(AClone), A, AClone),
                                  ACloneLast = lag(AClone),
                                  ACloneLast = if_else(is.na(ACloneLast), AClone, ACloneLast),
                                  event = if_else(ZClone == 2, 0,
                                                  if_else(ZClone == 0, 1, 2)),
                                  censor = if_else(ZClone == 2, 1, 0),
                                  death = if_else(ZClone == 1, 1, 0),
                                  IDClone = if_else(is.na(IDClone), paste0(ID, "A"), IDClone))
  
  #Inverse probability of censoring weighting analysis
  
  d_cloned$den_d <- glm((censor==0) ~ AClone + ACloneLast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=d_cloned)$fitted.values
  
  d_cloned$num_d <- glm((censor==0) ~ as.factor(Int), 
                        family=binomial(link="logit"), data=d_cloned)$fitted.values
  
  
  # Build time-specific weights
  # wt=0 when someone is censored
  d_cloned$wt_d <- ifelse(d_cloned$censor==1, 0, d_cloned$num_d/d_cloned$den_d)
  
  # Multiply censoring weights across time
  d_cloned <- d_cloned %>% 
    group_by(IDClone) %>%  
    mutate(cum_wt=cumprod(wt_d)) %>% 
    ungroup(IDClone) 

  d_cloned <- d_cloned %>% filter(cum_wt > 0)
  
  # Analysis applying censoring weights
  ipcw.cox <- coxph(Surv(Int0, Tv, death)  ~ AClone + ACloneLast, data=d_cloned, weights=cum_wt, ties="efron", timefix=FALSE)
  
  all_res[i, "CCW1"] <- coef(ipcw.cox)[1]
  
  ###################################################
  ## Clone-Censor-Weighting Estimator: Approach 2
  ###################################################
  # assumes that the immortal person time is fixed for everyone (i.e. all time prior to interval 1 is misclassified)
  d_cloned <- d %>% filter(Int == 1)
  d_cloned <- d_cloned %>% mutate(AClone = if_else(A == 0, 1,
                                             if_else(A == 1, 0, NA)),
                            ZClone = if_else(Z == 0, 2, Z),
                            IDclone = paste0(ID, "B"))
  
  d_cloned <- bind_rows(d, d_cloned)
  
  d_cloned <- d_cloned %>% mutate(AClone = if_else(is.na(AClone), A, AClone)) %>%
    group_by(IDclone) %>%
    mutate(ACloneLast = lag(AClone),
                                  ACloneLast = if_else(is.na(ACloneLast), AClone, ACloneLast),
                                  ZClone = if_else(is.na(ZClone), Z, ZClone),
                                  death = if_else(ZClone == 1, 1, 0),
                                  censor = if_else(ZClone == 2, 1, 0),
                                  IDclone = if_else(is.na(IDclone), paste0(ID, "A"), IDclone)) %>%
    ungroup()
  
  #T01 is the interval where treatment flips from 0 to 1
  d_clone01 <- d %>%
    group_by(ID) %>%
    mutate(Asum = cumsum(A),
           T01 = if_else(Asum == 1, Int, NA)) %>%
    filter(!is.na(T01)) %>%
    select(ID, T01) %>%
    ungroup() %>%
    right_join(d, by = "ID") %>%
    mutate(AClone = if_else(Int < T01, 1, A),
           IDClone = paste0(ID, "A"),
           ZClone = Z)

  #Generate clone & censored data
  d_clone00 <- d_clone01 %>%
    group_by(ID) %>%
    filter(Int <= T01) %>%
    mutate(AClone = 0,
           IDClone = paste0(ID, "B"),
           ZClone = if_else(row_number()==n(), 2, ZClone)) %>%
    ungroup()

  #Bind data together
  d_cloned <- bind_rows(d_clone01, d_clone00)

  d_cloned <- d_cloned %>% mutate(AClone = if_else(is.na(AClone), A, AClone)) %>%
  group_by(IDClone) %>%
                                  mutate(ACloneLast = lag(AClone),
                                  ACloneLast = if_else(is.na(ACloneLast), AClone, ACloneLast),
                                  event = if_else(ZClone == 2, 0,
                                                  if_else(ZClone == 0, 1, 2)),
                                  censor = if_else(ZClone == 2, 1, 0),
                                  death = if_else(ZClone == 1, 1, 0))

  #Inverse probability of censoring weighting analysis

  d_cloned$den_d <- glm((censor==0) ~ AClone + ACloneLast + L + LLast + as.factor(Int), family=binomial(link="logit"), data=d_cloned)$fitted.values
  
  d_cloned$num_d <- glm((censor==0) ~ as.factor(Int), 
                    family=binomial(link="logit"), data=d_cloned)$fitted.values
  
  
  # Build time-specific weights
  # wt=0 when someone is censored
  d_cloned$wt_d <- ifelse(d_cloned$censor==1, 0, d_cloned$num_d/d_cloned$den_d)
  
  # Multiply censoring weights across time
  d_cloned <- d_cloned %>% 
    group_by(IDClone) %>%  
    mutate(cum_wt=cumprod(wt_d)) %>% 
    ungroup(IDClone) 
  
  d_cloned <- d_cloned %>% filter(cum_wt > 0)
  # Analysis applying censoring weights
  iptw.pl <- suppressWarnings(glm(death ~ AClone + ACloneLast + as.factor(Int), data=d_cloned, family=binomial(link="logit"), weights=cum_wt))
  ipcw.cox <- coxph(Surv(Int0, Tv, death)  ~ AClone + ACloneLast, data=d_cloned, weights=cum_wt, ties="efron", timefix=FALSE)
  
  all_res[i, "CCW2"] <- coef(ipcw.cox)[1]
  
  ###################################################
  ## Time-Fixed 
  ## Scenario 1 (Everyone's exposure is measured at time 1)
  ## Immortal Person Time
  ###################################################
  #Generate Time-Fixed Data (where exposure is determined at Time 1)
  tte <- d %>% 
    mutate(first_exposure = as.numeric(A+ALast == 1),
           last_id = as.numeric(Int==N|Z==1))
  d0 <- tte %>% 
    filter(A==0) %>% 
    group_by(ID) %>% 
    mutate(max_time = max(Tv),
           min_time = min(Int0)) %>% 
    group_by(ID) %>% 
    mutate(last_id = !duplicated(ID, fromLast = T)) %>% 
    filter(last_id == 1) %>% 
    select(ID, A, Z, min_time, max_time)
  
  d1 <- tte %>% 
    filter(A==1) %>% 
    group_by(ID) %>% 
    mutate(max_time = max(Tv),
           min_time = min(Int0)) %>% 
    filter(last_id == 1) %>% 
    select(ID, A, Z, min_time, max_time)
  
  d_temp <- rbind(d0,d1) %>% 
    arrange(ID)  %>% 
    group_by(ID) %>% 
    mutate(ident = 1,
           Int = cumsum(ident)) %>% 
    mutate(last_id = !duplicated(ID, fromLast = T))
  
  d_a1 <- d_temp %>% filter(Int == 1) %>% select(ID,A) %>% rename(A1=A)
  
  d1_tf <- left_join(d_temp, d_a1) %>% select(-A) %>% rename(A=A1)
  
  max_time <- d1_tf %>% group_by(ID) %>% filter(row_number()==n()) %>% select(ID, max_time) %>% ungroup()
  
  d1_tf <- d1_tf %>% group_by(ID) %>% filter(row_number()==1) %>% select(-c(max_time)) %>% ungroup() %>% left_join(max_time, by = "ID")
    
  #Time-Fixed IPW (with exposure values from Time 1)
  d1_tf$ps <- glm(A ~ 1, data = d1_tf, family = binomial("logit"))$fitted.values 
  
  #Stabilized weights
  d1_tf$sw <- (mean(d1_tf$A)/d1_tf$ps)*d1_tf$A + ((1 - mean(d1_tf$A))/(1-d1_tf$ps))*(1-d1_tf$A)
  summary(d1_tf$sw)
  
  iptw.cox <- coxph(Surv(min_time, max_time, Z)  ~ A, data=d1_tf, weights=sw, ties="efron", timefix=FALSE)
  all_res[i, "IPTW.TF1"] <- coef(iptw.cox)[1]
  
  ###################################################
  ## Time-Fixed 
  ## Scenario 2 (Everyone's exposure is measured at time 2)
  ## Immortal Person Time
  ###################################################
  #Generate Time-Fixed Data (where exposure is determined at Time 2)
  d_a2 <- d_temp %>% filter(last_id == TRUE) %>% select(ID,A) %>% rename(A2=A)
  
  d2_tf <- left_join(d_temp, d_a2) %>% select(-A) %>% rename(A=A2)
  
  d2_tf <- d2_tf %>% filter(last_id == TRUE) %>% mutate(min_time = 0)
  
  #Time-Fixed IPW (with exposure values from Time 1)
  d2_tf$ps <- glm(A ~ 1, data = d2_tf, family = binomial("logit"))$fitted.values 
  
  #Stabilized weights
  d2_tf$sw <- (mean(d2_tf$A)/d2_tf$ps)*d2_tf$A + ((1 - mean(d2_tf$A))/(1-d2_tf$ps))*(1-d2_tf$A)
  summary(d2_tf$sw)
  
  iptw.cox <- coxph(Surv(min_time, max_time, Z)  ~ A, data=d2_tf, weights=sw, ties="efron", timefix=FALSE)
  all_res[i, "IPTW.TF2"] <- coef(iptw.cox)[1]
  #robust variance code for iptw.cox notes
  #add bootstrap function for the robust variance estimators
}

end <- Sys.time()

run_time <- end-start

#Summarize over simulations
res.sum <- data.frame(
  method = c("Truth", "Oracle", "Cox MSM", "CCW1", "CCW2", "Time-Fixed 1", "Time-Fixed 2"),
  average = rep(NA, 7), 
  bias = rep(NA, 7),
  mse = rep(NA, 7))

res.sum <- res.sum %>%
  mutate(
    average = case_when(
      method == "Truth" ~ log(1),
      method == "Oracle" ~ mean(all_res$Oracle),
      method == "Cox MSM" ~ mean(all_res$IPTW.COX),
      method == "CCW1" ~ mean(all_res$CCW1),
      method == "CCW2" ~ mean(all_res$CCW2),
      method == "Time-Fixed 1" ~ mean(all_res$IPTW.TF1),
      method == "Time-Fixed 2" ~ mean(all_res$IPTW.TF2),
      TRUE ~ NA_real_  
    ),
    average = round(average, digits = 3),
    bias = average - log(1),
    sd = case_when(
      method == "Oracle" ~ sd(all_res$Oracle),
      method == "Cox MSM" ~ sd(all_res$IPTW.COX),
      method == "CCW1" ~ sd(all_res$CCW1),
      method == "CCW2" ~ sd(all_res$CCW2),
      method == "Time-Fixed 1" ~ sd(all_res$IPTW.TF1),
      method == "Time-Fixed 2" ~ sd(all_res$IPTW.TF2),
      TRUE ~ NA_real_  
    ),
    mse = sd^2 + bias^2
  )

res.sum
