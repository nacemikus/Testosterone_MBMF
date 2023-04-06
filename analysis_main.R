# =============================================================================
#### Info #### 
# =============================================================================
# hierarchical model
run_model_fit <- function(modelfile, savemodelname, InfType = "Sampling") {
   # =============================================================================
   #### variables #### 
   # =============================================================================
   # modelfile = name of stan script
   # savemodelname = name to save
   # InfType: "Sampling" - mcmc, "VB", var bayes, "Optim" - point est, "try" - try with 2 mcmc samples
   # =============================================================================
   #### set priors #### 
   # =============================================================================
   # simulate_data = as.integer(0)
   # prior_mu = c(0,0,0)#,0)
   # prior_sd = c(5,5,5)#,5)
   # prior_sd_random = 1
    # =============================================================================
    #### Construct Data #### 
    # =============================================================================
    # clear workspace
    # rm(list=ls())
    library(rstan)
    library(ggplot2)
    # library(R.matlab)
    library(tidyr)
    library(dplyr)
  
   
    no_subjects = NA
    belief = 0
    
### load data
data_beh <-readRDS(file = "Group Data/data_beh.rds")
data_group <- readRDS(file = "Group Data/data_group.rds")

data_group$comt_z <- as.numeric(scale(data_group$comt))
data_group$cag_z <- as.numeric(scale(data_group$cag))

# View(data_beh)
data_beh$trials <- data_beh$trial +1
# dim(data_beh)
data_beh <- data_beh %>% filter(choice1 != -1)
### prepare data
subjList <- unique(data_beh$ID)
data_group <- data_group[data_group$ID %in% subjList,  ]

if (belief == 1) {
   ExcludeSubj <- subjList[is.na(data_group$belief_in_T)]
   # ExcludeSubj <- union(ExcludeSubj, subjList[is.na(ankk)])
   
   subjList <- setdiff(subjList, ExcludeSubj)
   
   ### remove subjects with no genetic data
   data_beh <- data_beh[data_beh$ID %in% subjList,  ]
   data_group <- data_group[data_group$ID %in% subjList,  ]
   print("removed sub with no belief data")
}

if (!is.na(no_subjects)) { 
subjList <- subjList[1:no_subjects]
data_beh <- data_beh[data_beh$ID %in% subjList,  ]

cat("running only ", no_subjects, " subjects\n")
}

# ankk <- data_beh$ankk[data_beh$trials==9]
# darpp <- data_beh$darpp[data_beh$trials==9]
# oprm <- data_beh$oprm[data_beh$trials==9]
# oprm <- as.numeric(oprm)

subjList <- unique(data_beh$ID)


cat("running on ", length(subjList), " subjects\n")
# sum(comt==2, na.rm = T)
# #
# dat<- data_beh$dat[data_beh$trials==9]
# dat <- (as.numeric(dat==3))
# dat <- (as.numeric(dat==0))
# comt <- data_beh$comt[data_beh$trials==9]
# comt <- (as.numeric(comt==3))
# comt <-(comt - mean(comt))/sd(comt)



# comt <- (as.numeric(comt==3))
# 
# subjList <- unique(data_beh$s)
numSub <- length(subjList)




Tsubj <- as.vector(rep(0, numSub))
for (ss in 1:numSub) {
  Tsubj[ss] <- length(data_beh$trials[data_beh$ID == subjList[ss]]);
}


maxTrials <- 200 # max(Tsubj)

level1_choice <- array(1,c(numSub, maxTrials))
stim_left <- array(1,c(numSub, maxTrials))
reward <- array(1,c(numSub, maxTrials))
state1 <- array(1,c(numSub, maxTrials))
state2 <-array(1,c(numSub, maxTrials))

for (i in 1:numSub) {
  tmp <- subset(data_beh, data_beh$ID==subjList[i])
  level1_choice[i,1:Tsubj[i]] = tmp$choice1
  stim_left[i,1:Tsubj[i]] = tmp$stim_left
  reward[i,1:Tsubj[i]] = tmp$points
  state1[i,1:Tsubj[i]] = tmp$state1
  state2[i,1:Tsubj[i]] = tmp$state2
 
}



# # drug <-factor(data_beh_drug$drug, level = c(1,2,3), labels = c("ami", "nal", "pla"))
testosterone <- data_group$admin == "Testosterone"
testosterone<- as.numeric(testosterone)
comt  <- data_group$comt_z
cag <- data_group$cag_z
dat1<- as.numeric(data_group$dat1 == "other") -1/2

belief_in_T <- data_group$belief_in_T
# belief_in_T <- as.numeric(belief_in_T)

# # length(amisulpride)
# amisulpride <- data_beh_drug$drug == 1 # amisul
# amisulpride<- as.numeric(amisulpride)
# serum <- data_beh$serum[data_beh$trials==1] # Nal
# serum[is.na(serum)] <- 0
# serum[amisulpride==1]
dataList <- list(N = numSub, T = maxTrials, Tsubj = Tsubj, level1_choice = level1_choice, stim_left = stim_left,
                 reward = reward, state1 = state1, state2 = state2,testosterone= testosterone, belief_in_T=belief_in_T,
                 comt = comt, dat1=dat1, cag =cag)
#                 simulate_data = simulate_data, prior_mu = prior_mu, prior_sd  = prior_sd, prior_sd_random=prior_sd_random)


    # =============================================================================
    #### Running Stan #### 
    # =============================================================================
    rstan_options(auto_write = TRUE)
    options(mc.cores = 4)
    
    if (InfType == "try") {
       nIter     <- 2
       nChains   <- 1
       nWarmup   <- 1 # floor(nIter/2)
       nThin     <- 1
       InfType = "Sampling"
    } else {
    nIter     <- 2000
    nChains   <- 4
    nWarmup   <- 800 # floor(nIter/2)
    nThin     <- 1
    }
    cat("Estimating", modelfile, "model... \n")
    
    startTime = Sys.time(); print(startTime)
  
    modelcode <-  stan_model(modelfile)
    
    if (InfType == "VB") {
       print('Using Variational Bayesian Inference')
      
    
      
        fit_rl <- vb(modelcode, 
                      data    = dataList, 
                      iter =  50000,
                      init    = "0",
                      seed    = 1450154637,
                     tol_rel_obj = 0.01
                      )
       
    } else if (InfType == "Optim") {
       print('Optimizing for point estimation only')
         fit_rl <- optimizing(modelcode, 
                              data    = dataList,
                              init    = "random",
                              verbose = TRUE
                              )
       
    } else if (InfType == "Sampling") {  # sampling (by  default) 
       cat("Calling", nChains, "simulations in Stan... \n")
       
       fit_rl <-  sampling(modelcode, 
                           data = dataList,
                           chains = nChains,
                           iter = nIter,
                           warmup = nWarmup,
                           thin = nThin,
                           init = "0",
                           seed = 1450154637,
                           control = list(adapt_delta = 0.9, max_treedepth=10)
       )
    }
    # fit_rl <- stan(modelfile, 
    #                data    = dataList, 
    #                chains  = nChains,
    #                iter    = nIter,
    #                warmup  = nWarmup,
    #                thin    = nThin,
    #                init    = "0", # "random",
    #                seed    = 1450154637,
    #                control = list(
    #                  adapt_delta = 0.90)
    #                )
    
   
    
    # 
    # stepsize = 2.0,
    # max_treedepth = 10
    cat("Finishing", modelfile, "model simulation ... \n")
    endTime = Sys.time(); print(endTime)
    cat("It took",as.character.Date(endTime - startTime), "\n")
    cat("Saving in ", savemodelname, "... \n")
    saveRDS(fit_rl, file = savemodelname)
    
    
}
