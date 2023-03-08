# Forage Fish Aggregation Models: Posterior Predictive Check -------------------

# Workspace Set-up -------------------------------------------------------------

#Set the directories I will be using

input_dir<-"P:/NYSERDA Ecosystem Dynamics/Forage Fish Modeling/Aerial Survey Forage Fish Aggregations/data"
fixed_grid_model<-"P:/NYSERDA Ecosystem Dynamics/Forage Fish Modeling/Aerial Survey Forage Fish Aggregations/output/updated fitted_models/Fitted Model Correct Pgrid"
kfold_cv_updt_grd_dir<-"P:/NYSERDA Ecosystem Dynamics/Forage Fish Modeling/Aerial Survey Forage Fish Aggregations/output/k_fold_cv_updt_grd"
ppc_dir<-"P:/NYSERDA Ecosystem Dynamics/Forage Fish Modeling/Aerial Survey Forage Fish Aggregations/output/ppc_updt_grd_model"

draft_fig_dir<-"P:/NYSERDA Ecosystem Dynamics/Forage Fish Modeling/Draft Figures"


# Load libraries

library(nimble)
library(sf)
library(rgdal)
library(tidyverse)
library(reshape2)
library(spdep)
library(tibble)
library(snow)
library(coda)
library(data.table)
library(narray)
library(nlist)
library(stringr)
library(scales)




# Load nimble model components --------------------------------------------------

# save(dat, yy, obs, eff, cgrid, adj.sparse, num, bou, C, M, cgrid.link, sea.key, sea1, sea2,
#      sea3, cov.sum, cov.set, fprob, fmean, params, const, inits, mod,
#      file = paste0(pcc_dir, "/data_for_pcc_updt_grd_fprob_model"))

#'[Note: data saved here is for the fprob nonjoint model]
#load(file= paste0(pcc_dir, "/data_for_pcc_updt_grd_fprob_model"))

load(file= paste0(fixed_grid_model, "/updt_grd_data_for_model_run.Rdata"))


## regenerate R objects const, dat, inits, and mod -----------------------------

# #take the grid and make a CAR object out of it for nimble
# #create adjacency list from polygon data
# 
# adj <- poly2nb(cgrid)
# adjm <- nb2mat(adj, style = 'B')
# adj.sparse <- unlist(adj)
# num <- unlist(lapply(adj, length))
# 
# #now create the tools I need to parameterize a CAR in nimble using spdep
# CM2<-as.carCM(adj.sparse, rep(1, length(adj.sparse)), num)
# C<- CM2$C
# M <- CM2$M
# bou<-carBounds(C, adj.sparse, num, M)
# # save(bou, file = 'ffa_bou_cgrid.Rdata')
# #load('ffa_bou_cgrid.Rdata')
# mumod<-rep(0,length(num))
# 
# select covariates for model 
selectcovs <- c('log_depth_4km', 'plcurv_2km',  'rugosity', 'bpi', 'sediment', 'slope_2km', 'sst', 'chl', 'ssha', 'fsle', 'salinity', 'mld', 'uo', 'vo', 'sst_fprob', 'sst_fmean', 'chl_fprob', 'chl_fmean', 'sst_anomaly')

#covs that we used, make sure the same vector is used in the const line (350)

#'[set which covariate set to use: "fprob" or "fmean"]
cov.set<-"fprob"

fprob<-c(1, 4, 5, 8, 10, 11, 12, 15, 17, 19)
fmean<-c(1, 4, 5, 8, 10, 11, 12, 16, 18, 19)

selectcovs[ switch(cov.set, "fprob" = fprob, "fmean" = fmean)]


#dat <- list(y = yy, size = obs$Area_sqm)
#to check logged data Gaussian size model use dat below
dat <- list(y = yy, size = obs$log_Area_sqm ) #switching to the log transformed data

const <- list(ngrid = nrow(eff), ncgrid = nrow(cgrid), nobs = nrow(obs), L = length(adj.sparse), wei = rep(1, length(adj.sparse)), 
              adj = adj.sparse, num = num, bou = bou, C= C, M = M, gridID = cgrid.link$ID,
              obsID = obs$fullgridID, obsgridID = obs$ID, obsurvID = obs$surveyID, obssurvtype = obs$survtype,
              e = eff, surveyID = cgrid.link$survey, nsurvtype = 2, seaID = sea.key, sea1 = sea1, sea2 = sea2, sea3 = sea3,
              envcov = cov.sum[,,switch(cov.set, "fprob" = fprob, "fmean" = fmean)], nenvcov = 10,
              nsurv = 15)

# params <- c('a.psi', 'u', 'tauu', 'gammau', 'tau', 'u1', 'u2', 'mu0', 'b.pois', 'b.size', 'z',
#             'b1.p', 'b2.p', 'b3.p', 'lambda.star', 'mu')

inits <- list(u = structure(.Data= runif(const$ncgrid * 2, -.1, .1),.Dim=c(const$ncgrid, 2)),
              sdu = runif(2, 2, 4),
              mu0 = c(runif(1, -1, 0), runif(1, 3, 5)),
              a.psi = runif(1, -1, 1),
              b1.p = runif(1, -1, 1),
              b2.p = runif(1, -1, 1),
              b3.p = runif(1, -1, 1),
              b.pois = runif(const$nenvcov, -1, 1), #old: runif(const$nenvcov)
              b.size = runif(const$nenvcov)
              #r = runif(1, 1, 2) #check this with evan
)

#'[Note: Must choose by comment/uncomment the count and size model within the call to mod]
#'
mod <- nimbleCode(
  {
    #priors
    a.psi ~ dnorm(0, 0.8) #for the negbin mixture - changing these from 0,0.01 to 0,0.5 to improve convergence
    tau <- 1/sig^2
    sig ~ dgamma(2, 0.5)
    b1.p ~ dnorm(0, 0.5)
    b2.p ~ dnorm(0, 0.5)
    b3.p ~ dnorm(0, 0.5)
    b1.s ~ dnorm(0, 0.01)
    b2.s ~ dnorm(0, 0.01)
    b3.s ~ dnorm(0, 0.01)
    b0.p ~ dnorm(0, 0.5)
    b0.s ~ dnorm(0, 0.01)
    
    #add prior for negative Poisson - pick one
    #my version
    r ~ dunif(0,50) #tightening the r prior
    #r ~ dunif(0,250)
    
    #Evan's version 
    #'[Note: giving r values stuck near zero, which seems extreme - was suggested we use prior like above]
    #r ~ dgamma(0.01, 0.01)
    
    #beta cov priors #bumping up these priors from 0, 0.1 to 0,0.5
    for(n in 1:nenvcov){
      b.pois[n] ~ dnorm(0, 0.5)
      b.size[n] ~ dnorm(0, 0.4)
      
    }
    
    #CAR priors
    for(n in 1:2){
      u[1:ncgrid, n] ~ dcar_proper(muu[1:ncgrid, n], C[1:L], adj[1:L], num[1:ncgrid], M[1:ncgrid], tauu[n], gammau[n])
      tauu[n] <- 1/sdu[n]^2
      sdu[n] ~ dgamma(2, 0.5)
      gammau[n] ~ dunif(bou[1], bou[2])
      
      mu0[n] ~ dnorm(0, 0.1) #use this if we don't want the CAR to be zero-centered
      
      
      for(m in 1:ncgrid){
        muu[m, n] <- mu0[n]
        #muu[m, n] <- 0
      }#n
      
    }#m
    
    u1 <- mean(u[1:ncgrid, 1])
    u2 <- mean(u[1:ncgrid, 2])
    
    #likelihood statement for number of aggregations - must pick a model and comment out the others
    #models: ZIP, ZINB as poisson-gamma mixture model, ZINB with dnegbin with code conversion
    for(k in 1:ngrid){
      for(j in 1:nsurv){
        
        # #Zero-inflated Negative binomial (ZINB) via dnegbin, modifying p to represent dispersion 
        # y[k, j] ~ dnegbin(p[k, j], r)
        # p[k, j] <- r/(r + lambda.star[k,j]) - 1e-10*z[k, j]
        # lambda.star[k, j] <- lambda[k, j] * (1-z[k, j]) #this might be inverted from what we want 
        # log(lambda[k, j]) <- log(theta[k, j]) + e[k, j]
        # z[k, j] ~ dbern(psi[k, j])
        # log(theta[k, j]) <- u[gridID[k], 1] + inprod(envcov[k, j, 1:nenvcov], b.pois[1:nenvcov]) # + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        # logit(psi[k, j]) <- a.psi + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        # 
        
        #Zero-inflated Negative Binomial (ZINB) via poisson-gamma mixture model
        y[k, j] ~ dpois(lambda.star2[k, j])
        lambda.star2[k, j] <- lambda.star1[k, j] * rho[k, j]
        lambda.star1[k,j]<- lambda[k, j] * z[k, j]
        log(lambda[k, j]) <- log(theta[k, j]) + e[k, j]
        z[k, j] ~ dbern(psi[k, j])
        log(theta[k, j]) <- u[gridID[k], 1] + inprod(envcov[k, j, 1:nenvcov], b.pois[1:nenvcov]) # + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        logit(psi[k, j]) <- a.psi + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        rho[k, j] ~ dgamma(r, r)
        
        
        #negative binomial (no zero-inflation)
        # y[k, j] ~ dpois(lambda.star[k, j])
        # lambda.star[k, j] <- lambda[k, j] * rho[k, j] 
        # log(lambda[k, j]) <- log(theta[k, j]) + e[k, j]
        # #z[k, j] ~ dbern(psi[k, j])
        # log(theta[k, j]) <- u[gridID[k], 1] + inprod(envcov[k, j, 1:nenvcov], b.pois[1:nenvcov]) + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        # #logit(psi[k, j]) <- a.psi + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        # rho[k, j] ~ dgamma(r, r)
        
        #Zero-inflated Poisson (ZIP)
        # y[k, j] ~ dpois(lambda.star[k, j])
        # lambda.star[k, j] <- lambda[k, j] * z[k, j]
        # log(lambda[k, j]) <- log(theta[k, j]) + e[k, j]
        # z[k, j] ~ dbern(psi[k, j])
        # log(theta[k, j]) <- u[gridID[k], 1] + inprod(envcov[k, j, 1:nenvcov], b.pois[1:nenvcov]) # + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        # logit(psi[k, j]) <- a.psi + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]

        

        
      }#j
    }#k
    
    #size observation model - have to select which model we want to run now: OG Gaussian, Gaussian w/ log data, Lognormal, or Gamma
    for(i in 1:nobs){
      
      #lognormal model
      # size[i] ~ dlnorm(mu[i], tau) # * z[k, j]
      # #non-joint model
      # mu[i] <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
      
      #Gaussian model
      size[i] ~ dnorm(mu[i], tau) # * z[k, j]
      
      # #non-joint model - log tranformed data
      mu[i] <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
      
      # #non-joint model - OG model
      # #log(mu[i]) <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
      # 
      # #joint model
      # #log(mu[i]) <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.pois[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
      # 
      

    }#i
    
  })



# Set up the model and mcmc object ---------------------------------------------
#'[Note: for the fprob nonjoint model]
#'[Note: this is for the ppc sampler]

fm <- nimbleModel(code = mod, data = dat, constants = const, inits = inits) #define our model and initialize
fmc <- compileNimble(fm) #compile the model

# create vectors of node names 
# #nodes needed to simulate new datasets
dataNodes <- fm$getNodeNames(dataOnly = TRUE)
parentNodes <- fm$getParents(dataNodes, stochOnly = TRUE, returnScalarComponents = TRUE) #
cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n") 
#Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- fm$getDependencies(parentNodes, self = FALSE, downstream = TRUE, returnScalarComponents = TRUE)# self = TRUE, downstream = TRUE, returnScalarComponents = TRUE

parentNodes2<- fm$getParents(dataNodes, stochOnly = TRUE)
#troubleshooting
# simNodes <- fm$getDependencies(parentNodes, self = FALSE)
# simNodes2 <- fm$getDependencies(parentNodes, self = FALSE, downstream = TRUE) #same as simNodes
# identical(simNodes,simNodes2)
# simNodes3 <- fm$getDependencies(parentNodes, self = TRUE, downstream = TRUE, returnScalarComponents = TRUE)
# View(as.data.frame(simNodes))
# View(as.data.frame(simNodes3))

#Get the mcmc object from the uncompiled model
fm.mcmc <- buildMCMC(fm, monitors = parentNodes) #build an MCMC object
fmc.mcmc <- compileNimble(fm.mcmc, project = fm)



# Rerun model in parallel with only the parent nodes monitored -----------------
nimbrun <- function(seed, ndata, nconst){
  
  require(nimble) 

  mod <- nimbleCode(
    {
      #priors
      a.psi ~ dnorm(0, 0.01)
      tau <- 1/sig^2
      sig ~ dgamma(2, 0.5)
      b1.p ~ dnorm(0, 0.01)
      b2.p ~ dnorm(0, 0.01)
      b3.p ~ dnorm(0, 0.01)
      b1.s ~ dnorm(0, 0.01)
      b2.s ~ dnorm(0, 0.01)
      b3.s ~ dnorm(0, 0.01)
      b0.p ~ dnorm(0, 0.01)
      b0.s ~ dnorm(0, 0.01)
      
      #beta cov priors
      for(n in 1:nenvcov){
        b.pois[n] ~ dnorm(0, 0.1)
        b.size[n] ~ dnorm(0, 0.1)
        
      }
      
      #CAR priors
      for(n in 1:2){
        u[1:ncgrid, n] ~ dcar_proper(muu[1:ncgrid, n], C[1:L], adj[1:L], num[1:ncgrid], M[1:ncgrid], tauu[n], gammau[n])
        tauu[n] <- 1/sdu[n]^2
        sdu[n] ~ dgamma(2, 0.5)
        gammau[n] ~ dunif(bou[1], bou[2])
        
        mu0[n] ~ dnorm(0, 0.1) #use this if we don't want the CAR to be zero-centered
        
        
        for(m in 1:ncgrid){
          muu[m, n] <- mu0[n]
          #muu[m, n] <- 0
        }#n
        
      }#m
      
      u1 <- mean(u[1:ncgrid, 1])
      u2 <- mean(u[1:ncgrid, 2])
      
      #likelihood statement for number of aggregations
      for(k in 1:ngrid){
        for(j in 1:nsurv){
          
          y[k, j] ~ dpois(lambda.star[k, j])
          lambda.star[k, j] <- lambda[k, j] * z[k, j]
          log(lambda[k, j]) <- log(theta[k, j]) + e[k, j]
          z[k, j] ~ dbern(psi[k, j])
          log(theta[k, j]) <- u[gridID[k], 1] + inprod(envcov[k, j, 1:nenvcov], b.pois[1:nenvcov]) # + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
          logit(psi[k, j]) <- a.psi + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
          
          #derive variables
          
          #tp[k, j] <- theta[k, j] * psi[k, j]
          #mu2[k, j] <- u[gridID[k], 2] # + b.size[1]*envcov[k, j, 1] + b.size[7]*envcov[k, j, 7] + b.size[8]*envcov[k, j, 8] + b.size[10]*envcov[k, j, 10] + b.size[11]*envcov[k, j, 11]
          
          #tp.mu[k, j] <- tp[k, j] * mu2[k, j]
          
          #model checking
          
          #ynew[k, j] ~ dpois(lambda.star[k, j])
          
        }#j
      }#k
      
      #size observation model - have to select which model we want to run now: OG Gaussian, Gaussian w/ log data, Lognormal, or Gamma
      for(i in 1:nobs){
        
        #gamma
        
        #lognormal model
        size[i] ~ dlnorm(mu[i], tau) # * z[k, j]
        #non-joint model
        mu[i] <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
        
        
        #Gaussian model
        # size[i] ~ dnorm(mu[i], tau) # * z[k, j]
        # #non-joint model - log tranformed data
        # mu[i] <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
        # 
        # #non-joint model - OG model
        # #log(mu[i]) <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
        # 
        # #joint model
        # #log(mu[i]) <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.pois[1:nenvcov]) # + b.size[1]*envcov[obsID[i], obsurvID[i], 1] + b.size[7]*envcov[obsID[i], obsurvID[i], 7] + b.size[8]*envcov[obsID[i], obsurvID[i], 8] + b.size[10]*envcov[obsID[i], obsurvID[i], 10]  + b.size[11]*envcov[obsID[i], obsurvID[i], 11] #+ inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov]) # + b1.s*sea1[obsurvID[i], obssurvtype[i]] + b2.s*sea2[obsurvID[i], obssurvtype[i]] + b3.s*sea3[obsurvID[i], obssurvtype[i]]
        # 
        
        #use z as an indicator variable?
        
        #model checking
        
        #sizenew[i] ~ dnorm(mu[i], tau)
        
      }#i
      
    })
  
  dat <- ndata
  const <- nconst
  
  inits <- list(u = structure(.Data= runif(const$ncgrid * 2, -.1, .1),.Dim=c(const$ncgrid, 2)),
                sdu = runif(2, 2, 4),
                mu0 = c(runif(1, -1, 0), runif(1, 3, 5)),
                a.psi = runif(1, -1, 1),
                b.pois = runif(const$nenvcov),
                b.size = runif(const$nenvcov))

  #set MCMC run parameters
  ni <- 100000
  na <- 10
  nb <- 50000
  nt <- 10
  nc <- 1
  
  fm <- nimbleModel(code = mod, data = dat, constants = const, inits = inits) #define our model and initialize
  # create vectors of node names 
  dataNodes <- fm$getNodeNames(dataOnly = TRUE)
  parentNodes <- fm$getParents(dataNodes, stochOnly = TRUE)
  simNodes <- fm$getDependencies(parentNodes, self = FALSE)
  
  fmc <- compileNimble(fm) #compile the model
  #Get the mcmc object from the uncompiled model
  fm.mcmc <- buildMCMC(fm, monitors = parentNodes) #build an MCMC object
  fmc.mcmc <- compileNimble(fm.mcmc, project = fm)
  
  samples <- runMCMC(fmc.mcmc, niter = ni, nburnin = nb, nchains = nc, thin = nt, 
                     setSeed = seed,
                     summary = FALSE, samplesAsCodaMCMC = TRUE)
  
  
  return(samples)
  
}

start.time <- Sys.time()
this_cluster <- snow::makeCluster(3)

chain_out <- snow::parLapply(cl = this_cluster, x = 1:3, fun = nimbrun, ndata = dat, nconst = const)
stopCluster(this_cluster)
stop.time <- Sys.time()
stop.time - start.time

save(chain_out, file= paste0(ppc_dir,'/sample_for_pcc_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_finalcovs_logdepth_fprob.Rdata'))

# Load model data --------------------------------------------------------------
#'[Note: model results are saved  with Robject name 'chain_out']
#'[Note: this is exactly the same as our fitted model except we only monitored the parent nodes, to match
#'[      when using the ppc sampler


#import samples output either from fitted model or model rerun to just get parentNodes

#OG non-joint fprob model 1/27/2022 - correct v3 grid, corrected covariates with aligned grid IDs, updated sst_anomaly
#load(file= paste0(ppc_dir,'/sample_for_pcc_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_finalcovs_logdepth_fprob.Rdata'))


#updated model 3/16/2022 with log transformed size data (count model the same)
# load(file= paste0(fixed_grid_model,'/logsizedat_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_finalcovs_logdepth_fprob_pred.Rdata'))
# chain_out <- as.mcmc.list(list(chain_out[[1]]$samples, chain_out[[2]]$samples, chain_out[[3]]$samples))

# # nonjoint negbin count model and lognormal size model
# load(file= paste0(fixed_grid_model,'/negbincount_lgnmlsize_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_fprob.Rdata'))
# #fixing the data for the size model
# load(file= paste0(fixed_grid_model,'/take2_negbincount_lgnmlsize_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_fprob.Rdata'))
# #doubling iterations and burn-in to converge the negbin count model
# load(file= paste0(fixed_grid_model,'/take3_negbincount_lgnmlsize_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_fprob.Rdata'))
# #model gamma-poisson mixture with uniform dist for r prior
# load(file= paste0(fixed_grid_model,'/take2_mixture_negbincount_lgnmlsize_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_fprob.Rdata'))
#final? model gamma-poisson mix long run, tighter priors
#load(file= paste0(fixed_grid_model,'/take5_mixture_negbincount_lgnmlsize_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_fprob.Rdata'))
load(file= paste0(fixed_grid_model,'/take6_mixture_negbincount_Gauslogdtasize_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_fprob.Rdata'))


chain_out <- as.mcmc.list(list(chain_out[[1]]$samples, chain_out[[2]]$samples, chain_out[[3]]$samples))

#get the nodes in the samples
# parentNodes_cols<-colnames(chain_out[[1]]) #saved this from model where only ParentNodes were monitored
View(as.data.frame(colnames(chain_out[[1]])))

#if did not rerun model to just get parentNodes: subset and sort samples for just parentNode columns
#if reran the model for just parent nodes skip to combine chains for ppc

#get indices of the parentNodes and subset chain_out
idx<-dimnames(chain_out[[2]])
idx<-idx[[2]]
idx.b.pois<-grep('^b.pois', idx) 
idx.b.size<-grep('^b.size', idx) 
#idx.r<-grep('^r', idx)
idx.rho<-grep('^rho', idx)
idx.sig<-grep('^sig', idx)
idx.u<-grep('^u', idx)
idx.z<-grep('^z', idx)
idx.parentNodes<-c(idx.b.pois,idx.b.size,idx.rho,idx.sig,idx.u[1:(length(idx.u)-2)],idx.z)
#idx.parentNodes<-c(idx.b.pois,idx.b.size,idx.sig,idx.u[1:(length(idx.u)-2)],idx.z)
chain_out<-chain_out[, idx.parentNodes]

# #thin the posterior to prevent erroring out - thin by 3
# chain_out<-as.mcmc.list(list(as.mcmc(chain_out[[1]][seq(1,nrow(chain_out[[1]]),3),]), 
#                             as.mcmc(chain_out[[2]][seq(1,nrow(chain_out[[2]]),3),]), 
#                             as.mcmc(chain_out[[3]][seq(1,nrow(chain_out[[3]]),3),])))
# 
#thin by 2
chain_out<-as.mcmc.list(list(as.mcmc(chain_out[[1]][seq(1,nrow(chain_out[[1]]),2),]),
                             as.mcmc(chain_out[[2]][seq(1,nrow(chain_out[[2]]),2),]),
                             as.mcmc(chain_out[[3]][seq(1,nrow(chain_out[[3]]),2),])))


#combine chains for ppc - 
#'[Note: samples needs to be a matrix
samples<-as.matrix(mcmc(rbind(chain_out[[1]],chain_out[[2]],chain_out[[3]])))

#check the nodes in the model
#View(as.data.frame(fm$expandNodeNames(fm.mcmc$mvSamples$getVarNames())))

#save(samples, file= paste0(ppc_dir,'/chain_out_data_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(samples, file = paste0(ppc_dir,'/parentNodes_samples_logsizedat_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(samples, file = paste0(ppc_dir, '/parentNodes_samples_negbincount_lnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(samples, file = paste0(ppc_dir, '/parentNodes_samples_take2_mixture_negbincount_lnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(samples, file = paste0(ppc_dir, '/parentNodes_samples_ZIPgausmix_lgnmlsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
# save(samples, file = paste0(ppc_dir, '/parentNodes_samples_take5_mixture_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
# save(samples, file = paste0(ppc_dir, '/parentNodes_thinned_samples_take5_mixture_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
# save(samples, file = paste0(ppc_dir, '/parentNodes_thinnedx2_samples_take5_mixture_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
save(samples, file = paste0(ppc_dir, '/parentNodes_samples_take6_mixture_negbincount_lgdtaGaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
save(samples, file = paste0(ppc_dir, '/parentNodes_thinnedX2_samples_take6_mixture_negbincount_lgdtaGaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))

load(file = paste0(ppc_dir, '/parentNodes_thinnedX2_samples_take6_mixture_negbincount_lgdtaGaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))


#save parentNodes column names from output of model run just with parentNodes monitored 
#save(parentNodes_cols, file = paste0(ppc_dir, '/parentNodes_colnames_modelrun_onlyParentNodesMonitored.Rdata'))

# Set up nimble function to sample the posterior and recreate the data nodes ----

#all querying of the model structure needs to happen in the setup code
#need to pass the mcmc object to the nimble function, to determine the names of 
# the variables we are copying from the posterior samples into the model

ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE, returnScalarComponents = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE, downstream = TRUE, returnScalarComponents = TRUE)
    #vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix - this is not working I think
    postNames <- colnames(samples)
    cat("Using posterior samples of:", paste(postNames, collapse = ','), ".\n") #changed vars to postNames here
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)
    for(i in 1:nSamp) {
      values(model, postNames) <<- samples[i, ] #changed vars to postNames here
      model$simulate(simNodes, includeData = TRUE)
      ppSamples[i, ] <- values(model, dataNodes)
    }
    returnType(double(2))
    return(ppSamples)
  })



# Create the sampler for this model and this MCMC ----
ppSampler <- ppSamplerNF(fm, fm.mcmc) 
cppSampler <- compileNimble(ppSampler, project = fm)


# Calculate the simulated data values ----
set.seed(1)
system.time(ppSamples_via_nf <- cppSampler$run(samples)) #here is where we put in the output from the fitted model

# And more data wrangling . . . 

#set column names as datanodes names
colnames(ppSamples_via_nf)<-dataNodes
ppSamples_via_nf<-as.data.frame(ppSamples_via_nf)

#get indices of the y and size columns
idx.y<-grep('^y', dataNodes)
idx.size<-grep('^size', dataNodes)
ppc.y<-ppSamples_via_nf[, min(idx.y):(max(idx.y))]
ppc.size<-ppSamples_via_nf[, 1:length(idx.size)]

#save ppc output
#save(ppc.y, ppc.size, file= paste0(ppc_dir,'/ppsamples_simulated_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(ppc.size, ppSamples_via_nf, file = paste0(ppc_dir, '/ppsamples_size_logdata_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(ppc.y, ppc.size, file = paste0(ppc_dir, '/ppsamples_negbincount_lgnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(ppc.y, ppc.size, file = paste0(ppc_dir, '/ppsamples_take2_mixture_negbincount_lgnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#save(ppc.y, ppc.size, file = paste0(ppc_dir, '/ppsamples_ZIPgausmix_lgnmlsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))

save(ppc.y, ppc.size, file = paste0(ppc_dir, '/ppsamples_take5_mixture_negbincount_lgnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
save(ppc.y, ppc.size, file = paste0(ppc_dir, '/thinby2_ppsamples_take5_mixture_negbincount_lgnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))


save(ppc.y, ppc.size, file = paste0(ppc_dir, '/thinby2_ppsamples_take6_mixture_negbincount_lggddtasize_updt_grid_ffa_car_nonjoint_finalcovs_fprob.Rdata'))


# Load and organize data for ppc ----

#'[Note: if satisfied that ppSampler data is correct, can delete samples, fm, fmc, fm.mcmc, 
#'[fmc.mcmc, ppSampler, and cppSampler from environment to free up memory]


##load data for ppc ----

# 1.load lambda.star and mu from fitted model
#load(file = paste0(fixed_grid_model,'/updt_grid_ffa_car_nonjoint_fprob_postdata.Rdata'))

# #take5 mixture model
# load(file = paste0(fixed_grid_model,'/outlambdastar2_take5_mix_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
# load(file = paste0(fixed_grid_model,'/outmu_take5_mix_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))

#take6 - mixture abundance model, logged data Gaussian size model
load( file = paste0(fixed_grid_model,'/outlambdastar2_take6_mix_negbincount_lggdgaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
load( file = paste0(fixed_grid_model,'/outmu_take6_mix_negbincount_lggdgaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))


#'[only if posterior wasn't thinned, otherwise we have to recalculate below]
#2. load lambda.star.fitsum and mu.fitsum - if using thinned ppc posterior skip loading fitsum

#load(file = paste0(fixed_grid_model,'/updt_grid_ffa_car_nonjoint_fprob_fit_sum_fullmodel.Rdata'))
# load( file = paste0(fixed_grid_model,'/lambdastar2_fitsum_take5_mix_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
# load( file = paste0(fixed_grid_model,'/mu_fitsum_take5_mix_negbincount_lgnmlsize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))

# #take6 - mixture abundance model, logged data Gaussian size model
# load( file = paste0(fixed_grid_model,'/lambdastar2_fitsum_take6_mix_negbincount_lggdgaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))
# load( file = paste0(fixed_grid_model,'/mu_fitsum_take6_mix_negbincount_lggdgaussize_updt_grid_ffa_car_nonjoint_fprob.Rdata'))


#3. load ppc data - if not in environment

#load(file = paste0(ppc_dir,'/ppsamples_simulated_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
#load(file = paste0(ppc_dir, '/ppsamples_take5_mixture_negbincount_lgnormsize_updt_grid_ffa_car_nonjoint_finalcovs_logdepth_fprob.Rdata'))
load(file = paste0(ppc_dir, '/thinby2_ppsamples_take6_mixture_negbincount_lggddtasize_updt_grid_ffa_car_nonjoint_finalcovs_fprob.Rdata'))


## thin the posterior if we did that for the sampler ----

#subset the samples in outlambda.star2 
#thin the posterior by 3 - same as we thinned the sampler by
# out_lambda.star2<-as.mcmc.list(list(as.mcmc(out_lambda.star2[[1]][seq(1,nrow(out_lambda.star2[[1]]),3),]), 
#                              as.mcmc(out_lambda.star2[[2]][seq(1,nrow(out_lambda.star2[[2]]),3),]), 
#                              as.mcmc(out_lambda.star2[[3]][seq(1,nrow(out_lambda.star2[[3]]),3),])))
# 

# #thin the posterior by 3 - same as we thinned the sampler by
# out_mu<-as.mcmc.list(list(as.mcmc(out_mu[[1]][seq(1,nrow(out_mu[[1]]),3),]), 
#                                     as.mcmc(out_mu[[2]][seq(1,nrow(out_mu[[2]]),3),]), 
#                                     as.mcmc(out_mu[[3]][seq(1,nrow(out_mu[[3]]),3),])))

#thin the posterior by 2 - same as we thinned the sampler by

out_lambda.star2<-as.mcmc.list(list(as.mcmc(out_lambda.star2[[1]][seq(1,nrow(out_lambda.star2[[1]]),2),]),
                                    as.mcmc(out_lambda.star2[[2]][seq(1,nrow(out_lambda.star2[[2]]),2),]),
                                    as.mcmc(out_lambda.star2[[3]][seq(1,nrow(out_lambda.star2[[3]]),2),])))


out_mu<-as.mcmc.list(list(as.mcmc(out_mu[[1]][seq(1,nrow(out_mu[[1]]),2),]), 
                          as.mcmc(out_mu[[2]][seq(1,nrow(out_mu[[2]]),2),]), 
                          as.mcmc(out_mu[[3]][seq(1,nrow(out_mu[[3]]),2),])))

## combine chains into a matrix ----

#lambda.star.matrix
# lambda.star.mat<-as.data.frame(rbind(out_lambda.star[[1]],out_lambda.star[[2]],out_lambda.star[[3]]))
# lambda.star.mat<-as.data.frame(rbind(out_lambda.star1[[1]],out_lambda.star1[[2]],out_lambda.star1[[3]]))

lambda.star.mat<-as.data.frame(rbind(out_lambda.star2[[1]],out_lambda.star2[[2]],out_lambda.star2[[3]]))

mu.mat<-as.data.frame(rbind(out_mu[[1]],out_mu[[2]],out_mu[[3]]))


## actual observations, convert to vector ----
y_act<-c(yy)
#size_act<-obs$Area_sqm
size_act<-obs$log_Area_sqm

#lambda.star.mat and y_act are already the same size and in the same order
#mu.mat and size_act are already the same size and in the same order

#remove out_lambda.star and out_mu to free up memory if necessary

## if using a thinned posterior - redo fitsum ----
lambda.star2.fitsum<-summary(out_lambda.star2)
mu.fitsum<-summary(out_mu)

#Wrangling count data ----

## fitsum to dataframe ----

#lambda.star.fit

# lambdastar.fit<-as.data.frame(lambda.star.fitsum$statistics[,1])
# colnames(lambdastar.fit)<-"lambda.star"

lambdastar.fit<-as.data.frame(lambda.star2.fitsum$statistics[,1])
colnames(lambdastar.fit)<-"lambda.star"

lambdastar.fit<-lambdastar.fit %>%
  tibble::rownames_to_column(var = "nodeName")

#concatenate y_act with lambda.star means from fit.sum, remove cases with y_act = NA
ppc.dat1<-cbind(lambdastar.fit,y_act)
ppc.dat1<-ppc.dat1[complete.cases(ppc.dat1),]
ppc.dat1$dataNode<-str_sub(ppc.dat1$nodeName, start = 12)

#need to do some data gymnastics here to get the lambda.star and ppc.y matrices to match:
# 1) the same size - need to eliminate lambda.star columns where y_act was NA
# 2) the same order - the ppc.y output is ordered topologically which is different from the model output

# names_ppc_y<-as.character(colnames(ppc.y)) #this is the correct # of cols in topological order
# names_ppc_y<-str_sub(names_ppc_y, start = 2) #remove the y from the names
# names_use<-paste0("lambda.star", names_ppc_y) #get colnames in correct order for lambda.star.mat

#for lambda.star2
names_ppc_y<-as.character(colnames(ppc.y)) #this is the correct # of cols in topological order
names_ppc_y<-str_sub(names_ppc_y, start = 2) #remove the y from the names
names_use<-paste0("lambda.star2", names_ppc_y) #get colnames in correct order for lambda.star.mat

#subset AND sort the columns of lambda.star.mat 
lambda.star.mat.2<-lambda.star.mat %>%
  dplyr::select(all_of(names_use))

# #make a vector of y_act  to match ppc.y/lambda.star.mat
# #(can use all of these values and they are already in the correct order)
# y_act.mat<-as.data.frame(t(y_act[complete.cases(y_act)]))
# 
# #add column names to y_act.mat
# y_act_cols<-paste0("y.act", ppc.dat1$dataNode)
# colnames(y_act.mat)<-y_act_cols
# 
# #make column names
# y_names_use<-paste0("y.act", names_ppc_y)
# #subset AND sort the columns of y_act.mat 
# 
# y_act.mat<-y_act.mat %>%
#   dplyr::select(all_of(y_names_use))
# 
# #replicate to match size of ppc_y
# y_act.mat<-bind_rows(y_act.mat, y_act.mat[rep(1, 14999), ])

# Wrangling the size data ----

#Data gymnastics for the size data
# 1) the same size - need to eliminate mu columns where size_act was NA
# 2) the same order - the ppc.size output is ordered topologically which is different from the model output

names_ppc_size<-as.character(colnames(ppc.size)) #this is the correct # of cols in topological order
names_ppc_size<-str_sub(names_ppc_size, start = 5) #remove the 'size' from the names
names_use<-paste0("mu", names_ppc_size) #get colnames in correct order for mu.mat

#subset AND sort the columns of mu.mat 
mu.mat.2<-mu.mat %>%
  dplyr::select(all_of(names_use))

mu.mat.2<-log(mu.mat.2)

#mu.fit
mu.fit<-as.data.frame(mu.fitsum$statistics[,1])
colnames(mu.fit)<-"mu"

mu.fit<-mu.fit %>%
  tibble::rownames_to_column(var = "nodeName")

#concatenate size_act with mu means from fit.sum, remove cases with size_act = NA
ppc.dat2<-cbind(mu.fit,size_act)
ppc.dat2<-ppc.dat2[complete.cases(ppc.dat2),]
ppc.dat2$dataNode<-str_sub(ppc.dat2$nodeName, start = 3)

# Calculate summary stats for the data ----

##ppSampler y values
hist(c(as.matrix(ppc.y)))
max(ppc.y)
min(ppc.y)
median(as.matrix(ppc.y))
mean(as.matrix(ppc.y))
var(c(as.matrix(ppc.y)))
sd(c(as.matrix(ppc.y)))

# actual y observations
hist(y_act)
max(y_act, na.rm = TRUE)
min(y_act, na.rm = TRUE)
median(y_act, na.rm = TRUE)
mean(y_act, na.rm = TRUE)
var(y_act, na.rm = TRUE)
sd(y_act, na.rm = TRUE)

#lambda.star values
hist(c(as.matrix(lambda.star.mat)))
max(lambda.star.mat)
min(lambda.star.mat)
median(as.matrix(lambda.star.mat))
mean(as.matrix(lambda.star.mat))
var(c(as.matrix(lambda.star.mat)))
sd(c(as.matrix(lambda.star.mat)))

#ppc.size values
#'[if using the logged data Gaussian size model, the values are logged]
max(ppc.size)
min(ppc.size)
median(as.matrix(ppc.size))
mean(as.matrix(ppc.size))
hist(as.matrix(ppc.size))
#plot exponentiated values
hist(as.matrix(exp(ppc.size)))

#actual size observations - these are logged values if using take6 model run
max(size_act, na.rm = TRUE)
min(size_act, na.rm = TRUE)
median(size_act, na.rm = TRUE)
mean(size_act, na.rm = TRUE)
var(size_act, na.rm = TRUE)
sd(size_act, na.rm = TRUE)
hist(size_act)
#plot exponentiated values
hist(exp(size_act))

hist(obs$log_Area_sqm)

#mu values - I exponentiated these values 

# max(log(mu.mat))
# min(log(mu.mat))
# median(as.matrix(log(mu.mat)))
# mean(as.matrix(log(mu.mat)))
# var(c(as.matrix(log(mu.mat))))
# sd(c(as.matrix(log(mu.mat))))
# hist(c(as.matrix(log(mu.mat))))


max(mu.mat.2)
min(mu.mat.2)
median(as.matrix(mu.mat.2))
mean(as.matrix(mu.mat.2))
var(c(as.matrix(mu.mat.2)))
sd(c(as.matrix(mu.mat.2)))
hist(c(as.matrix(mu.mat.2)))




#' #plot y_act.mat vs. ppc.y
#' #'[Note: this takes a huge amount of RAM]
#' plot.data<-cbind(c(as.matrix(y_act.mat)),c(as.matrix(ppc.y)))
#' colnames(plot.data)<-c("y.act","ppc.y")
#' plot.data<-as.data.frame(plot.data)
#' head(plot.data)
#' 
#' #let's just plot part of the data first: 1 chain is 1:153205 rows
#' #plot(plot.data[1:153205,1],plot.data[1:153205,2], pch = 20, cex = 1)
#' 
#' 
#' # create a color palette to use in smoothed scatterplot
#' library(RColorBrewer)
#' palette = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF",
#'            "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
#' # myColRamp = colorRampPalette(c(buylrd))
#' 
#' palette <- hcl.colors(30, palette = "viridis")
#' palette <- topo.colors(30)
#' 
#' #or smoothScatter
#' smoothScatter(plot.data$y.act, plot.data$ppc.y, colramp = colorRampPalette(palette), bandwidth = 0.1, nrpoints = 1000, pch = 20, cex = 0.5, col = "white")
#' 

   
# PPC checks --------------------------------------------------------------------

## Chi-square PPC ----

### Count model ----

# Calculate Chi-square test statistics
#'[Note: in Kery, he calculates the chi-squared by square-rooting the denominator lambda.star (expected value)]

##method 1: use vector of lambda.star means
#chi_y<- (ppc.dat1$y_act - ppc.dat1$lambda.star)^2/(ppc.dat1$lambda.star + 0.5)

#method 2: subtract a matrix of y_act values with the matrix of lambda.stars
    #sweep is subtracting the vector from each row in lambda.star.mat 
chi_y2<-(sweep(lambda.star.mat, 2, y_act, '-'))^2/(lambda.star.mat + 0.5)

#spot check that this is doing what we want
# test<-(lambda.star.mat[1,]-y_act)^2/(lambda.star.mat[1,] + 0.5)
# identical(chi_y2[1,], test)

#calculate chi-square of simulated data from ppSampler
chi_ynew<-(ppc.y - lambda.star.mat.2)^2/(lambda.star.mat.2 + 0.5)

#spotcheck
# test<-(ppc.y[3,] - lambda.star.mat.2[3,])^2/(lambda.star.mat.2[3,] + 0.5)
# identical(chi_ynew[3,], test)
# test.vect<-c(3,4,5,6,7)
# test.vect2<-c(1,2,3,4,5)
# test.vect-test.vect2

#calculate the omnibus test statistic for each - sum discrepancies across columns (a sum for each iteration of the model)
fit.chi_y<-rowSums(chi_y2, na.rm = TRUE)
fit.chi_ynew<-rowSums(chi_ynew, na.rm = TRUE)

#Summarize the fit statistics
summary(fit.chi_y)
summary(fit.chi_ynew)

#calculate the ratio
mean(fit.chi_y)/mean(fit.chi_ynew)

fname<-paste0(draft_fig_dir, "/ffa_ppc_thinnedpost_negbin_count_take6.tif")
tiff(filename = fname, width = 5.5, height = 5, units = "in", compression = "lzw", 
     res = 300)

#plot the posteriors of the chi^2 fit statistics
#par(mfrow = c(1,2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(fit.chi_y, fit.chi_ynew, xlim = c(0,1000), ylim = c(0,1000),
     main = "", xlab = "Discrepancy observed data", ylab = "Discrepancy expected data", 
     frame.plot = F, cex = 1, pch = 20, col = alpha("blue", 0.4))
abline(0,1, lwd = 2, lty = "dotted", col = "red")
text(600,200, "bpv = 0.31", cex = 0.8, pos = 4, col = "black")


dev.off()

#calculate the Bayesian p-value
(bpv<-mean(fit.chi_ynew >fit.chi_y))

#we could also do some graphical checks, or look at some other test statistics besides chi^2 that may
#give more info on how our model is doing badly: like min, max, sd, mean, etc.

#can look at the bayesian residuals p158 to see how much error we are seeing

### Size model ----

# Calculate Chi-square test statistics
#'[Note: in Kery, he calculates the chi-squared by square-rooting the denominator lambda.star (expected value)]

##method 1: use vector of mu means
#chi_size<- (ppc.dat1$size_act - ppc.dat1$mu)^2/(ppc.dat1$mu + 0.5)

#method 2: subtract a matrix of size_act values with the matrix of mus
#sweep is subtracting the vector from each row in lambda.star.mat 
#for logged data gaussian model - I need to re-log the mu values
chi_size2<-(sweep(log(mu.mat), 2, size_act, '-'))^2/(log(mu.mat) + 0.5)

# #logged - for when we were modeling the non-logged data
# chi_size2_log<-(sweep(log(mu.mat), 2, log(size_act), '-'))^2/(log(mu.mat) + 0.5)

#calculate chi-square of simulated data from ppSampler
chi_sizenew<-(ppc.size - mu.mat.2)^2/(mu.mat.2 + 0.5)

# #logged
# chi_sizenew_log<-(log(ppc.size) - log(mu.mat.2))^2/(log(mu.mat.2) + 0.5)

#calculate the omnibus test statistic for each - sum discrepancies across columns (a sum for each iteration of the model)
fit.chi_size<-rowSums(chi_size2, na.rm = TRUE)
fit.chi_sizenew<-rowSums(chi_sizenew, na.rm = TRUE)

#logged
fit.chi_size<-rowSums(chi_size2_log, na.rm = TRUE)
fit.chi_sizenew<-rowSums(chi_sizenew_log, na.rm = TRUE)

#Summarize the fit statistics
summary(fit.chi_size)
summary(fit.chi_sizenew)


#calculate the ratio
mean(fit.chi_size)/mean(fit.chi_sizenew)

#plot the posteriors of the chi^2 fit statistics
#par(mfrow = c(1,2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(fit.chi_size, fit.chi_sizenew, xlim = c(1000000,51569600 ), ylim = c(1000000,90800182 ),
     main = "", xlab = "Discrepancy observed data", ylab = "Discrepancy expected data", 
     frame.plot = F, cex = 1, pch = 20, col = alpha("blue", 0.4))
abline(0,1, lwd = 2, lty = "dotted", col = "red")


#logged data plot - gaussian model

fname<-paste0(draft_fig_dir, "/ffa_ppc_thinnedpost_logged_dta_gaus_size.tif")
tiff(filename = fname, width = 5.5, height = 5, units = "in", compression = "lzw", 
     res = 300)

plot(fit.chi_size, fit.chi_sizenew, xlim = c(-36000,21000 ), ylim = c(-36000,21000 ),
     main = "", xlab = "Discrepancy observed data", ylab = "Discrepancy expected data", 
     frame.plot = F, cex = 1, pch = 20, col = alpha("blue", 0.4))
abline(0,1, lwd = 2, lty = "dotted", col = "red")
text(0,-20000, "bpv = 0.35", cex = 0.8, pos = 4, col = "black")

dev.off()

#calculate the Bayesian p-value
(bpv<-mean(fit.chi_sizenew >fit.chi_size))



## Other PPC discrepancy measures ----
#'[Note: These are from the speed of light example in the Gelman book, where we are choosing a test stat
#'[      and comparing it to posterior predicted distribution; I think this is the same as the parametric
#'[      bootstrap example in Kery book p196, he says it does not incorporate the parameter estimation uncertainty]]]

# Minimum observations - example from Gelman et al. (3rd edition), Figure 6.3 and nimble

obsMin<-min(y_act, na.rm = TRUE)
ppMin<-apply(ppc.y, 1, min, na.rm = TRUE)
summary(ppMin)

hist(ppMin, xlim = c(-10, 10),
     main = "Discrepancy = min(y)",
     xlab = "min(y_rep)")
abline(v = obsMin, col = 'red')

obsMin<-min(size_act, na.rm = TRUE)
ppMin<-apply(ppc.size, 1, min, na.rm = TRUE)
summary(ppMin)

hist(ppMin, xlim = c(-5, 0),
     main = "Discrepancy = min(size)",
     xlab = "min(size_rep)")
abline(v = obsMin, col = 'red')

#code from Gelman p 596 to calculate the p-value
# p-value definition (Gelman) - the proportion of the simulations where 
# the test quantity meets or exceeds its realized value
cat("\nPr (T(ppc.y) > T(y) =", round(mean(ppMin>obsMin),2), "\n")

#an extreme p-value (near 0 or 1) implies that the model cannot be expected to capture
#this aspect of the data

# Max observations 

obsMax<-max(y_act, na.rm = TRUE)
ppMax<-apply(ppc.y, 1, max, na.rm = TRUE)
summary(ppMax)

hist(ppMax, xlim = c(700, 1300),
     main = "Discrepancy = max(y)",
     xlab = "max(y_rep)")
abline(v = obsMax, col = 'red')

obsMax<-max(size_act, na.rm = TRUE)
ppMax<-apply(ppc.size, 1, max, na.rm = TRUE)
summary(ppMax)

hist(ppMax, xlim = c(4000, 713590),
     main = "Discrepancy = max(size)",
     xlab = "max(size_rep)")
abline(v = obsMax, col = 'red')

# Median 
obsMed<-median(y_act, na.rm = TRUE)
ppMed<-apply(ppc.y, 1, median, na.rm = TRUE)
obsMed
summary(ppMed)

#no point in plotting - all zeros
# hist(ppMed, xlim = c(700, 1300),
#      main = "Discrepancy = median(y)",
#      xlab = "median(y_rep)")
# abline(v = obsMed, col = 'red')

obsMed<-median(size_act, na.rm = TRUE)
ppMed<-apply(ppc.size, 1, median, na.rm = TRUE)
obsMed
summary(ppMed)

hist(ppMed, xlim = c(25, 33),
     main = "Discrepancy = median(size)",
     xlab = "median(size_rep)")
abline(v = obsMed, col = 'red')

#Mean
obsMean<-mean(y_act, na.rm = TRUE)
ppMean<-apply(ppc.y, 1, mean, na.rm = TRUE)
obsMean
summary(ppMean)

hist(ppMean, xlim = c(0.65, 0.8),
     main = "Discrepancy = mean(y)",
     xlab = "mean(y_rep)")
abline(v = obsMean, col = 'red')

obsMean<-mean(size_act, na.rm = TRUE)
ppMean<-apply(ppc.size, 1, mean, na.rm = TRUE)
obsMean
summary(ppMean)

hist(ppMean, xlim = c(80, 135),
     main = "Discrepancy = mean(size)",
     xlab = "mean(size_rep)",
     breaks = 25)
abline(v = obsMean, col = 'red')

#Variance
obsVar<-var(y_act, na.rm = TRUE)
ppVar<-apply(ppc.y, 1, var, na.rm = TRUE)
obsVar
summary(ppVar)

hist(ppVar, xlim = c(100, 250),
     main = "Discrepancy = var(y)",
     xlab = "var(y_rep)")
abline(v = obsVar, col = 'red')

obsVar<-var(size_act, na.rm = TRUE)
ppVar<-apply(ppc.size, 1, var, na.rm = TRUE)
obsVar
summary(ppVar)

hist(ppVar, xlim = c(2, 2.5),
     main = "Discrepancy = var(size)",
     xlab = "var(size_rep)")
abline(v = obsVar, col = 'red')

# standard deviation
obsSd<-sd(y_act, na.rm = TRUE)
ppSd<-apply(ppc.y, 1, sd, na.rm = TRUE)
obsSd
summary(ppSd)

hist(ppSd, xlim = c(11, 16),
     main = "Discrepancy = sd(y)",
     xlab = "sd(y_rep)")
abline(v = obsSd, col = 'red')

obsSd<-sd(size_act, na.rm = TRUE)
ppSd<-apply(ppc.size, 1, sd, na.rm = TRUE)
obsSd
summary(ppSd)

hist(ppSd, xlim = c(1.5, 1.6),
     main = "Discrepancy = sd(size)",
     xlab = "sd(size_rep)")
abline(v = obsSd, col = 'red')

#Chi-square 
obsChi2<-sum((ppc.dat1$y_act-ppc.dat1$lambdastar)^2/(ppc.dat1$lambdastar + 0.5))
ppChi2<-rowSums((ppc.y - lambda.star.mat.2)^2/(lambda.star.mat.2 + 0.5))
obsChi2
summary(ppChi2)

hist(ppChi2, xlim = c(850, 112200),
     main = "Discrepancy = chi2(y)",
     xlab = "chi2(y_rep)")
abline(v = obsChi2, col = 'red')

obsChi2<-sum((ppc.dat2$size_act-ppc.dat2$mu)^2/(ppc.dat2$mu + 0.5))
ppChi2<-rowSums((ppc.size - mu.mat.2)^2/(mu.mat.2 + 0.5))
obsChi2
summary(ppChi2)

hist(ppChi2, xlim = c(9000000, 150000000), breaks = 75,
     main = "Discrepancy = chi2(size)",
     xlab = "chi2(size_rep)")
abline(v = obsChi2, col = 'red')

#range
obsRange<-diff(range(y_act, na.rm = TRUE))
ppRange<-apply(ppc.y, 1, range, na.rm = TRUE)
ppRange<-apply(t(ppRange), 1, diff, na.rm = TRUE)
obsRange
summary(ppRange)

obsRange<-diff(range(size_act, na.rm = TRUE))
ppRange<-apply(ppc.size, 1, range, na.rm = TRUE)
ppRange<-apply(t(ppRange), 1, diff, na.rm = TRUE)
obsRange
summary(ppRange)

hist(ppRange, xlim = c(4000, 720000), 
     main = "Discrepancy = range(size)",
     xlab = "range(size_rep)")
abline(v = obsRange, col = 'red')

##Graphical posterior checke ----

hist(y_act)

#Take a random sample of the posterior predictions
rand_ppc.y<-ppc.y %>%
  do(sample_n(.,30)) 

par(mfrow=c(5,6))
for (i in 1:30) {
  hist(as.matrix(rand_ppc.y[i,]),
       xlab = paste0("sim.ydata.",i),
       main = "")
}
  
summary(t(rand_ppc.y))

hist(as.matrix(rand_ppc.y[i, which(rand_ppc.y<20)]),
     xlab = paste0("sim.ydata.",i),
     main = "")

# Checking new data values generated from rerunning model ------------------------

#'[Reran fprob nonjoint model code with ynew and sizenew turned on and monitored]

#load data
load(file= paste0(ppc_dir,'/data_sim_updt_grid_ffa_car_spatiotemp_envcov_32grid_seapsi_nonjoint_finalcovs_logdepth_fprob_pred.Rdata'))

#convert to mcmc.list
out_mcmc <- as.mcmc.list(list(chain_out[[1]], chain_out[[2]], chain_out[[3]]))

#or

out_mcmc <- as.mcmc.list(list(chain_out[[1]]$samples, chain_out[[2]]$samples, chain_out[[3]]$samples))


#let's check out the ynew and sizenew output generated from rerunning the model

#get indices of the ynew and sizenew
idx<-dimnames(out_mcmc[[2]])
idx<-idx[[2]]
idx.ynew<-grep('^ynew', idx) 
idx.sizenew<-grep('^sizenew', idx) #

#Pull out the new data generated from model

idx.ynew<-out_mcmc[, min(idx.ynew):(max(idx.ynew))]
out_sizenew<-out_mcmc[, min(idx.sizenew):max((idx.sizenew))]

save(out_ynew,out_sizenew, file = paste0(ppc_dir,'/simulated_dataNodes_updt_grid_ffa_car_nonjoint_fprob_postdata.Rdata'))

load(file = paste0(ppc_dir,'/simulated_dataNodes_updt_grid_ffa_car_nonjoint_fprob_postdata.Rdata'))
#ynew
out_ynew<-as.matrix(mcmc(rbind(out_ynew[[1]],out_ynew[[2]],out_ynew[[3]])))
max(out_ynew)
min(out_ynew)
median(out_ynew)
mean(out_ynew)
hist(out_ynew)

#sizenew
out_sizenew<-as.matrix(mcmc(rbind(out_sizenew[[1]],out_sizenew[[2]],out_sizenew[[3]])))
max(out_sizenew)
min(out_sizenew)
median(out_sizenew)
mean(out_sizenew)
hist(out_sizenew)



