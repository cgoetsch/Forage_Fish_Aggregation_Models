## An example specification of the FFA models using the Nimble package in R ----

### Forage fish aggregation models (abundance and size) 

### As described in Goetsch et al.:
### Surface and subsurface features drive forage fish distributions and aggregations in the US Northeast Shelf ecosystem

### packages
library(spdep)
library(nimble)  

### Create CAR object from spatial grid

## cgrid = 32 km grid overlaying the FFA survey area

## create adjacency table
adj <- poly2nb(cgrid)
adjm <- nb2mat(adj, style = 'B')
adj.sparse <- unlist(adj)
num <- unlist(lapply(adj, length))

## create tools to parameterize the CAR in Nimble
CM2<-as.carCM(adj.sparse, rep(1, length(adj.sparse)), num)
C<- CM2$C
M <- CM2$M
bou<-carBounds(C, adj.sparse, num, M)
mumod<-rep(0,length(num))

### Data required for the models:

## ngrid = number of 4 km grid cells with effort in FFA study area
## ncgrid = number of 32 km grid cells in cgrid
## nobs = number of aggregations observed across all surveys
## nsurvtype = number of survey projects (i.e., New York project + Mid-Atlantic project)
## nsurv = number of survey transects (use number from the survey project with the most transects)
## nenvcov = number of environmental covariates 
## L = length of adj.sparse
## wei = vector of ones the length of L
## adj = adj.sparse, the unlisted adjacency table
## gridID = vector the length of ngrid, linking each 4 km grid cell to its corresponding 32 km grid cell
## obsID = vector the length of nobs, linking each aggregation to its corresponding 4 km grid cell
## obsgridID = vector linking each aggregation to its corresponding 32 km grid cell
## obssurvID = vector linking each aggregation to its survey transect number
## obssurvtype = vector linking each aggregation to its survey project, for season association
## surveyID = vector linking each 4 km grid cell to its survey project
## e = matrix of effort (km2) per 4 km grid cell for each transect
## seaID = a matrix linking season to each transect by survey project
## sea1, sea2, sea3 = seasonal fixed effect variables; logical vectors specifying whether a transect is in a given season
## y = matrix of the number of aggregations per 4 km grid cell for each transect
## size = vector of logged size for each aggregation observation
## envcov = environmental covariate array (ncgrid X survey number X nenvcov)

### Model code

model <- nimbleCode(
  {
    # ecological process priors
    a.psi ~ dnorm(0, 0.8) 
    tau <- 1/sig^2
    sig ~ dgamma(2, 0.5)
    b1.p ~ dnorm(0, 0.5)
    b2.p ~ dnorm(0, 0.5)
    b3.p ~ dnorm(0, 0.5)
    b0.p ~ dnorm(0, 0.5)
    b0.s ~ dnorm(0, 0.01)
    r ~ dunif(0,50) 
    
    # beta covariate priors 
    for(n in 1:nenvcov){
      b.pois[n] ~ dnorm(0, 0.5)
      b.size[n] ~ dnorm(0, 0.4)
    }
    
    # CAR priors
    for(n in 1:2){
      u[1:ncgrid, n] ~ dcar_proper(muu[1:ncgrid, n], C[1:L], adj[1:L], num[1:ncgrid], M[1:ncgrid], tauu[n], gammau[n])
      tauu[n] <- 1/sdu[n]^2
      sdu[n] ~ dgamma(2, 0.5)
      gammau[n] ~ dunif(bou[1], bou[2])
      mu0[n] ~ dnorm(0, 0.1)
      
      for(m in 1:ncgrid){
        muu[m, n] <- mu0[n]
      }#n
    }#m
    
    u1 <- mean(u[1:ncgrid, 1])
    u2 <- mean(u[1:ncgrid, 2])
    
    ## ecological process models
    
    # likelihood statement for number of aggregations
    for(k in 1:ngrid){
      for(j in 1:nsurv){
        #Zero-inflated Negative Binomial (ZINB): Poisson-Gamma mixture model
        y[k, j] ~ dpois(lambda.star2[k, j])
        lambda.star2[k, j] <- lambda.star1[k, j] * rho[k, j]
        lambda.star1[k,j]<- lambda[k, j] * z[k, j]
        log(lambda[k, j]) <- log(theta[k, j]) + e[k, j]
        z[k, j] ~ dbern(psi[k, j])
        log(theta[k, j]) <- u[gridID[k], 1] + inprod(envcov[k, j, 1:nenvcov], b.pois[1:nenvcov]) 
        logit(psi[k, j]) <- a.psi + b1.p*sea1[j, surveyID[k]] + b2.p*sea2[j, surveyID[k]] + b3.p*sea3[j, surveyID[k]]
        rho[k, j] ~ dgamma(r, r)
      }#j
    }#k
    
    # likelihood statement for size of aggregations
    for(i in 1:nobs){
      #Gaussian process
      size[i] ~ dnorm(mu[i], tau) 
      mu[i] <- u[obsgridID[i], 2] + inprod(envcov[obsID[i], obsurvID[i], 1:nenvcov], b.size[1:nenvcov])
    }#i
    
  })
