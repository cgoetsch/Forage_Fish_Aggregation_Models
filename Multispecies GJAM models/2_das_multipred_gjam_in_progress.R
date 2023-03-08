# Multispecies Predator Models: GJAM Model -------------------------------------

# Workspace Set-up -------------------------------------------------------------

#install.packages("gjam")

library(gjam)
library(knitr)
library(ggpubr)
library(ggplot2)
library(cowplot)
#library(sf)
#library(rgdal)
library(tidyverse)
library(reshape2)
#library(spdep)
library(tibble)
library(coda)
library(tictoc)
library(matrixStats)


#Directories
input_dir<-"P:/NYSERDA Ecosystem Dynamics/Species Interactions Modeling/Aerial Survey Multispecies Predator Analysis/data"
# test_output_dir<-"P:/NYSERDA Ecosystem Dynamics/Species Interactions Modeling/Aerial Survey Multispecies Predator Analysis/output/gjam models/test models"
fit_mods_dir<-"P:/NYSERDA Ecosystem Dynamics/Species Interactions Modeling/Aerial Survey Multispecies Predator Analysis/output/gjam models/fitted models"


# Load multipredator and prey observation data ----------------------------------

#'[Choose scale "16km" or "32km"
sp_scale <- "32km"

#FFA with no 0 size data 
#'[1/12/2023 I think this is our final data set
load(file = switch(sp_scale, "16km" = paste0(input_dir, "/ms_predprey_alldata_pts_covs_16kmgrid_no0sizeFFA.Rdata"),
     "32km" = paste0(input_dir, "/ms_predprey_alldata_pts_covs_32kmgrid_no0sizeFFA.Rdata"))
     )

# #4km data
# # load(file = paste0(input_dir, "/ms_predprey_alldata_pts.Rdata")) #this is not the no 0 size data
# # 
# # View(ms.y.pts.df)
# 
# #16km data
# #all FFA
# # load(file = paste0(input_dir, "/ms_predprey_alldata_pts_covs_16kmgrid.Rdata"))
# 
# #FFA with no 0 size data
# load(file = paste0(input_dir, "/ms_predprey_alldata_pts_covs_16kmgrid_no0sizeFFA.Rdata"))
# 
# #32km data
# #all FFA
# #load(file = paste0(input_dir, "/ms_predprey_alldata_pts_covs_32kmgrid.Rdata"))
# 
# #FFA with no 0 size data
# load(file = paste0(input_dir, "/ms_predprey_alldata_pts_covs_32kmgrid_no0sizeFFA.Rdata"))
# 

#'[Note: See code in 1_das_multispecies_data_prep_step2.R]
#'


# Create base for plots --------------------------------------------------------

base <- ggplot() +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))



# Full Model -------------------------------------------------------------------

## data wrangling --------

#'[Notes from example: https://rstudio-pubs-static.s3.amazonaws.com/790824_b60052ce8c38436895f64120b2effb4a.html#the-data]
# ydata = dataframe of observed counts, each row is a site (grid cell)-year(survey.visit), 
#         each column is a species - I don't think this allows NAs
# xdata = dataframe of env. data at each site-year
# edata = vector containing effort data for each site-year


ms.data<-ms.y.pts.df[complete.cases(ms.y.pts.df),] # Remove incomplete cases


#shorten names for some categories
ms.data<- ms.data %>% 
  dplyr::rename(Pisc.Fish = 'Piscivorous Fishes',
                Shear.Pets = Shearwaters.Petrels,
                FF = "Forage Fish", NOGA = 'Northern Gannet')

#observation data

#option 1: FF as abundance
# y.ms.dat<-ms.data %>% 
#   dplyr::select(c(Delphinids, Pisc.Fish, Alcids, Shear.Pets, 
#                   NOGA, Gulls.Lge, Gulls.M.Sm, Terns, Loons, FF))

#Balaenids, 

#option 2: FF as cum size
y.ms.dat.2<-ms.data %>% 
  dplyr::select(c(Delphinids, Pisc.Fish, Alcids, Shear.Pets, 
                  NOGA, Gulls.Lge, Gulls.M.Sm, Terns, Loons, FFA.cum.size))

#Gulls.all,

#option 3: switch to PA data
# y.ms.dat.3<- ifelse(as.matrix(y.ms.dat > 0), 1, as.matrix(y.ms.dat)) 
# #reduced
# y.ms.dat.3.red<- as.data.frame(y.ms.dat.3) %>% 
#   dplyr::select(c(Pisc.Fish, Alcids, 
#                   NOGA, Gulls.all, Terns, Loons))


# #option 4: reduced set - 6 species
# y.ms.dat.red<-ms.data %>% 
#   dplyr::select(c(Pisc.Fish, Alcids, 
#                   NOGA, Gulls.all, Terns, Loons))
# 


#covariate data
#with year and season.2 (this is recoding the NY survey that was "winter" but
# actually conducted in a spring month)
# x.ms.dat<-ms.data %>% 
#   dplyr::select(c(space.grid.id, year, season.2, log_depth, rugosity, bpi, sediment, 
#                   sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
#   dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
#                 chl.fprob = chl_fprob)

#without year
x.ms.dat<-ms.data %>%
  dplyr::select(c(space.grid.id, season.2, log_depth, rugosity, bpi, sediment,
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>%
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob)


x.ms.dat$season<-as.factor(x.ms.dat$season)
x.ms.dat$space.grid.id<-as.factor(x.ms.dat$space.grid.id)
x.ms.dat$season.2<-as.factor(x.ms.dat$season.2)
# x.ms.dat$year<-as.factor(x.ms.dat$year)

#with only season
# x.ms.dat<-ms.data %>% 
#   dplyr::select(c(year, season.2)) 

#effort data
e.ms.data<-ms.data$effort.km2

#format effort data
# elist <- list(columns = 1:ncol(y.ms.dat.2),
              # values = e.ms.data)

# convert observations to SPUE -sightings per unit effort - continuous response

# y.ms.dat.cpue <- sqrt(y.ms.dat/e.ms.data)

#with ffa.cum.size per unit effort
y.ms.dat.spue <- sqrt(y.ms.dat.2/e.ms.data)

# #and reduced species to test
# y.ms.dat.cpue.red<- y.ms.dat.cpue%>% 
#   dplyr::select(c(Pisc.Fish, Alcids, 
#                   NOGA, Gulls.all, Terns, Loons))


## Data distribution ---------

# base +
#   geom_histogram(aes(x = y.ms.dat$Delphinids)) +
#   labs(title = "Distribution of Delphinid data",
#        y = "Frequency", x = "Counts of Delphinids")
# 
# base +
#   geom_histogram(aes(x = y.ms.dat.2$FFA.cum.size)) +
#   labs(title = "Distribution of FFA size data",
#        y = "Frequency", x = "Cum Size of FFA")
# 
# base +
#   geom_histogram(aes(x = e.ms.data)) +
#   labs(title = "Distribution of effort toy data",
#        y = "Frequency", x = "Effort (km2)")

## Model specification ----------

#formulas

# formula.full<-as.formula(~ season.2 + log.depth + bpi + sediment + 
#                         sst + chl + ssha + fsle + salinity + 
#                         mld + sst.fprob + chl.fprob)

#interactions with season
formula.full.inter<-as.formula(~ season.2 * (log.depth + bpi + sediment + 
                           sst + chl + ssha + fsle + salinity + 
                           mld + sst.fprob + chl.fprob))

#model information for MCMC resampling and model specs
#ng = number of iterations of Gibbs sampler, burnin = discard; typeNames - datatypes

# mlist <- list(ng = 30000, burnin = 20000, typeNames = "DA", effort = elist)


#spue response data
mlist <- list(ng = 4000, burnin = 1000, 
              typeNames = "CA",
              #effort = elist, 
              random = "space.grid.id",
              PREDICTX = F #speed up computation
)


## Run the model -----------

#32km models - spue and ffa abundance * size term
#Run the model - no interactions
# tic()
# ms.mod <- gjam(formula.full, xdata = x.ms.dat, ydata = y.ms.dat.spue, modelList = mlist)
# toc()
# 
# save(ms.mod, file = paste0(fit_mods_dir, "/ms_gjam_fullmod1_no0sizeFFA_32km_10spp_SPUE_allcov_randeff.Rdata"))
# 


#Run season interactions model - 
tic()
ms.mod <- gjam(formula.full.inter, xdata = x.ms.dat, ydata = y.ms.dat.spue, modelList = mlist)
toc()


#save model inputs for cross-validation

save(ms.mod, 
     file = switch(sp_scale, 
                   "16km" = paste0(fit_mods_dir, "/ms_gjam_interaction_fullmod_no0sizeFFA_16km_10spp_SPUE_allcov_randeff_4000iter.Rdata"),
                   "32km" = paste0(fit_mods_dir, "/ms_gjam_interaction_fullmod_no0sizeFFA_32km_10spp_SPUE_allcov_randeff_4000iter.Rdata"))
     )


save(rand_levels, season_levels, 
     file = switch(sp_scale, 
                   "16km" = paste0(fit_mods_dir, "/factor_levs_ms_gjam_interaction_fullmod_no0sizeFFA_16km_10spp_SPUE_allcov_randeff_4000iter.Rdata"),
                   "32km" = paste0(fit_mods_dir, "/factor_levs_ms_gjam_interaction_fullmod_no0sizeFFA_32km_10spp_SPUE_allcov_randeff_4000iter.Rdata"))
     )


## Initial output and chain stabilization -----------

gjamPlot(ms.mod, plotPars = list(PLOTALLY = T, GRIDPLOTS = T, CLUSTERPLOTS = T))


## Evaluate convergence  -------------------------------------------------

#Put gjam chains in mcmc list format

# PSRF/R hat/Gelman diagnostic

#'[Note: for gelman.diag function run the model 3 times and compare
#'[      the chains from each model run

# gelman.diag(x, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
#             multivariate=TRUE)
# 
# #Effective samples size
# lapply(x,effectiveSize)


# Seasonal Models --------------------------------------------------------------

## data wrangling ---------------

ms.data<-ms.y.pts.df[complete.cases(ms.y.pts.df),] # Remove incomplete cases

#shorten names for some categories
ms.data<- ms.data %>% 
  dplyr::rename(Pisc.Fish = 'Piscivorous Fishes',
                Shear.Pets = Shearwaters.Petrels,
                FF = "Forage Fish", NOGA = 'Northern Gannet')

#by season.2 (reclassifying one of Nyserda survey transects from winter to Spring)
fall.ms.data<-ms.data %>% 
  dplyr::filter(season.2 == "Fall")

win.ms.data<-ms.data %>% 
  dplyr::filter(season.2 == "Winter")

spr.ms.data<-ms.data %>% 
  dplyr::filter(season.2 == "Spring")

sum.ms.data<-ms.data %>% 
  dplyr::filter(season.2 == "Summer")

#observation data
fall.y.dat<-fall.ms.data %>% 
  dplyr::select(c(Delphinids, Pisc.Fish, Alcids, Shear.Pets, 
                  NOGA, Gulls.Lge, Gulls.M.Sm, Terns, Loons, FFA.cum.size))

win.y.dat<-win.ms.data %>% 
  dplyr::select(c(Delphinids, Pisc.Fish, Alcids, Shear.Pets, 
                  NOGA, Gulls.Lge, Gulls.M.Sm, Terns, Loons, FFA.cum.size))

spr.y.dat<-spr.ms.data %>% 
  dplyr::select(c(Delphinids, Pisc.Fish, Alcids, Shear.Pets, 
                  NOGA, Gulls.Lge, Gulls.M.Sm, Terns, Loons, FFA.cum.size))

sum.y.dat<-sum.ms.data %>% 
  dplyr::select(c(Delphinids, Pisc.Fish, Alcids, Shear.Pets, 
                  NOGA, Gulls.Lge, Gulls.M.Sm, Terns, Loons, FFA.cum.size))

#check which species have observations in each season
# colSums(fall.y.dat)
# colSums(win.y.dat)
# colSums(spr.y.dat)
# colSums(sum.y.dat)

#Cull the species from each season without enough observations

#Winter - remove Shear.Pets, Terns, and FFA
win.y.dat<-win.y.dat %>% 
  dplyr::select(-c(Shear.Pets, Terns, FFA.cum.size))

#Summer - remove Alcids, Noga, Loons
sum.y.dat<-sum.y.dat %>% 
  dplyr::select(-c(Alcids, NOGA, Loons))

#Covariate data
#'[Note: 1/23/2023 taking out year,season.2 and rugosity so I don't have to make dummy 
#'[      variables for predictions; also tested chl_fmean and sst_fmean on fall but
#'[      decided to remove because of high pairwise correlation and vif]]]
#'[Note: 2/16/2023 - updated model replaces the 64km2 (space.grid.id) RE 
#'[                  with the 128 km2 RE (space.grid.id.2)]

fall.x.dat<-fall.ms.data %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob )#,
                 # sst_fmean, chl_fmean)
                ) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob) #, sst.fmean = sst_fmean,
                # chl.fmean = chl_fmean)


# fall.x.dat$season<-as.factor(fall.x.dat$season)
# fall.x.dat$space.grid.id<-as.factor(fall.x.dat$space.grid.id)
# fall.x.dat$season.2<-as.factor(fall.x.dat$season.2)
# fall.x.dat$year<-as.factor(fall.x.dat$year)
fall.x.dat$space.grid.id.2<-as.factor(fall.x.dat$space.grid.id.2)

win.x.dat<-win.ms.data %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob)


# win.x.dat$season<-as.factor(win.x.dat$season)
# win.x.dat$space.grid.id<-as.factor(win.x.dat$space.grid.id)
# win.x.dat$season.2<-as.factor(win.x.dat$season.2)
# win.x.dat$year<-as.factor(win.x.dat$year)
win.x.dat$space.grid.id.2<-as.factor(win.x.dat$space.grid.id.2) 


spr.x.dat<-spr.ms.data %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob)

# spr.x.dat$season<-as.factor(spr.x.dat$season)
# spr.x.dat$space.grid.id<-as.factor(spr.x.dat$space.grid.id)
# spr.x.dat$season.2<-as.factor(spr.x.dat$season.2)
# spr.x.dat$year<-as.factor(spr.x.dat$year)
spr.x.dat$space.grid.id.2<-as.factor(spr.x.dat$space.grid.id.2)

sum.x.dat<-sum.ms.data %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob)

# sum.x.dat$season<-as.factor(sum.x.dat$season)
# sum.x.dat$space.grid.id<-as.factor(sum.x.dat$space.grid.id)
# sum.x.dat$season.2<-as.factor(sum.x.dat$season.2)
# sum.x.dat$year<-as.factor(sum.x.dat$year)
sum.x.dat$space.grid.id.2<-as.factor(sum.x.dat$space.grid.id.2)

#effort data by season
fall.e.data<-fall.ms.data$effort.km2
win.e.data<-win.ms.data$effort.km2
spr.e.data<-spr.ms.data$effort.km2
sum.e.data<-sum.ms.data$effort.km2

#format effort data- only if running abundance and NOT spue
# fall.elist <- list(columns = 1:ncol(fall.y.dat),
#               values = fall.e.data)
# win.elist <- list(columns = 1:ncol(win.y.dat),
#                    values = win.e.data)
# spr.elist <- list(columns = 1:ncol(spr.y.dat),
#                    values = spr.e.data)
# sum.elist <- list(columns = 1:ncol(sum.y.dat),
#                    values = sum.e.data)

# convert observations to SPUE 

#with ffa.cum.size per unit effort
fall.y.dat.spue <- sqrt(fall.y.dat/fall.e.data)
win.y.dat.spue <- sqrt(win.y.dat/win.e.data)
spr.y.dat.spue <- sqrt(spr.y.dat/spr.e.data)
sum.y.dat.spue <- sqrt(sum.y.dat/sum.e.data)

#look at the number of the random levels
# fall.x.dat %>% 
#   dplyr::count(space.grid.id.2)
# 
# win.x.dat %>% 
#   dplyr::count(space.grid.id.2)
# 
# spr.x.dat %>% 
#   dplyr::count(space.grid.id.2)
# 
# sum.x.dat %>% 
#   dplyr::count(space.grid.id.2)

## Model specification ----------

#formula

# formula.season<-as.formula(~ log.depth + bpi + sediment + 
#                            sst + chl + ssha + fsle + salinity + 
#                            mld + sst.fprob + chl.fprob)
# 
# formula.season.ext<-as.formula(~ log.depth + plcurv + bpi + sediment + 
#                                  sst + chl + ssha + fsle + salinity + 
#                                  mld + uo + vo + sst.fprob + chl.fprob +
#                                  sst.fmean + chl.fmean)
# 
# formula.season.fin<-as.formula(~ log.depth + bpi + sediment + 
#                                  sst + chl + ssha + fsle + salinity + 
#                                  mld + sst.fprob + chl.fprob+
#                                  sst.fmean + chl.fmean
#                                )

formula.season.fin2<-as.formula(~ log.depth + bpi + sediment + 
                                 sst + chl + ssha + fsle + salinity + 
                                 mld + sst.fprob + chl.fprob)

# formula.season.falltest<-as.formula(~ log.depth +bpi + sediment + 
#                                   chl + ssha + fsle + salinity + 
#                                   mld + sst.fprob + chl.fprob)




form<-formula.season.fin2
iterations<-30000
burnin <-10000
#fall needs a longer run #I'm actually thinking that they all may need longer runs tbh
mlist <- list(ng = iterations, burnin = burnin, 
              typeNames = "CA",
              # effort = elist, 
              random = "space.grid.id.2",
              FULL = TRUE,
              PREDICTX = TRUE #speed up computation make F
)

#same model inputs (with random effect) for cross validation
save(ms.data, fall.ms.data, win.ms.data, spr.ms.data, sum.ms.data, fall.y.dat, 
     fall.y.dat.spue, win.y.dat, win.y.dat.spue, spr.y.dat, spr.y.dat.spue,
     sum.y.dat, sum.y.dat.spue, fall.x.dat, win.x.dat, spr.x.dat, sum.x.dat,
     fall.e.data, win.e.data, spr.e.data, sum.e.data, formula.season.fin2, 
     form, iterations, burnin, mlist, 
     file = paste0(fit_mods_dir, "/seasonal_data_gjamfin_no0sizeFFA_32km_SPUE_11cov_RE128km_30kiter.Rdata"))

#test with no random effect
# mlist_norand <- list(ng = 20000, burnin = 10000, 
#                     typeNames = "CA",
#                     #effort = elist, 
#                     # random = "space.grid.id",
#                     FULL = TRUE,
#                     PREDICTX = T #speed up computation make F
# )



## Run the seasonal models -----------

### With random effect ----

runs<-1:3
# seeds<-c(1,1.5,0.5) #for testing convergence issues 

fall.ms.mod <- list()
win.ms.mod <- list()
spr.ms.mod <- list()
sum.ms.mod <- list()
  
for (i in 1:length(runs)) {
  print(i)
  set.seed(i)
  # set.seed(seeds[i])
 
#final mods with spue and 20k samples
tic("fall")
fall.ms.mod[[i]] <- gjam(form, xdata = fall.x.dat, ydata = fall.y.dat.spue, modelList = mlist)
toc()
tic("winter")
win.ms.mod[[i]] <- gjam(form, xdata = win.x.dat, ydata = win.y.dat.spue, modelList = mlist)
toc()
tic("spring")
spr.ms.mod[[i]] <- gjam(form, xdata = spr.x.dat, ydata = spr.y.dat.spue, modelList = mlist)
toc()
tic("summer")
sum.ms.mod[[i]] <- gjam(form, xdata = sum.x.dat, ydata = sum.y.dat.spue, modelList = mlist)
toc()

}


save(fall.ms.mod, file = paste0(fit_mods_dir, "/fall_ms_gjam_no0sizeFFA_32km_10spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))
save(win.ms.mod, file = paste0(fit_mods_dir, "/winter_ms_gjam_no0sizeFFA_32km_7spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))
save(spr.ms.mod, file = paste0(fit_mods_dir, "/spring_ms_gjam_no0sizeFFA_32km_10spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))
save(sum.ms.mod, file = paste0(fit_mods_dir, "/summer_ms_gjam_no0sizeFFA_32km_7spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))

#more fall model testing
# save(fall.ms.mod, file = paste0(fit_mods_dir, "/test_fall_ms_gjam_drop_no0sizeFFA_32km_10spp_SPUE_11cov_randeff_20k_iter_3runs_setseed1_3.Rdata"))
# save(fall.ms.mod, file = paste0(fit_mods_dir, "/test_fall_ms_gjam_drop_no0sizeFFA_32km_10spp_SPUE_11cov_randeff_20k_iter_3runs_setseed1_15.Rdata"))
# save(fall.ms.mod, file = paste0(fit_mods_dir, "/test_fall_ms_gjam_drop_no0sizeFFA_32km_10spp_SPUE_11cov_NORE_20k_iter_3runs_setseed1_1_3.Rdata"))
# save(fall.ms.mod, file = paste0(fit_mods_dir, "/test_fall_ms_gjam_drop_no0sizeFFA_32km_10spp_SPUE_11cov_RE128km_20k_iter_3runs_setseed1_1_3.Rdata"))


## Initial output and chain stabilization -----------
# fall_test_mod_dir<-paste0(fit_mods_dir, '/fall_test_converge_gjam_output')

for (i in 1:length(runs)) {
gjamPlot(fall.ms.mod[[i]], plotPars = list(PLOTALLY = T, GRIDPLOTS = T, CLUSTERPLOTS = T, 
                                           SMALLPLOTS = F, SAVEPLOTS = T,
                                           outFolder = paste0(fit_mods_dir, '/final_gjamPlot_output/gjam_fall_plots_run_', i)))
                                           # outFolder = paste0(fall_test_mod_dir, '/test_20kiter_30kburn_withRE_seed1_3_new/test_fall_plots_run_', i)))
gjamPlot(win.ms.mod[[i]], plotPars = list(PLOTALLY = T, GRIDPLOTS = T, CLUSTERPLOTS = T,
                                          SMALLPLOTS = F, SAVEPLOTS = T,
                                          outFolder = paste0(fit_mods_dir, '/final_gjamPlot_output/gjam_winter_plots_run_', i)))
gjamPlot(spr.ms.mod[[i]], plotPars = list(PLOTALLY = T, GRIDPLOTS = T, CLUSTERPLOTS = T,
                                          SMALLPLOTS = F, SAVEPLOTS = T,
                                          outFolder = paste0(fit_mods_dir, '/final_gjamPlot_output/gjam_spring_plots_run_', i)))
gjamPlot(sum.ms.mod[[i]], plotPars = list(PLOTALLY = T, GRIDPLOTS = T, CLUSTERPLOTS = T,
                                          SMALLPLOTS = F, SAVEPLOTS = T,
                                          outFolder = paste0(fit_mods_dir, '/final_gjamPlot_output/gjam_summer_plots_run_', i)))
}

## Evaluate convergence  -------------------------------------------------

#Extract the chains from the gjam model objects and format as mcmc list

#'[Note: only have to do this once, unless of course we have to rerun the model]  


season= c("fall", "winter", "spring","summer")
season_mods  <- list(fall.ms.mod, win.ms.mod, spr.ms.mod, sum.ms.mod)
names(season_mods) <- season

# season_mods  <- list(fall.ms.mod)
chain_out <- list()

for (i in 1:length(season)) {

  for (j in 1:length(runs)) {
#get the standardized beta chains from the models
bchains <- season_mods[[i]][[j]]$chains$bgibbs
#remove the burn-in
bchains <- bchains[(burnin+1):iterations,]

#get the correlation chains
corchains <- season_mods[[i]][[j]]$chains$sgibbs
corchains <- corchains[(burnin+1):iterations,]

chain_out[[j]] <- as.mcmc(cbind(bchains, corchains))

}

  chain_out <- as.mcmc.list(list(chain_out[[1]], chain_out[[2]], chain_out[[3]]))
  fname<-paste0(fit_mods_dir, "/", season[[i]], "_beta_cor_chains_noburnin_notThinned.Rdata")
  save(chain_out, file = fname)

}


#check convergence 

#fall
load(file = paste0(fit_mods_dir, "/fall_beta_cor_chains_noburnin_notThinned.Rdata"))


list.names <- colnames(chain_out[[1]])
#fall_test_converge_3chain

for (i in 1:length(list.names)) {
  
  #plot chains
  fname <- paste0(fit_mods_dir, "/fall_32km_3chain_diagnostics/fall_", list.names[i], "_traceplot_noburnin_notThinned.png")
  png(filename = fname, units = "in", width = 8, height = 4.5, res = 75)
  plot(chain_out[, list.names[i]], main = list.names[i])
  dev.off()
}  
  #thin the posterior - for testing
  # out_mcmc<-as.mcmc.list(list(as.mcmc(chain_out[[1]][seq(1,nrow(chain_out[[1]]),20),]), 
  #                             as.mcmc(chain_out[[2]][seq(1,nrow(chain_out[[2]]),20),]), 
  #                             as.mcmc(chain_out[[3]][seq(1,nrow(chain_out[[3]]),20),])))
  # 
  # out_mcmc<-as.mcmc.list(list(as.mcmc(chain_out[[1]][1:10000,]), 
  #                             as.mcmc(chain_out[[2]][1:10000,]), 
  #                             as.mcmc(chain_out[[3]][1:10000,])))
  # 
  # plot(out_mcmc[, list.names[97]], main = list.names[97])
  
  #effective sample size
  ess <- effectiveSize(chain_out)
  ess[which(ess<=1000)]
  fname <- paste0(fit_mods_dir, "/fall_32km_3chain_diagnostics/fall_ess_noburnin_notThinned_30iter.csv")
  write.csv(ess, file = fname, row.names = TRUE)
  
  # PSRF/R hat/Gelman diagnostic  
  psrf<-gelman.diag(chain_out, multivariate=FALSE)$psrf
  psrf[which(psrf[,1]>=1.1),]
  psrf[which(psrf[,1]>=1.01),]
  fname <- paste0(fit_mods_dir, "/fall_32km_3chain_diagnostics/fall_psrf_noburnin_notThinned_30kiter.csv")
  write.csv(psrf, file = fname, row.names = TRUE)



#winter
load(file = paste0(fit_mods_dir, "/winter_beta_cor_chains_noburnin_notThinned.Rdata"))


list.names <- colnames(chain_out[[1]])

for (i in 1:length(list.names)) {
  
  #plot chains
  fname <- paste0(fit_mods_dir, "/winter_32km_3chain_diagnostics/winter_", list.names[i], "_traceplot_noburnin_notThinned.png")
  png(filename = fname, units = "in", width = 8, height = 4.5, res = 75)
  plot(chain_out[, list.names[i]], main = list.names[i])
  dev.off()
}  


#effective sample size
ess <- effectiveSize(chain_out)
ess[which(ess<=1000)]
fname <- paste0(fit_mods_dir, "/winter_32km_3chain_diagnostics/winter_ess_noburnin_notThinned.csv")
write.csv(ess, file = fname, row.names = TRUE)

# PSRF/R hat/Gelman diagnostic  
psrf<-gelman.diag(chain_out, multivariate=FALSE)$psrf
psrf[which(psrf[,1]>=1.1),]
psrf[which(psrf[,1]>=1.01),]
fname <- paste0(fit_mods_dir, "/winter_32km_3chain_diagnostics/winter_psrf_noburnin_notThinned.csv")
write.csv(psrf, file = fname, row.names = TRUE)

#spring
load(file = paste0(fit_mods_dir, "/spring_beta_cor_chains_noburnin_notThinned.Rdata"))

list.names <- colnames(chain_out[[1]])

for (i in 1:length(list.names)) {
  
  #plot chains
  fname <- paste0(fit_mods_dir, "/spring_32km_3chain_diagnostics/spring_", list.names[i], "_traceplot_noburnin_notThinned.png")
  png(filename = fname, units = "in", width = 8, height = 4.5, res = 75)
  plot(chain_out[, list.names[i]], main = list.names[i])
  dev.off()
}  


#effective sample size
ess <- effectiveSize(chain_out)
ess[which(ess<=1000)]
fname <- paste0(fit_mods_dir, "/spring_32km_3chain_diagnostics/spring_ess_noburnin_notThinned.csv")
write.csv(ess, file = fname, row.names = TRUE)

# PSRF/R hat/Gelman diagnostic  
psrf<-gelman.diag(chain_out, multivariate=FALSE)$psrf
psrf[which(psrf[,1]>=1.1),]
psrf[which(psrf[,1]>=1.01),]
fname <- paste0(fit_mods_dir, "/spring_32km_3chain_diagnostics/spring_psrf_noburnin_notThinned.csv")
write.csv(psrf, file = fname, row.names = TRUE)


#summer
load(file = paste0(fit_mods_dir, "/summer_beta_cor_chains_noburnin_notThinned.Rdata"))

list.names <- colnames(chain_out[[1]])

for (i in 1:length(list.names)) {
  
  #plot chains
  fname <- paste0(fit_mods_dir, "/summer_32km_3chain_diagnostics/summer_", list.names[i], "_traceplot_noburnin_notThinned.png")
  png(filename = fname, units = "in", width = 8, height = 4.5, res = 75)
  plot(chain_out[, list.names[i]], main = list.names[i])
  dev.off()
}  


#effective sample size
ess <- effectiveSize(chain_out)
ess[which(ess<=1000)]
fname <- paste0(fit_mods_dir, "/summer_32km_3chain_diagnostics/summer_ess_noburnin_notThinned.csv")
write.csv(ess, file = fname, row.names = TRUE)

# PSRF/R hat/Gelman diagnostic  
psrf<-gelman.diag(chain_out, multivariate=FALSE)$psrf
psrf[which(psrf[,1]>=1.1),]
psrf[which(psrf[,1]>=1.01),]
fname <- paste0(fit_mods_dir, "/summer_32km_3chain_diagnostics/summer_psrf_noburnin_notThinned.csv")
write.csv(psrf, file = fname, row.names = TRUE)

## Full Model Evaluation  -------------------------------------------------

### with random effect ----
runs<-1:3

# Pearson correlation for observed vs. predicted

#fall
spp.list<-colnames(fall.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = fall.ms.mod[[k]]$inputs$y[,i], y = fall.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("fall_32km_SPUE_11cov_RE128km_30kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/fall_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/fall_full_species_correlations/fall_correlations_32km_SPUE_11cov_RE128km_30kiter.csv"),
          row.names = TRUE)

#Winter

spp.list<-colnames(win.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = win.ms.mod[[k]]$inputs$y[,i], y = win.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("winter_32km_SPUE_11cov_RE128km_30kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/winter_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/winter_full_species_correlations/winter_correlations_32km_SPUE_11cov_RE128km_30kiter.csv"),
          row.names = TRUE)

#Spring

spp.list<-colnames(spr.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = spr.ms.mod[[k]]$inputs$y[,i], y = spr.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    # spec<-ggscatter(data=df, x = "x", y = "y",
    #                 size = 0.5,
    #                 add = "reg.line",  # Add regressin line
    #                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    #                 conf.int = TRUE, # Add confidence interval
    #                 title = spp.list[i],
    #                 xlab = "Observed", ylab = "Predicted"
    # )
    # 
    # spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))
    # 
    # ggsave(filename = paste0("spring_32km_SPUE_11cov_RE128km_30kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
    #        path = paste0(eval_mods_dir, "/spring_full_species_correlations"),
    #        width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/spring_full_species_correlations/spring_correlations_32km_SPUE_11cov_RE128km_30kiter.csv"),
          row.names = TRUE)

#Summer

spp.list<-colnames(sum.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = sum.ms.mod[[k]]$inputs$y[,i], y = sum.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("summer_32km_SPUE_11cov_RE128km_30kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/summer_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/summer_full_species_correlations/summer_correlations_32km_SPUE_11cov_RE128km_30kiter.csv"),
          row.names = TRUE)


#Check DIC and covariate VIF/correlations
View(fall.ms.mod)
View(fall.ms.mod$inputs$designTable)
View(win.ms.mod)
View(win.ms.mod$inputs$designTable)
View(spr.ms.mod)
View(spr.ms.mod$inputs$designTable)
View(sum.ms.mod)
View(sum.ms.mod$inputs$designTable)


### No random effect ----
# runs<-1:3

# Pearson correlation for observed vs. predicted

#fall
spp.list<-colnames(fall.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = fall.ms.mod[[k]]$inputs$y[,i], y = fall.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("fallmod_NOrandeff_32km_SPUE_11cov_20kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/fall_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/fall_full_species_correlations/fall_NOrandeff_correlations_32km_SPUE_11cov_20kiter.csv"),
          row.names = TRUE)

#Winter

spp.list<-colnames(win.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = win.ms.mod[[k]]$inputs$y[,i], y = win.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("wintermod_NOrandeff_32km_SPUE_11cov_20kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/winter_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/winter_full_species_correlations/winter_NOrandeff_correlations_32km_SPUE_11cov_20kiter.csv"),
          row.names = TRUE)

#Spring

spp.list<-colnames(spr.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = spr.ms.mod[[k]]$inputs$y[,i], y = spr.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("springmod_NOrandeff_32km_SPUE_11cov_20kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/spring_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/spring_full_species_correlations/spring_NOrandeff_correlations_32km_SPUE_11cov_20kiter.csv"),
          row.names = TRUE)

#Summer

spp.list<-colnames(sum.y.dat.spue)
full.mod.corr<-matrix(NA, nrow = length(spp.list), ncol = length(runs))

for (k in 1:length(runs)) {
  
  for (i in 1:length(spp.list)) {
    
    df<-as.data.frame(cbind(x = sum.ms.mod[[k]]$inputs$y[,i], y = sum.ms.mod[[k]]$prediction$ypredMu[,i]))
    
    spec<-ggscatter(data=df, x = "x", y = "y",
                    size = 0.5,
                    add = "reg.line",  # Add regressin line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    title = spp.list[i],
                    xlab = "Observed", ylab = "Predicted"
    )

    spec + stat_cor(method = "pearson", label.x = 0, label.y = max(spec$data$y))

    ggsave(filename = paste0("summermod_NOrandeff_32km_SPUE_11cov_20kiter_run_", k, "_", spp.list[i], ".png"), device = "png",
           path = paste0(eval_mods_dir, "/summer_full_species_correlations"),
           width = 3, height = 3, units = "in")

    val<-cor(df, method = "pearson", use = "pairwise")
    full.mod.corr[i,k]<-val[1,2]
    
  }
}

rownames(full.mod.corr)<-spp.list

full.mod.corr <- as.data.frame(full.mod.corr) %>% 
  mutate(mean.cor = rowMeans(full.mod.corr),
         sd.cor = rowSds(as.matrix(full.mod.corr)))

full.mod.corr <- full.mod.corr %>% 
  mutate(CV = sd.cor/mean.cor)

write.csv(full.mod.corr, 
          file = paste0(eval_mods_dir, "/summer_full_species_correlations/summer_NOrandeff_correlations_32km_SPUE_11cov_20kiter.csv"),
          row.names = TRUE)


#Check DIC and covariate VIF/correlations ----
View(fall.ms.mod)
View(fall.ms.mod$inputs$designTable)
View(win.ms.mod)
View(win.ms.mod$inputs$designTable)
View(spr.ms.mod)
View(spr.ms.mod$inputs$designTable)
View(sum.ms.mod)
View(sum.ms.mod$inputs$designTable)








