# Multispecies Predator Models: GJAM Model Predictions -------------------------

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
library(tictoc)
library(raster)
library(viridis)
library(vegan)

#Directories
input_dir<-"P:/NYSERDA Ecosystem Dynamics/Species Interactions Modeling/Aerial Survey Multispecies Predator Analysis/data"
fit_mods_dir<-"P:/NYSERDA Ecosystem Dynamics/Species Interactions Modeling/Aerial Survey Multispecies Predator Analysis/output/gjam models/fitted models"
mod_pred_dir<-"P:/NYSERDA Ecosystem Dynamics/Species Interactions Modeling/Aerial Survey Multispecies Predator Analysis/output/gjam models/model predictions"



#Load fitted model objects----


# # 32km full model with seasonal interactions
# load(file = paste0(fit_mods_dir, "/ms_gjam_interaction_fullmod_no0sizeFFA_32km_10spp_SPUE_allcov_randeff_4000iter.Rdata"))


# 32km seasonal models - RE 128km
load(file = paste0(fit_mods_dir, "/fall_ms_gjam_no0sizeFFA_32km_10spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))
load(file = paste0(fit_mods_dir, "/winter_ms_gjam_no0sizeFFA_32km_7spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))
load(file = paste0(fit_mods_dir, "/spring_ms_gjam_no0sizeFFA_32km_10spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))
load(file = paste0(fit_mods_dir, "/summer_ms_gjam_no0sizeFFA_32km_7spp_SPUE_11cov_RE128km_30kiter_3runs.Rdata"))

#Load polygon with bays and estuaries clipped ----------------------------------
clip1 <- readOGR(input_dir, 'ffa_bound_updt2_finalextent')
clip1 <- raster::aggregate(clip1)

#Load prediction grid data -----------------------------------------------------

#'[Choose scale "16km" or "32km"
sp_scale <- "32km"

#load the polygon pgrid
load(file = paste0(input_dir, 'FFA_parea_v4_32km_RE128km.Rdata'))
plot(pgrid)

#fall
load(file= paste0(input_dir,"/fall_ms_pts_pgrid_cov_", sp_scale,".Rdata")) #this is the correct grid

#winter
load(file= paste0(input_dir,"/winter_ms_pts_pgrid_cov_", sp_scale,".Rdata")) #this is the correct grid

#spring
load(file= paste0(input_dir,"/spring_ms_pts_pgrid_cov_", sp_scale,".Rdata")) #this is the correct grid

#summer
load(file= paste0(input_dir,"/summer_ms_pts_pgrid_cov_", sp_scale,".Rdata")) #this is the correct grid


#load the updated pts pgrid with space.grid.id.2 - to add to seasonal pgrids
load(file=paste0(input_dir,"/ffa_pts_pgrid_", sp_scale, "_RE128km.Rdata"))

#Isolate the space.grid.id.2 column and drop geometry
pt_pgrid<-st_as_sf(pt_pgrid)
plot(pt_pgrid)
tmp<-pt_pgrid %>%
  dplyr::select(c(surv.grid.id, space.grid.id.2))
plot(tmp)
tmp<-st_drop_geometry(tmp)

#Data wrangling ----------------------------------------------------------------

## Prediction Effort ----
#get prediction grid cell size in km2 #this would be perfect effort
pgrid <- pgrid %>%
  mutate(area = (st_area(.) %>% as.numeric())/1E6)

#get the mean effort per grid cell per season
fall.mean.eff <- fall.ms.data %>% 
  group_by(surv.grid.id) %>% 
  summarize(mean.eff = mean(effort.km2))

# summary(fall.mean.eff$mean.eff)
max(round(fall.mean.eff$mean.eff))
#'[100]

winter.mean.eff <- win.ms.data %>% 
  group_by(surv.grid.id) %>% 
  summarize(mean.eff = mean(effort.km2))

# summary(winter.mean.eff$mean.eff)
max(round(winter.mean.eff$mean.eff))
#'[91]


spring.mean.eff <- spr.ms.data %>% 
  group_by(surv.grid.id) %>% 
  summarize(mean.eff = mean(effort.km2))
max(round(spring.mean.eff$mean.eff))
#'[113]

summer.mean.eff <- sum.ms.data %>% 
  group_by(surv.grid.id) %>% 
  summarize(mean.eff = mean(effort.km2))
max(round(summer.mean.eff$mean.eff))
#'[97]

#seasonal mean of the max mean
mean(c(100,91,113,97))
#'[100.25]

## Fall ----
fall_ms_pts_pgrid<-st_as_sf(fall_ms_pts_pgrid)
fall_ms_pts_pgrid<-left_join(fall_ms_pts_pgrid, tmp, by = 'surv.grid.id')

fall_ms_pts_pgrid<-fall_ms_pts_pgrid %>% 
  dplyr::mutate(xcoord = sf::st_coordinates(.)[,1],
                ycoord = sf::st_coordinates(.)[,2],
                Lon = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,1],
                Lat = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,2],
  )

fall_ms_pts_pgrid<-st_drop_geometry(fall_ms_pts_pgrid)

#'[only for the full model with season effect]
# #add season  - named same as fitted model 
# fall_ms_pts_pgrid$season.2<-"Fall"

#add effort
#perfect effort
fall_ms_pts_pgrid$effort.km2<-st_drop_geometry(pgrid$area)
#max mean effort/grid cell
fall_ms_pts_pgrid$max.eff.km2 <- max(round(fall.mean.eff$mean.eff))

## Winter ----

winter_ms_pts_pgrid<-st_as_sf(winter_ms_pts_pgrid)
winter_ms_pts_pgrid<-left_join(winter_ms_pts_pgrid, tmp, by = 'surv.grid.id')

winter_ms_pts_pgrid<-winter_ms_pts_pgrid %>% 
  dplyr::mutate(xcoord = sf::st_coordinates(.)[,1],
                ycoord = sf::st_coordinates(.)[,2],
                Lon = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,1],
                Lat = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,2],
  )

winter_ms_pts_pgrid<-st_drop_geometry(winter_ms_pts_pgrid)

#add effort
#perfect effort
winter_ms_pts_pgrid$effort.km2<-st_drop_geometry(pgrid$area)
#max mean effort/grid cell
winter_ms_pts_pgrid$max.eff.km2 <- max(round(fall.mean.eff$mean.eff))

## Spring ----
spring_ms_pts_pgrid<-st_as_sf(spring_ms_pts_pgrid)
spring_ms_pts_pgrid<-left_join(spring_ms_pts_pgrid, tmp, by = 'surv.grid.id')

spring_ms_pts_pgrid<-spring_ms_pts_pgrid %>% 
  dplyr::mutate(xcoord = sf::st_coordinates(.)[,1],
                ycoord = sf::st_coordinates(.)[,2],
                Lon = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,1],
                Lat = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,2],
  )

spring_ms_pts_pgrid<-st_drop_geometry(spring_ms_pts_pgrid)

#add effort
#perfect effort
spring_ms_pts_pgrid$effort.km2<-st_drop_geometry(pgrid$area)
#max mean effort/grid cell
spring_ms_pts_pgrid$max.eff.km2 <- max(round(fall.mean.eff$mean.eff))


## Summer ----
summer_ms_pts_pgrid<-st_as_sf(summer_ms_pts_pgrid)
summer_ms_pts_pgrid<-left_join(summer_ms_pts_pgrid, tmp, by = 'surv.grid.id')

summer_ms_pts_pgrid<-summer_ms_pts_pgrid %>% 
  dplyr::mutate(xcoord = sf::st_coordinates(.)[,1],
                ycoord = sf::st_coordinates(.)[,2],
                Lon = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,1],
                Lat = sf::st_coordinates(st_transform(., crs = '+proj=longlat +datum=WGS84'))[,2],
  )

summer_ms_pts_pgrid<-st_drop_geometry(summer_ms_pts_pgrid)

#add effort
#perfect effort
summer_ms_pts_pgrid$effort.km2<-st_drop_geometry(pgrid$area)
#max mean effort/grid cell
summer_ms_pts_pgrid$max.eff.km2 <- max(round(fall.mean.eff$mean.eff))

#'[can attach these back to the polygon pgrid later for plotting




#Set up prediction runs --------------------------------------------------------

#for predictions to new data
#xdata = the prediction grid xdata, effort = original effort or standard effort (1),
# for the random effect I think I just need to include the random effect column in 
# xdata and it will automatically pull the random variable

#'[new_xdata has to perfectly match the OG model xdata (order and variables) or it messes up]


## fall ----
fall_new_xdata<-fall_ms_pts_pgrid %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob) 

fall_new_xdata$space.grid.id.2<-as.factor(fall_new_xdata$space.grid.id.2)

## winter ----
winter_new_xdata<-winter_ms_pts_pgrid %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob) 

winter_new_xdata$space.grid.id.2<-as.factor(winter_new_xdata$space.grid.id.2)

## spring ----
spring_new_xdata<-spring_ms_pts_pgrid %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob) 

spring_new_xdata$space.grid.id.2<-as.factor(spring_new_xdata$space.grid.id.2)

## summer ----
summer_new_xdata<-summer_ms_pts_pgrid %>% 
  dplyr::select(c(space.grid.id.2, log_depth, bpi, sediment, 
                  sst, chl, ssha, fsle, salinity, mld, sst_fprob, chl_fprob)) %>% 
  dplyr::rename(log.depth = log_depth, sst.fprob = sst_fprob,
                chl.fprob = chl_fprob) 

summer_new_xdata$space.grid.id.2<-as.factor(summer_new_xdata$space.grid.id.2)

# #factor levels must match the fitted model
# new_xdata$space.grid.id.2<-factor(new_xdata$space.grid.id.2, levels = rand_levels)
# new_xdata$year<-factor(new_xdata$year, levels = levels(x.ms.dat$year))
# new_xdata$season<-factor(new_xdata$season, levels = levels(x.ms.dat$season))

str(new_xdata)

#test<-new_xdata[complete.cases(new_xdata),] # Remove incomplete cases

#'[Note: for the spue models we do not add effort because it is already incorporated
#'[      into the ydata - but I will need the full area later to convert back to abundance] ]
# #new effort
# new_effort<-fall_ms_pts_pgrid$effort.km2
# new_effort_list<-list(columns = 1:ncol(ms.mod$inputs$y),
#                       values = new_effort)




#Make predictions --------------------------------------------------------------

#fall
newdata <- list(xdata = fall_new_xdata)
ms.mod <- fall.ms.mod[[1]]
fall_pred <- gjamPredict(output = ms.mod, newdata = newdata, FULL = FALSE)

#winter
newdata = list(xdata = winter_new_xdata)
ms.mod <- win.ms.mod[[1]]
winter_pred<-gjamPredict(output = ms.mod, newdata = newdata, FULL = FALSE)

#spring
newdata = list(xdata = spring_new_xdata)
ms.mod <- spr.ms.mod[[1]]
spring_pred<-gjamPredict(output = ms.mod, newdata = newdata, FULL = FALSE)

#summer
newdata = list(xdata = summer_new_xdata)
ms.mod <- sum.ms.mod[[1]]
summer_pred<-gjamPredict(output = ms.mod, newdata = newdata, FULL = FALSE)


#Format predictions --------------------------------------------------------------

##fall ----

#convert spue predictions back to abundance
#perfect effort 
# fall_pred_abund<-(as.data.frame(fall_pred$sdList$yMu)*fall_ms_pts_pgrid$effort.km2)^2

#max mean effort/grid cell
fall_pred_abund<-(as.data.frame(fall_pred$sdList$yMu)*fall_ms_pts_pgrid$max.eff.km2)^2

#species richness
fall_pred_abund <- fall_pred_abund %>% 
  mutate(Cum.abund.all = rowSums(round(fall_pred_abund[,1:10])),
         Cum.abund.predator = rowSums(round(fall_pred_abund[,1:9])),
         S.all = rowSums(ifelse(fall_pred_abund[,1:10]>1, 1, 0)),
         S.predator = rowSums(ifelse(fall_pred_abund[,1:9]>1, 1, 0)),
         H.all = diversity(round(fall_pred_abund[,1:10])),#Shannon diversity
         H.predator = diversity(round(fall_pred_abund[,1:9])) 
         )

fall_pred_df<-cbind(fall_pred_abund,
                    surv.grid.id = fall_ms_pts_pgrid$surv.grid.id,
                    # perf.effort.km2 = fall_ms_pts_pgrid$effort.km2,
                    max.eff.km2 = fall_ms_pts_pgrid$max.eff.km2,
                    xcoord = fall_ms_pts_pgrid$xcoord,
                    ycoord = fall_ms_pts_pgrid$ycoord,
                    Lon = fall_ms_pts_pgrid$Lon, 
                    Lat = fall_ms_pts_pgrid$Lon)


fall_pred_sf<-left_join(pgrid, fall_pred_df, by = "surv.grid.id", .keep_all = TRUE)


#clip to final prediction grid and save
clip1<- st_as_sf(clip1)
plot(clip1)

fall_pred_sf<-st_intersection(fall_pred_sf, clip1)

#save raw predictions and converted predictions
save(fall_pred, fall_pred_df, fall_pred_sf, 
     file = paste0(mod_pred_dir, "/fall_32km_RE128km_newdata_predictions.Rdata"))

##winter ----

#convert spue predictions back to abundance
#perfect effort 
# winter_pred_abund<-(as.data.frame(winter_pred$sdList$yMu)*winter_ms_pts_pgrid$effort.km2)^2

#max mean effort/grid cell
winter_pred_abund<-(as.data.frame(winter_pred$sdList$yMu)*winter_ms_pts_pgrid$max.eff.km2)^2
#no ffa in winter
winter_pred_abund <- winter_pred_abund %>% 
  mutate(Cum.abund.all = rowSums(round(fall_pred_abund[,1:7])),
         S.all = rowSums(ifelse(fall_pred_abund[,1:7]>1, 1, 0)),
         H.all = diversity(round(fall_pred_abund[,1:7])),#Shannon diversity
  )

winter_pred_df<-cbind(winter_pred_abund,
                    surv.grid.id = winter_ms_pts_pgrid$surv.grid.id,
                    # perf.effort.km2 = winter_ms_pts_pgrid$effort.km2,
                    max.eff.km2 = winter_ms_pts_pgrid$max.eff.km2,
                    xcoord = winter_ms_pts_pgrid$xcoord,
                    ycoord = winter_ms_pts_pgrid$ycoord,
                    Lon = winter_ms_pts_pgrid$Lon, 
                    Lat = winter_ms_pts_pgrid$Lon)

winter_pred_sf<-left_join(pgrid, winter_pred_df, by = "surv.grid.id", .keep_all = TRUE)

#clip to final prediction grid and save
# clip1<- st_as_sf(clip1)
# plot(clip1)
winter_pred_sf<-st_intersection(winter_pred_sf, clip1)

#save raw predictions and converted predictions
save(winter_pred, winter_pred_df, winter_pred_sf, 
     file = paste0(mod_pred_dir, "/winter_32km_RE128km_newdata_predictions.Rdata"))


##spring ----

#convert spue predictions back to abundance

#perfect effort 
# spring_pred_abund<-(as.data.frame(spring_pred$sdList$yMu)*spring_ms_pts_pgrid$effort.km2)^2

#max mean effort/grid cell
spring_pred_abund<-(as.data.frame(spring_pred$sdList$yMu)*spring_ms_pts_pgrid$max.eff.km2)^2

spring_pred_df<-cbind(spring_pred_abund,
                    surv.grid.id = spring_ms_pts_pgrid$surv.grid.id,
                    # perf.effort.km2 = spring_ms_pts_pgrid$effort.km2,
                    max.eff.km2 = spring_ms_pts_pgrid$max.eff.km2,
                    xcoord = spring_ms_pts_pgrid$xcoord,
                    ycoord = spring_ms_pts_pgrid$ycoord,
                    Lon = spring_ms_pts_pgrid$Lon, 
                    Lat = spring_ms_pts_pgrid$Lon)


spring_pred_sf<-left_join(pgrid, spring_pred_df, by = "surv.grid.id", .keep_all = TRUE)


#clip to final prediction grid and save
clip1<- st_as_sf(clip1)
plot(clip1)

spring_pred_sf<-st_intersection(spring_pred_sf, clip1)

#save raw predictions and converted predictions
save(spring_pred, spring_pred_df, spring_pred_sf, 
     file = paste0(mod_pred_dir, "/spring_32km_RE128km_newdata_predictions.Rdata"))

##summer ----

#convert spue predictions back to abundance
#perfect effort 
# summer_pred_abund<-(as.data.frame(summer_pred$sdList$yMu)*summer_ms_pts_pgrid$effort.km2)^2

#max mean effort/grid cell
summer_pred_abund<-(as.data.frame(summer_pred$sdList$yMu)*summer_ms_pts_pgrid$max.eff.km2)^2

summer_pred_df<-cbind(summer_pred_abund,
                    surv.grid.id = summer_ms_pts_pgrid$surv.grid.id,
                    # perf.effort.km2 = summer_ms_pts_pgrid$effort.km2,
                    max.eff.km2 = summer_ms_pts_pgrid$max.eff.km2,
                    xcoord = summer_ms_pts_pgrid$xcoord,
                    ycoord = summer_ms_pts_pgrid$ycoord,
                    Lon = summer_ms_pts_pgrid$Lon, 
                    Lat = summer_ms_pts_pgrid$Lon)


summer_pred_sf<-left_join(pgrid, summer_pred_df, by = "surv.grid.id", .keep_all = TRUE)


# Clip to final prediction grid and save

# clip1<- st_as_sf(clip1)
# plot(clip1)

summer_pred_sf<-st_intersection(summer_pred_sf, clip1)

#save raw predictions and converted predictions
save(summer_pred, summer_pred_df, summer_pred_sf, 
     file = paste0(mod_pred_dir, "/summer_32km_RE128km_newdata_predictions.Rdata"))


#plot predictions -----------------------------------------------------------


spp.list<-colnames(fall_pred_abund)
seas = c("Summer", "Fall", "Spring", "Winter")

spp.list<-colnames(winter_pred_abund)
    
    
#predictions
  for (i in 1:length(spp.list)) {
    
    ggplot() +
      geom_sf(data = winter_pred_sf, aes_string(fill = spp.list[i]), color = NA)  +
      ggtitle(paste0("Fall_", spp.list[i])) +
      scale_fill_gradientn(colours = viridis(200)) #+
      # geom_sf(data = clip1, color = "black", size = 1, fill = NA)
    ggsave(filename = paste0("inter_mod_pred_32km_", seas[j], "_", spp.list[i], ".png"), device = "png",
           path = mod_pred_dir, width = 3, height = 3.5, units = "in")

}




























































































