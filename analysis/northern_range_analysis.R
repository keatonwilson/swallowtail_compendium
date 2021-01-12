#Leroux et al. style Northern Range Limit things
#Keaton Wilson
#keatonwilson@me.com
#2019-12-30

#packages
#packages and libraries
library(tidyverse)
library(dismo)
library(lubridate)
library(rgdal)
library(sp)
library(raster)
library(maptools)
library(ggmap)
library(viridis)
library(ggthemes)
library(rgeos)
library(maps)
library(ggpubr)
library(blockCV)
library(ENMeval)

#Loading in raw occurence data
swallowtail = read_csv("./data/raw_data/swallowtail_data.csv")
swallowtail = swallowtail[,-1] %>%
  dplyr::select(longitude, latitude, date, year, time_frame)

#hostplant
hostplant = read_csv("./data/raw_data/hostplant_data.csv")
hostplant = hostplant[,-1]

#Environmental Data 
bv_t1 = raster::brick("./data/raw_data/biovar_avg_t1.grd")
bv_t2 = raster::brick("./data/raw_data/biovar_avg_t2.grd")

#renaming
names_seq = paste("Bio",seq(1:19), sep = "")
names(bv_t1) = names_seq
names(bv_t2) = names_seq

# Determine geographic extent of our data
max_lat_swallowtail <- ceiling(max(swallowtail$latitude))
min_lat_swallowtail <- floor(min(swallowtail$latitude))
max_lon_swallowtail <- ceiling(max(swallowtail$longitude))
min_lon_swallowtail <- floor(min(swallowtail$longitude))
geographic.extent <- extent(x = c(min_lon_swallowtail, max_lon_swallowtail, min_lat_swallowtail, max_lat_swallowtail))

#Loading threshold objects and models
mx_best_st_t1 = readRDS("./data/swallowtail_t1.rds")
mx_best_st_t2 = readRDS("./data/swallowtail_t2.rds")
mx_best_hp_1_t2 = readRDS("./data/hostplant_1_t2.rds")
mx_best_hp_1_t1 = readRDS("./data/hostplant_1_t1.rds")
mx_best_hp_2_t2 = readRDS("./data/hostplant_2_t2.rds")
mx_best_hp_2_t1 = readRDS("./data/hostplant_2_t1.rds")
mx_best_hp_3_t2 = readRDS("./data/hostplant_3_t2.rds")
mx_best_hp_3_t1 = readRDS("./data/hostplant_3_t1.rds")

# evaluation objs and thresholds
ev_st_t1 = evaluations[[1]]
ev_st_t2 = evaluations[[2]]
ev_hp_1_t1 = evaluations[[3]]
ev_hp_1_t2 = evaluations[[4]]
ev_hp_2_t1 = evaluations[[5]]
ev_hp_2_t2 = evaluations[[6]]
ev_hp_3_t1 = evaluations[[7]]
ev_hp_3_t2 = evaluations[[8]]

#finding the threshold for presence/absence for each model
st_t1_threshold = threshold(ev_st_t1, 'spec_sens')
st_t2_threshold = threshold(ev_st_t2, 'spec_sens')
hp_1_t1_threshold = threshold(ev_hp_1_t1, 'spec_sens')
hp_1_t2_threshold = threshold(ev_hp_1_t2, 'spec_sens')
hp_2_t1_threshold = threshold(ev_hp_2_t1, 'spec_sens')
hp_2_t2_threshold = threshold(ev_hp_2_t2, 'spec_sens')
hp_3_t1_threshold = threshold(ev_hp_3_t1, 'spec_sens')
hp_3_t2_threshold = threshold(ev_hp_3_t2, 'spec_sens')


#Predictions
#Predictions from full model (Hostplant 3 T1)
predict_presence_hp_3_t1 = dismo::predict(object = mx_best_hp_3_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_3_t1 <- as(predict_presence_hp_3_t1, "SpatialPixelsDataFrame")
pred_sp_df_hp_3_t1 <- as.data.frame(pred_sp_hp_3_t1)
colnames(pred_sp_df_hp_3_t1) <- c("value", "x", "y")

#Predictions from full model (Hostplant 3 T2)
predict_presence_hp_3_t2 = dismo::predict(object = mx_best_hp_3_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_3_t2 <- as(predict_presence_hp_3_t2, "SpatialPixelsDataFrame")
pred_sp_df_hp_3_t2 <- as.data.frame(pred_sp_hp_3_t2)
colnames(pred_sp_df_hp_3_t2) <- c("value", "x", "y")

#Predictions from full model (Hostplant 2 T1)
predict_presence_hp_2_t1 = dismo::predict(object = mx_best_hp_2_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_2_t1 <- as(predict_presence_hp_2_t1, "SpatialPixelsDataFrame")
pred_sp_df_hp_2_t1 <- as.data.frame(pred_sp_hp_2_t1)
colnames(pred_sp_df_hp_2_t1) <- c("value", "x", "y")

#Predictions from full model (Hostplant 2 T2)
predict_presence_hp_2_t2 = dismo::predict(object = mx_best_hp_2_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_2_t2 <- as(predict_presence_hp_2_t2, "SpatialPixelsDataFrame")
pred_sp_df_hp_2_t2 <- as.data.frame(pred_sp_hp_2_t2)
colnames(pred_sp_df_hp_2_t2) <- c("value", "x", "y")

predict_presence_hp_1_t1 = dismo::predict(object = mx_best_hp_1_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_1_t1 <- as(predict_presence_hp_1_t1, "SpatialPixelsDataFrame")
pred_sp_df_hp_1_t1 <- as.data.frame(pred_sp_hp_1_t1)
colnames(pred_sp_df_hp_1_t1) <- c("value", "x", "y")

#Predictions from full model (Hostplant T2)
predict_presence_hp_1_t2 = dismo::predict(object = mx_best_hp_1_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_1_t2 <- as(predict_presence_hp_1_t2, "SpatialPixelsDataFrame")
pred_sp_df_hp_1_t2 <- as.data.frame(pred_sp_hp_1_t2)
colnames(pred_sp_df_hp_1_t2) <- c("value", "x", "y")

#Predictions from full model (Swallowtail T1)
predict_presence_st_t1 = dismo::predict(object = mx_best_st_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_st_t1 <- as(predict_presence_st_t1, "SpatialPixelsDataFrame")
pred_sp_df_st_t1 <- as.data.frame(pred_sp_st_t1)
colnames(pred_sp_df_st_t1) <- c("value", "x", "y")


#Predictions from full model (Swallowtail T2)
predict_presence_st_t2 = dismo::predict(object = mx_best_st_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_st_t2 <- as(predict_presence_st_t2, "SpatialPixelsDataFrame")
pred_sp_df_st_t2 <- as.data.frame(pred_sp_st_t2)
colnames(pred_sp_df_st_t2) <- c("value", "x", "y")

#building filtered dataframes of predictions
st_t1_threshold = pred_sp_df_st_t1 %>%
  filter(value > st_t1_threshold)

st_t2_threshold = pred_sp_df_st_t2 %>%
  filter(value > st_t2_threshold)

hp_1_t1_threshold = pred_sp_df_hp_1_t1 %>%
  filter(value > hp_1_t1_threshold)

hp_1_t2_threshold = pred_sp_df_hp_1_t2 %>%
  filter(value > hp_1_t2_threshold)

hp_2_t1_threshold = pred_sp_df_hp_2_t1 %>%
  filter(value > hp_2_t1_threshold)

hp_2_t2_threshold = pred_sp_df_hp_2_t2 %>%
  filter(value > hp_2_t2_threshold)

hp_3_t1_threshold = pred_sp_df_hp_3_t1 %>%
  filter(value > hp_3_t1_threshold)

hp_3_t2_threshold = pred_sp_df_hp_3_t2 %>%
  filter(value > hp_3_t2_threshold)

#binding
threshold_df_st = bind_rows("t1" = st_t1_threshold, "t2" = st_t2_threshold, .id = "timeframe")
threshold_df_hp_1 = bind_rows("t1" = hp_1_t1_threshold, "t2" = hp_1_t2_threshold, .id = "timeframe")
threshold_df_hp_2 = bind_rows("t1" = hp_2_t1_threshold, "t2" = hp_2_t2_threshold, .id = "timeframe")
threshold_df_hp_3 = bind_rows("t1" = hp_3_t1_threshold, "t2" = hp_3_t2_threshold, .id = "timeframe")

#Trying to generate northern range for swallowtails t1
threshold_df_st %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y)) %>%
  ggplot(aes(x = max_y, fill = timeframe)) +
  geom_density(alpha = 0.8) +
  theme_classic() +
  ylab("Kernel Density Estimate") +
  xlab("Maximum Latitude (º)") +
  scale_fill_discrete(name = "Time Frame", 
                      labels = c("t1", "t2"))

max_y = threshold_df_st %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y))

max_y_t1 = max_y %>% 
  ungroup() %>%
  filter(timeframe == "t1") %>%
  pull(max_y)

max_y_t2 = max_y %>% 
               ungroup() %>%
               filter(timeframe == "t2") %>%
               pull(max_y)

st_test = t.test(max_y_t1, max_y_t2, paired = TRUE)

#So the problem with a paired t-test is that there are different numbers, because the there are chunks of the northern border that may not be present in either data set. 
threshold_df_st %>%
  group_by(timeframe) %>%
  summarize(n())
#finding the threshold for presence/absence for each model
st_t1_threshold = threshold(ev_st_t1, 'spec_sens')
st_t2_threshold = threshold(ev_st_t2, 'spec_sens')
hp_1_t1_threshold = threshold(ev_hp_1_t1, 'spec_sens')
hp_1_t2_threshold = threshold(ev_hp_1_t2, 'spec_sens')
hp_2_t1_threshold = threshold(ev_hp_2_t1, 'spec_sens')
hp_2_t2_threshold = threshold(ev_hp_2_t2, 'spec_sens')
hp_3_t1_threshold = threshold(ev_hp_3_t1, 'spec_sens')
hp_3_t2_threshold = threshold(ev_hp_3_t2, 'spec_sens')
#So, we need to do the predictions first... then change threshold values below that to NA, instead of just removing

# limits
n_limit_st = pred_sp_df_st_t1 %>%
  bind_rows(pred_sp_df_st_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= st_t1_threshold & timeframe == 1, 1,
                        ifelse(value >= st_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()


n_limit_hp_1 = pred_sp_df_hp_1_t1 %>%
  bind_rows(pred_sp_df_hp_1_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= hp_1_t1_threshold & timeframe == 1, 1,
                      ifelse(value >= hp_2_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()

limits_combined = bind_rows(list(st = n_limit_st, hp_1 = n_limit_hp_1), 
                            .id = "species")

# paired t.tests
st_test = t.test(limits_combined %>%
                   filter(species == "st" & timeframe == 1) %>%
                   pull(max_y), 
                 limits_combined %>%
                   filter(species == "st" & timeframe == 2) %>%
                   pull(max_y), 
                 paired = TRUE)

# graphical validation
ggplot(limits_combined %>%
         filter(species == "st"), 
       aes(x = max_y, fill = timeframe)) +
  geom_density()

hp_1_test = t.test(limits_combined %>%
                     filter(species == "hp_1" & timeframe == 1) %>%
                     pull(max_y), 
                   limits_combined %>%
                     filter(species == "hp_1" & timeframe == 2) %>%
                     pull(max_y), 
                   paired = TRUE)

# graphical validation
ggplot(limits_combined %>%
         filter(species == "hp_1"), 
       aes(x = max_y, fill = timeframe)) +
  geom_density()
