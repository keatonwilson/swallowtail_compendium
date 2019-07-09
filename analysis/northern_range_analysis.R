#Leroux et al. style Northern Range Limit things
#Keaton Wilson
#keatonwilson@me.com
#2019-05-09

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
swallowtail = read_csv("./data/swallowtail_data.csv")
swallowtail = swallowtail[,-1] %>%
  dplyr::select(longitude, latitude, date, year, time_frame)

#hostplant
hostplant = read_csv("./data/hostplant_data.csv")
hostplant = hostplant[,-1]

#Environmental Data 
bv_t1 = raster::brick("./data/terraclim/biovar_avg_t1.grd")
bv_t2 = raster::brick("./data/terraclim/biovar_avg_t2.grd")

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
mx_best_st_t1 = readRDS("./models/full_best_st_t1.rds")
mx_best_st_t2 = readRDS("./models/full_best_st_t2.rds")
mx_best_hp_1_t2 = readRDS("./models/full_best_hp_1_t2.rds")
mx_best_hp_1_t1 = readRDS("./models/full_best_hp_1_t1.rds")
mx_best_hp_2_t2 = readRDS("./models/full_best_hp_2_t2.rds")
mx_best_hp_2_t1 = readRDS("./models/full_best_hp_2_t1.rds")
mx_best_hp_3_t2 = readRDS("./models/full_best_hp_3_t2.rds")
mx_best_hp_3_t1 = readRDS("./models/full_best_hp_3_t1.rds")

#Evaluate objects
ev_st_t1 = readRDS("./data/ev_st_t1.RDS")
ev_st_t2 = readRDS("./data/ev_st_t2.RDS")
ev_hp_1_t1 = readRDS("./data/ev_hp_1_t1.RDS")
ev_hp_1_t2 = readRDS("./data/ev_hp_1_t2.RDS")
ev_hp_2_t1 = readRDS("./data/ev_hp_2_t1.RDS")
ev_hp_2_t2 = readRDS("./data/ev_hp_2_t2.RDS")
ev_hp_3_t1 = readRDS("./data/ev_hp_3_t1.RDS")
ev_hp_3_t2 = readRDS("./data/ev_hp_3_t2.RDS")

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

#
n_limit_st = pred_sp_df_st_t1 %>%
  bind_rows(pred_sp_df_st_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= st_t1_threshold & timeframe == 1, 1,
                        ifelse(value >= st_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()


pred_sp_df_st_t1 %>%
  filter(x > -70) %>%
  mutate(occ = ifelse(value >= st_t1_threshold, 1, 0)) %>%
  filter(occ == 1)


#Density plots
fig_5_a = ggplot(data = n_limit_st, aes(x = max_y, fill = timeframe)) +
  geom_density(alpha = 0.7) +
  theme_classic() +
  geom_vline(data = n_limit_st[n_limit_st$timeframe == "1",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  geom_vline(data = n_limit_st[n_limit_st$timeframe == "2",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  xlab("Max Northern Occurence (ºC)") +
  ylab("Kernel Density Estimate") +
  scale_fill_discrete(name = "Time Frame", 
                      labels = c("t1 - 1959-1999", "t2 - 2000-2018")) + 
  labs(title = expression(italic("P. cresphontes"))) +
  xlim(c(38, 50)) +
  ylim(c(0, 0.90))

#Mapping to see if it makes sense
ggplot(n_limit_st, aes(x = x, y = max_y)) +
  geom_point(aes(color = timeframe)) +
  theme_classic() 

#Are there different numbers?
n_limit_st %>%
  group_by(timeframe) %>%
  summarize(n_na = sum(is.na(max_y)))

#Breaking apart to feed into the paired t-test
st_max_y_t1 = n_limit_st %>%
  filter(timeframe == "1") %>%
  pull(max_y)

st_max_y_t2 = n_limit_st %>%
  filter(timeframe == "2") %>%
  pull(max_y)

test = t.test(st_max_y_t1, st_max_y_t2, paired = TRUE)
sd(st_max_y_t1, na.rm = TRUE)
median(st_max_y_t2, na.rm = TRUE) - median(st_max_y_t1, na.rm = TRUE)
median(st_max_y_t1, na.rm = TRUE)

#Same thing for Z.americanum
n_limit_hp_1 = pred_sp_df_hp_1_t1 %>%
  bind_rows(pred_sp_df_hp_1_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= hp_1_t1_threshold & timeframe == 1, 1,
                      ifelse(value >= hp_1_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()

#Breaking apart to feed into the paired t-test
hp_1_max_y_t1 = n_limit_hp_1 %>%
  filter(timeframe == "1") %>%
  pull(max_y)

hp_1_max_y_t2 = n_limit_hp_1 %>%
  filter(timeframe == "2") %>%
  pull(max_y)

#Tests and calculations
test = t.test(hp_1_max_y_t1, hp_1_max_y_t2, paired = TRUE)
sd(hp_1_max_y_t2, na.rm = TRUE)
median(hp_1_max_y_t2, na.rm = TRUE) - median(hp_1_max_y_t1, na.rm = TRUE)
median(st_max_y_t1, na.rm = TRUE)

#Can we see if there is a difference between the the butterfly and host plant 1?
n_limit_st_hp1_t2 = pred_sp_df_st_t2 %>%
  bind_rows(pred_sp_df_hp_1_t2, .id = "species") %>%
  mutate(occ = ifelse(value >= st_t2_threshold & species == 1, 1,
                      ifelse(value >= hp_1_t2_threshold & species == 2, 1, 0))) %>%
  group_by(x, species) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()

#Breaking apart to feed into the paired t-test
st_max_y_t2 = n_limit_st_hp1_t2 %>%
  filter(species == "1") %>%
  pull(max_y)

hp_1_max_y_t2 = n_limit_st_hp1_t2 %>%
  filter(species == "2") %>%
  pull(max_y)

#Tests and calculations
test = t.test(st_max_y_t2, hp_1_max_y_t2, paired = TRUE)
sd(hp_1_max_y_t2, na.rm = TRUE)
median(hp_1_max_y_t2, na.rm = TRUE) - median(hp_1_max_y_t1, na.rm = TRUE)
median(st_max_y_t1, na.rm = TRUE)

#Can we see if there is a difference between the the butterfly and host plant 1?
n_limit_st_hp1_t1 = pred_sp_df_st_t1 %>%
  bind_rows(pred_sp_df_hp_1_t1, .id = "species") %>%
  mutate(occ = ifelse(value >= st_t2_threshold & species == 1, 1,
                      ifelse(value >= hp_1_t2_threshold & species == 2, 1, 0))) %>%
  group_by(x, species) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()

#Breaking apart to feed into the paired t-test
st_max_y_t1 = n_limit_st_hp1_t1 %>%
  filter(species == "1") %>%
  pull(max_y)

hp_1_max_y_t1 = n_limit_st_hp1_t1 %>%
  filter(species == "2") %>%
  pull(max_y)

#Tests and calculations
test = t.test(st_max_y_t1, hp_1_max_y_t1, paired = TRUE)
sd(hp_1_max_y_t2, na.rm = TRUE)
median(hp_1_max_y_t2, na.rm = TRUE) - median(hp_1_max_y_t1, na.rm = TRUE)
median(st_max_y_t1, na.rm = TRUE)

median(hp_1_max_y_t1, na.rm = TRUE) - median(st_max_y_t1, na.rm = TRUE)
median(hp_1_max_y_t2, na.rm = TRUE) - median(st_max_y_t2, na.rm = TRUE)

#Density plots
fig_5_d = ggplot(data = n_limit_st_hp1_t2, aes(x = max_y, fill = species)) +
  geom_density(alpha = 0.7) +
  theme_classic() +
  geom_vline(data = n_limit_st_hp1_t2[n_limit_st_hp1_t2$species == "1",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  geom_vline(data = n_limit_st_hp1_t2[n_limit_st_hp1_t2$species == "2",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  xlab("Max Northern Occurence (ºC)") +
  ylab("Kernel Density Estimate") +
  scale_fill_manual(name = "Species", 
                    labels = c("P. cresphontes", "Z. americanum"),
                    values = c("#E69F00", "#56B4E9")) +
  xlim(c(38, 50)) +
  ylim(c(0, 0.90)) +
  ggtitle("t2")

#T1
n_limit_st_hp1_t1 = pred_sp_df_st_t1 %>%
  bind_rows(pred_sp_df_hp_1_t1, .id = "species") %>%
  mutate(occ = ifelse(value >= st_t1_threshold & species == 1, 1,
                      ifelse(value >= hp_1_t1_threshold & species == 2, 1, 0))) %>%
  group_by(x, species) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()


#Density plots
fig_5_c = ggplot(data = n_limit_st_hp1_t1, aes(x = max_y, fill = species)) +
  geom_density(alpha = 0.7) +
  theme_classic() +
  geom_vline(data = n_limit_st_hp1_t1[n_limit_st_hp1_t1$species == "1",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  geom_vline(data = n_limit_st_hp1_t1[n_limit_st_hp1_t1$species == "2",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  xlab("Max Northern Occurence (ºC)") +
  ylab("Kernel Density Estimate") +
  scale_fill_manual(name = "Species", 
                      labels = c("P. cresphontes", "Z. americanum"),
                    values = c("#E69F00", "#56B4E9")) +
  xlim(c(38, 50)) +
  ylim(c(0, 0.90)) +
  ggtitle("t1")

#T1
n_limit_st_hp1 = pred_sp_df_hp_1_t1 %>%
  bind_rows(pred_sp_df_hp_1_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= hp_1_t1_threshold & timeframe == 1, 1,
                      ifelse(value >= hp_1_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()


#Density plots
fig_5_b = ggplot(data = n_limit_st_hp1, aes(x = max_y, fill = timeframe)) +
  geom_density(alpha = 0.7) +
  theme_classic() +
  geom_vline(data = n_limit_st_hp1[n_limit_st_hp1$timeframe == "1",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  geom_vline(data = n_limit_st_hp1[n_limit_st_hp1$timeframe == "2",], 
             aes(xintercept = median(max_y, na.rm = TRUE)),
             lty = 2) +
  xlab("Max Northern Occurence (ºC)") +
  ylab("Kernel Density Estimate") +
  scale_fill_discrete(name = "Time Frame", 
                      labels = c("t1 - 1959-1999", "t2 - 2000-2018")) + 
  labs(title = expression(italic("Z. americanum")))
  xlim(c(38, 50)) +
  ylim(c(0, 0.90))


fig_5 = ggarrange(fig_5_a, fig_5_b, fig_5_c, fig_5_d, labels = "AUTO")
ggsave(plot = fig_5, filename = "./output/fig_5.png", device = "png")


#Same thing for P. trifoliata
n_limit_hp_3 = pred_sp_df_hp_3_t1 %>%
  bind_rows(pred_sp_df_hp_3_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= hp_3_t1_threshold & timeframe == 1, 1,
                      ifelse(value >= hp_3_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()

#Breaking apart to feed into the paired t-test
hp_3_max_y_t1 = n_limit_hp_3 %>%
  filter(timeframe == "1") %>%
  pull(max_y)

hp_3_max_y_t2 = n_limit_hp_3 %>%
  filter(timeframe == "2") %>%
  pull(max_y)

#Tests and calculations
test = t.test(hp_3_max_y_t1, hp_3_max_y_t2, paired = TRUE)
sd(hp_3_max_y_t2, na.rm = TRUE)
median(hp_3_max_y_t2, na.rm = TRUE) - median(hp_3_max_y_t1, na.rm = TRUE)
median(st_max_y_t1, na.rm = TRUE)


threshold_df_hp_1 %>%
  filter(timeframe == "t1") %>%
  filter(value > hp_1_t1_threshold) %>%
  summarize(median_lat = median(y))

threshold_df_hp_1 %>%
  filter(timeframe == "t2") %>%
  filter(value > hp_1_t2_threshold) %>%
  summarize(median_lat = median(y))

threshold_df_hp_3 %>%
  filter(timeframe == "t1") %>%
  filter(value > hp_3_t1_threshold) %>%
  summarize(median_lat = median(y))

threshold_df_hp_3 %>%
  filter(timeframe == "t2") %>%
  filter(value > hp_3_t2_threshold) %>%
  summarize(median_lat = median(y))

