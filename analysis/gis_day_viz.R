# GIS Day Lightning Talk Scripts
# Keaton Wilson
# keatonwilson@me.com
# 2019-11-4


# Packages ----------------------------------------------------------------

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
library(ggridges)
library(tidyverse)
library(gganimate)
library(dismo)
library(ENMeval)

# Loading Data ------------------------------------------------------------

#Loading in raw occurence data
swallowtail = read_csv("./data/raw_data/swallowtail_data.csv")
swallowtail = swallowtail[,-1] %>%
  dplyr::select(longitude, latitude, date, year, time_frame)

swallowtail_t1 = swallowtail %>%
  filter(time_frame == "T1")

swallowtail_t2 = swallowtail %>%
  filter(time_frame == "T2")

#hostplant
hostplant = read_csv("./data/raw_data/hostplant_data.csv")
hostplant = hostplant[,-1]

# #bioclim environmental variables
# bioclim.data <- raster::getData(name = "worldclim",
#                                 var = "bio",
#                                 res = 2.5,
#                                 path = "./data/")

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


#Loading in model objects
mx_best_st_t1 = readRDS("./data/swallowtail_t1.rds")
mx_best_st_t2 = readRDS("./data/swallowtail_t2.rds")
mx_best_hp_1_t2 = readRDS("./data/hostplant_1_t1.rds")
mx_best_hp_1_t1 = readRDS("./data/hostplant_1_t2.rds")
mx_best_hp_2_t2 = readRDS("./data/hostplant_2_t1.rds")
mx_best_hp_2_t1 = readRDS("./data/hostplant_2_t2.rds")
mx_best_hp_3_t2 = readRDS("./data/hostplant_3_t1.rds")
mx_best_hp_3_t1 = readRDS("./data/hostplant_3_t2.rds")


# Geographic Mapping Data ---------------------------------------

#Pulling in polygons for states and provinces
#Getting map data
usa = getData(country = 'USA', level = 1)

#extract states (need to uppercase everything)
to_remove = c("Alaska", "Hawaii", "North Dakota", "South Dakota", "Montana", 
              "Wyoming", "Idaho", "Washington", "Oregon", "Nevada", "California", 
              "Arizona", "Utah", "New Mexico", "Colorado", "Nebraska", "Texas", 
              "Oklahoma", "Kansas")

#filtering
mapping = usa[-match(toupper(to_remove), toupper(usa$NAME_1)),]

#simplying polygons
simple_map_US = gSimplify(mapping, tol = 0.01, topologyPreserve = TRUE)

#Pulling Canada Province data
can = getData(country = 'CAN', level = 1)
province = c("Ontario")
can_mapping = can[match(toupper(c("Ontario", "Québec", "New Brunswick", "Prince Edward Island", "Nova Scotia")), toupper(can$NAME_1)),]
simple_map_can = gSimplify(can_mapping, tol = 0.01, topologyPreserve = TRUE)

#Great lakes issues
lakes <- rgdal::readOGR("./data/raw_data/ne_10m_lakes/ne_10m_lakes.shp")
lakes = lakes[lakes$scalerank==0,]
lakes = crop(lakes, geographic.extent)
# Predictions -------------------------------------------------------------

#Swallowtail
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

#Z. americanum
#Predictions from full model (Hostplant 1 T1)
predict_presence_hp_1_t1 = dismo::predict(object = mx_best_hp_1_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_1_t1 <- as(predict_presence_hp_1_t1, "SpatialPixelsDataFrame")
pred_sp_df_hp_1_t1 <- as.data.frame(pred_sp_hp_1_t1)
colnames(pred_sp_df_hp_1_t1) <- c("value", "x", "y")

#Predictions from full model (Hostplant T2)
predict_presence_hp_1_t2 = dismo::predict(object = mx_best_hp_1_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_1_t2 <- as(predict_presence_hp_1_t2, "SpatialPixelsDataFrame")
pred_sp_df_hp_1_t2 <- as.data.frame(pred_sp_hp_1_t2)
colnames(pred_sp_df_hp_1_t2) <- c("value", "x", "y")

#Z. clava-herculis
predict_presence_hp_2_t1 = dismo::predict(object = mx_best_hp_2_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_2_t1 <- as(predict_presence_hp_2_t1, "SpatialPixelsDataFrame")
pred_sp_df_hp_2_t1 <- as.data.frame(pred_sp_hp_2_t1)
colnames(pred_sp_df_hp_2_t1) <- c("value", "x", "y")

#Predictions from full model (Hostplant 2 T2)
predict_presence_hp_2_t2 = dismo::predict(object = mx_best_hp_2_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_2_t2 <- as(predict_presence_hp_2_t2, "SpatialPixelsDataFrame")
pred_sp_df_hp_2_t2 <- as.data.frame(pred_sp_hp_2_t2)
colnames(pred_sp_df_hp_2_t2) <- c("value", "x", "y")

#P. trifoliata
predict_presence_hp_3_t1 = dismo::predict(object = mx_best_hp_3_t1, x = bv_t1, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_3_t1 <- as(predict_presence_hp_3_t1, "SpatialPixelsDataFrame")
pred_sp_df_hp_3_t1 <- as.data.frame(pred_sp_hp_3_t1)
colnames(pred_sp_df_hp_3_t1) <- c("value", "x", "y")

#Predictions from full model (Hostplant 3 T2)
predict_presence_hp_3_t2 = dismo::predict(object = mx_best_hp_3_t2, x = bv_t2, ext = geographic.extent, args = 'outputformat=cloglog')

pred_sp_hp_3_t2 <- as(predict_presence_hp_3_t2, "SpatialPixelsDataFrame")
pred_sp_df_hp_3_t2 <- as.data.frame(pred_sp_hp_3_t2)
colnames(pred_sp_df_hp_3_t2) <- c("value", "x", "y")

# Thresholds --------------------------------------------------------------
#Threshold maps

#Loading the evaluate objects from the model building script
evaluations = readRDS("./data/evaluations.rds")

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

#Threshold dfs
hostplant_clean = 
  hostplant %>%
  mutate(name = ifelse(str_detect(name, "Zanthoxylum americanum"), "Zanthoxylum americanum",
                       ifelse(str_detect(name, "Zanthoxylum clava-herculis"), "Zanthoxylum claca-herculis",
                              "Ptelea trifoliata")))
#hp_df_t1
hp_thresholds_df_t1 = 
  bind_rows(hp_1_t1_threshold, hp_2_t1_threshold, hp_3_t1_threshold, .id = "id") %>%
  mutate(name = ifelse(id == 1, "Zanthoxylum americanum",
                       ifelse(id == 2, "Zanthoxylum clava-herculis", "Ptelea trifoliata"))) %>%
  dplyr::select(-id)

#hp_t2
hp_thresholds_df_t2 = 
  bind_rows(hp_1_t2_threshold, hp_2_t2_threshold, hp_3_t2_threshold, .id = "id") %>%
  mutate(name = ifelse(id == 1, "Zanthoxylum americanum",
                       ifelse(id == 2, "Zanthoxylum clava-herculis", "Ptelea trifoliata"))) %>%
  dplyr::select(-id)


# All Occurences ----------------------------------------------------------

anim_map = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_point(data = swallowtail, aes(x = longitude, y = latitude), col = "yellow", alpha = 0.5, size = 5) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.position="right") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 60)) +
  #coord_equal(ylim = c(22, 50), xlim = c(-100, -65)) +
  theme_nothing(legend = TRUE) +
  coord_quickmap() +
  transition_states(year(swallowtail$date), transition_length = 0, state_length = 50) +
  shadow_mark() +
  labs(title = 'Year: {closest_state}')

?gifski_renderer
animation_1 = animate(anim_map, duration = 20, 
                      renderer = gifski_renderer(loop=FALSE),
                      width = 2500, height = 2500)
save_animation(animation_1, "./output/st_occurences.gif")


# Swallowtail Threshold Maps -----------------------------------------
g13 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "grey10") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "grey10") +
  geom_tile(data = st_t1_threshold, aes(x=x, y=y), fill = "lightgrey") + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_point(data = swallowtail_t1, aes(x = longitude, y = latitude), alpha = 0.5, color = "yellow", shape = 3, size = 0.5) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("1960-1999") +
  coord_quickmap()

#T2
g14 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "grey10") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "grey10") +
  geom_tile(data=st_t2_threshold, aes(x=x, y=y), fill = "lightgrey") + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_point(data = swallowtail_t2, aes(x = longitude, y = latitude), alpha = 0.2, color = "yellow", shape = 3, size = 0.5) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("2000-2019") +
  coord_quickmap()

fig_2_b = ggarrange(g13, g14, common.legend = TRUE)
ggsave("./output/gis_threshold_maps.png", fig_2_b)



# Plant threshold maps ----------------------------------------------------


cols = c("yellow", "red", "blue")

g13 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "white") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "white") +
  geom_tile(data = hp_thresholds_df_t1, aes(x = x, y = y, fill = name), alpha = 0.6) +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  # geom_point(data = swallowtail_t2, aes(x = longitude, y = latitude), alpha = 0.2, color = "yellow", shape = 3) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25, color = "grey50") +
  # geom_point(data = hostplant %>%
  #              filter(time_frame == "T1") %>%
  #              filter(str_detect(name, "Zanthoxylum americanum")),
  #            aes(x = longitude, y = latitude), color ="#F8766D",  alpha = 0.5, shape = 3) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
        text = element_text(size = 18)) +
  theme_nothing(legend = TRUE) +
  scale_fill_manual(values = cols,
                    name = NULL, 
                    breaks = c("Zanthoxylum americanum", 
                               "Zanthoxylum clava-herculis", 
                               "Ptelea trifoliata"),
                    labels = c(expression(italic("Zanthoxylum americanum")), 
                               expression(italic("Zanthoxylum clava-herculis")), 
                               expression(italic("Ptelea trifoliata")))) +
  ggtitle("1959 - 1999") +
  coord_quickmap() +
  scale_color_manual(guide = FALSE, 
                     values = cols)

#T2
g14 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "white") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "white") +
  geom_tile(data = hp_thresholds_df_t2, aes(x = x, y = y, fill = name), alpha = 0.6) +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  # geom_point(data = swallowtail_t2, aes(x = longitude, y = latitude), alpha = 0.2, color = "yellow", shape = 3) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25, color = "grey50") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"), 
        text = element_text(size = 18)) +
  theme_nothing(legend = TRUE) +
  scale_fill_manual(values = cols, 
                    name = NULL, 
                    breaks = c("Zanthoxylum americanum", 
                               "Zanthoxylum clava-herculis", 
                               "Ptelea trifoliata"),
                    labels = c(expression(italic("Zanthoxylum americanum")), 
                               expression(italic("Zanthoxylum clava-herculis")), 
                               expression(italic("Ptelea trifoliata")))) +
  ggtitle("2000 - 2018") +
  coord_quickmap() +
  scale_color_manual(guide = FALSE, 
                     values = cols)

fig_3_b = ggarrange(g13, g14, common.legend = TRUE, legend = "bottom")


# Northern limit figure ---------------------------------------------------

n_limit_st = pred_sp_df_st_t1 %>%
  bind_rows(pred_sp_df_st_t2, .id = "timeframe") %>%
  mutate(occ = ifelse(value >= st_t1_threshold & timeframe == 1, 1,
                      ifelse(value >= st_t2_threshold & timeframe == 2, 1, 0))) %>%
  group_by(x, timeframe) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()

ggplot(n_limit_st, aes(x = x, y = max_y)) +
  geom_point(aes(color = timeframe)) +
  theme_classic() 



ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "grey10") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "grey10") +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_point(data = n_limit_st, aes(x = x, y = max_y, color = timeframe), alpha = 0.5) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("2000-2019") +
  coord_quickmap(ylim = c(37, 65)) +
  scale_color_manual(values = cols)


# ggmap approach 
library(ggmap)
register_google(key = "AIzaSyDyAqUc4o9p_DOBSF_JOXH5c_JXPqoU4Yw")
#Pulling in polygons for states and provinces
#Getting map data
usa = getData(country = 'USA', level = 1)

#extract states (need to uppercase everything)
to_remove = c("Alaska", "Hawaii", "North Dakota", "South Dakota", "Montana", 
              "Wyoming", "Idaho", "Washington", "Oregon", "Nevada", "California", 
              "Arizona", "Utah", "New Mexico", "Colorado", "Nebraska", "Texas", 
              "Oklahoma", "Kansas")

#filtering
mapping = usa

#simplying polygons
simple_map_US = gSimplify(mapping, tol = 0.01, topologyPreserve = TRUE)

#Pulling Canada Province data
can = getData(country = 'CAN', level = 1)
province = c("Ontario")
can_mapping = can[match(toupper(c("Ontario", "Québec", "New Brunswick", "Prince Edward Island", "Nova Scotia")), toupper(can$NAME_1)),]
simple_map_can = gSimplify(can, tol = 0.01, topologyPreserve = TRUE)

#Great lakes issues
lakes <- rgdal::readOGR("./data/raw_data/ne_10m_lakes/ne_10m_lakes.shp")
lakes = lakes[lakes$scalerank==0,]
lakes = crop(lakes, geographic.extent)
map = get_map(location = c(-100, 37, -65, 50), maptype = "toner-lines")  
ggmap(map) +
  geom_point(data = n_limit_st, aes(x = x, y = max_y, color = timeframe), alpha = 0.7, size = 3, shape = 3)

cols = c('purple', 'darkgreen')

north_range = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "white") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "white") +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_point(data = n_limit_st, aes(x = x, y = max_y, color = timeframe), alpha = 0.7, size = 6, shape = 15) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25, color = "grey50") +

  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), size = 0.25) +
  theme_nothing(legend = TRUE) +
  coord_quickmap(xlim = c(-100, -60), ylim = c(37, 50)) +
  scale_color_manual(labels = c("1960-1999", "2000-2018"), 
                     values = cols, 
                     name = NULL) +
  theme(panel.background = element_rect(fill = "black"), 
        legend.text = element_text(size = 15), 
        legend.key = element_rect(colour = NA, fill = NA))

ggsave("./output/northern_range_gis.png", north_range, width = 38, height = 17, units = "in")
