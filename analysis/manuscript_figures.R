#Manuscript Figures
#Keaton Wilson
#keatonwilson@me.com
#2019-06-19

# Packages ----------------------------------------------------------------
#packages
library(tidyverse)
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

# Figure 1 -  Evidence of northward shift from raw occurence data------------------------------------------------- 
#Panel A

#Maximum latitude by year figure
swallowtail_inset = swallowtail %>%
  group_by(year) %>%
  summarize(max_lat = max(latitude), 
            median_lat = median(latitude),
            n = n()) %>%
  filter(max_lat > 35) #Filtering - there are some weird years that only have a few records at really low lattitudes.

fig_1_a = ggplot(data = swallowtail_inset, aes(x = year, y = max_lat, size = n)) +
  geom_point(alpha = 0.8) +
  # geom_smooth(data = swallowtail_inset %>%
  #               filter(year < 2000), aes(x = year, y = max_lat), method = "lm", show.legend = FALSE) +
  # geom_smooth(data = swallowtail_inset %>%
  #               filter(year >= 2000), aes(x = year, y = max_lat), method = "lm", show.legend = FALSE) +
  geom_smooth(data = swallowtail_inset, show.legend = FALSE) +
  theme_classic() +
  scale_size_continuous(name = "Number of Observations") +
  labs(x = "Year", y = "Maximum Latitude (º)") +
  geom_vline(xintercept = 2000, lty = 2) +
  annotate(geom = "text", label = "Timeframe Break Point", x = 1992, y = 47.5) +
  theme(axis.title = element_text(size = 18),
        legend.position = "top")


#Panel B
#Ridgeplot
years = swallowtail %>%
  group_by(year) %>%
  summarize(n = n()) %>%
  filter(n > 5) %>%
  filter(year != 2019) %>%
  dplyr::select(year) %>%
  pull()

fig_1_b = swallowtail %>%
  filter(year %in% years) %>%
  ggplot(aes(y = factor(year), x = latitude)) +
  geom_density_ridges(scale = 4) +
  geom_linerange(x = 46.8139, ymin = 1, ymax = 38, lty = 2, size = 0.25, alpha = 0.6) +
  geom_linerange(x = 45.4215, ymin = 1, ymax = 39, lty = 2, size = 0.25, alpha = 0.6) +
  geom_linerange(x= 43.6532, ymin = 1, ymax = 40, lty = 2, size = 0.25, alpha = 0.6) +
  geom_linerange(x = 39.7684, ymin = 1, ymax = 38, lty = 2, size = 0.25, alpha = 0.6) +
  theme_classic() +
  annotate(geom = "text", 
           label = "Indianapolis", 
           x = 39.7684,
           y = 39.5) +
  annotate(geom = "text", 
           label = "Toronto", 
           x = 43.6532,
           y = 41) +
  annotate(geom = "text", 
           label = "Ottawa", 
           x = 45.4215,
           y = 39.9) +
  annotate(geom = "text", 
           label = "Quebec City", 
           x = 48.5,
           y = 39) +
  coord_cartesian(ylim = c(0,41)) +
  xlab("Latitude (º)") +
  ylab("Year") +
  theme(axis.title = element_text(size = 18))

fig_1 = ggarrange(fig_1_a, fig_1_b, labels = "AUTO")
ggsave(plot = fig_1, filename = "./output/fig_1.png", width = 16, height = 8.5, units = "in")

# Figure 2 - Swallowtail cloglog and Threshold----------------------------------------------------------------
#Panel A
g1 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_st_t1, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="right") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  #coord_equal(ylim = c(22, 50), xlim = c(-100, -65)) +
  theme_nothing(legend = TRUE) +
  ggtitle("1960 - 1999") +
  coord_quickmap()


#Plotting Swallowtail T2
g2 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_st_t2, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="right") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  theme_nothing(legend = TRUE) +
  ggtitle("2000 - 2019") +
  coord_quickmap()

fig_2_a = ggarrange(g1, g2, common.legend = TRUE, legend = "top")
annotate_figure(fig_2_a,
                top = text_grob("Papilio cresphontes", face = "italic", size = 22))

#Figure 2 Panel B
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

fig_2 = ggarrange(fig_2_a, fig_2_b, nrow = 2, align = "hv", 
                  heights = c(1.25, 1))
ggsave(plot = fig_2, filename = "./output/fig_2.png", device = "png",
       height = 11, width = 8.5, units = "in")

# Figure 3 - Hostplant cloglog and Threshold ------------------------------

g3 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_hp_1_t1, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  ggtitle("1959 - 1999") +
  coord_quickmap()



#Plotting hp 1 T2
g4 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_hp_1_t2, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  ggtitle("2000 - 2018") +
  coord_quickmap()

maxent_raw_hp_1 = ggarrange(g3, g4, common.legend = TRUE, legend = "top")

ggsave(plot = maxent_raw_hp_1, filename = "./output/hostplant_1_maxent_raw.png", device = "png")


#Plotting hostplant 2 t1
g5 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_hp_2_t1, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("1959 - 1999") +
  coord_quickmap()



#Plotting hostplant 2 T2
g6 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_hp_2_t2, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.4,0.5,0.5,0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("2000 - 2019") +
  coord_quickmap()

maxent_raw_hp_2 = ggarrange(g5, g6, common.legend = TRUE, legend = "none")
maxent_raw_hp_2

ggsave(plot = maxent_raw_hp_2, filename = "./output/hostplant_2_maxent_raw.png", device = "png")



#Plotting hp 3 t1
g7 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_hp_3_t1, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("1959 - 1999") +
  coord_quickmap() +
  annotate(geom = "text", label = "test")



#Plotting Hp 3 T2
g8 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "#440154FF") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "#440154FF") +
  geom_tile(data=pred_sp_df_hp_3_t2, aes(x=x, y=y, fill=value)) + 
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey50", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  scale_fill_viridis(name = "Probability of Occurence") +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  # ggtitle("2000 - 2019") +
  coord_quickmap()

maxent_raw_hp_3 = ggarrange(g7, g8, common.legend = TRUE, legend = "none")
maxent_raw_hp_3 

ggsave(plot = maxent_raw_hp_3, filename = "./output/hostplant_3_maxent_raw.png", device = "png")

#Big plot
library(gridExtra)

text.1 = ggparagraph(text = "Zanthoxylum americanum", face = "italic", size = 12)
text.2 = ggparagraph(text = "Zanthoxylum clava-herculis", face = "italic", size = 12)
text.3 = ggparagraph(text = "Ptelea trifoliata", face = "italic", size = 12)
fig_3_a = ggarrange(maxent_raw_hp_1, 
                    maxent_raw_hp_2,
                    maxent_raw_hp_3,
                    nrow = 3, ncol = 1, 
                    legend = "top", 
                    heights = c(1.4,1,1))
#Threholds Panel B

g13 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "grey10") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "grey10") +
  geom_tile(data = hp_thresholds_df_t1, aes(x = x, y = y, fill = name), alpha = 0.6) +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  # geom_point(data = hostplant %>%
  #              filter(time_frame == "T1") %>%
  #              filter(str_detect(name, "Zanthoxylum americanum")),
  #            aes(x = longitude, y = latitude), color ="#F8766D",  alpha = 0.5, shape = 3) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  scale_fill_discrete(name = "Species", 
                      breaks = c("Zanthoxylum americanum", 
                                 "Zanthoxylum clava-herculis", 
                                 "Ptelea trifoliata")) +
  ggtitle("1959 - 1999") +
  coord_quickmap() +
  scale_color_discrete(guide = FALSE)

#T2
g14 = ggplot() +  
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color=NA, size=0.25, fill = "grey10") +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = NA, size = 0.25, fill = "grey10") +
  geom_tile(data = hp_thresholds_df_t2, aes(x = x, y = y, fill = name), alpha = 0.6) +
  geom_polygon(data=simple_map_US, aes(x=long, y=lat, group=group), 
               color="grey75", size=0.25, fill = NA) +
  geom_polygon(data = simple_map_can, aes(x = long, y = lat, group = group), color = "grey50", size = 0.25, fill = NA) +
  # geom_point(data = swallowtail_t2, aes(x = longitude, y = latitude), alpha = 0.2, color = "yellow", shape = 3) +
  geom_polygon(data = lakes, aes(x = long, y = lat, group = group), fill = "white", size = 0.25) +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"),
        plot.title = element_text(hjust = 0.5, size = 24),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines")) +
  theme_nothing(legend = TRUE) +
  scale_fill_discrete(name = "Species", 
                      breaks = c("Zanthoxylum americanum", 
                                 "Zanthoxylum clava-herculis", 
                                 "Ptelea trifoliata")) +
  ggtitle("2000 - 2018") +
  coord_quickmap() +
  scale_color_discrete(guide = FALSE)

fig_3_b = ggarrange(g13, g14, common.legend = TRUE, legend = "bottom")

fig_3 = ggarrange(fig_3_a, fig_3_b, ncol = 1, nrow = 2)
ggsave(plot = fig_3, filename = "./output/fig_3.png", device = "png",
       width = 10, height = 16, units = "in")


# Figure 4 - Density Estimates - Whole range --------------------------------

#plotting st
g9 = ggplot(threshold_df_st, aes(x = y, fill = timeframe)) +
  geom_density(alpha = 0.8) +
  theme_classic() +
  labs(x = "Latitude", y = "Kernel Density Estimate") +
  scale_fill_discrete(name = "Time Frame", labels = c("1959-1999", "2000-2018")) +
  xlim(c(25,50)) +
  labs(title = expression(italic("P. cresphontes")))

#plotting hp 1
g10 = ggplot(threshold_df_hp_1, aes(x = y, fill = timeframe)) +
  geom_density(alpha = 0.8) +
  theme_classic() +
  labs(x = "Latitude", y = "Kernel Density Estimate") +
  scale_fill_discrete(name = "Time Frame", labels = c("1959-1999", "2000-2018")) +
  xlim(c(25,50)) +
  labs(title = expression(italic("Z. americanum")))

g11 = ggplot(threshold_df_hp_2, aes(x = y, fill = timeframe)) +
  geom_density(alpha = 0.8) +
  theme_classic() +
  labs(x = "Latitude", y = "Kernel Density Estimate") +
  scale_fill_discrete(name = "Time Frame", labels = c("1959-1999", "2000-2018")) +
  xlim(c(25,50)) +
  labs(title = expression(italic("Z. clava-herculis")))

g12 = ggplot(threshold_df_hp_3, aes(x = y, fill = timeframe)) +
  geom_density(alpha = 0.8) +
  theme_classic() +
  labs(x = "Latitude", y = "Kernel Density Estimate") +
  scale_fill_discrete(name = "Time Frame", labels = c("1959-1999", "2000-2018")) +
  xlim(c(25,50)) +
  labs(title = expression(italic("P. trifoliata")))


histograms_plot = ggarrange(g9,g10,g11,g12, common.legend = TRUE, nrow = 4)
histograms_plot

ggsave(plot = histograms_plot, filename = "./output/fig_4.png", device = "png",
       height = 8.5, width = 11, units = "in")


# Figure 5 - Northern Range Viz -------------------------------------------
#finding the threshold for presence/absence for each model
st_t1_threshold = threshold(ev_st_t1, 'spec_sens')
st_t2_threshold = threshold(ev_st_t2, 'spec_sens')
hp_1_t1_threshold = threshold(ev_hp_1_t1, 'spec_sens')
hp_1_t2_threshold = threshold(ev_hp_1_t2, 'spec_sens')
hp_2_t1_threshold = threshold(ev_hp_2_t1, 'spec_sens')
hp_2_t2_threshold = threshold(ev_hp_2_t2, 'spec_sens')
hp_3_t1_threshold = threshold(ev_hp_3_t1, 'spec_sens')
hp_3_t2_threshold = threshold(ev_hp_3_t2, 'spec_sens')

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
                      labels = c("T1 - 1959-1999", "T2 - 2000-2018")) + 
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

median(st_max_y_t2, na.rm = TRUE)
sd(st_max_y_t2, na.rm = TRUE)

median(st_max_y_t1, na.rm = TRUE)
sd(st_max_y_t1, na.rm = TRUE)

test = t.test(st_max_y_t1, st_max_y_t2, paired = TRUE)

#Can we see if there is a difference between the the butterfly and host plant 1?
n_limit_st_hp1_t2 = pred_sp_df_st_t2 %>%
  bind_rows(pred_sp_df_hp_1_t2, .id = "species") %>%
  mutate(occ = ifelse(value >= st_t2_threshold & species == 1, 1,
                      ifelse(value >= hp_1_t2_threshold & species == 2, 1, 0))) %>%
  group_by(x, species) %>%
  summarize(max_y = max(y[occ == 1])) %>%
  mutate(max_y = ifelse(max_y == -Inf, NA, max_y)) %>%
  ungroup()


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
                    labels = c(expression(italic("P. cresphontes")), expression(italic("Z. americanum"))),
                    values = c("#E69F00", "#56B4E9")) +
  xlim(c(38, 50)) +
  ylim(c(0, 0.90)) +
  ggtitle("T2 (2000-2018)")

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
                    labels = c(expression(italic("P. cresphontes")), expression(italic("Z. americanum"))),
                    values = c("#E69F00", "#56B4E9")) +
  xlim(c(38, 50)) +
  ylim(c(0, 0.90)) +
  ggtitle("T1 (1959-1999)")

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
                      labels = c("T1 - 1959-1999", "T2 - 2000-2018")) + 
  labs(title = expression(italic("Z. americanum"))) +
xlim(c(38, 50)) +
  ylim(c(0, 0.90))


fig_5 = ggarrange(fig_5_a, fig_5_b, fig_5_c, fig_5_d, labels = "AUTO")
ggsave(plot = fig_5, filename = "./output/fig_5.png", device = "png", 
       width = 11, height = 8.5, units = "in")


# Figure 6 - Environmental Feature Contribution ---------------------------

#Environmental Variable Importance

#Swallowtail time-frame 1
df = var.importance(mx_best_st_t1)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_1 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('P. cresphontes')*' 1959-1999')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Swallowtail time-frame 2
df = var.importance(mx_best_st_t2)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_2 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('P. cresphontes')*' 2000-2018')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#hostplant 1 time-frame 1
df = var.importance(mx_best_hp_1_t1)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_3 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('Z. americanum')*' 1959-1999')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#hostplant 1 time-frame 2
df = var.importance(mx_best_hp_1_t2)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_4 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('Z. americanum')*' 2000-2018')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#hostplant 2 time-frame 1 
df = var.importance(mx_best_hp_2_t1)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_5 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('Z. clava-herculis')*' 1959-1999')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#hostplant 2 time-frame 2
df = var.importance(mx_best_hp_2_t2)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_6 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('Z. clava-herculis')*' 2000-2018')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#hostplant 3 time-frame 1 
df = var.importance(mx_best_hp_3_t1)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_7 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('P. trifoliata')*' 1959-1999')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#hostplant 3 time-frame 2
df = var.importance(mx_best_hp_3_t2)
df$variable = factor(df$variable, levels = c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7",
                                             "Bio8", "Bio9", "Bio10", "Bio11", "Bio12", "Bio13", 
                                             "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19"))
env_plot_8 = ggplot(df, aes(x = variable, y = percent.contribution)) +
  geom_col() +
  theme_classic() +
  labs(x = "Environmental Variable", 
       y = "Percent Contribution") +
  labs(title = expression(''*italic('P. trifoliata')*' 2000-2018')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



env_plot = ggarrange(env_plot_1, env_plot_2, env_plot_3, env_plot_4,
                     env_plot_5, env_plot_6, env_plot_7, env_plot_8,
                     common.legend = TRUE, 
                     ncol = 2, nrow = 4)
env_plot

ggsave(plot = env_plot, filename = "./output/fig_6.png", device = "png")

