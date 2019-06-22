#Script for getting appropriate swallowtail data
#Keaton Wilson
#keatonwilson@me.com
#2019-02-15
#
#
#Updated on 2019-04-08 to include additional host plants

#Loading appropriate packages
library(tidyverse)
library(spocc)
library(mapr)
library(ggmap)
library(scrubr)
library(lubridate)
library(readxl)
library(daymetr)
library(tidync)
library(FedData)
library(ncdf4)
library(dismo)

#register google api for mapping stuff
register_google(key = "AIzaSyDyAqUc4o9p_DOBSF_JOXH5c_JXPqoU4Yw")

#Let's query inat and GBIF for swallowtail records
swallowtail = occ(query = "Papilio cresphontes*", from = c("gbif","inat"),
                  has_coords = TRUE, limit = 20000)


#Filtering out anything from inat that isn't research grade
swallowtail$inat$data$`Papilio_cresphontes*` = swallowtail$inat$data$`Papilio_cresphontes*` %>%
  filter(quality_grade == "research")
  
#Aggregating records
swallowtail_df = occ2df(swallowtail)

#Initial mapping
swallowtail_df %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>%
  map_ggmap()

#Lots of naming errors (let's do some filtering)
swallowtail_df = swallowtail_df %>%
  filter(str_detect(name, 'Papilio cresphontes'))

#let's check out time periods
swallowtail_df %>%
  mutate(year = year(date)) %>%
  filter(year >= 1960) %>%
  group_by(year) %>%
  summarize(n()) 

#Let's pull in Kent's data from ebutterfly, maine atlas, maritime atlas and MA butterfly club
ebutterfly = read_xlsx(path = "./data/e_butterfly.xlsx")
maine = read_xlsx(path = "./data/maine_butterfly_atlas.xlsx")
maritime = read_xlsx(path = "./data/maritime_atlas.xlsx")
ma_club = read_xlsx(path = "./data/ma_butterfly_club.xlsx")
bamona = read_csv("./data/bamona_data.csv")

#cleaning and merging
ebutterfly = ebutterfly %>%
  dplyr::select(OccuranceID, 'Date Observed', Latitude, Longitude) %>%
  dplyr::select(Latitude, Longitude, Date = 'Date Observed') %>%
  mutate(date = as.Date(Date)) %>%
  dplyr::select(-Date)

maine = maine %>%
  dplyr::select(Latitude, Longitude, Year, Month, Day) %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  mutate(date = date(paste(Year, Month, Day, sep = "-"))) %>%
  dplyr::select(-c(Year, Month, Day))

maritime = maritime %>%
  dplyr::select(Latitude, Longitude, Year, Month, Day) %>%
  filter(Day != "XX") %>%
  mutate(date = date(paste(Year, Month, Day, sep = "-"))) %>%
  dplyr::select(-c(Year, Month, Day))

ma_club = ma_club %>%
  dplyr::select(Latitude, Longitude, Date) %>%
  mutate(date = date(Date)) %>%
  dplyr::select(-Date)

bamona = bamona %>%
  dplyr::select(Latitude = `Lat/Long`, Longitude = Longitude, date =`Observation Date`) %>%
  mutate(date = as.Date(date, "%m/%d/%Y"))

swallowtail_df = swallowtail_df %>%
  select(Latitude = latitude, Longitude = longitude, date) %>%
  mutate(Latitude = as.numeric(Latitude), 
         Longitude = as.numeric(Longitude))

#binding together
swallowtail_master = bind_rows("inat_gbif" = swallowtail_df, 
                               "ebutterfly" = ebutterfly,
                               "maine" = maine, 
                               "maritime" = maritime,
                               "ma_club" = ma_club,
                               "bamona" = bamona,
                               .id = "data_source")


#So, lots more records as time has progressed - seems like probably an inat phenomenom. In the original manuscript, only went up to 2010. Would be nice to include more recent data

#Can we just build a simple faceted plot by decade (doesn't include new data from Kent)

#removing duplicates, filtering older data and restricting data to the chunk of the NE we're interested in
swallowtail_master = swallowtail_master %>%
  filter(year(date) > 1959) %>%
  distinct() %>%
  filter(Latitude < 51 & Latitude > 22) %>%
  filter(Longitude < -52 & Longitude > -94)

northeast = get_map("New York", zoom = 3)
ggmap(northeast) +
  geom_point(data = swallowtail_master, aes(x = Longitude, y = Latitude))


#-----------------------------------------------------------------------------------
#HOST PLANT DATA

#Pull in host plant data
host_plant_1 = occ(query = "Zanthoxylum americanum", has_coords = TRUE, from = c("inat", "gbif"), limit = 10000)
host_plant_2 = occ(query = "Zanthoxylum clava-herculis", has_coords = TRUE, from = c("inat", "gbif"), limit = 10000)
host_plant_3 = occ(query = "Ptelea trifoliata", has_coords = TRUE, from = c("inat", "gbif"), limit = 10000)

#Filtering out anything from inat that isn't research grade
host_plant_1$inat$data$`Zanthoxylum` = host_plant_1$inat$data$`Zanthoxylum` %>%
  filter(quality_grade == "research")
host_plant_2$inat$data$`Zanthoxylum` = host_plant_2$inat$data$`Zanthoxylum` %>%
  filter(quality_grade == "research")
host_plant_3$inat$data$`Ptelea` = host_plant_3$inat$data$`Ptelea` %>%
  filter(quality_grade == "research")

#Aggregating records
host_plant_df_1 = occ2df(host_plant_1)
host_plant_df_2 = occ2df(host_plant_2)
host_plant_df_3 = occ2df(host_plant_3)

#filtering names
host_plant_df_1 = host_plant_df_1 %>%
  filter(str_detect(name, "Zanthoxylum americanum")) %>%
  mutate(Longitude = as.numeric(longitude), 
         Latitude = as.numeric(latitude))

host_plant_df_2 = host_plant_df_2 %>%
  filter(str_detect(name,"Zanthoxylum clava-herculis")) %>%
  mutate(Longitude = as.numeric(longitude), 
         Latitude = as.numeric(latitude))

host_plant_df_3 = host_plant_df_3 %>%
  filter(str_detect(name,"Ptelea trifoliata")) %>%
  mutate(Longitude = as.numeric(longitude), 
         Latitude = as.numeric(latitude))

#Initial mapping
host_plant_df_1 %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>%
  map_ggmap()



#combining and filtering
hostplant_master = host_plant_df_1 %>%
  bind_rows(host_plant_df_2) %>%
  bind_rows(host_plant_df_3) %>%
  filter(year(date) > 1959) %>%
  distinct() %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>%
  filter(latitude < 50 & latitude > 22) %>%
  filter(longitude < -50 & longitude > -94)

hostplant_master %>%
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>%
  map_ggmap()

#quick map
ggmap(northeast) +
  geom_point(data = hostplant_master, aes(x = Longitude, y = Latitude, color = name))

#Final steps 

# Importing Data and Cleaning ---------------------------------------------
#Pulling in data that has been cleaned in previous script

#Generating some features that will be useful later on
swallowtail_master = swallowtail_master %>%
  mutate(year = year(date), 
         time_frame = ifelse(year >= 2000, "T2", "T1")) %>%
  rename(latitude = Latitude, longitude = Longitude) %>%
  select(-data_source)

hostplant_master = hostplant_master %>%
  mutate(year = year(date), 
         time_frame = ifelse(year >= 2000, "T2", "T1"),
         longitude = as.numeric(longitude), 
         latitude = as.numeric(latitude)) %>%
  select(-Longitude, -Latitude, -key)

#Filtering the data to include stuff east of texas (94º), and in the US, Canada. Should be done from data_import_clean script, but good to double check
lon_min = -94
lon_max = -65
lat_min = 25
lat_max = 55

swallowtail_master = swallowtail_master %>%
  filter(latitude >= lat_min & latitude <= lat_max) %>%
  filter(longitude >= lon_min & longitude <= lon_max)

hostplant_master = hostplant_master %>%
  filter(latitude >= lat_min & latitude <= lat_max) %>%
  filter(longitude >= lon_min & longitude <= lon_max)

#mapping to check
st = get_map("North Carolina", zoom = 4, maptype = "toner-background")
ggmap(st, maptype = "toner-background", extent = "panel") +
  geom_point(data = swallowtail_master, aes(x = longitude, y = latitude, color = time_frame), alpha = 0.5) +
  scale_color_discrete(name = "Time Frame", labels = c("Pre-2000", "Post-2000")) +
  labs(x = "Longitude (º) ", y = "Latitude (º)") +
  facet_wrap(~ time_frame)


# Importing Bioclim Data and Cropping -------------------------------------
#We can do this from the dismo package - Interesting point here, this is different from original methods. These climate data are representative of "current" conditions - averaged between 1970 and 2000 (https://www.researchgate.net/publication/316999789_WorldClim_2_New_1-km_spatial_resolution_climate_surfaces_for_global_land_areas). 

bioclim.data <- raster::getData(name = "worldclim",
                                var = "bio",
                                res = 2.5,
                                path = "./data/")

# Determine geographic extent of our data
max_lat_swallowtail <- ceiling(max(swallowtail_master$latitude))
min_lat_swallowtail <- floor(min(swallowtail_master$latitude))
max_lon_swallowtail <- ceiling(max(swallowtail_master$longitude))
min_lon_swallowtail <- floor(min(swallowtail_master$longitude))
geographic.extent <- extent(x = c(min_lon_swallowtail, max_lon_swallowtail, min_lat_swallowtail, max_lat_swallowtail))

# Crop bioclim data to geographic extent of swallowtails
bioclim.data <- crop(x = bioclim.data, y = geographic.extent)

#IMporting CRU files as raster stacks
precip_raster = raster::stack("./data/cru_ts4.02.1901.2017.pre.dat.nc")
tmin_raster = raster::stack("./data/cru_ts4.02.1901.2017.tmn.dat.nc")
tmax_raster = raster::stack("./data/cru_ts4.02.1901.2017.tmx.dat.nc")

#This doesn't work because it only does one year at a time
bio_vars_test = biovars(prec = precip_raster,
                        tmin = tmin_raster,
                        tmax = tmax_raster)

#So, first step is to split all of our rasters into the two groups we want (1959-1999 & 2000-2017)

#1901-01-01 - 1959-01-01
(1959-1901)*12 - 1 #number of months. 

#So the starting point for T1 is month 695
(1999-1901)*12 - 1 #Last month in T1

#QC 
(1176-696)/12 #Divisible by 12, we're good

precip_raster_t1 = precip_raster[[696:1175]]
tmin_raster_t1 = tmin_raster[[696:1175]]
tmax_raster_t1 = tmax_raster[[696:1175]]

#Starting point for T2
(2000-1901)

#QC
(1404-1188)/12

precip_raster_t2 = precip_raster[[1189:1404]]
tmin_raster_t2 = tmin_raster[[1189:1404]]
tmax_raster_t2 = tmax_raster[[1189:1404]]

#So essentially we now need to divide each set of data into 1 year's worth of stuff, run through biovars, compute our variables, and then average all of that over each time period. 
#
#Test on one iteration

precip_test = precip_raster_t1[[1:12]]
tmin_test = tmin_raster_t1[[1:12]]
tmax_test = tmax_raster_t1[[1:12]]

biovars_test = biovars(prec = precip_test,
                       tmin = tmin_test,
                       tmax = tmax_test)

#Great. This works - but it outputs a rasterbrick, not sure how that will work with stuff feeding into a raster stack below. 
biovar_list = list()
length = dim(precip_raster_t1)[3]/12
seq = 1:12

for(i in 1:length) {
  precip_sub = precip_raster_t1[[seq]]
  tmin_sub = tmin_raster_t1[[seq]]
  tmax_sub = tmax_raster_t1[[seq]]
  
  biovar_list[[i]] = biovars(prec = precip_sub,
                            tmin = tmin_sub,
                            tmax = tmax_sub)
  seq = seq + 12
  print(seq)
}

#Workflow 
#pull out a list of raster layers for each bioclim variable
#turn those layers into a rasterStack
#Compute average
#Recombine

biovar_avg_combined_t1 = raster::brick(nrows = 360, ncols = 720)
for(i in 1:19) {
  biovar_sublist = lapply(biovar_list, '[[', i) #pulls out each bioclim variable iteratively
  biovar_substack = stack(biovar_sublist) #combines all years into a raster stack
  biovar_avg = calc(biovar_substack, fun = mean) #Calculates the average for each var
  biovar_avg_combined_t1[[i]] = biovar_avg
}

#comparing structure of calculated bioclim versus what the function outputs
bioclim.data
biovar_avg_combined
plot(biovar_avg_combined[[1]])

#T2
biovar_list = list()
length = dim(precip_raster_t2)[3]/12
seq = 1:12

for(i in 1:length) {
  precip_sub = precip_raster_t2[[seq]]
  tmin_sub = tmin_raster_t2[[seq]]
  tmax_sub = tmax_raster_t2[[seq]]
  
  biovar_list[[i]] = biovars(prec = precip_sub,
                             tmin = tmin_sub,
                             tmax = tmax_sub)
  seq = seq + 12
  print(seq)
}

biovar_avg_combined_t2 = raster::brick(nrows = 360, ncols = 720)
for(i in 1:19) {
  biovar_sublist = lapply(biovar_list, '[[', i) #pulls out each bioclim variable iteratively
  biovar_substack = stack(biovar_sublist) #combines all years into a raster stack
  biovar_avg = calc(biovar_substack, fun = mean) #Calculates the average for each var
  biovar_avg_combined_t2[[i]] = biovar_avg #binding each averaged layer back into a brick
}

#comparing structure of calculated bioclim versus what the function outputs
bioclim.data
biovar_avg_combined_t2
plot(biovar_avg_combined_t2[[1]])
plot(biovar_avg_combined_t1[[1]])

#Saving these objects
writeRaster(biovar_avg_combined_t1, "./data/biovar_avg_combined_t1")
writeRaster(biovar_avg_combined_t2, "./data/biovar_avg_combined_t2")

#writing bioclim data
saveRDS(bioclim.data, "./data/bioclim.rds")

#Writing hostplant records
write.csv(hostplant_master, "./data/hostplant_data.csv")

#Writing butterfly records
write.csv(swallowtail_master, "./data/swallowtail_data.csv")

hostplant_master %>%
  filter(str_detect(name, "Zanthoxylum clava-herculis")) 
