#Function to bring in dataframes and prepare data for modeling
#Keaton Wilson
#keatonwilson@me.com
#2019-06-14

library(dismo)
library(raster)
#Need these bioclim data to run successfully
##Environmental Data 

swallowtail = read_csv("./data/raw_data/swallowtail_data.csv")
swallowtail = swallowtail[,-1] %>%
  dplyr::select(longitude, latitude, date, year, time_frame)

bv_t1 = raster::brick("./data/raw_data/biovar_avg_t1.gri")
bv_t2 = raster::brick("./data/raw_data/biovar_avg_t2.gri")

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

# Crop bioclim data to geographic extent of swallowtails
bv_t1 <- crop(x = bv_t1, y = geographic.extent)
bv_t2 <- crop(x = bv_t2, y = geographic.extent)

prep_data = function(df = NULL) {
  #Step 1. Split by time period into two data frames
  df_t1 = df %>% 
    filter(time_frame == "T1")
  df_t2 = df %>%
    filter(time_frame == "T2")
  
  #Step 2. Inspect the dataframe
  glimpse(df_t1)
  glimpse(df_t2)
  
  #Step 3. Generate 10k background points for each one. 
  bg_t1 = dismo::randomPoints(bv_t1, 10000)
  colnames(bg_t1) = c("longitude", "latitude")
  
  bg_t2 = randomPoints(bv_t2, 10000)
  colnames(bg_t2) = c("longitude", "latitude")
  
  #Step. 4 Merging background data and occurence data
  df_comb_t1 = data.frame(df_t1) %>%
    mutate(pb = 1) %>%
    dplyr::select(pb, longitude, latitude) %>%
    bind_rows(data.frame(bg_t1) %>% 
                mutate(pb = 0))  %>%
    mutate(Species = as.integer(pb)) %>%
    dplyr::select(-pb)
  
  df_comb_t2 = data.frame(df_t2) %>%
    mutate(pb = 1) %>%
    dplyr::select(pb, longitude, latitude) %>%
    bind_rows(data.frame(bg_t2) %>% 
                mutate(pb = 0)) %>%
    mutate(Species = as.integer(pb)) %>%
    dplyr::select(-pb)
  
  #Step 5. Changing to a spatial points data frame
  df_sp_t1 = SpatialPointsDataFrame(df_comb_t1[,c("longitude","latitude")], 
                                    df_comb_t1, 
                                    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  df_sp_t2 = SpatialPointsDataFrame(df_comb_t2[,c("longitude","latitude")], 
                                   df_comb_t2, 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  #Converting to a list with the two dataframes
  prepared_data_list = list(df_sp_t1, df_sp_t2)
  prepared_data_list
}

# #Test
# prep_st = prep_data(swallowtail)
