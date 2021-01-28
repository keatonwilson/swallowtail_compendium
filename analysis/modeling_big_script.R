#blockCV and Modeling Script for Swallowtail and hostplant data
#Keaton Wilson
#keatonwilson@me.com
#2019-06-19

#libraries
library(blockCV)
library(tidyverse)
library(raster)
library(maxnet)
library(dismo)
library(ENMeval)

# Data Preparation --------------------------------------------------------
#setting seed for reproducibility down the line
set.seed(223)

#importing swallowtail, hostplant and environmental data
#butterfly
swallowtail = read_csv("./data/raw_data/swallowtail_data.csv")
swallowtail = swallowtail[,-1] %>%
  dplyr::select(longitude, latitude, date, year, time_frame)

hp = read_csv("./data/raw_data/hostplant_data.csv") %>%
  dplyr::select(-1, longitude, latitude, date, year, time_frame)

hostplant_1 = hp %>%
  filter(str_detect(name, "Zanthoxylum americanum"))

hostplant_2 = hp %>%
  filter(str_detect(name, "Zanthoxylum clava-herculis"))

hostplant_3 = hp %>%
  filter(str_detect(name, "Ptelea trifoliata"))

bv_t1 = raster::brick("./data/raw_data/biovar_avg_t1.gri")
bv_t2 = raster::brick("./data/raw_data/biovar_avg_t2.gri")

#renaming
names_seq = paste("Bio",seq(1:19), sep = "")
names(bv_t1) = names_seq
names(bv_t2) = names_seq

# Determine geographic extent of our data - specific extents for each species
# Function to get extents

get_extent = function(species_data = NULL) {
  max_lat = ceiling(max(species_data$latitude))
  min_lat = floor(min(species_data$latitude))
  max_lon = ceiling(max(species_data$longitude))
  min_lon = floor(min(species_data$longitude))
  geographic_extent = extent(x = c(min_lon, 
                                   max_lon, 
                                   min_lat, 
                                   max_lat))
  return(geographic_extent)
}

geographic_extent_st = get_extent(swallowtail)
geographic_extent_hp_1 = get_extent(hostplant_1)
geographic_extent_hp_2 = get_extent(hostplant_2)
geographic_extent_hp_3 = get_extent(hostplant_3)

# Crop bioclim data to geographic extent of swallowtails
bv_t1_st <- crop(x = bv_t1, y = geographic_extent_st)
bv_t2_st <- crop(x = bv_t2, y = geographic_extent_st)

# Crop bioclim data to geographic extent of swallowtails
bv_t1_hp_1 <- crop(x = bv_t1, y = geographic_extent_hp_1)
bv_t2_hp_1 <- crop(x = bv_t2, y = geographic_extent_hp_1)

# Crop bioclim data to geographic extent of swallowtails
bv_t1_hp_2 <- crop(x = bv_t1, y = geographic_extent_hp_2)
bv_t2_hp_2 <- crop(x = bv_t2, y = geographic_extent_hp_2)

# Crop bioclim data to geographic extent of swallowtails
bv_t1_hp_3 <- crop(x = bv_t1, y = geographic_extent_hp_3)
bv_t2_hp_3 <- crop(x = bv_t2, y = geographic_extent_hp_3)


#using the custom prep_data function to ready the data for blockCV
source("./R/prep_data_func.R")

swallowtail_prepared_data = prep_data(swallowtail, bv_t1_st, bv_t2_st)
hostplant1_prepared_data = prep_data(hostplant_1, bv_t1_hp_1, bv_t2_hp_1)
hostplant2_prepared_data = prep_data(hostplant_2, bv_t1_hp_2, bv_t2_hp_2)
hostplant3_prepared_data = prep_data(hostplant_3, bv_t1_hp_3, bv_t2_hp_3)

#Merging all the mini-lists into a large list of dataframes
prepared_data_master = c(swallowtail_prepared_data, 
                         hostplant1_prepared_data, 
                         hostplant2_prepared_data, 
                         hostplant3_prepared_data)

names(prepared_data_master) = c("st_t1", "st_t2", "hp1_t1", "hp1_t2",
                                "hp2_t1", "hp2_t2", "hp3_t1", "hp3_t2")


# blockCV Train-Test Split for all 4 models ------------------------------------------------
#Running spatialBlock function over every dataframe in the list
#
block_list = list()
for (i in 1:length(prepared_data_master)) {
  if (str_detect(names(prepared_data_master[i]), "t1") == TRUE) {
    raster = bv_t1 
  } else {
    raster = bv_t2
  }
  
  block_list[[i]] = spatialBlock(speciesData = prepared_data_master[[i]],
                               species = "Species",
                               rasterLayer = raster,
                               theRange = 400000,
                               k = 5, 
                               selection = "random", 
                               iteration = 250, 
                               biomod2Format = TRUE, 
                               xOffset = 0, 
                               yOffset = 0, 
                               progress = T)
}

#Saving Spatial CV splits - these actually take a surprising amount of time to run, and are necessary building blocks for threshold maps in the figures script
saveRDS(block_list, "./data/block_list.rds")

#Getting dataframes to feed into the model (dropping NAs)
#Swallowtail

model_data_list = list()
for (i in 1:length(prepared_data_master)) {
  if (str_detect(names(prepared_data_master[i]), "t1") == TRUE) {
    raster = bv_t1 
  } else {
    raster = bv_t2
  }
  
  model_data_list[[i]] = raster::extract(raster, prepared_data_master[[i]][,-3], df = TRUE) %>%
    bind_cols(as.data.frame(prepared_data_master[[i]])) %>%
    drop_na() %>%
    dplyr::select(-ID, Species, longitude, latitude, Bio1:Bio19)
}

#vectors of presence-background
pb_list = lapply(prepared_data_master,  function(x) '['(x, 3))

#folds for each model
fold_list = lapply(block_list, function(x) '[['(x, 1))

#Writing a function that unlists and combines all training and test indices for each species-time period combination
extract_index = function(list_of_folds = NULL) {
for(k in 1:length(list_of_folds)){
  train_index <- unlist(list_of_folds[[k]][1]) # extract the training set indices
  test_index <- unlist(list_of_folds[[k]][2])# extract the test set indices
}
  mini_list = list(train_index, test_index)
  mini_list
}

#Applying the function to the list of folds
train_test_index_list = lapply(fold_list, extract_index)

train_test_data_list = list()
for (i in 1:length(model_data_list)) {
    train_index = train_test_index_list[[i]][[1]]
    test_index = train_test_index_list[[i]][[2]]
  
    train_data = model_data_list[[i]][train_index,]
    test_data = model_data_list[[i]][test_index,]
    
    mini_list = list(train_data, test_data)
    train_test_data_list[[i]] = mini_list
}

str(train_test_data_list)

#Adding on T1, T2 designations
for (i in 1:length(train_test_data_list)) {
  for (j in 1:length(train_test_data_list[[i]])) {
    if (j == 1) {
      train_test_data_list[[i]][[j]]$time_period = 1
    } else {
      train_test_data_list[[i]][[j]]$time_period = 2
    }
  }
}

# Modeling ----------------------------------------------------------------
if (Sys.getenv("JAVA_HOME")!="")
  Sys.setenv(JAVA_HOME="")
library(rJava)

model_func = function(data = NULL) {
  data_occ = data[[1]] %>%  #Generating occurence lat long
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude)
  
  data_bg = data[[1]] %>% #Generating background lat long
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude)
  
  if (data[[1]]$time_period[1] == 1) { #Setting the appropriate environmental layer for the time period
    env_data = bv_t1
  } else {
    env_data = bv_t2
  }
  
  #Running the model
  eval = ENMevaluate(occ = data_occ, 
                     bg.coords = data_bg,
                     env = env_data,
                     method = 'randomkfold', 
                     kfolds = 5, 
                     algorithm = 'maxent.jar', 
                     fc = c("L", "LQ", "H", "LQH", "QH", "LH"))
  eval
}

start = Sys.time()
#Running the model function over the list of data
big_model_list = lapply(train_test_data_list, model_func)

#Saving this bad boy
saveRDS(big_model_list, "./data/big_model_list.rds")

end = Sys.time()
# Model Evaluation --------------------------------------------------------

#Function to build set of evaluation plots - just plug in the appropriate eval model object from above

eval_plots = function(eval_object = NULL) {
  par(mfrow=c(2,3))
  eval.plot(eval_object@results)
  eval.plot(eval_object@results, 'avg.test.AUC', legend = F)
  eval.plot(eval_object@results, 'avg.diff.AUC', legend = F)
  eval.plot(eval_object@results, 'avg.test.or10pct', legend = F)
  eval.plot(eval_object@results, 'avg.test.orMTP', legend = F)
  plot(eval_object@results$avg.test.AUC, eval_object@results$delta.AICc, bg=factor(eval_object@results$features), pch=21, cex= eval_object@results$rm/2, xlab = "avg.test.AUC", ylab = 'delta.AICc', cex.lab = 1.5)
  legend("topright", legend=unique(eval_object@results$features), pt.bg=factor(eval_object@results$features), pch=21)
  mtext("Circle size proportional to regularization multiplier value", cex = 0.6)
}


#Evaluation plots
#Feed in a list of models and it outputs and saves the evaluation plots for each one
for (i in 1:length(big_model_list)) {
  name = paste("eval_plot", i, sep = "_")
  png(filename = paste0("./output/", name, ".png"), 
      width = 1080, height = 720)
  plot = eval_plots(big_model_list[[i]])
  dev.off()
}

#Picking the best model based on highest AUC for each set
#Pulling out indices of the "best" model based on AUC scores - if there are two models that are equal, it pulls the first.

model_selection_index_list = list()

for (i in 1:length(big_model_list)) {
  model_selection_index_list[[i]] = as.numeric(row.names(big_model_list[[i]]@results[which(big_model_list[[i]]@results$avg.test.AUC== max(big_model_list[[i]]@results$avg.test.AUC)),]))[1]
}

best_model_list = list()
for (i in 1:length(big_model_list)) {
  index = model_selection_index_list[[i]]
  model = big_model_list[[i]]@models[[index]]
  best_model_list[[i]] = model
}

#combining models and data into a master list

master_list = list()
for (i in 1:length(train_test_data_list)) {
  master_list[[i]] = append(train_test_data_list[[i]], best_model_list[[i]])
}
#Generating evaluate objects on test data
# 
# evaluate_models = function(master_list_sub = NULL) {
#   test_data_occ = master_list_sub[[2]] %>%
#     filter(Species == 1) %>%
#     dplyr::select(longitude, latitude)
#   
#   bg_data = master_list_sub[[2]] %>%
#     filter(Species == 0) %>%
#     dplyr::select(longitude, latitude)
#   
#   if (master_list_sub[[2]]$time_period[1] == 1) {
#     env_data = bv_t1
#   } else {
#     env_data = bv_t2
#   }
#   
#   model_sub = master_list_sub[[3]]
#   
#   ev = evaluate(test_data_occ, a = bg_data, model = model_sub, x = env_data)
#   ev
# }


names(master_list) = c("st_t1", 
                           "st_t2", 
                           "hp_1_t1", 
                           "hp_1_t2", 
                           "hp_2_t1", 
                           "hp_2_t2", 
                           "hp_3_t1", 
                           "hp_3_t2")

names_list = c("bv_t1_st", 
               "bv_t2_st", 
               "bv_t1_hp_1", 
               "bv_t2_hp_1", 
               "bv_t1_hp_2", 
               "bv_t2_hp_2", 
               "bv_t1_hp_3", 
               "bv_t2_hp_3")


evaluations = list()
for(i in 1:length(master_list)) {
  test_data_occ = master_list[[i]][[2]] %>%
    filter(Species == 1) %>%
    dplyr::select(longitude, latitude)
  
  bg_data = master_list[[i]][[2]] %>%
    filter(Species == 0) %>%
    dplyr::select(longitude, latitude)
  
  env_data = get(names_list[i])
  
  model_sub = master_list[[i]][[3]]
  
  ev = evaluate(test_data_occ, a = bg_data, model = model_sub, x = env_data)
  evaluations[[i]] = ev
  print(paste("Finished", names_list[i]))
}



#Saving evaluations
saveRDS(evaluations, file = "./data/evaluations.rds")

# Selecting Final Models and Running on All Data --------------------------
#Let's build final models

full_model = function(models = NULL, full_data = NULL, best_model_index = NULL, name = NULL, env_data_t1 = NULL, env_data_t2 = NULL) {
  auc_mod = models@results[best_model_index,]
  FC_best = as.character(auc_mod$features[1])
  rm_best = auc_mod$rm
  maxent.args = ENMeval::make.args(RMvalues = rm_best, fc = FC_best)
  
  if (full_data$time_frame[1] == "T1") {
    env_data = env_data_t1
  } else {
    env_data = env_data_t2
  }
  
  full_mod = maxent(env_data, as.matrix(full_data[,1:2]), args = maxent.args[[1]])
  saveRDS(full_mod, paste0("./data/", name, ".rds" ))
  
  full_mod
}

# As of 2020-12-21, this includes cropped versions of env vars
#Creating each master model
swallowtail_t1 = full_model(models = big_model_list[[1]], best_model_index = model_selection_index_list[[1]], 
           full_data = swallowtail %>% filter(time_frame == "T1"), name = "swallowtail_t1", 
           env_data_t1 = bv_t1_st, env_data_t2 = bv_t2_st)

swallowtail_t2 = full_model(models = big_model_list[[2]], best_model_index = model_selection_index_list[[2]], 
                            full_data = swallowtail %>% filter(time_frame == "T2"), name = "swallowtail_t2", 
                            env_data_t1 = bv_t1_st, env_data_t2 = bv_t2_st)

hostplant_1_t1 = full_model(models = big_model_list[[3]], best_model_index = model_selection_index_list[[3]], 
                            full_data = hostplant_1[,-1] %>% filter(time_frame == "T1"), name = "hostplant_1_t1", 
                            env_data_t1 = bv_t1_hp_1, env_data_t2 = bv_t2_hp_1)

hostplant_1_t2 = full_model(models = big_model_list[[4]], best_model_index = model_selection_index_list[[4]], 
                            full_data = hostplant_1[,-1] %>% filter(time_frame == "T2"), name = "hostplant_1_t2", 
                            env_data_t1 = bv_t1_hp_1, env_data_t2 = bv_t2_hp_1)

hostplant_2_t1 = full_model(models = big_model_list[[5]], best_model_index = model_selection_index_list[[5]], 
                            full_data = hostplant_2[,-1] %>% filter(time_frame == "T1"), name = "hostplant_2_t1", 
                            env_data_t1 = bv_t1_hp_2, env_data_t2 = bv_t2_hp_2)

hostplant_2_t2 = full_model(models = big_model_list[[6]], best_model_index = model_selection_index_list[[6]], 
                            full_data = hostplant_2[,-1] %>% filter(time_frame == "T2"), name = "hostplant_2_t2", 
                            env_data_t1 = bv_t1_hp_2, env_data_t2 = bv_t2_hp_2)

hostplant_3_t1 = full_model(models = big_model_list[[7]], best_model_index = model_selection_index_list[[7]], 
                            full_data = hostplant_3[,-1] %>% filter(time_frame == "T1"), name = "hostplant_3_t1", 
                            env_data_t1 = bv_t1_hp_3, env_data_t2 = bv_t2_hp_3)

hostplant_3_t2 = full_model(models = big_model_list[[8]], best_model_index = model_selection_index_list[[8]], 
                            full_data = hostplant_3[,-1] %>% filter(time_frame == "T2"), name = "hostplant_3_t2", 
                            env_data_t1 = bv_t1_hp_3, env_data_t2 = bv_t2_hp_3)

