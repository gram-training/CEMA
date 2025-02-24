# 0. load needed libraries
library(caret)
library(rnaturalearth)
library(readxl)
library(spdep)
library(mgcv)
library(broom)
library(tidyverse)
library(glm)
library(matrixStats)
filepath <- "c:/users/giselara/Documents/workshops/cema/"

clsi_factors <- read.csv(paste0(filepath, "guideline_factor.csv"))
clsi_factors <- clsi_factors[clsi_factors$CLSI_2023_adjustment_factor<1,]
tertiary_factors <- read.csv(paste0(filepath, "tertiary_factor.csv")) %>% 
  mutate(tert_group = paste0(tert_group," ", super_region)) 
tertiary_factors$super_region <- NULL

raw_data <- read.csv("C:/Users/giselara/Documents/workshops/cema/raw_data.csv") %>% 
  mutate(resistant = ifelse(interpretation %in% c('resistant'), 1, 0)) %>%
  group_by(country,super_region_name,year,guideline_group,tert_group,hosp_type,
           combo,combo2,species,antibiotic,abx_class) %>% 
  summarise(total_test = n(), resistant = sum(resistant)) %>%
  mutate(prop = resistant / total_test)

###########################
# CLSI adjustment
############################

raw_data <- raw_data %>% left_join(clsi_factors)
raw_data$CLSI_2023_adjustment_factor[is.na(raw_data$CLSI_2023_adjustment_factor)] = 1

raw_data$resistant_adjusted <- raw_data$resistant * raw_data$CLSI_2023_adjustment_factor
raw_data <- raw_data %>% group_by(country, super_region_name, year,
                                  tert_group, hosp_type, combo,combo2, 
                                  species, abx_class) %>% 
  summarise(total_test = sum(total_test), resistant = sum(resistant_adjusted)) 


###########################
# Tertiary adjustments
############################
tertiary <- raw_data[(raw_data$hosp_type == "tertiary") & (raw_data$resistant > 0),] %>% left_join(tertiary_factors)
raw_data <- raw_data[(raw_data$hosp_type == "non-tertiary") | (raw_data$resistant == 0),]

tertiary$adj[is.na(tertiary$adj)] = 0
#raw_data$adj[raw_data$resistant <= 0] = 0
tertiary$prop <- tertiary$resistant / tertiary$total_test
tertiary$prop_adjusted <- exp(log(tertiary$prop) - tertiary$adj)
tertiary$resistant_adjusted <- tertiary$total_test * tertiary$prop_adjusted

summary(tertiary[,c('prop','prop_adjusted')])
tertiary$prop <- tertiary$prop_adjusted
tertiary$resistant <- tertiary$resistant_adjusted

tertiary[,c('resistant_adjusted','prop_adjusted','prop','adj','beta')] <- NULL

raw_data <- rbind(raw_data, tertiary)

raw_data <- raw_data %>% group_by(country, year, combo2) %>% 
  summarise(total = sum(total_test), resistant = sum(resistant)) %>%
  mutate(prop = resistant / total) 

##########################
# Linear model and prediction
##########################

abxcovs <- read.csv(paste0(filepath, 'abxcovs.csv'))# %>% filter(!is.na(location_name)) #%>% mutate(Country = location_name) %>% select(!location_name)
raw_data <- raw_data %>% left_join(abxcovs, by = c(c('country' = 'location_name'), c('year' = 'year_id')))

###########################
# Linear model for one combo and one country
############################
i='acinetobacter baumannii aminoglycoside'
combo_data <- raw_data[raw_data$combo2 == paste0(i) & raw_data$country == "United States of America",]

model1 <- lm(prop ~ year , data = combo_data)

print(model1)

combo_data$linear_full_pred <- predict(model1, type = 'response', allow.new.levels = TRUE)

print(cor(combo_data$prop,combo_data$linear_full_pred))

###########################
# Loop for linear and other modelling frameworks 
############################

metric_table <- tibble()
table_cor <- c()
for(i in unique(raw_data$combo2)) {
  combo_data <- raw_data[raw_data$combo2 == paste0(i),]
  
  model1 <- lm(prop ~ year , data = combo_data)
  
  combo_data$linear_full_pred <- predict(model1, type = 'response', allow.new.levels = TRUE)
  
  response = cbind(successes = combo_data$resistant,
                   failures = (combo_data$total - combo_data$resistant))
  
  xg_grid <- expand.grid(nrounds = c(50, 100, 200),
                         max_depth = c(4, 6, 8),
                         eta = (3:6) / 100, #learning rate 
                         colsample_bytree = .5,
                         min_child_weight = 1,
                         subsample = 1,
                         gamma = 0)
  
  xg_fit <- caret::train( prop ~ year + ddd_per_1000_fitted ,
                          data = combo_data,
                          #                       trControl = train_control,
                          verbose = T,
                          tuneGrid = xg_grid,
                          metric = "RMSE",
                          method = "xgbTree",
                          objective = "reg:logistic",iteration_range=200,
                          #trControl = cctrl,
                          # nrounds = 200000)
                          #n_estimators=1000, learning_rate=0.1, subsample=0.8, colsample_bytree=0.8,
                          #weights = combo_data$w
  )
  combo_data$xgboost_full_pred <- predict(xg_fit, combo_data, iteration_range=(best_iteration))
  
  ##########GAM
  
  covs_to_include = c('year','ddd_per_1000_fitted')
  gam_formula <- paste0('response ~ 1+ s(', paste(covs_to_include, collapse = ", bs = 'ts', k = 3) + s("), ", bs = 'ts', k = 3)")
  
  full_gam = mgcv::gam(response ~ 1 + s(year, bs = 'ts', k = 3), 
                       data = combo_data, 
                       family = 'quasibinomial',
                       #                     weights = mydata$w, 
                       control = list(nthreads = 2))
  full_gam$model_name = 'GAM'
  
  combo_data$gam_full_pred <- predict(full_gam, combo_data, type = 'response')
  
  weight_table <-tibble()
  for(j in c("linear_full_pred"      ,    "gam_full_pred"  ,    "xgboost_full_pred" )) {
    mt <- c(model  = j )
    mt$cor <- cor(combo_data$prop,combo_data[,j])
    mt$combo <- i
    metric_table <- rbind(metric_table,mt)
    weight_table <- rbind(weight_table, mt)
  }
  
  # for(x in  1:length(child_models)){
  #   for(j in 1:length(child_models)){
  #     if(x==j){
  #     } else{
  #         cor1 <- c(child_models[x],child_models[j],formatC(cor(combo_data[,(paste0(child_models[x], '_full_pred'))],combo_data[,(paste0(child_models[j], '_full_pred'))])^2, format = 'f', digits = 6))
  #                 row1 <- c(cor1,combo = paste0(i))
  #         table_cor <- rbind(table_cor,row1)
  #         rm(cor1, row1)
  #        }
  #       }
  #     }
  #   }
  denominator <- weight_table %>% summarise(cor = sum(cor))
  weights <- weight_table$cor / as.numeric(denominator)
  X = as.matrix(combo_data[c("linear_full_pred"      ,    "gam_full_pred"  ,    "xgboost_full_pred")])
  combo_data$stacked_full_pred <- rowWeightedMeans(X, w = weights)   
  stacked_pred <- combo_data[c("country", "year", "combo2", "stacked_full_pred")]
  
}
write.csv(combo_data, paste0(filepath,'combos.csv'), row.names = F)
write.csv(stacked_pred, paste0(filepath,'stacked.csv'), row.names = F)
write.csv(metric_table, paste0(filepath,'metric_table.csv'), row.names = F)
write.csv(table_cor, paste0(filepath,'table_cor.csv'), row.names = F)
