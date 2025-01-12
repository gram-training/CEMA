# get spneu out
library(ggplot2)
library(ggforce)
library(ggrepel)
library(raster)
library(gridExtra)
library(latticeExtra)
library(RColorBrewer)
source('/ihme/code/st_gpr/central/src/stgpr/api/public.R')
source("/ihme/cc_resources/libraries/current/r/get_location_metadata.R")
source("/ihme/cc_resources/libraries/current/r/get_covariate_estimates.R")
library(tidyverse)
amr_repo <- '~/amr/' 
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
location_meta <- subset(get_location_metadata(location_set_id=35, release_id = 16),select = c(location_id,
level,location_name,location_name_short,super_region_id,super_region_name, region_name,region_id,ihme_loc_id,parent_id), 
                                      droplevels = TRUE) 
spmeid <-  read.csv(paste0(amr_repo,'/maps/bug_drug_combos.csv'), stringsAsFactors = FALSE) %>% 
  filter(pathogen=="streptococcus_pneumoniae" & abx_class == "penicillin") 

output_path <- '/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/'
outputdir <-  paste0(output_path, '/stackers_7adj/streptococcus_pneumoniae/penicillin/') 

datapath <- '/ihme/limited_use/LIMITED_USE/LU_AMR/Oxford_antibiotic/'
abxcovs <- read.csv(paste0(datapath,"Gram_abx_covs_1990to2021.csv")) %>%
  dplyr::mutate(as.vector(scale(ddd_per_1000_fitted, center = TRUE, scale = TRUE))) %>%
  dplyr::mutate(as.vector(scale(prop_j01c_final, center = TRUE, scale = TRUE))) 
abxcovs <- cbind(abxcovs$location_id,abxcovs$year_id,abxcovs[,c(4,7,13:14)])
colnames(abxcovs) <-c('location_id','year_id','ddd_per_1000_est','j01c_prop','centered_ddd','centered_j01c')
abxcovs <- abxcovs %>% left_join(location_meta[,c('location_id','level')]) %>% dplyr::filter(level > 2)
forecast <- abxcovs[abxcovs$year_id==2021,]
for(i in c(2022,2023,2024)){
  forecast$year_id <- as.numeric(paste0(i))
  abxcovs <- rbind(abxcovs,forecast)
}
forecast <- abxcovs[abxcovs$location_id == 44855,]
forecast$location_id <- 60908L
abxcovs <- rbind(abxcovs,forecast)
forecast$location_id <- 94364L #95069
abxcovs <- rbind(abxcovs,forecast)
forecast$location_id <- 95069L
abxcovs <- rbind(abxcovs,forecast)

covs_list <- c(1093,50,463,1995,210,2314,2437,2447)

covs_dict <- read.csv('/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/longcovariates.csv', stringsAsFactors = F) %>%
  filter(!cov %in% c(2436,2437,319,320,160,142))

covbind <- tibble()
for(i in covs_list) { 
  j <- as.numeric(i) #[is.numeric(as.numeric(covs_list))]) {
  if(is.numeric(j) & !is.na(j)){
    covtest <- get_covariate_estimates(eval(parse(text = paste0("covariate_id =", j))), release_id = 16, year_id = 1990:2024) %>% 
      dplyr::group_by(location_id,year_id,covariate_name_short) %>% dplyr::summarise(mean_value = mean(mean_value,na.rm=T)) 
    z <- scale(covtest$mean_value, center = TRUE, scale = TRUE)
    covtest$val <- z[,1] 
    covtest <- covtest[,c('location_id','year_id','covariate_name_short','val')]
    covs_dict$covariate_name_short[covs_dict$cov == j] <- paste0(unique(covtest$covariate_name_short))[1]
    covbind <- rbind(covbind, covtest) #%>% left_join(covtest, on = c('location_id','year_id'))
  }
}
covs_dict$covariate_name_short[is.na(covs_dict$covariate_name_short)] <- covs_dict$cov[is.na(covs_dict$covariate_name_short)]
covbind <- pivot_wider(covbind, id_cols = c('location_id','year_id'), names_from = 'covariate_name_short', values_from = 'val')
covbind <- covbind %>% dplyr::select(colnames(covbind)[!colnames(covbind) %in% c('LDI_pc','level',"hib_indirect")])
covbind <- data.frame(covbind)  

run_id <- spmeid$resistance_run_id
lm <- read.csv(paste0(outputdir,'/custom_stage1_df.csv'))
st <- get_estimates(version_id = run_id, entity = 'spacetime') 
  setnames(st, c("val"), c("st"))
  st$st <- logit2prob(st$st) 
final <- get_estimates(run_id, entity = "final")
  setnames(final, c("val", "upper", "lower"), c("gpr_v", "gpr_u", "gpr_l"))
  data_used <- get_input_data(run_id, data_stage_name = "original")
  setnames(data_used, c("val", "upper", "lower"), c("data_val", "data_upper", "data_lower"))
  mydata <- merge(final, st, all.x = TRUE, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  mydata <- merge(mydata, lm, all.x = TRUE, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  mydata <- merge(mydata, data_used, all.x = TRUE, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  mydata <- merge(mydata, covbind, all.x = TRUE, by = c('location_id', 'year_id'))
  mydata <- mydata %>% left_join(location_meta[,c("location_id",'level','super_region_id','super_region_name','region_name','ihme_loc_id','location_name_short')])
  mydata$super_region_name[mydata$super_region_name == ''] <- 'Global'
  mydata <- mydata[mydata$level <= 3,]
  mydata$pathogen<- 'streptococcus_pneumoniae'
  mydata$abx_class<- 'penicillin'
  min_year = 1990
  max_year = 2024
    gpr<-list()
    for(i in 1:length(unique(location_meta$super_region_id))){
      subset <- mydata[mydata$super_region_id == unique(mydata$super_region_id)[i],]
      if(length(unique(subset$location_id))>1) {
        gpr[[i]] <- ggplot(subset)+
          geom_pointrange(aes(x=year_id, y = data_val, ymin = data_lower, ymax = data_upper, shape = as.factor(is_outlier)))+ 
          geom_line(aes(x=year_id, y = gpr_v,colour = 'gpr'))+
          geom_ribbon(aes(ymin = gpr_l, ymax=gpr_u, x = year_id), alpha = 0.5, fill = 'red', colour = 'red') +
          geom_line(aes(x=year_id, y = cv_custom_stage_1, colour = 'stage1_ensemble'))+
          geom_line(aes(x=year_id, y = st, colour = 'spatio_temporal'))+
          scale_colour_manual("", 
                              breaks = c("spatio_temporal","stage1_ensemble", "gpr"),
                              values = c("red", "green", "blue")) +
          theme(legend.position = "top")+
          scale_x_continuous("Year", 
                             breaks = seq(min_year, max_year, 5),
                             labels = seq(min_year, max_year, 5))+
          ylim(0,1)+
          ylab('Proportion of resistance')+
          theme_bw()+
          theme(legend.position = "bottom")+
          ggtitle(paste0(str_to_title(gsub('[-]|[_]',' ', z))))+
          facet_wrap_paginate(~location_name_short, nrow = ceiling(sqrt(length(unique(subset$location_id)))),page=i)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
          theme(plot.title = element_text(hjust = 0.5))
      }
    }
    ggsave(paste0(outputdir,"plot.pdf"), marrangeGrob(grobs = gpr, nrow = 1, ncol = 1),#newpage = T),
           height = 20, width = 40, units = 'cm')
write.csv(mydata, paste0(outputdir,'estimates_and_data.csv'), row.names = F)
