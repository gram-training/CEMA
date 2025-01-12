library(tidyverse)
library(dplyr)
library(foreign)
library(ini)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(raster)
library(gridExtra)
library(latticeExtra)
library(RColorBrewer)
library(sf)
library(viridis)
source("/ihme/cc_resources/libraries/current/r/get_population.R")
source("/ihme/cc_resources/libraries/current/r/get_location_metadata.R")
###########################
source('/ihme/code/st_gpr/central/src/stgpr/api/public.R')
#source('/share/code/st_gpr/central/stgpr/r_functions/utilities/utility.r')
username <- Sys.getenv("USER")
amr_repo <- sprintf("/ihme/homes/%s/amr/", username)
pmap <- read.csv('/mnt/share/homes/giselaa/amr/maps/pathogen_map.csv') #%>% select(pathogen,pathogen_name_long) 
amap <- read.csv('/mnt/share/homes/giselaa/amr/maps/antibiotic_map.csv') #%>% select(abx_class,abx_class_name_short) 

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#appending gpr estimates and input data to plot
#1. load files
location_md2 <- subset(get_location_metadata(location_set_id=35, release_id=16),
                       select = c(location_id,level,location_name,location_name_short,super_region_name,region_name,super_region_id, region_id,ihme_loc_id,parent_id), droplevels = TRUE)

combos<-  read.csv(paste0(amr_repo,'/maps/bug_drug_combos.csv'), stringsAsFactors = FALSE) %>% filter(kevin_catrin_exclude==0)# %>% dplyr::select(pathogen, abx_class,modelable_entity_id) 
combos$pathogen[combos$pathogen == "acinetobacter_baumanii"] <- 'acinetobacter_baumannii'
combos$pathogen[grepl('gonorrheae', combos$pathogen)] <- 'neisseria_gonorrhoeae'
combos$abx_class[combos$abx_class=='rifampicin_new']<-'mdr'
combos$abx_class[combos$abx_class=='isoniazid_new']<-'xdr'
combos <- combos[!grepl('prop',combos$abx_class),]
combos$abx_class <- sub('_retreated',"",combos$abx_class)
combos$combo<-paste0(combos$pathogen,"-",combos$abx_class)
combos <- combos[!grepl("tuberculosis",combos$combo),]

#file_runs_complete <- as.data.frame(read.csv("/ihme/cod/prep/amr/04_resistance/stgpr/configfinal2.csv", stringsAsFactors = F) %>% filter(complete_status == 1))
#file_runs_complete$run_id <- as.numeric(file_runs_complete$run_id)
ucovpath <- "/share/covariates/ubcov/model/output/"
playpath <- "/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/stgpr_amp5/"
pathconfig<-'/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/'

# I read the file that contains runIDs AND a column of which are complete
files = c('ab_long_config_table.csv','secondl_config_table.csv','third_config_table.csv')
completed <- tibble()
for(i in files){
  f <- read.csv(paste0(pathconfig, i))
  completed <- rbind(completed, f)
}
#write.csv(completed,paste0(pathconfig,'flagconfig.csv'),row.names=F)
#completed <- read.csv(paste0(playpath,'flagcomplete.csv')) %>% filter(complete==1)
#completed$combo <- gsub('raw-data test for ','', completed$description)
completed <- read.csv(paste0(output_path,'stackers_4adj/interim_runids.csv')) %>% filter(latestat=='success')
completed$dup <- ave(completed$description, completed$description, FUN = length) > 1L
completed <- completed[completed$dup==F | ((completed$dup == T) & (completed$latestat=="success")),] #global_run_id == 211467)), ]
completed<-combos
gpr<-list()
fileform<-tibble()
i<-0
for(combination in unique(completed$combo)) { #[48:length(completed$run_id)]) { #j<-unique(completed$run_id[grepl('aureus-methicillin',completed$combo)])
  #get some locals 
  #combination <- completed$description[completed$latest == j]
  y <- substring(combination,1,stringr::str_locate(combination,"-")-1)[1]
  z <- substring(combination,stringr::str_locate(combination,"-")+1,)[1]
  #input <- new_data[new_data$combo == ]  
  
  ##### Time series PDFs saved
  max_year<-2024
#  savepdf<-paste0(playpath,"global/")
  #    national <- read.csv(paste0(playpath,name),stringsAsFactors = F)
  #   x <- substring(name,1,stringr::str_locate(name,"-")-1)[1]
  #  z <- substring(name,stringr::str_locate(name,"-")+1,(nchar(name)-4))[1]
  run_id <- completed$resistance_run_id[completed$combo==combination]#[i]
  if(exists(paste0(output_path,'stackers_5adj/',y,'/', z,'/custom_stage1_df.csv'))){
    lm <-  read.csv(paste0(output_path,'stackers_5adj/',y,'/', z,'/custom_stage1_df.csv')) %>% rename('lm'='cv_custom_stage_1')
  }else{
    lm <-  read.csv(paste0(output_path,'stackers_4adj/',y,'/', z,'/custom_stage1_df.csv')) %>% rename('lm'='cv_custom_stage_1')
  }
  #read in model parameters
#  params <- get_parameters(run_id)
  st <- get_estimates(version_id = run_id, entity = 'spacetime') 
  if(length(unique(st$val))>1) {
    setnames(st, c("val"), c("st"))
  #lm <- get_estimates(version_id = run_id, entity = 'stage1') 
  #setnames(lm, c("val"), c("lm"))
  # back transform the estimates
  st$st <- logit2prob(st$st) #inv.logit(st$st)
  #lm$stage1 <- logit2prob(lm$lm) #inv.logit(lm$lm)
  # Get final (raked) results and merge on the original data
  final <- get_estimates(run_id, entity = "final")
  setnames(final, c("val", "upper", "lower"), c("gpr_v", "gpr_u", "gpr_l"))
  data_used <- get_input_data(run_id, data_stage_name = "original")
  setnames(data_used, c("val", "upper", "lower"), c("data_val", "data_upper", "data_lower"))
  #data_used <- read.csv("/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/stackers_raw/staphylococcus_aureus/methicillin/input.csv") %>% dplyr::mutate(data_val = val, data_upper = (1.96*se+val),data_lower = (1.96*se-val)) 
  mydata <- merge(final, st, all.x = TRUE, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  mydata <- merge(mydata, lm, all.x = TRUE, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  mydata <- merge(mydata, data_used, all.x = TRUE, by = c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  #mydata <- left_join(mydata, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
#  rm(lm, st,gpr, data_used)

  # mydata$location_id <- ifelse((mydata$level == 5), mydata$parent_id, mydata$location_id)
  # mydata$level<-NULL
  # mydata$parent_id<-NULL
  # mydata <- left_join(mydata, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
  # mydata$location_id <- ifelse((mydata$level == 4), mydata$parent_id, mydata$location_id)
  # mydata$level<-NULL
  # mydata$parent_id<-NULL
#  mydata <- left_join(mydata, location_md2[,c('location_id','level','parent_id')], by = c('location_id'))
  mydata <- mydata %>% left_join(location_md2[,c("location_id",'level','super_region_id','super_region_name','region_name','ihme_loc_id','location_name_short')])
  mydata$super_region_name[mydata$super_region_name == ''] <- 'Global'
    mydata <- mydata[mydata$level <= 3,]
 mydata$combo<-combination   
 mydata$pathogen<-y
 mydata$abx_class<-z
    fileform <- rbind(fileform,mydata)
  
      #rm(subset,name, x, z,gpr,national)
  }
}
write.csv(fileform,paste0(output_path,'stgpr_amp5/fileform.csv'), row.names = F)
fileform$label = str_to_title(gsub('[_]',' ', fileform$combo))
gpr<-list()
i<-0
for(z in unique(fileform$abx_class)) {
  mydata<-fileform[(fileform$abx_class==z)& (fileform$level == 0),]
i <- i+1
gpr[[i]] <-
ggplot()+
  geom_line(data = mydata, aes(x = year_id, y = gpr_v, colour = super_region_name), size = 1.5)+
  geom_ribbon(data = mydata, aes(x = year_id, ymin = gpr_l, ymax = gpr_u, fill = super_region_name), alpha = 0.1)+
  # geom_line(data = income[income == 'High income countries',], aes(x = year, y = ddd_per_1000_per_day, group = income, linetype = income), colour = 'black', size = 1)+
  #geom_line(data = income[income != 'High income countries',], aes(x = year, y = ddd_per_1000_per_day, group = income, linetype = income), colour = 'black', size = 1)+
  #  geom_ribbon(data = income[income != 'High income countries',], aes(x = year, ymin = ddd_per_1000_per_day_lower, ymax = ddd_per_1000_per_day_upper), alpha = 0.1, fill = 'black')+
  scale_colour_manual(values = c(
    '#099999', #'#000000', #WB lmic
    "#fe9929", #central europe
                                 "#8c2d04", #HI
                                 # '#000000', #WB HIC
                                 "#e78ac3", #latin America
                                 "#984ea3", #NAME
                                 "#4daf4a", #South Asia
                                 "#377eb8", # SE Asia
                                 "#e41a1c"))+ #SSA+
  scale_fill_manual(values = c(
    '#099999', #'#000000', #WB lmic
    "#fe9929", #central europe
                               "#8c2d04", #HI
                               #'#000000', #WB HIC
                               "#e78ac3", #latin America
                               "#984ea3", #NAME
                               "#4daf4a", #South Asi
                               "#377eb8", # SE Asia
                               "#e41a1c"), guide = 'none')+
  scale_linetype_manual(values = c('longdash', 'dotdash'))+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle(paste0(str_to_title(mydata$abx_class))) +
  labs(x = 'Year', y = 'Proportion of resistance', colour = 'GBD Super Region') + #, linetype = 'World Bank Income Group')+
  ylim(0,1)+
  facet_wrap(~label)
}
ggsave(paste0(output_path,"stgpr_amp5/g_abx.pdf"), marrangeGrob(grobs = gpr, nrow = 1, ncol = 1),#newpage = T),
       height = 20, width = 40, units = 'cm')

  

  }
}


# trace down logs for failures
for(i in c(211800,211801, 211802,211803)){#},211416,211417,211418)){
logs <- c()
  login = list.files('/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/logs/', pattern = paste0(i), recursive = T)
  logs = c(logs,login)
tracer <- substring(logs,1,8) #str_locate("/",logs))#gsub('[/*.*$]','',logs)
print(tail(tracer))
a<-tracer[length(tracer)]
b<-list.files(paste0('/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/logs/',a), pattern = '_0.e')
#z<-fread(b)
#print(ifelse(grepl()))
}

b<-tibble()
for(i in c(211641,
           211648,
           211649,
           211670,
           211746,
           211750)){
#  a<-as.data.frame(get_parameters(i))
#  b<-bind_rows(b,a)
 # path_to_data <-  a$pathpaste0(inputdir,'input.csv') 
  #path_to_custom_stage_1 <-  paste0(inputdir,'custom_stage1_df.csv')
    combination<-unique(filer3$description[filer3$latest == i]) 
  y <- substring(combination,1,stringr::str_locate(combination,"-")-1)[1]
  z <- substring(combination,stringr::str_locate(combination,"-")+1,)[1]
  y <- ifelse(grepl('anii',y), gsub("acinetobacter_baumanii", 'acinetobacter_baumannii',y),y)
  y <- ifelse(grepl('orrheae',y), gsub("neisseria_gonorrheae", 'neisseria_gonorrhoeae',y),y)
  inputdir <-  paste0(output_path, 'stackers_adjusted_fxid2/',y,'/', z,'/') 
  over1fix <- read.csv(paste0(inputdir,'input.csv'))
  over1fix$nid[is.na(over1fix$nid)] <-   over1fix$source[is.na(over1fix$nid)]
write.csv(over1fix,paste0(inputdir,'input.csv'), row.names = F)
  stgpr_sendoff(i, project = "proj_amr", nparallel = 50, log_path = paste0(output_path,'logs/'))
  }

  # Read in input data, I have no access to original files for MRSA nor SPNEU
  # if(combination == 'staphylococcus_aureus-methicillin') {
  #   input_data <- read.csv("/ihme/homes/annieb6/AMR/staph/input_data/staph_data_plus.csv")
  # } else {
  #   input_data <- read.csv(params$path_to_data)
  # }
  # 
  # #create a summary table to export
  # if(i >= 91) {
  #   input_data2 <-input_data %>% left_join(location_md2,by='location_id')
  #   line_input <- input_data2 %>% dplyr::summarize (data_points = n(), regions = n_distinct(unique(region_id), rm.na = T)) 
  #   rm(input_data2)
  #   line_input$pathogen <- y
  #   line_input$abx_class <- z
  # } else {
  #   line_input <- input_data %>% dplyr::summarize (data_points = n(), regions = n_distinct(unique(region_id), rm.na = T)) 
  #   line_input$pathogen <- y
  #   line_input$abx_class <- z
  # }
  # summary <- rbind(summary,line_input)
  # 
  # # Calculate confidence intervals for the input data and bound to sensible values
  # input_data$upper_ci <- input_data$val+(1.96*sqrt(input_data$var))
  # input_data$lower_ci <- input_data$val-(1.96*sqrt(input_data$var))
  # input_data$lower_ci[input_data$lower_ci<0] <- 0
  # input_data$upper_ci[input_data$upper_ci>1] <- 1
  # #we reasign subnational points to national
  # input_data$location_id[input_data$level == 4] <- input_data$parent_id[input_data$level == 4]
  # input_data$location_id[input_data$location_id == 44767] <- 95
  # input_data <- input_data[,c('location_id','year_id','nid','age_group_id','sex_id','measure_id','val','lower_ci','upper_ci','is_outlier')]
  # #merge input data onto output data
  # mydata <- merge(mydata,input_data, by = c('location_id','year_id','age_group_id', 'sex_id'), all.x = T, all.y = T)
  # #Subset to national data
  # mydata <- mydata[mydata$level == 3,]
  # mydata <- mydata[!is.na(mydata$gpr_mean),]
  # data_to_export <- mydata %>% dplyr::select('location_id','ihme_loc_id','year_id','super_region_id','gpr_mean','gpr_lower','gpr_upper','st','stage1','val','lower_ci','upper_ci','nid','is_outlier')
  # 
  # #order data by alphabetical country and region
  # write.csv(data_to_export, paste0(playpath,y,"-",z,".csv"), row.names = F)
  # rm(data_to_export,mydata,input_data,combination,y,z,line_input)
}

write.csv(summary, paste0(playpath,"all/summary.csv"), row.names = F)
summary<-read.csv(paste0(playpath,"all/summary.csv"), stringsAsFactors = F)
summary$studies<-NULL
summary$locations<-NULL
names(summary)<-c('dp','r','pathogen','abx_class') #summary$pivot1 <- str_locate(summary$combination,"-")-1 #summary$pivot2 <- str_locate(summary$combination,"-")+1 #summary$pathogen <- substr(summary$combination,1,summary$pivot1) #summary$abx_class <- substr(summary$combination,summary$pivot2,nchar(summary$combination)) #summary$combination<-NULL #summary$pivot1<-NULL #summary$pivot2<-NULL
summary <- summary %>% pivot_wider(id_cols = c(pathogen), names_from = abx_class, values_from = c(dp,r),names_sep = "_")
summary<-as.data.table(summary)
pabs <- read.csv("/ihme/homes/giselaa/res/pathogens_assessed_by_syndrome.csv", stringsAsFactors = F)
pabs$pathogen<-trim(pabs$pathogen)
pabs <- left_join(pabs,summary,by = 'pathogen')
write.csv(pabs, "/ihme/homes/giselaa/res/pathogens_by_syndrome.csv", row.names = F, na = "")

##################################################################################
#####GLOBAL AND SUPEREGION, 1PAGE PER COMBO = 88 PAGES
i<-i+1
gpr[[i]] <-
  ggplot()+
  geom_line(data = mydata, aes(x = year_id, y = gpr_v, colour = super_region_name), size = 1.5)+
  geom_ribbon(data = mydata, aes(x = year_id, ymin = gpr_l, ymax = gpr_u, fill = super_region_name), alpha = 0.1)+
  # geom_line(data = income[income == 'High income countries',], aes(x = year, y = ddd_per_1000_per_day, group = income, linetype = income), colour = 'black', size = 1)+
  #geom_line(data = income[income != 'High income countries',], aes(x = year, y = ddd_per_1000_per_day, group = income, linetype = income), colour = 'black', size = 1)+
  #  geom_ribbon(data = income[income != 'High income countries',], aes(x = year, ymin = ddd_per_1000_per_day_lower, ymax = ddd_per_1000_per_day_upper), alpha = 0.1, fill = 'black')+
  scale_colour_manual(values = c("#fe9929", #central europe
                                 "#8c2d04", #HI
                                 # '#000000', #WB HIC
                                 "#e78ac3", #latin America
                                 '#099999', #'#000000', #WB lmic
                                 "#984ea3", #NAME
                                 "#4daf4a", #South Asia
                                 "#377eb8", # SE Asia
                                 "#e41a1c"))+ #SSA+
  scale_fill_manual(values = c("#fe9929", #central europe
                               "#8c2d04", #HI
                               #'#000000', #WB HIC
                               "#e78ac3", #latin America
                               '#099999', #'#000000', #WB lmic
                               "#984ea3", #NAME
                               "#4daf4a", #South Asi
                               "#377eb8", # SE Asia
                               "#e41a1c"), guide = 'none')+
  scale_linetype_manual(values = c('longdash', 'dotdash'))+
  theme_bw()+
  ggtitle(paste0(str_to_title(gsub('[-]|[_]',' ', combination))))+
  labs(x = 'Year', y = 'Proportion of resistance', colour = 'GBD Super Region') + #, linetype = 'World Bank Income Group')+
  ylim(0,1)
}
}
ggsave(paste0(savepdf,"/nosmart.pdf"), marrangeGrob(grobs = gpr, nrow = 1, ncol = 1),#newpage = T),
       height = 20, width = 40, units = 'cm')

############################################################################################
#########################################################
# RAW+TREND DATA FOR ALL LOCATIONS, 1PAGE PER SUPEREGION = 8 PAGES,1 FILE PER COMBO
fileform$stage1_ensemble <- fileform$lm
fileform$stage1 <- fileform$lm
for(z in unique(fileform$combo)) {
  mydata<-fileform[(fileform$combo==z)& (fileform$level <= 3),]
  min_year <- 2000#min(mydata$year_id)
max_year<-2023#max(mydata$year_id)
gpr<-list()
for(i in 1:length(unique(location_md2$super_region_id))){
  subset <- mydata[mydata$super_region_id == unique(mydata$super_region_id)[i],]
  if(length(unique(subset$location_id))>1) {
    gpr[[i]] <- ggplot(mydata)+
      geom_pointrange(aes(x=year_id, y = data_val, ymin = data_lower, ymax = data_upper, shape = as.factor(is_outlier)))+ 
      geom_line(aes(x=year_id, y = gpr_v,colour = 'gpr'))+
      geom_ribbon(aes(ymin = gpr_l, ymax=gpr_u, x = year_id), alpha = 0.5, fill = 'red', colour = 'red') +
      geom_line(aes(x=year_id, y = stage1, colour = 'stage1_ensemble'))+
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
ggsave(paste0(playpath,"combos/",z,".pdf"), marrangeGrob(grobs = gpr, nrow = 1, ncol = 1),#newpage = T),
       height = 20, width = 40, units = 'cm')
}
#rm(subset,name, x, z,gpr,national)
###############################################################################################


#A loop that compiles all csv in one large document
files = list.files(playpath,pattern="*.csv$")
logs <- c()
for(i in c(209715, 209704, 209696, 209673, 209653, 209650, 209629, 209627, 209626)){
login = list.files('/mnt/team/amr/priv/intermediate_files/04_prevalence_resistance/time_series/logs/', pattern = paste0(i), recursive = T)
logs = c(logs,login)
}
tracer <- substring(logs,1,8) #str_locate("/",logs))#gsub('[/*.*$]','',logs)
all_data<-c()
for(j in 1:length(files)){
  name <- paste0(files[j])
  temp <- read.csv(paste0(playpath,name),stringsAsFactors = F)
  y <- substring(name,1,stringr::str_locate(name,"-")-1)[1]
  z <- substring(name,stringr::str_locate(name,"-")+1,(nchar(name)-4))[1]
  temp$pathogen <- paste0(y)
  temp$abx_class <- paste0(z)
  all_data <- rbind(all_data,temp)
  rm(temp,name,y,z)
}

write.csv(all_data,paste0(playpath,'summary/all_data.csv'),row.names = F)


##### Time series PDFs saved
max_year<-2018
savepdf<-paste0(playpath,"/pdf/")
for(j in 1:length(files)){
  name <- paste0(files[j])
  national <- read.csv(paste0(playpath,name),stringsAsFactors = F)
  x <- substring(name,1,stringr::str_locate(name,"-")-1)[1]
  z <- substring(name,stringr::str_locate(name,"-")+1,(nchar(name)-4))[1]
  min_year <- min(national$year_id)
  gpr<-list()
 for(i in 1:length(unique(national$super_region_id))){
    subset <- national[national$super_region_id == unique(national$super_region_id)[i],]
    gpr[[i]] <- ggplot(subset)+
        geom_pointrange(aes(x=year_id, y = val, ymin = lower_ci, ymax = upper_ci, shape = as.factor(is_outlier)))+ 
      geom_line(aes(x=year_id, y = gpr_mean,colour = 'gpr'))+
      geom_ribbon(aes(ymin = gpr_lower, ymax=gpr_upper, x = year_id), alpha = 0.5, fill = 'red', colour = 'red') +
        geom_line(aes(x=year_id, y = stage1, colour = 'stage1_ensemble'))+
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
        facet_wrap_paginate(~ihme_loc_id, nrow = ceiling(sqrt(length(unique(subset$location_id)))),page=i)+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
        theme(plot.title = element_text(hjust = 0.5))
 }
  ggsave(paste0(savepdf,"/",x,"-",z,".pdf"), marrangeGrob(grobs = gpr, nrow = 1, ncol = 1,newpage = T),
         height = 20, width = 40, units = 'cm')
    rm(subset,name, x, z,gpr,national)
}

