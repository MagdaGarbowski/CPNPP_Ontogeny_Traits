setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
library(dplyr)
library(ggplot2)

TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")

# Growth Rates
SpeciesData_Growthrate<-read.csv("Data_Generated/TraitData_GrowthRates_2017.csv")
SpeciesData_Growthrate$H_num<-as.character(SpeciesData_Growthrate$H_num)
SpeciesData_Growthrate<-SpeciesData_Growthrate[!(SpeciesData_Growthrate$H_num) %in% c("H1"),] # Remove since no values at H_num1
SpeciesData_Growthrate<-SpeciesData_Growthrate[c("H_num","POP_ID","SPECIES","GrowthForm","RER","RGR_Tot")]
SpeciesData_Growthrate_long<-reshape(SpeciesData_Growthrate, 
                                     direction = "long",
                                     varying = c("RER", "RGR_Tot"),
                                     v.names = "value",
                                     timevar = "trait",
                                     times = c("RER", "RGR_Tot"))


TraitData_Pop_Avg_2017$value<-ifelse((TraitData_Pop_Avg_2017$GrowthForm %in% c("FORB","SHRUB") & TraitData_Pop_Avg_2017$H_num == "H1" & TraitData_Pop_Avg_2017$trait == "SLA"), NA, TraitData_Pop_Avg_2017$value) # Add NA for H1 SLA for forb species 
TraitData_Pop_Avg_2017$value<-ifelse((TraitData_Pop_Avg_2017$GrowthForm %in% c("FORB","SHRUB") & TraitData_Pop_Avg_2017$H_num == "H1" & TraitData_Pop_Avg_2017$trait == "LDMC"), NA, TraitData_Pop_Avg_2017$value) # Add NA for H1 SLA for forb species 

TraitData_Pop_FullHarvests<-as.data.frame(TraitData_Pop_Avg_2017[TraitData_Pop_Avg_2017$trait %in% c("RMR", "SRL", "RDMC", "RTD"),])
TraitData_Pop_noH1<-as.data.frame(TraitData_Pop_Avg_2017[TraitData_Pop_Avg_2017$trait %in% c("SLA", "LDMC"),])
TraitData_Pop_noH1<-TraitData_Pop_noH1[!TraitData_Pop_noH1$H_num %in% c("H1"),]

TraitData_Pop_FullHarvests_splits<-split(TraitData_Pop_FullHarvests, TraitData_Pop_FullHarvests$trait)
TraitData_Pop_noH1_splits<-split(TraitData_Pop_noH1, TraitData_Pop_noH1$trait)


summarise_function<-function(df){
  ddply(df, c("GrowthForm"), summarise,
        N = length(value),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        se = sd/sqrt(N))
}  

TraitData_Pop_sum_full<-lapply(TraitData_Pop_FullHarvests_splits, summarise_function)
TraitData_Pop_sum_leaf<-lapply(TraitData_Pop_noH1_splits, summarise_function)
TraitData_Pop_sum_RER<-summarise_function(SpeciesData_Growthrate_long[SpeciesData_Growthrate_long$trait == "RER",])
TraitData_Pop_sum_RGR<-summarise_function(SpeciesData_Growthrate_long[SpeciesData_Growthrate_long$trait == "RGR_Tot",])

TraitData_Pop_sum_full$RMR
TraitData_Pop_sum_full$RTD
TraitData_Pop_sum_full$SRL
TraitData_Pop_sum_full$RDMC
TraitData_Pop_sum_leaf$LDMC
TraitData_Pop_sum_leaf$SLA

