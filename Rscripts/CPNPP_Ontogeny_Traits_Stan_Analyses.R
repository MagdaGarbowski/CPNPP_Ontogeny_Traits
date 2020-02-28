if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

SpeciesData<-read.csv("Data_Generated/TraitData_2017.csv")
TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(bayesplot)
bayesplot::color_scheme_set("brightblue")


SpeciesData$H_num<-as.character(as.factor(SpeciesData$H_num))
SpeciesData$SLA<-ifelse((SpeciesData$GrowthForm %in% c("FORB","SHRUB") & SpeciesData$H_num == "H1"), NA, SpeciesData$SLA) # Add NA for H1 SLA for forb species 
SpeciesData$SLA<-ifelse((SpeciesData$SLA == 0), NA, SpeciesData$SLA)
Species_splits = split(SpeciesData, paste(SpeciesData$SPECIES))

TraitData_Pop_Avg_2017$H_num<-as.character(as.factor(TraitData_Pop_Avg_2017$H_num))
TraitData_Pop_Avg_2017$value<-ifelse((TraitData_Pop_Avg_2017$GrowthForm %in% c("FORB","SHRUB", "GRASS") & TraitData_Pop_Avg_2017$H_num == "H1" & TraitData_Pop_Avg_2017$trait %in% c("SLA", "LDMC")), NA, TraitData_Pop_Avg_2017$value)

# Data for variance estimates of traits
na.omit_traits<-na.omit.fun(TraitData_Pop_Avg_2017, c("POP_ID", "value","trait","SPECIES"))
traits_stan_df<-mk_data_traitvar_function(na.omit_traits)

# Data for full model - mean estimates of traits by species and time 
na.omit_species<-na.omit.fun(TraitData_Pop_Avg_2017, c("POP_ID", "value","trait","SPECIES", "H_num"))
Traits_all<-na.omit_species[na.omit_species$trait %in% c ("RDMC","RTD","RMR","SRL"),]
Traits_leaf<-na.omit_species[na.omit_species$trait %in% c ("SLA","LDMC"),]

Trait_splits_full = split(Traits_all, paste(Traits_all$trait))
Trait_splits_noH1<-split(Traits_leaf, paste(Traits_leaf$trait))


################################################################################
data_full<-lapply(Trait_splits_full, prep_data) 
data_noH1<-lapply(Trait_splits_noH1, prep_data) 



# Full model (Species and H_num together) --------------------------------------#
# Will still need to add random effect of population 
mod = stan_model("stan_models/All_Traits_Random_Effects_HS.stan")
  
mods_function_all <- function(df,
                          mod_file = "stan_models/All_Traits_Random_Effects_HS.stan",
                          iter = 1000,
                          cores = 2,
                          mod = stan_model(mod_file), ...){
  
  sampling(mod, df, iter = iter, cores = cores, ...) 
}

# The new model does not sample as efficiently as the matrix
# representation, so higher adapt delta and iterations are needed

all_mods_full = lapply(data_full, mods_function_all, mod = mod,
                       pars = c("beta_Hnum_raw", "beta_sp_raw",
                                "beta_pop_raw", "beta_sp_Hnum_raw"),
                       include = FALSE,
                       warmup = 1000, iter = 1500,
                       control = list(adapt_delta = 0.95))  # Bulk Effective Sample Size too low 

all_mods_noH1 = lapply(data_noH1, mods_function_all, mod = mod,
                       pars = c("beta_Hnum_raw", "beta_sp_raw",
                                "beta_pop_raw", "beta_sp_Hnum_raw"),
                       include = FALSE,
                       warmup = 1000, iter = 1500,
                       control = list(adapt_delta = 0.95)) # Divergent transitions and Bulk Effective Sample Size too low 

# RMR test ------------------------------------------------------------# 

out_names_full<-colnames(mk_data_full$RDMC$M)

RMR_test<-mods_function_all(data_full$RMR, mod = mod, 
                            pars = c("beta_Hnum_raw","beta_sp_raw", "beta_pop_raw","beta_sp_Hnum_raw"),
                            include = FALSE, 
                            warmup = 1000, iter = 1500, 
                            control = list(adapt_delta = 0.95)) # Bulk effective sample size too low 
View(summary(RMR_test)[["summary"]])
plot(RMR_test, pars = "beta_Hnum")

# All mods  ----------------------------------------------------------# 

View(summary(all_mods_full$RMR)[["summary"]])
plot(all_mods_full$RMR, pars = "beta_Hnum")

View(summary(all_mods_full$RDMC)[["summary"]])
plot(all_mods_full$RDMC, pars = "beta_Hnum")

View(summary(all_mods_full$SRL)[["summary"]])
plot(all_mods_full$SRL, pars = "beta_Hnum")

View(summary(all_mods_full$RTD)[["summary"]])
plot(all_mods_full$RTD, pars = "beta_Hnum")

View(summary(all_mods_noH1$SLA)[["summary"]])
plot(all_mods_noH1$SLA, pars = "beta_Hnum")

View(summary(all_mods_noH1$LDMC)[["summary"]])
plot(all_mods_noH1$LDMC, pars = "beta_Hnum")


