if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

SpeciesData<-read.csv("Data_Generated/TraitData_2017.csv")
TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(bayesplot)

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
Traits_all<-na.omit_species[na.omit_species$trait %in% c ("RDMC","RLD","RMR","SRL"),]
Traits_leaf<-na.omit_species[na.omit_species$trait %in% c ("SLA","LDMC"),]

Trait_splits_full = split(Traits_all, paste(Traits_all$trait))
Trait_splits_noH1<-split(Traits_leaf, paste(Traits_leaf$trait))

mk_data_full<-lapply(Trait_splits_full, make_matrix_function) # matrix - will need to add pop_id random effect
mk_data_noH1<-lapply(Trait_splits_noH1, make_matrix_function) 

# Data for running individual models by Species and trait 
dat_na.omit_SLA<-lapply(Species_splits, na.omit.fun, c("SAMPLE_ID","POP_ID","H_num","SLA"))
mk_data_SLA<-lapply(dat_na.omit_SLA, mk_data_function, "SLA")

dat_na.omit_RMR<-lapply(Species_splits, na.omit.fun, c("SAMPLE_ID","POP_ID","H_num","RMR"))
mk_data_RMR<-lapply(dat_na.omit_RMR, mk_data_function, "RMR")

# --------------------------------- Full model (Species and H_num together) --------------------------------------#
mod = stan_model("stan_models/All_TraitsbySpecies_matrix.stan")
  
mods_function_all <- function(df,
                          mod_file = "stan_models/All_TraitsbySpecies_matrix.stan",
                          iter = 1000,
                          cores = 2,
                          mod = stan_model(mod_file), ...){
  
  sampling(mod, df, iter = iter, cores = cores, ...) 
}

all_mods_full = lapply(mk_data_full, mods_function_all, mod = mod)
all_mods_noH1 = lapply(mk_data_noH1, mods_function_all, mod = mod) # Bulk effective sample size too low 

# -------------------------------- Plotting ------------------------------# 

out_names_full<-c("overall_mean", "ARTR_H1","ELTR_H1", "HEAN_H1","HECO_H1","HEVI_H1","MACA_H1","MUPO_H1","PAMU_H1","PLPA_H1","VUOC_H1",
                  "ACMI_H2","ARTR_H2","ELTR_H2", "HEAN_H2","HECO_H2","HEVI_H2","MACA_H2","MUPO_H2","PAMU_H2","PLPA_H2","VUOC_H2",
                  "ACMI_H3","ARTR_H3","ELTR_H3", "HEAN_H3","HECO_H3","HEVI_H3","MACA_H3","MUPO_H3","PAMU_H3","PLPA_H3","VUOC_H3",
                  "ACMI_H4","ARTR_H4","ELTR_H4", "HEAN_H4","HECO_H4","HEVI_H4","MACA_H4","MUPO_H4","PAMU_H4","PLPA_H4","VUOC_H4")

out_names_noH1<-c("overall_mean","ARTR_H2","ELTR_H2", "HEAN_H2","HECO_H2","HEVI_H2","MACA_H2","MUPO_H2","PAMU_H2","PLPA_H2","VUOC_H2",
                  "ACMI_H2", "ARTR_H3","ELTR_H3", "HEAN_H3","HECO_H3","HEVI_H3","MACA_H3","MUPO_H3","PAMU_H3","PLPA_H3","VUOC_H3",
                  "ACMI_H3", "ARTR_H4","ELTR_H4", "HEAN_H4","HECO_H4","HEVI_H4","MACA_H4","MUPO_H4","PAMU_H4","PLPA_H4","VUOC_H4")

names(all_mods_noH1$SLA)[grep("beta",names(all_mods_noH1$SLA))]<-out_names_noH1
plot(all_mods_noH1$SLA, pars = "beta")

#
#
#
#
#
#
#
#
#
# ------------------------------------------------------------- Variance model----------------------------------------------------------------------------
# Questions: Which traits are most variable through time and does this depend on GrowthForm?
# To "match" y ~ H_num + GrowthForm + H_num * GrowthForm + (1|POP_ID) I want: 
# CV ~ GrowthForm + (1|POP_ID) - H_num within a POP_ID is what the CV is calculated from 
# CV for every trait x POP_ID? 

fit_trait_var<- stan(file = "stan_models/example.stan", data = traits_stan_df, warmup = 500, iter = 2000, chains = 4, cores = 2, control = list(adapt_delta = 0.90))  

summary(fit_trait_var)$summary
aggregate(.~H_num, dat_na.omit_RMR$ACMI, FUN = mean)

# ----------------------------- Run SLA models one by one to see if they are working. Divergent transitions for most -------------------------------------------------------

mods_function <- function(df,
                              mod_file = "stan_models/IndividualTraitbySpecies.stan",
                              iter = 1000,
                              cores = 2,
                              mod = stan_model(mod_file), ...){

    sampling(mod, df, iter = iter, cores = cores, ...) 
}


# Why is the second function better? Because it makes many of the
# parameters user-changeable, you can do things like this

# precompile the model once
mod = stan_model("stan_models/IndividualTraitbySpecies.stan")

# then use that model for every element in the list
all_SLA_models = lapply(mk_data_SLA, mods_function, mod = mod)
all_RMR_models = lapply(mk_data_RMR, mods_function, mod = mod, warmup = 800, options(adapt_delta = 0.9))



# -----------------------------------------------

summary(all_SLA_models$ACMI)$summary
summary(all_RMR_models$ACMI)$summary
summary(all_RMR_models$ELTR)$summary


summary(fit_SLA_ARTR) 

summary(fit_SLA_ELTR) 

summary(fit_SLA_HECO) 

summary(fit_SLA_HEAN) 

summary(fit_SLA_HEVI) 

summary(fit_SLA_MUPO) 

summary(fit_SLA_MACA) 

summary(fit_SLA_PAMU) 

summary(fit_SLA_PLPA) 

summary(fit_SLA_VUOC) 


# ----------------------------- Run RMR models one by one to see if they are working. Divergent transitions for most -------------------------------------------------------

fit_RMR_ACMI <- stan(file = "stan_models/example.stan", data = mk_data_RMR$ACMI, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99))  
summary(fit_RMR_ACMI) 
aggregate(.~H_num, dat_na.omit_RMR$ACMI, FUN = mean)

fit_RMR_ARTR <- stan(file = "stan_models/example.stan", data = mk_data_RMR$ARTR, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_ARTR) 
aggregate(.~H_num, dat_na.omit_RMR$ARTR, FUN = mean)

fit_RMR_ELTR <- stan(file = "stan_models/example.stan", data = mk_data_RMR$ELTR, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99))  
summary(fit_RMR_ELTR) 
aggregate(.~H_num, dat_na.omit_RMR$ELTR, FUN = mean)

fit_RMR_HECO <- stan(file = "stan_models/example.stan", data = mk_data_RMR$HECO, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_HECO) 
aggregate(.~H_num, dat_na.omit_RMR$HECO, FUN = mean)

fit_RMR_HEAN <- stan(file = "stan_models/example.stan", data = mk_data_RMR$HEAN, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_HEAN) 
aggregate(.~H_num, dat_na.omit_RMR$HEAN, FUN = mean)

fit_RMR_HEVI <- stan(file = "stan_models/example.stan", data = mk_data_RMR$HEVI, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_HEVI) 
aggregate(.~H_num, dat_na.omit_RMR$HEVI, FUN = mean)

fit_RMR_MUPO <- stan(file = "stan_models/example.stan", data = mk_data_RMR$MUPO, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_MUPO) 
aggregate(.~H_num, dat_na.omit_RMR$MUPO, FUN = mean)

fit_RMR_MACA <- stan(file = "stan_models/example.stan", data = mk_data_RMR$MACA, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_MACA) 
aggregate(.~H_num, dat_na.omit_RMR$MACA, FUN = mean)

fit_RMR_PAMU <- stan(file = "stan_models/example.stan", data = mk_data_RMR$PAMU, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_PAMU) 
aggregate(.~H_num, dat_na.omit_RMR$PAMU, FUN = mean)

fit_RMR_PLPA <- stan(file = "stan_models/example.stan", data = mk_data_RMR$PLPA, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_PLPA) 
aggregate(.~H_num, dat_na.omit_RMR$PLPA, FUN = mean)

fit_RMR_VUOC <- stan(file = "stan_models/example.stan", data = mk_data_RMR$VUOC, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_VUOC) 
aggregate(.~H_num, dat_na.omit_RMR$VUOC, FUN = mean)
