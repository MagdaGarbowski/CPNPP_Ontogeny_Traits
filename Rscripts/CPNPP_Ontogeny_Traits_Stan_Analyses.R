if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

SpeciesData<-read.csv("Data_Generated/TraitData_2017.csv")
TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(bayesplot)
library(multcomp)
library(broom)

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
# Will still need to add random effect of population 
mod = stan_model("stan_models/All_TraitsbySpecies_matrix.stan")
  
mods_function_all <- function(df,
                          mod_file = "stan_models/All_TraitsbySpecies_matrix.stan",
                          iter = 1000,
                          cores = 2,
                          mod = stan_model(mod_file), ...){
  
  sampling(mod, df, iter = iter, cores = cores, ...) 
}

all_mods_full = lapply(mk_data_full, mods_function_all, mod = mod) # Bulk effective sample size too low 
all_mods_noH1 = lapply(mk_data_noH1, mods_function_all, mod = mod) # Bulk effective sample size too low 

#------------------------------- Contrasts -------------------------------# 
#--------------------- Species comparisons at Hnum_1----------------------------# 

RMR_test<-mods_function_all(mk_data_full$RMR, mod = mod) # Bulk Effective Samples Size (ESS) is too low
mcmc_RMR = RMR_test

Sps_comp_H1<-rep(1/11,11) #Define contribution from each species 
coefs_Sps_H1<-as.matrix(mcmc_RMR)[,1:11] # pick coeffecients to compare 
dat_Sps_H1<-data.frame(x = colnames(coefs_Sps_H1)) # get colnames to make Tukey matrix
tuk.mat_Sps_H1<-contrMat(n = table(dat_Sps_H1$x), type = "Tukey") # make Tukey matrix
Xmat_Sps_H1<-model.matrix(~x, dat_Sps_H1) 
pairwise.mat_Sps_H1<-tuk.mat_Sps_H1 %*% t(Xmat_Sps_H1)
comps_Sps_H1 = tidyMCMC(coefs_Sps_H1 %*% t(pairwise.mat_Sps_H1), conf.int = TRUE, conf.method = "HPDinterval") # get effect sie CI - is HPD appropriate? 

ggplot(comps_Sps_H1, aes(y = estimate, x = term)) + 
  geom_pointrange(aes(ymin = conf.low,
                      ymax = conf.high)) + geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous("Effect size") + scale_x_discrete("") + coord_flip() +
  theme_classic()

# Compare to raw_data boxplot to see if results make sense 
dat_raw_plotted<-boxplot(value ~ SPECIES, data = Trait_splits_full$RMR[Trait_splits_full$RMR$H_num == "H1",])

# Questions for Matt: 
# (1) Am I doing this right? 
# (2) Interpretation: If for a given comparison the CI does not intersect 0, the two are different. 
# (2) e.g. SPECIESVUOC is "higher" than all other species besides intercept (ACMI)? 
# (3) Is HPD (Highest posterior density interval) appropriate to use? 


#--------------------- H_num main effects----------------------------# 
# STUCK HERE
# Questions for Matt: 
# (1) How do I pull comparisons from ALL species at each H_num for comparisons? 
# There are 44  values (every species at every H_num) but I want to end up with 4 (i.e. "main effect" of H_num) 

comp_H<-rep(1/44, 44) # Define contribution from each H_num
coefs_comp_H<-as.matrix(mcmc_RMR)[,1:44] # All coefs should go into comparisons for H_num? 
dat_comp_H<-data.frame(x = colnames(coefs_comp_H)) # get colnames to make Tukey matrix
Xmat_comp_H<-model.matrix(~x, dat_comp_H) 
pairwise.mat_comp_H<-comp_H %*% t(Xmat_comp_H)

comps_comp_H = tidyMCMC(coefs_comp_H %*% t(pairwise.mat_comp_H), conf.int = TRUE, conf.method = "HPDinterval") # get effect sie CI - is HPD appropriate? 
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
