if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")
Pop_avg_data_wrates<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(bayesplot)
library(HDInterval)
bayesplot::color_scheme_set("brightblue")

Pop_avg_data_wrates$H_num<-as.character(as.factor(Pop_avg_data_wrates$H_num))
Pop_avg_data_wrates$value<-ifelse((Pop_avg_data_wrates$GrowthForm %in% c("FORB","SHRUB", "GRASS") & Pop_avg_data_wrates$H_num == "H1" & Pop_avg_data_wrates$trait %in% c("ln.SLA", "ln.LDMC", "value.SLA","value.LDMC")), NA, Pop_avg_data_wrates$value)

# Data for variance estimates of traits
# na.omit_traits<-na.omit.fun(TraitData_Pop_Avg_2017, c("POP_ID", "value","trait","SPECIES"))
# traits_stan_df<-mk_data_traitvar_function(na.omit_traits)

# Data for full model - mean estimates of traits by species and time 
na.omit_species<-na.omit.fun(Pop_avg_data_wrates, c("POP_ID", "value","trait","SPECIES", "H_num"))
Traits_all<-na.omit_species[na.omit_species$trait %in% c ("ln.RDMC","ln.RTD","ln.RMR","ln.SRL"),]
Traits_noH1<-na.omit_species[na.omit_species$trait %in% c ("ln.SLA","ln.LDMC", "RER_ln","RGR_Tot_ln"),]

Trait_splits_full = split(Traits_all, paste(Traits_all$trait))
Trait_splits_noH1<-split(Traits_noH1, paste(Traits_leaf$trait))


################################################################################
data_full<-lapply(Trait_splits_full, prep_data) 
data_noH1<-lapply(Trait_splits_noH1, prep_data) 



# Full model (Species and H_num together) --------------------------------------#
# Will still need to add random effect of population 
mod = stan_model("stan_models/All_Traits_Random_Effects.stan")
  
mods_function_all <- function(df,
                          mod_file = "stan_models/All_Traits_Random_Effects.stan",
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

# Estimate differences between H_nums for all traits ----------------------------

H_num_diff_RMR<-Hnum_difference_function(all_mods_full$ln.RMR)
H_num_diff_RDMC<-Hnum_difference_function(all_mods_full$ln.RDMC)
H_num_diff_SRL<-Hnum_difference_function(all_mods_full$ln.SRL)
H_num_diff_RTD<-Hnum_difference_function(all_mods_full$ln.RTD)
H_num_diff_SLA<-Hnum_difference_function_noH1(all_mods_noH1$ln.SLA)
H_num_diff_LDMC<-Hnum_difference_function_noH1(all_mods_noH1$ln.LDMC)
H_num_diff_RER<-Hnum_difference_function_noH1(all_mods_noH1$RER_ln)
H_num_diff_RGR<-Hnum_difference_function_noH1(all_mods_noH1$RGR_Tot_ln)

# Back transformed values for plotting --------------------------------------------
RMR_bt_Hnum<-back_trans_function_H_num(all_mods_full$ln.RMR)
RDMC_bt_Hnum<-back_trans_function_H_num(all_mods_full$ln.RDMC)
SRL_bt_Hnum<-back_trans_function_H_num(all_mods_full$ln.SRL)
RTD_bt_Hnum<-back_trans_function_H_num(all_mods_full$ln.RTD)

#
#
#
#
#
#
#
#
#
#
#
#





# Get names  ----------------------------------------------------------# 
beta_names<-unique(Traits_all[,c("SPECIES", "H_num")])
beta_names$names<-paste(beta_names$SPECIES, beta_names$H_num, sep = "_")
beta_names<- beta_names[order(beta_names$names),]
beta_names<-beta_names$names

beta_names_noh1<-unique(Traits_noH1[,c("SPECIES", "H_num")])
beta_names_noh1$names_noh1<-paste(beta_names_noh1$SPECIES, beta_names_noh1$H_num, sep = "_")
beta_names_noh1<- beta_names_noh1[order(beta_names_noh1$names_noh1),]
beta_names_noh1<-beta_names_noh1$names_noh1

# Look at results ---------------------------------------------------# 

names(all_mods_full$ln.RMR)[grep("beta_sp_Hnum",names(all_mods_full$ln.RMR))]<-beta_names
View(summary(all_mods_full$ln.RMR)[["summary"]])
plot(all_mods_full$ln.RMR, pars = "beta_Hnum")
plot(all_mods_full$ln.RMR, pars = "beta_sp_Hnum")

names(all_mods_full$ln.RDMC)[grep("beta_sp_Hnum",names(all_mods_full$ln.RDMC))]<-beta_names
View(summary(all_mods_full$ln.RDMC)[["summary"]])
plot(all_mods_full$ln.RDMC, pars = "beta_Hnum")
plot(all_mods_full$ln.RDMC, pars = "beta_sp_Hnum")

names(all_mods_full$ln.SRL)[grep("beta_sp_Hnum",names(all_mods_full$ln.SRL))]<-beta_names
View(summary(all_mods_full$ln.SRL)[["summary"]])
plot(all_mods_full$ln.SRL, pars = "beta_Hnum")
plot(all_mods_full$ln.SRL, pars = "beta_sp_Hnum")

names(all_mods_full$ln.RTD)[grep("beta_sp_Hnum",names(all_mods_full$ln.RTD))]<-beta_names
View(summary(all_mods_full$ln.RTD)[["summary"]])
plot(all_mods_full$ln.RTD, pars = "beta_Hnum")
plot(all_mods_full$ln.RTD, pars = "beta_sp_Hnum")

names(all_mods_noH1$ln.SLA)[grep("beta_sp_Hnum",names(all_mods_noH1$ln.SLA))]<-beta_names_noh1
View(summary(all_mods_noH1$ln.SLA)[["summary"]])
plot(all_mods_noH1$ln.SLA, pars = "beta_Hnum")
plot(all_mods_noH1$ln.SLA, pars = "beta_sp_Hnum")

names(all_mods_noH1$ln.LDMC)[grep("beta_sp_Hnum",names(all_mods_noH1$ln.LDMC))]<-beta_names_noh1
View(summary(all_mods_noH1$ln.LDMC)[["summary"]])
plot(all_mods_noH1$ln.LDMC, pars = "beta_Hnum")
plot(all_mods_noH1$ln.LDMC, pars = "beta_sp_Hnum")

names(all_mods_noH1$RER_ln)[grep("beta_sp_Hnum",names(all_mods_noH1$RER_ln))]<-beta_names_noh1
View(summary(all_mods_noH1$RER_ln)[["summary"]])
plot(all_mods_noH1$RER_ln, pars = "beta_Hnum")
plot(all_mods_noH1$RER_ln, pars = "beta_sp_Hnum")

names(all_mods_noH1$RGR_Tot_ln)[grep("beta_sp_Hnum",names(all_mods_noH1$RGR_Tot_ln))]<-beta_names_noh1
View(summary(all_mods_noH1$RGR_Tot_ln)[["summary"]])
plot(all_mods_noH1$RGR_Tot_ln, pars = "beta_Hnum")
plot(all_mods_noH1$RGR_Tot_ln, pars = "beta_sp_Hnum")

