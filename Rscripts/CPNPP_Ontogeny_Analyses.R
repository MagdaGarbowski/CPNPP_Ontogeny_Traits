
# 2017 Trait data - Traits through time   

# (1) Clean up code - make more functions? 
# (2) Remove SLA and LDMC at H1 from graphs and possibly from models as well? 
# (3) Figure out divergence issues in models - which ones are they? How can they be fixed? 
# (4) Same model is being used for all traits - is this "right"? 

setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/") 
source("Rscripts/Functions/Functions_Analyses.R")


library(lmerTest)
library(rstanarm)
library(rstan)
library(plyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidybayes)
library(ggridges)
library(ggmcmc)
library(bayesplot)
library(cowplot)

# Data - from "Traits_CommitteeMeeting_2019_DataMgmt.R"

# Full Data
SpeciesData<-read.csv("Data_Generated/TraitData_2017.csv")  
SpeciesData$H_num<-as.factor(SpeciesData$H_num)

# Growth Rates
SpeciesData_Growthrate<-read.csv("Data_Generated/TraitData_GrowthRates_2017.csv")
SpeciesData_Growthrate<-SpeciesData_Growthrate[!(SpeciesData_Growthrate$H_num) %in% c("1"),] # Remove since no values at H_num1
SpeciesData_Growthrate$H_num<-as.factor(SpeciesData_Growthrate$H_num)

# Plasticity 
SpeciesData_Plasticity<-read.csv("Data_Generated/TraitData_Plasticity_2017.csv")
Trait_splits_PI = split(SpeciesData_Plasticity, paste(SpeciesData_Plasticity$Trait))

# Plasticity - no h1 values 
SpeciesData_Plasticity_10day<-read.csv("Data_Generated/TraitData_Plasticity_2017_10day.csv")
Trait_splits_PI_10day = split(SpeciesData_Plasticity_10day, paste(SpeciesData_Plasticity_10day$Trait))

# Relative Differences
Relative_diff<- read.csv("Data_Generated/TraitData_RelativeDiff_2017.csv")
Relative_diff<-Relative_diff[!(Relative_diff$H_num) %in% c("H1"),] # Remove since no values at H_num1
Relative_diff$H_num<-as.factor(Relative_diff$H_num)

# Population avgs. to check models with POP_ID as random effects and those with pop avgs 
Trait_avg_pop<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")
Trait_avg_pop_splits<-split(Trait_avg_pop, paste(Trait_avg_pop$trait))
#
#
# ---------------------------------- Rstanarm models -----------------------------
#
# Traits, growth rates, relative differences 
Trait_rstanarm <- function (df, var) {
  stan_lmer(df[[var]] ~ 0 + SPECIES + H_num + SPECIES*H_num + (1|SPECIES/POP_ID), 
            data = df,
            prior_intercept = student_t(3,0,30), 
            adapt_delta = 0.999)
}
# Plasticity 
Plasticity_rstanarm <- function (df) {
  stan_aov(dist ~ 0 + Species, 
           data = df,
           prior = R2(0.5))
}

# ---------------------------Run models, get datasets------------------ 
# ---------------------------- Traits through time --------------------

Trait_rstanarm_SLA<-Trait_rstanarm(SpeciesData, "SLA")
Trait_rstanarm_LDMC<-Trait_rstanarm(SpeciesData, "LDMC")
Trait_rstanarm_RMR<-Trait_rstanarm(SpeciesData, "RMR")
Trait_rstanarm_RDMC<-Trait_rstanarm(SpeciesData, "RDMC")
Trait_rstanarm_RTD<-Trait_rstanarm(SpeciesData, "RTD")
Trait_rstanarm_SRL<-Trait_rstanarm(SpeciesData, "SRL")

Trait_rstanarm_mods<-list(Trait_rstanarm_SLA = Trait_rstanarm_SLA, Trait_rstanarm_LDMC = Trait_rstanarm_LDMC, Trait_rstanarm_RMR = Trait_rstanarm_RMR, 
                          Trait_rstanarm_RDMC = Trait_rstanarm_RDMC, Trait_rstanarm_RTD = Trait_rstanarm_RTD, Trait_rstanarm_SRL=Trait_rstanarm_SRL) 

mods_out<-lapply(Trait_rstanarm_mods, df_function)
median_out<-lapply(mods_out, Median_CI_calc_function)

# COME BACK - How to make this a function that passes the name of the list to a column? 

SLA<-as.data.frame(median_out[[1]])
SLA$trait<-"SLA"

LDMC<-as.data.frame(median_out[[2]])
LDMC$trait<-"LDMC"

RMR<-as.data.frame(median_out[[3]])
RMR$trait<-"RMR"

RDMC<-as.data.frame(median_out[[4]])
RDMC$trait<-"RDMC"

RTD<-as.data.frame(median_out[[5]])
RTD$trait<-"RTD"

SRL<-as.data.frame(median_out[[6]])
SRL$trait<-"SRL"

median_out_df<-rbind(SLA, LDMC,RMR, RDMC, RTD, SRL)

# ---------------------- Growth rates-------------------------

Trait_rstanarm_RER<-Trait_rstanarm(SpeciesData_Growthrate, "RER_ln")
Trait_rstanarm_RGR_tot<-Trait_rstanarm(SpeciesData_Growthrate, "RGR_Tot_ln")
Trait_rstanarm_RGR_below<-Trait_rstanarm(SpeciesData_Growthrate, "RGR_Below_ln")
Trait_rstanarm_RGR_above<-Trait_rstanarm(SpeciesData_Growthrate, "RGR_Above_ln") # Divergent transition 

RGR_rstanarm_mods<-list(Trait_rstanarm_RER = Trait_rstanarm_RER, Trait_rstanarm_RGR_tot = Trait_rstanarm_RGR_tot, Trait_rstanarm_RGR_below = Trait_rstanarm_RGR_below, 
                        Trait_rstanarm_RGR_above = Trait_rstanarm_RGR_above) 

RGR_mods_out<-lapply(RGR_rstanarm_mods, df_function)
RGR_median_out<-lapply(RGR_mods_out, Median_CI_calc_function)

# COME BACK - How to make this a function that passes the name of the list to a column? 

RER<-as.data.frame(RGR_median_out[[1]])
RER$trait<-"RER"

RGR_Tot<-as.data.frame(RGR_median_out[[2]])
RGR_Tot$trait<-"RGR_Tot"

RGR_Below<-as.data.frame(RGR_median_out[[3]])
RGR_Below$trait<-"RGR_Below"

RGR_Above<-as.data.frame(RGR_median_out[[4]])
RGR_Above$trait<-"RGR_Above"

median_out_RGR_df<-rbind(RER, RGR_Tot,RGR_Below, RGR_Above)
levels(median_out_RGR_df$Harvest)[levels(median_out_RGR_df$Harvest) %in% c("h1")]<-c("h2")


# ------------------------------ Relative differences -----------------------------

Trait_rstanarm_rel_SLA2<-Trait_rstanarm(Relative_diff, "rel_SLA2")
Trait_rstanarm_rel_LDMC2<-Trait_rstanarm(Relative_diff, "rel_LDMC2")
Trait_rstanarm_rel_RMR2<-Trait_rstanarm(Relative_diff, "rel_RMR2")
Trait_rstanarm_rel_RDMC2<-Trait_rstanarm(Relative_diff, "rel_RDMC2")
Trait_rstanarm_rel_RTD2<-Trait_rstanarm(Relative_diff, "rel_RTD2")
Trait_rstanarm_rel_SRL2<-Trait_rstanarm(Relative_diff, "rel_SRL2")

Trait_rstanarm_mods_RD<-list(Trait_rstanarm_rel_SLA2 = Trait_rstanarm_rel_SLA2, Trait_rstanarm_rel_LDMC2 = Trait_rstanarm_rel_LDMC2, Trait_rstanarm_rel_RMR2 = Trait_rstanarm_rel_RMR2, 
                             Trait_rstanarm_rel_RDMC2 = Trait_rstanarm_rel_RDMC2, Trait_rstanarm_rel_RTD2 = Trait_rstanarm_rel_RTD2, Trait_rstanarm_rel_SRL2=Trait_rstanarm_rel_SRL2) 

RD_mods_out<-lapply(Trait_rstanarm_mods_RD, df_function)
RGR_median_out<-lapply(RD_mods_out, Median_CI_calc_function)

# COME BACK - How to make this a function that passes the name of the list to a column? 

rel_SLA2<-as.data.frame(RGR_median_out[[1]])
rel_SLA2$trait<-"rel_SLA2"

rel_LDMC2<-as.data.frame(RGR_median_out[[2]])
rel_LDMC2$trait<-"rel_LDMC2"

rel_RMR2<-as.data.frame(RGR_median_out[[3]])
rel_RMR2$trait<-"rel_RMR2"

rel_RDMC2<-as.data.frame(RGR_median_out[[4]])
rel_RDMC2$trait<-"rel_RDMC2"

rel_RTD2<-as.data.frame(RGR_median_out[[5]])
rel_RTD2$trait<-"rel_RTD2"

rel_SRL2<-as.data.frame(RGR_median_out[[6]])
rel_SRL2$trait<-"rel_SRL2"

median_out_RD_df<-rbind(rel_SLA2, rel_LDMC2,rel_RMR2, rel_RDMC2, rel_RTD2, rel_SRL2)
levels(median_out_RD_df$Harvest)[levels(median_out_RD_df$Harvest) %in% c("H_num1")]<-c("H_num2")

# --------------------------- Plasticity -------------------------
# Where is Species = VUOC? 

PI_mods<-lapply(Trait_splits_PI, Plasticity_rstanarm)
PI_mods_out<-lapply(PI_mods, dat_frame_function_PI)
PI_median_out<-PI_mods_out

LDMC<-as.data.frame(PI_median_out[[1]])
LDMC$trait<-"LDMC"

RDMC<-as.data.frame(PI_median_out[[2]])
RDMC$trait<-"RDMC"

RER_ln<-as.data.frame(PI_median_out[[3]])
RER_ln$trait<-"RER"

RGR_Above_ln<-as.data.frame(PI_median_out[[4]])
RGR_Above_ln$trait<-"RGR_Above_ln"

RGR_Below_ln<-as.data.frame(PI_median_out[[5]])
RGR_Below_ln$trait<-"RGR_Below_ln"

RGR_Tot_ln<-as.data.frame(PI_median_out[[6]])
RGR_Tot_ln$trait<-"RGR_Tot_ln"

RMR<-as.data.frame(PI_median_out[[7]])
RMR$trait<-"RMR"

RTD<-as.data.frame(PI_median_out[[8]])
RTD$trait<-"RTD"

SLA<-as.data.frame(PI_median_out[[9]])
SLA$trait<-"SLA"

SRL<-as.data.frame(PI_median_out[[10]])
SRL$trait<-"SRL"


plasticity_df<-rbind(LDMC, RDMC, RER_ln, RGR_Above_ln, RGR_Below_ln, RGR_Tot_ln, RMR, RTD, SLA, SRL)
plasticity_df<-plasticity_df[!plasticity_df$Species %in% c("(Intercept)"),]

### Plasticity of traits - without H1
# PI_mods_10day<-lapply(Trait_splits_PI_10day, PI_rstanarm)
# posterior_medians_intervals_PI_df_10day<-lapply(PI_mods_10day, posterior_medians_intervals_PI)
# Medians_intervals_PI_10day<-dat_frame_function_PI(posterior_medians_intervals_PI_df_10day)
# plasticity_df_10day<-Medians_intervals_PI_10day


### Write .csv

write.csv(median_out_df, "Data_Generated/SpeciesxTimeMedians.csv")
write.csv(median_out_RGR_df, "Data_Generated/RelativegrowthrateMedians.csv")
write.csv(median_out_RD_df, "Data_Generated/RelativeDifferencesMedians.csv")
write.csv(plasticity_df, "Data_Generated/PlasticityMedians.csv")





