setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
library(lme4)
library(emmeans)

# Doing regular old stats, like I know how. All data together. Functional group as a predictor.

TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")

# Growth Rates
SpeciesData_Growthrate<-read.csv("Data_Generated/TraitData_GrowthRates_2017.csv")
SpeciesData_Growthrate$H_num<-as.character(SpeciesData_Growthrate$H_num)
SpeciesData_Growthrate<-SpeciesData_Growthrate[!(SpeciesData_Growthrate$H_num) %in% c("1"),] # Remove since no values at H_num1

TraitData_Pop_Avg_2017$H_num<-as.character(as.factor(TraitData_Pop_Avg_2017$H_num))
TraitData_Pop_Avg_2017$trait<-as.character(as.factor(TraitData_Pop_Avg_2017$trait))

TraitData_Pop_Avg_2017$value<-ifelse((TraitData_Pop_Avg_2017$GrowthForm %in% c("FORB","SHRUB") & TraitData_Pop_Avg_2017$H_num == "H1" & TraitData_Pop_Avg_2017$trait == "SLA"), NA, TraitData_Pop_Avg_2017$value) # Add NA for H1 SLA for forb species 
TraitData_Pop_Avg_2017$value<-ifelse((TraitData_Pop_Avg_2017$GrowthForm %in% c("FORB","SHRUB") & TraitData_Pop_Avg_2017$H_num == "H1" & TraitData_Pop_Avg_2017$trait == "LDMC"), NA, TraitData_Pop_Avg_2017$value) # Add NA for H1 SLA for forb species 

TraitData_Pop_FullHarvests<-as.data.frame(TraitData_Pop_Avg_2017[TraitData_Pop_Avg_2017$trait %in% c("RMR", "SRL", "RDMC", "RTD"),])
TraitData_Pop_noH1<-as.data.frame(TraitData_Pop_Avg_2017[TraitData_Pop_Avg_2017$trait %in% c("SLA", "LDMC"),])
TraitData_Pop_noH1<-TraitData_Pop_noH1[!TraitData_Pop_noH1$H_num %in% c("H1"),]

TraitData_Pop_FullHarvests_splits<-split(TraitData_Pop_FullHarvests, TraitData_Pop_FullHarvests$trait)
TraitData_Pop_noH1_splits<-split(TraitData_Pop_noH1, TraitData_Pop_noH1$trait)
                      
lmer_function<-function(df){
  mod<-lmer(log(value) ~ H_num + GrowthForm + H_num*GrowthForm + (1|SPECIES), data = df)
  summary<-anova(mod)
}

contrasts_lmer_function<-function(df){
  mod<-lmer(log(value) ~ H_num + GrowthForm + H_num*GrowthForm + (1|SPECIES), data = df)
  means<-emmeans(mod, pairwise ~ GrowthForm)
  contrasts<-as.data.frame(CLD(means, Letters=letters, alpha=0.05))
}

lmer_out_RDMC_RMR_RTD_SRL<-lapply(TraitData_Pop_FullHarvests_splits, lmer_function)
lmer_out_SLA_LDMC<-lapply(TraitData_Pop_noH1_splits, lmer_function)
lmer_out_RGR<-anova(lmer(RGR_Tot_ln ~ H_num + GrowthForm + H_num*GrowthForm + (1|SPECIES), data = SpeciesData_Growthrate))
lmer_out_RER<-anova(lmer(RER_ln ~ H_num + GrowthForm + H_num*GrowthForm + (1|SPECIES), data = SpeciesData_Growthrate))

get_p_f_df<-function(df){
  p = df$"Pr(>F)"
  f = df$"F value"
  p = format(round(p,4), nsmall = 3)
  f = format(round(f, 2), nsmall = 2)
  p_f = paste(f,p, sep = ", ")
  df_num = df$NumDF
  df_den = df$DenDF
  df_num = format(round(df_num, 0), nsmall = 0)
  df_den = format(round(df_den, 0), nsmall = 0)
  df_tot = paste(df_num, df_den, sep = ", ")
  dat = as.data.frame(rbind(p_f, df_tot))
  dat
}

Root_ts<-lapply(lmer_out_RDMC_RMR_RTD_SRL, get_p_f_df)
Leaf_ts<-lapply(lmer_out_SLA_LDMC, get_p_f_df)
P_F_out<-do.call("rbind", c(Root_ts,Leaf_ts))

RGR<-get_p_f_df(lmer_out_RGR)
RER<-get_p_f_df(lmer_out_RER)

PF_out_growth<-rbind(RGR,RER)
P_F_out_2<-rbind(P_F_out,PF_out_growth)

contrasts_out_all<-lapply(TraitData_Pop_FullHarvests_splits, contrasts_lmer_function)
contrasts_out_SLA_LDMC<-lapply(TraitData_Pop_noH1_splits, contrasts_lmer_function)

write.csv(P_F_out_2, "~/Downloads/P_Fvalues_out.csv")
write.csv(contrasts_out, "~/Downloads/SEM-EDS_contrasts_out.csv")



