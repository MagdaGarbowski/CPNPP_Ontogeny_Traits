if(Sys.info()["login"] == "MagdaGarbowski")
  setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(ggplot2)
library(rstan)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

TraitData_2017<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")
TRY_data<-read.csv("Data_Raw_TRY_Seedweights/TRY_data/TRY.csv")
#FRED_data<-read.csv(file.choose())

# Trait data into format for stan - SLA and SRL only 
TraitData_2017_SLA<-TraitData_2017[TraitData_2017$trait %in% c("value.SLA_w_cots"),]
TraitData_2017_SLA<-TraitData_2017_SLA[c("H_num", "POP_ID","SPECIES","value")]
TraitData_2017_SLA<-na.omit.fun(TraitData_2017_SLA, c("POP_ID","SPECIES", "H_num", "value"))
TraitData_2017_SLA<-TraitData_2017_SLA[!TraitData_2017_SLA$H_num == "H4",]   # Drop H4 

# TraitData_2017_SRL<-TraitData_2017[c("SAMPLE_ID","H_num", "POP_ID","SPECIES","SRL")]
# TraitData_2017_SRL<-na.omit.fun(TraitData_2017_SRL, c("SAMPLE_ID","POP_ID", "SRL","SPECIES", "H_num"))
# TraitData_2017_SRL<-TraitData_2017_SRL[!TraitData_2017_SRL$H_num == "H4",]   # Drop H4 
# names(TraitData_2017_SRL)[names(TraitData_2017_SRL) == "SRL"] <- "value"

# TRY data into format for stan 
TRY_data$TraitName<-as.character(TRY_data$TraitName)
SLA_TRY<-TRY_data[grepl("leaf",TRY_data$TraitName),]
SLA_TRY$SPECIES<-as.factor(NA)

SLA_TRY$SPECIES<-ifelse(SLA_TRY$SpeciesName %in% c("Achillea  millefolium", "Achillea millefolia", "Achillea millefolium","Achillea millefolium L.","Achillea millefolium L. (s.str.)"), "ACMI", NA)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Agropyron subsecundum","Agropyron trachycaulum", "Elymus trachycaulus") & is.na(SLA_TRY$SPECIES)), "ELTR", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Artemisia tridentata") & is.na(SLA_TRY$SPECIES)), "ARTR", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Helianthus annuus", "Helianthus annuus  cv. Simple  Giant", "Helianthus annuus L." ) & is.na(SLA_TRY$SPECIES)), "HEAN", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Hesperostipa comata", "Stipa comata" ) & is.na(SLA_TRY$SPECIES)), "HECO", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Muhlenbergia porterii " ) & is.na(SLA_TRY$SPECIES)), "MUPO", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Packera multilobata" ) & is.na(SLA_TRY$SPECIES)), "PAMU", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("Plantago patagonica","Plantago purshii" ) & is.na(SLA_TRY$SPECIES)), "PLPA", SLA_TRY$SPECIES)
SLA_TRY$SPECIES<-ifelse((SLA_TRY$SpeciesName %in% c("vulpia octoflora" ) & is.na(SLA_TRY$SPECIES)), "VUOC", SLA_TRY$SPECIES)

# Make TRY data similar to CPNPP data to go into Stan model 
SLA_TRY<-SLA_TRY[SLA_TRY$UnitName %in% c("cm/g","mm2 mg-1"),]
SLA_TRY<-SLA_TRY[c("StdValue","SPECIES")]
SLA_TRY$POP_ID <- "TRY"
SLA_TRY$H_num <- "TRY" 
SLA_TRY<- SLA_TRY[c("POP_ID","SPECIES","H_num","StdValue")]
names(SLA_TRY)[names(SLA_TRY) == "StdValue"] <- "value"

SLA_wTRY<-rbind(TraitData_2017_SLA, SLA_TRY)
# SLA_TRY_splits<-split(SLA_wTRY, paste(SLA_wTRY$SPECIES))

# FRED data into format for stan
# FRED_data<-FRED_data[c("POP_ID","SPECIES","SRL")]
# FRED_data$H_num<-"FRED"
# FRED_data$POP_ID<-"FRED"
# FRED_data$SAMPLE_ID<-as.factor(seq.int(nrow(FRED_data)))
# FRED_data<- FRED_data[c("SAMPLE_ID","POP_ID","SRL","SPECIES","H_num")]
# names(FRED_data)[names(FRED_data) == "SRL"] <- "value"
# FRED_data<-na.omit.fun(FRED_data, c("SAMPLE_ID","POP_ID", "value","SPECIES", "H_num"))

# SRL_with_FRED<-rbind(TraitData_2017_SRL, FRED_data)
# SRL_with_FRED<-SRL_with_FRED[SRL_with_FRED$SPECIES %in% c("ACMI","ARTR","HEAN"),]
# SRL_with_FRED_splits<-split(SRL_with_FRED, paste(SRL_with_FRED$SPECIES))

# ---------------------------- Make data ----------------------------------------------
SLA_wTRY$value[SLA_wTRY$value < 1] <- NA 
SLA_wTRY_ln<-SLA_wTRY
SLA_wTRY_ln$value<-log(SLA_wTRY_ln$value)
mk_SLA_TRY_ln<-prep_data(SLA_wTRY_ln)

SLA_no_TRY_ln<-SLA_wTRY_ln[!SLA_wTRY_ln$POP_ID %in% c("TRY"),]
mk_SLA_no_TRY_ln<-prep_data(SLA_no_TRY_ln)

# ---------------------------- Model ----------------------------------------------
mod = stan_model("stan_models/All_Traits_Random_Effects.stan")

#SLA_TRY_mod<-mods_function_all(mk_SLA_TRY_ln, mods_function_all, mod = mod,
#                               pars = c("beta_Hnum_raw", "beta_sp_raw",
#                                  "beta_pop_raw", "beta_sp_Hnum_raw"),
#                               include = FALSE,
#                               warmup = 1000, iter = 1500,
#                               control = list(adapt_delta = 0.95))

SLA_noTRY_mod<-mods_function_all(mk_SLA_no_TRY_ln, mods_function_all, mod = mod,
                                 pars = c("beta_Hnum_raw", "beta_sp_raw",
                                          "beta_pop_raw", "beta_sp_Hnum_raw"),
                                 include = FALSE,
                                 warmup = 1000, iter = 1500,
                                 control = list(adapt_delta = 0.95))
#### Stan model with TRY 
#summary(SLA_TRY_mod)
#inter_names = levels(factor(paste(SLA_wTRY$SPECIES, SLA_wTRY$H_num, sep = "_")))
#tmp = names(SLA_TRY_mod)[grep("beta_sp_Hnum", names(SLA_TRY_mod))] = inter_names
#SLA_TRY_mod<-as.data.frame(SLA_TRY_mod)
#interactions = SLA_TRY_mod[,grep("_H[1-9]{1}|TRY", colnames(SLA_TRY_mod))]
#diff_TRY<-diff_by_group(interactions)
#diff_TRY_out<-summarize_diffs(diff_TRY, p=c(0.05, 0.5, 0.95))
#lapply(diff_TRY_out, function(y) y [apply(y,1, includes_zero),])

### Stan model without TRY 
summary(SLA_noTRY_mod)
inter_names_noTRY = levels(factor(paste(SLA_no_TRY_ln$SPECIES, SLA_no_TRY_ln$H_num, sep = "_")))
tmp_noTRY = names(SLA_noTRY_mod)[grep("beta_sp_Hnum", names(SLA_noTRY_mod))] = inter_names_noTRY
SLA_noTRY_mod<-as.data.frame(SLA_noTRY_mod)
interactions_noTRY = SLA_noTRY_mod[,grep("_H[1-9]{1}", colnames(SLA_noTRY_mod))]
diff_noTRY<-diff_by_group(interactions_noTRY)
diff_noTRY_out<-summarize_diffs(diff_noTRY, p=c(0.05, 0.5, 0.95))
SLA_NoTRY_stan<-lapply(diff_noTRY_out, function(y) y [apply(y,1, includes_zero),])

#
#
#
#### Check against lmer model - with TRY 
#lmer_SLA<-lmer((value)~ SPECIES + H_num + SPECIES * H_num + (1|POP_ID), data = SLA_wTRY_ln)
#plot(lmer_SLA)
#anova(lmer_SLA)
#mm<-emmeans::emmeans(lmer_SLA, ~SPECIES|H_num)
#emmeans::CLD(mm, Letters=letters, alpha=0.05)

### Check against lmer model - without TRY 
lmer_SLA_noTRY<-lmer((value)~ SPECIES + H_num + SPECIES * H_num + (1|POP_ID), data = SLA_no_TRY_ln)
plot(lmer_SLA_noTRY)
anova(lmer_SLA_noTRY)
mm_noTRY<-emmeans::emmeans(lmer_SLA_noTRY, ~SPECIES|H_num) 
SLA_noTRY_lmer<-emmeans::CLD(mm_noTRY, Letters=letters, alpha=0.05)

# Different results when comparing lmer and stan without TRY 
SLA_NoTRY_stan
SLA_noTRY_lmer



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
#
#
#
#
#

################# Plotting #############################

plot_sla_try<-function(df, Species){
  ggplot(df, aes(x=H_num, y = CI50)) + 
    geom_errorbar(aes(ymin = CI05, ymax = CI95), size = .3, width = 0.3)+
    geom_point()+
    labs(x = " ", y = "SRL")+
    scale_x_discrete(labels = c("10","24","42","FRED"))+
    ylim(-150, 500)+
    ggtitle(Species)+
    theme_bw()+
    theme(axis.text = element_text(size = 12))
}

plot_ACMI_SLA<-plot_sla_try(ACMI_posts, "ACMI")
plot_ARTR_SLA<-plot_sla_try(ARTR_posts, "ARTR")
plot_ELTR_SLA<-plot_sla_try(ELTR_posts, "ELTR")
plot_HEAN_SLA<-plot_sla_try(HEAN_posts, "HEAN")
plot_HECO_SLA<-plot_sla_try(HECO_posts, "HECO")
plot_PLPA_SLA<-plot_sla_try(PLPA_posts, "PLPA")
plot_PAMU_SLA<-plot_sla_try(PAMU_posts, "PAMU")
plot_VUOC_SLA<-plot_sla_try(VUOC_posts, "VUOC")

plot_ACMI_SRL<-plot_sla_try(ACMI_SRL_posts, "ACMI")
plot_ARTR_SRL<-plot_sla_try(ARTR_SRL_posts, "ARTR")
plot_HEAN_SRL<-plot_sla_try(HEAN_SRL_posts, "HEAN")

pdf("Output/Figures/SLA_TRY.pdf", height = 8, width = 6)
grid.arrange( plot_ACMI_SLA,plot_ARTR_SLA,plot_PAMU_SLA, plot_HEAN_SLA, plot_PLPA_SLA, plot_ELTR_SLA, plot_HECO_SLA, plot_VUOC_SLA, ncol =2 )
dev.off()

pdf("Output/Figures/SRL_FRED.pdf", height = 6, width = 3)
grid.arrange( plot_ACMI_SRL, plot_ARTR_SRL, plot_HEAN_SRL )
dev.off()



