if(Sys.info()["login"] == "MagdaGarbowski")
  setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(ggplot2)
library(rstan)
library(rstanarm)
library(plyr)
library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidybayes)
library(ggridges)
library(ggmcmc)

TraitData_2017<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")
TRY_data<-read.csv("Data_Raw_TRY_Seedweights/TRY_data/TRY.csv")
FRED_data<-read.csv(file.choose())

# Trait data into format for stan - SLA only 
TraitData_2017_SLA<-TraitData_2017[TraitData_2017$trait %in% c("value.SLA_w_cots"),]
TraitData_2017_SLA<-TraitData_2017_SLA[c("H_num", "POP_ID","SPECIES","value")]
TraitData_2017_SLA<-na.omit.fun(TraitData_2017_SLA, c("POP_ID","SPECIES", "H_num", "value"))
TraitData_2017_SLA<-TraitData_2017_SLA[!TraitData_2017_SLA$H_num == "H4",]   # Drop H4 

TraitData_2017_SRL<-TraitData_2017[c("SAMPLE_ID","H_num", "POP_ID","SPECIES","SRL")]
TraitData_2017_SRL<-na.omit.fun(TraitData_2017_SRL, c("SAMPLE_ID","POP_ID", "SRL","SPECIES", "H_num"))
TraitData_2017_SRL<-TraitData_2017_SRL[!TraitData_2017_SRL$H_num == "H4",]   # Drop H4 
names(TraitData_2017_SRL)[names(TraitData_2017_SRL) == "SRL"] <- "value"

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
SLA_TRY_splits<-split(SLA_wTRY, paste(SLA_wTRY$SPECIES))

# FRED data into format for stan
FRED_data<-FRED_data[c("POP_ID","SPECIES","SRL")]
FRED_data$H_num<-"FRED"
FRED_data$POP_ID<-"FRED"
FRED_data$SAMPLE_ID<-as.factor(seq.int(nrow(FRED_data)))
FRED_data<- FRED_data[c("SAMPLE_ID","POP_ID","SRL","SPECIES","H_num")]
names(FRED_data)[names(FRED_data) == "SRL"] <- "value"
FRED_data<-na.omit.fun(FRED_data, c("SAMPLE_ID","POP_ID", "value","SPECIES", "H_num"))

SRL_with_FRED<-rbind(TraitData_2017_SRL, FRED_data)
SRL_with_FRED<-SRL_with_FRED[SRL_with_FRED$SPECIES %in% c("ACMI","ARTR","HEAN"),]
SRL_with_FRED_splits<-split(SRL_with_FRED, paste(SRL_with_FRED$SPECIES))

# ---------------------------- Model ----------------------------------------------
SLA_wTRY$value<-log(SLA_wTRY$value)
mk_SLA_TRY<-prep_data(SLA_wTRY)

SLA_TRY_mod<-mods_function_all(mk_SLA_TRY, mods_function_all, mod = mod,
                         pars = c("beta_Hnum_raw", "beta_sp_raw",
                                  "beta_pop_raw", "beta_sp_Hnum_raw"),
                         include = FALSE,
                         warmup = 1000, iter = 1500,
                         control = list(adapt_delta = 0.95))

summary(SLA_TRY_mod)
plot(SLA_TRY_mod, pars = "beta_Hnum")
inter_names = levels(factor(paste(SLA_wTRY$SPECIES, SLA_wTRY$H_num, sep = "_")))




Trait_rstanarm_ind<-function (df){
  stan_lmer(value ~ 0 + H_num + (1|POP_ID),
  data=df, 
  adapt_delta = 0.95)
}

# -----------------------------run mods ---------------------------

all_mods_SLA_TRY<-lapply(SLA_TRY_splits, Trait_rstanarm_ind)
all_mods_SRL_FRED<-lapply(SRL_with_FRED_splits, Trait_rstanarm_ind)


# ------------------------- Get values ----------------------------------------------

posts<-function(df, sps){
  postsdf<-as.data.frame(summary(df,
                                 pars = c("H_numH1","H_numH2","H_numH3","H_numFRED"), 
                                 probs = c(.05,.5,.95)))
  postsdf$species<-sps
  postsdf$H_num<-rownames(postsdf)
  colnames(postsdf)<-c("mean", "mcse","sd","CI05","CI50","CI95", "n_eff", "Rhat", "species","H_num")
  return(postsdf)
}

ACMI_posts<-posts(all_mods_SLA_TRY$ACMI, "ACMI")
ARTR_posts<-posts(all_mods_SLA_TRY$ARTR, "ARTR")
ELTR_posts<-posts(all_mods_SLA_TRY$ELTR, "ELTR")
HEAN_posts<-posts(all_mods_SLA_TRY$HEAN, "HEAN")
HECO_posts<-posts(all_mods_SLA_TRY$HECO, "HECO")
PLPA_posts<-posts(all_mods_SLA_TRY$PLPA, "PLPA")
PAMU_posts<-posts(all_mods_SLA_TRY$PAMU, "PAMU")
VUOC_posts<-posts(all_mods_SLA_TRY$VUOC, "VUOC")

ACMI_SRL_posts<-posts(all_mods_SRL_FRED$ACMI, "ACMI")
ARTR_SRL_posts<-posts(all_mods_SRL_FRED$ARTR, "ARTR")
HEAN_SRL_posts<-posts(all_mods_SRL_FRED$HEAN, "HEAN")

# Plotting 

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

# Plotting 

# Specific leaf area 
SLA_plot_ACMI<-density_plot(posts_ACMI, c(-10, 110), c(0,2.6), "Days","") + annotate("text", x = 90, y = 11, label = "ACMI (F)", size = 4)
SLA_plot_ARTR<-density_plot(posts_ARTR, c(-10, 110), c(0,2.6), "Days", "") + annotate("text", x = 90, y = 11, label = "ARTR (S)", size = 4)
SLA_plot_ELTR<-density_plot(posts_ELTR, c(-10, 150), c(0,2.2), "Days", "") + annotate("text", x = 110, y = 11, label = "ELTR (G)", size = 4)
SLA_plot_HEAN<-density_plot(posts_HEAN, c(-10, 110), c(0,2.6), "Days", "") + annotate("text", x = 90, y = 11, label = "HEAN (AF)", size = 4)
SLA_plot_HECO<-density_plot(posts_HECO, c(-10, 150), c(0,2.2), "Days", "") + annotate("text", x = 110, y = 11, label = "HECO (G)", size = 4)
SLA_plot_PAMU<-density_plot(posts_PAMU, c(-10, 110), c(0,2.6), "Days", "") + annotate("text", x = 90, y = 11, label = "PAMU (F)", size = 4)
SLA_plot_PLPA<-density_plot(posts_PLPA, c(-10, 180), c(0,2.6), "Days", "") + annotate("text", x = 150, y = 11, label = "PLPA (AF)", size = 4)
SLA_plot_VUOC<-density_plot(posts_VUOC, c(-10, 150), c(0,2.2), "Days", "") + annotate("text", x = 110, y = 11, label = "VUOC (AG)", size = 4)
SLA_plot_MACA<-density_plot(posts_MACA, c(-10, 110), c(0,5.5), "Days", "") + annotate("text", x = 90, y = 11, label = "MACA (BF)", size = 4)
SLA_plot_HEVI<-density_plot(posts_HEVI, c(-10, 110), c(0,5.5), "Days", "") + annotate("text", x = 90, y = 11, label = "HEVI (AF)", size = 4)
SLA_plot_MUPO<-density_plot(posts_MUPO, c(-10, 150), c(0,3.5), "Days", "") + annotate("text", x = 90, y = 11, label = "MUPO (G)", size = 4)

RMR_plot_ACMI<-density_plot(posts_RMR$ACMI, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ACMI (F)", size = 4)
RMR_plot_ARTR<-density_plot(posts_RMR$ARTR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ARTR (S)", size = 4)
RMR_plot_ELTR<-density_plot(posts_RMR$ELTR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ELTR (G)", size = 4)
RMR_plot_HEAN<-density_plot(posts_RMR$HEAN, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEAN (AF)", size = 4)
RMR_plot_HECO<-density_plot(posts_RMR$HECO, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HECO (G)", size = 4)
RMR_plot_PAMU<-density_plot(posts_RMR$PAMU, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PAMU (F)", size = 4)
RMR_plot_PLPA<-density_plot(posts_RMR$PLPA, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PLPA (AF)", size = 4)
RMR_plot_VUOC<-density_plot(posts_RMR$VUOC, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "VUOC (AG)", size = 4)
RMR_plot_MACA<-density_plot(posts_RMR$MACA, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MACA (BF)", size = 4)
RMR_plot_HEVI<-density_plot(posts_RMR$HEVI, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEVI (F)", size = 4)
RMR_plot_MUPO<-density_plot(posts_RMR$MUPO, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MUPO (G)", size = 4)

# LDMC from different functions - how to pull them out? 
LDMC_plot_ACMI<-density_plot(posts_ACMI_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "ACMI (F)", size = 4)
LDMC_plot_ARTR<-density_plot(posts_ARTR_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "ARTR (S)", size = 4)
LDMC_plot_ELTR<-density_plot(posts_ELTR_LDMC, c(-10, 50), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ELTR (G)", size = 4)
LDMC_plot_HEAN<-density_plot(posts_HEAN_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "HEAN (AF)", size = 4)
LDMC_plot_HECO<-density_plot(posts_HECO_LDMC, c(-10, 50), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HECO (G)", size = 4)
LDMC_plot_PAMU<-density_plot(posts_PAMU_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "PAMU (F)", size = 4)
LDMC_plot_PLPA<-density_plot(posts_PLPA_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "PLPA (AF)", size = 4)
LDMC_plot_VUOC<-density_plot(posts_VUOC_LDMC, c(-10, 50), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "VUOC (AG)", size = 4)
LDMC_plot_MACA<-density_plot(posts_MACA_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "MACA (BF)", size = 4)
LDMC_plot_HEVI<-density_plot(posts_HEVI_LDMC, c(-10, 50), c(0,5.5), "", "") + annotate("text", x = 0, y = 12, label = "HEVI (F)", size = 4)
LDMC_plot_MUPO<-density_plot(posts_MUPO_LDMC, c(-10, 50), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MUPO (G)", size = 4)

RDMC_plot_ACMI<-density_plot(posts_RDMC$ACMI, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ACMI (F)", size = 4)
RDMC_plot_ARTR<-density_plot(posts_RDMC$ARTR, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ARTR (S)", size = 4)
RDMC_plot_ELTR<-density_plot(posts_RDMC$ELTR, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ELTR (G)", size = 4)
RDMC_plot_HEAN<-density_plot(posts_RDMC$HEAN, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEAN (AF)", size = 4)
RDMC_plot_HECO<-density_plot(posts_RDMC$HECO, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HECO (G)", size = 4)
RDMC_plot_PAMU<-density_plot(posts_RDMC$PAMU, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PAMU (F)", size = 4)
RDMC_plot_PLPA<-density_plot(posts_RDMC$PLPA, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PLPA (AF)", size = 4)
RDMC_plot_VUOC<-density_plot(posts_RDMC$VUOC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "VUOC (AG)", size = 4)
RDMC_plot_MACA<-density_plot(posts_RDMC$MACA, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MACA (BF)", size = 4)
RDMC_plot_HEVI<-density_plot(posts_RDMC$HEVI, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEVI (F)", size = 4)
RDMC_plot_MUPO<-density_plot(posts_RDMC$MUPO, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MUPO (G)", size = 4)

RTD_plot_ACMI<-density_plot(posts_RTD$ACMI, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "ACMI (F)", size = 4)
RTD_plot_ARTR<-density_plot(posts_RTD$ARTR, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "ARTR (S)", size = 4)
RTD_plot_ELTR<-density_plot(posts_RTD$ELTR, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "ELTR (G)", size = 4)
RTD_plot_HEAN<-density_plot(posts_RTD$HEAN, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "HEAN (AF)", size = 4)
RTD_plot_HECO<-density_plot(posts_RTD$HECO, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "HECO (G)", size = 4)
RTD_plot_PAMU<-density_plot(posts_RTD$PAMU, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "PAMU (F)", size = 4)
RTD_plot_PLPA<-density_plot(posts_RTD$PLPA, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "PLPA (AF)", size = 4)
RTD_plot_VUOC<-density_plot(posts_RTD$VUOC, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "VUOC (AG)", size = 4)
RTD_plot_MACA<-density_plot(posts_RTD$MACA, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "MACA (BF)", size = 4)
RTD_plot_HEVI<-density_plot(posts_RTD$HEVI, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "HEVI (F)", size = 4)
RTD_plot_MUPO<-density_plot(posts_RTD$MUPO, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "MUPO (G)", size = 4)

SRL_plot_ACMI<-density_plot(posts_SRL$ACMI, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "ACMI (F)", size = 4)
SRL_plot_ARTR<-density_plot(posts_SRL$ARTR, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "ARTR (S)", size = 4)
SRL_plot_ELTR<-density_plot(posts_SRL$ELTR, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "ELTR (G)", size = 4)
SRL_plot_HEAN<-density_plot(posts_SRL$HEAN, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "HEAN (AF)", size = 4)
SRL_plot_HECO<-density_plot(posts_SRL$HECO, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "HECO (G)", size = 4)
SRL_plot_PAMU<-density_plot(posts_SRL$PAMU, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "PAMU (F)", size = 4)
SRL_plot_PLPA<-density_plot(posts_SRL$PLPA, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "PLPA (AF)", size = 4)
SRL_plot_VUOC<-density_plot(posts_SRL$VUOC, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "VUOC (AG)", size = 4)
SRL_plot_MACA<-density_plot(posts_SRL$MACA, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "MACA (BF)", size = 4)
SRL_plot_HEVI<-density_plot(posts_SRL$HEVI, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "HEVI (F)", size = 4)
SRL_plot_MUPO<-density_plot(posts_SRL$MUPO, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "MUPO (G)", size = 4)


pdf("Output/Figures/Density_SLA_TRY.pdf", height = 8, width = 12)
sla_plot<-grid.arrange(SLA_plot_ELTR, SLA_plot_ACMI,SLA_plot_HEAN,SLA_plot_HECO,SLA_plot_ARTR,SLA_plot_MACA,SLA_plot_MUPO, SLA_plot_HEVI, SLA_plot_PLPA,SLA_plot_VUOC, SLA_plot_PAMU, ncol = 3)
dev.off()

pdf("Output/Figures/Density_RMR_TRY.pdf", height = 8, width = 12)
RMR_plot<-grid.arrange(RMR_plot_ELTR, RMR_plot_ACMI,RMR_plot_HEAN,RMR_plot_HECO,RMR_plot_ARTR,RMR_plot_MACA,RMR_plot_MUPO, RMR_plot_HEVI, RMR_plot_PLPA,RMR_plot_VUOC, RMR_plot_PAMU, ncol = 3)
dev.off()

pdf("Output/Figures/Density_LDMC_TRY.pdf", height = 8, width = 12)
LDMC_plot<-grid.arrange(LDMC_plot_ELTR, LDMC_plot_ACMI,LDMC_plot_HEAN,LDMC_plot_HECO,LDMC_plot_ARTR,LDMC_plot_MACA,LDMC_plot_MUPO, LDMC_plot_HEVI, LDMC_plot_PLPA,LDMC_plot_VUOC, LDMC_plot_PAMU, ncol = 3)
dev.off()

pdf("Output/Figures/Density_RDMC_TRY.pdf", height = 8, width = 12)
RDMC_plot<-grid.arrange(RDMC_plot_ELTR, RDMC_plot_ACMI,RDMC_plot_HEAN,RDMC_plot_HECO,RDMC_plot_ARTR,RDMC_plot_MACA,RDMC_plot_MUPO, RDMC_plot_HEVI, RDMC_plot_PLPA,RDMC_plot_VUOC, RDMC_plot_PAMU, ncol = 3)
dev.off()

pdf("Output/Figures/Density_RTD_TRY.pdf", height = 8, width = 12)
RTD_plot<-grid.arrange(RTD_plot_ELTR, RTD_plot_ACMI,RTD_plot_HEAN,RTD_plot_HECO,RTD_plot_ARTR,RTD_plot_MACA,RTD_plot_MUPO, RTD_plot_HEVI, RTD_plot_PLPA,RTD_plot_VUOC, RTD_plot_PAMU, ncol = 3)
dev.off()

pdf("Output/Figures/Density_SRL_TRY.pdf", height = 8, width = 12)
SRL_plot<-grid.arrange(SRL_plot_ELTR, SRL_plot_ACMI,SRL_plot_HEAN,SRL_plot_HECO,SRL_plot_ARTR,SRL_plot_MACA,SRL_plot_MUPO, SRL_plot_HEVI, SRL_plot_PLPA,SRL_plot_VUOC, SRL_plot_PAMU, ncol = 3)
dev.off()


pdf("Output/Figures/Density_Grasses.pdf", height = 10, width = 20)
grid.arrange(arrangeGrob(textGrob("SLA"), textGrob("LDMC"),textGrob("RMR"), textGrob("SRL"),textGrob("RDMC"), textGrob("RTD"),
             SLA_plot_ELTR, LDMC_plot_ELTR, RMR_plot_ELTR, SRL_plot_ELTR, RDMC_plot_ELTR, RTD_plot_ELTR,
             SLA_plot_HECO, LDMC_plot_HECO, RMR_plot_HECO, SRL_plot_HECO, RDMC_plot_HECO, RTD_plot_HECO,
             SLA_plot_MUPO, LDMC_plot_MUPO, RMR_plot_MUPO, SRL_plot_MUPO, RDMC_plot_MUPO, RTD_plot_MUPO,
             SLA_plot_VUOC, LDMC_plot_VUOC, RMR_plot_VUOC, SRL_plot_VUOC, RDMC_plot_VUOC, RTD_plot_VUOC, ncol = 6,
             heights = c(1,4,4,4,4)))
dev.off()



pdf("Output/Figures/Density_forbs.pdf", height = 16, width = 20)
grid.arrange(arrangeGrob(textGrob("SLA"), textGrob("LDMC"),textGrob("RMR"), textGrob("SRL"),textGrob("RDMC"), textGrob("RTD"),
            SLA_plot_ACMI, LDMC_plot_ACMI, RMR_plot_ACMI, SRL_plot_ACMI, RDMC_plot_ACMI, RTD_plot_ACMI,
             SLA_plot_ARTR, LDMC_plot_ARTR, RMR_plot_ARTR, SRL_plot_ARTR, RDMC_plot_ARTR, RTD_plot_ARTR,
             SLA_plot_HEVI, LDMC_plot_HEVI, RMR_plot_HEVI, SRL_plot_HEVI, RDMC_plot_HEVI, RTD_plot_HEVI,
             SLA_plot_PAMU, LDMC_plot_PAMU, RMR_plot_PAMU, SRL_plot_PAMU, RDMC_plot_PAMU, RTD_plot_PAMU,
             SLA_plot_HEAN, LDMC_plot_HEAN, RMR_plot_HEAN, SRL_plot_HEAN, RDMC_plot_HEAN, RTD_plot_HEAN,
             SLA_plot_MACA, LDMC_plot_MACA, RMR_plot_MACA, SRL_plot_MACA, RDMC_plot_MACA, RTD_plot_MACA,
             SLA_plot_PLPA, LDMC_plot_PLPA, RMR_plot_PLPA, SRL_plot_PLPA, RDMC_plot_PLPA, RTD_plot_PLPA, ncol = 6,
            heights = c(1,4,4,4,4,4,4,4)))
dev.off()

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
# RMR model has issues. All others okay? 
# lmer models to see if would be different from 
lmer_SLA<-lmer(sqrt(value)~ SPECIES + H_num + SPECIES * H_num + (1|POP_ID), data = SLA_wTRY)
anova(lmer_SLA)
mm<-emmeans::emmeans(lmer_SLA, ~SPECIES|H_num)
emmeans::CLD(mm, Letters=letters, alpha=0.05)

TraitData_2017_LDMC<-TraitData_2017[c("SAMPLE_ID","H_num", "POP_ID","SPECIES","LDMC")]
TraitData_2017_LDMC<-na.omit.fun(TraitData_2017_LDMC, c("SAMPLE_ID","POP_ID", "LDMC","SPECIES", "H_num"))
TraitData_2017_LDMC<-TraitData_2017_LDMC[!TraitData_2017_LDMC$H_num == "H4",]   # Drop H4 
names(TraitData_2017_LDMC)[names(TraitData_2017_LDMC) == "LDMC"] <- "value"
