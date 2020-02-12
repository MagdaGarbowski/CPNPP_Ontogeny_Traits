setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
source("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Rscripts/Functions/Functions_TRY_Density.R")
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

SpeciesData<-read.csv("Data_Generated/TraitData_2017.csv")

SpeciesData$RMR_100<-SpeciesData$RMR*100
SpeciesData$RDMC_100<-SpeciesData$RDMC*100
SpeciesData$LDMC_100<-SpeciesData$LDMC*100
SpeciesData$RTD_100<-SpeciesData$RTD*100
SpeciesData$H_num<-as.factor(SpeciesData$H_num)

# Add NA for H1 SLA for forb species 
SpeciesData$SLA<-ifelse((SpeciesData$GrowthForm %in% c("FORB","SHRUB") & SpeciesData$H_num == "H1"), NA, SpeciesData$SLA)
SpeciesData$SLA<-ifelse((SpeciesData$SLA == 0), NA, SpeciesData$SLA)

Species_splits = split(SpeciesData, paste(SpeciesData$SPECIES))


# Missing HEVI and MACA? 
TRY_data<-read.csv(file.choose())

TRY_data$TraitName<-as.character(TRY_data$TraitName)
SLA_TRY<-TRY_data[grepl("leaf",TRY_data$TraitName),]
levels(SLA_TRY$TraitName)
levels(SLA_TRY$SpeciesName)
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

# For now keep only the UnitNames "cm/g" and "mm2 mg-1"
SLA_TRY<-SLA_TRY[SLA_TRY$UnitName %in% c("cm/g","mm2 mg-1"),]
SLA_TRY_splits<-split(SLA_TRY, paste(SLA_TRY$SPECIES))


# ---------------------------- Models 
# To get TRY distributions
Trait_rstanarm_TRY <- function (df) {
  stan_glm(StdValue ~ 1 , 
            data = df,
            adapt_delta = 0.95)
}


# To get distributions
Trait_rstanarm_indv <- function (df, var) {
  stan_lmer(df[[var]] ~ 0 + H_num + (1|POP_ID), 
            data = df,
            prior_intercept = student_t(3,0,30), 
            adapt_delta = 0.95)
}

# -----------------------------run mods ---------------------------

SLA_TRY_mod_out<-lapply(SLA_TRY_splits, Trait_rstanarm_TRY) # Error - Constant variable(s) found: StdValue 

SLA_TRY_ACMI<-Trait_rstanarm_TRY(SLA_TRY_splits$ACMI) # These all work 
SLA_TRY_ARTR<-Trait_rstanarm_TRY(SLA_TRY_splits$ARTR)
SLA_TRY_ELTR<-Trait_rstanarm_TRY(SLA_TRY_splits$ELTR)
SLA_TRY_HEAN<-Trait_rstanarm_TRY(SLA_TRY_splits$HEAN)
SLA_TRY_HECO<-Trait_rstanarm_TRY(SLA_TRY_splits$HECO)
SLA_TRY_PLPA<-Trait_rstanarm_TRY(SLA_TRY_splits$PLPA)
SLA_TRY_PAMU<-Trait_rstanarm_TRY(SLA_TRY_splits$PAMU)
SLA_TRY_VUOC<-Trait_rstanarm_TRY(SLA_TRY_splits$VUOC)

SLA_Species_Ind_SLA_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "SLA")

Species_SLA_ACMI<-Trait_rstanarm_indv(Species_splits$ACMI,"SLA") # Divergent transitions
Species_SLA_ARTR<-Trait_rstanarm_indv(Species_splits$ARTR,"SLA") # Divergent transitions
Species_SLA_ELTR<-Trait_rstanarm_indv(Species_splits$ELTR,"SLA") 
Species_SLA_HEAN<-Trait_rstanarm_indv(Species_splits$HEAN,"SLA")
Species_SLA_HECO<-Trait_rstanarm_indv(Species_splits$HECO,"SLA") # Divergent transitions
Species_SLA_PLPA<-Trait_rstanarm_indv(Species_splits$PLPA,"SLA") # Divergent transitions
Species_SLA_PAMU<-Trait_rstanarm_indv(Species_splits$PAMU,"SLA")
Species_SLA_VUOC<-Trait_rstanarm_indv(Species_splits$VUOC,"SLA") # Divergent transitions
Species_SLA_MACA<-Trait_rstanarm_indv(Species_splits$MACA,"SLA")
Species_SLA_HEVI<-Trait_rstanarm_indv(Species_splits$HEVI,"SLA") # Divergent transitions
Species_SLA_MUPO<-Trait_rstanarm_indv(Species_splits$MUPO,"SLA")

RMR_Species_Ind_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "RMR_100")
LDMC_Species_Ind_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "LDMC_100")
RDMC_Species_Ind_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "RDMC_100")
RTD_Species_Ind_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "RTD_100")
SRL_Species_Ind_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "SRL")

# ------------------------- Get values ----------------------------------------------

posts_ACMI<-get_posts_forbs(SLA_Species_Ind_SLA_mod_out$ACMI, SLA_TRY_ACMI)
posts_ARTR<-get_posts_forbs(SLA_Species_Ind_SLA_mod_out$ARTR, SLA_TRY_ARTR)
posts_ELTR<-get_posts(SLA_Species_Ind_SLA_mod_out$ELTR, SLA_TRY_ELTR)
posts_HEAN<-get_posts_forbs(SLA_Species_Ind_SLA_mod_out$HEAN, SLA_TRY_HEAN)
posts_HECO<-get_posts(SLA_Species_Ind_SLA_mod_out$HECO, SLA_TRY_HECO)
posts_PLPA<-get_posts_forbs(SLA_Species_Ind_SLA_mod_out$PLPA, SLA_TRY_PLPA)
posts_PAMU<-get_posts_forbs(SLA_Species_Ind_SLA_mod_out$PAMU, SLA_TRY_PAMU)
posts_VUOC<-get_posts(SLA_Species_Ind_SLA_mod_out$VUOC, SLA_TRY_VUOC)
posts_MACA<-get_posts_noTRY(Species_SLA_MACA)
posts_MUPO<-get_posts_noTRY_grasses(Species_SLA_MUPO)
posts_HEVI<-get_posts_noTRY(Species_SLA_HEVI)

# Why dont these keep their name?
posts_LDMC_grasses<-lapply(list(LDMC_Species_Ind_mod_out$ELTR,LDMC_Species_Ind_mod_out$HECO, LDMC_Species_Ind_mod_out$MUPO, LDMC_Species_Ind_mod_out$VUOC),  get_posts_noTRY_grasses)
posts_LDMC_forbs<-lapply(list(LDMC_Species_Ind_mod_out$ACMI, LDMC_Species_Ind_mod_out$ARTR, LDMC_Species_Ind_mod_out$HEAN,
                           LDMC_Species_Ind_mod_out$HEVI, LDMC_Species_Ind_mod_out$HEVI, LDMC_Species_Ind_mod_out$MACA,
                           LDMC_Species_Ind_mod_out$PAMU, LDMC_Species_Ind_mod_out$PLPA), get_posts_noTRY)
posts_RDMC<-lapply(RDMC_Species_Ind_mod_out, get_posts_noTRY_grasses)
posts_RMR<-lapply(RMR_Species_Ind_mod_out, get_posts_noTRY_grasses)
posts_RTD<-lapply(RTD_Species_Ind_mod_out, get_posts_noTRY_grasses)
posts_SRL<-lapply(SRL_Species_Ind_mod_out, get_posts_noTRY_grasses)

# - Plotting 

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
