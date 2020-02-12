setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
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

SLA_Species_Ind_RMR_mod_out<-lapply(Species_splits, Trait_rstanarm_indv, "RMR_100")


Species_RMR_ACMI<-Trait_rstanarm_indv(Species_splits$ACMI,"RMR_100")
Species_RMR_ARTR<-Trait_rstanarm_indv(Species_splits$ARTR,"RMR_100")
Species_RMR_ELTR<-Trait_rstanarm_indv(Species_splits$ELTR,"RMR_100")
Species_RMR_HEAN<-Trait_rstanarm_indv(Species_splits$HEAN,"RMR_100")
Species_RMR_HECO<-Trait_rstanarm_indv(Species_splits$HECO,"RMR_100")
Species_RMR_PLPA<-Trait_rstanarm_indv(Species_splits$PLPA,"RMR_100")
Species_RMR_PAMU<-Trait_rstanarm_indv(Species_splits$PAMU,"RMR_100")
Species_RMR_VUOC<-Trait_rstanarm_indv(Species_splits$VUOC,"RMR_100")
Species_RMR_MACA<-Trait_rstanarm_indv(Species_splits$MACA,"RMR_100")
Species_RMR_HEVI<-Trait_rstanarm_indv(Species_splits$HEVI,"RMR_100")
Species_RMR_MUPO<-Trait_rstanarm_indv(Species_splits$MUPO,"RMR_100")

Species_LDMC_ACMI<-Trait_rstanarm_indv(Species_splits$ACMI,"LDMC_100")
Species_LDMC_ARTR<-Trait_rstanarm_indv(Species_splits$ARTR,"LDMC_100")
Species_LDMC_ELTR<-Trait_rstanarm_indv(Species_splits$ELTR,"LDMC_100")
Species_LDMC_HEAN<-Trait_rstanarm_indv(Species_splits$HEAN,"LDMC_100")
Species_LDMC_HECO<-Trait_rstanarm_indv(Species_splits$HECO,"LDMC_100")
Species_LDMC_PLPA<-Trait_rstanarm_indv(Species_splits$PLPA,"LDMC_100")
Species_LDMC_PAMU<-Trait_rstanarm_indv(Species_splits$PAMU,"LDMC_100")
Species_LDMC_VUOC<-Trait_rstanarm_indv(Species_splits$VUOC,"LDMC_100")
Species_LDMC_MACA<-Trait_rstanarm_indv(Species_splits$MACA,"LDMC_100")
Species_LDMC_HEVI<-Trait_rstanarm_indv(Species_splits$HEVI,"LDMC_100")
Species_LDMC_MUPO<-Trait_rstanarm_indv(Species_splits$MUPO,"LDMC_100")


Species_RDMC_ACMI<-Trait_rstanarm_indv(Species_splits$ACMI,"RDMC_100")
Species_RDMC_ARTR<-Trait_rstanarm_indv(Species_splits$ARTR,"RDMC_100")
Species_RDMC_ELTR<-Trait_rstanarm_indv(Species_splits$ELTR,"RDMC_100")
Species_RDMC_HEAN<-Trait_rstanarm_indv(Species_splits$HEAN,"RDMC_100")
Species_RDMC_HECO<-Trait_rstanarm_indv(Species_splits$HECO,"RDMC_100")
Species_RDMC_PLPA<-Trait_rstanarm_indv(Species_splits$PLPA,"RDMC_100")
Species_RDMC_PAMU<-Trait_rstanarm_indv(Species_splits$PAMU,"RDMC_100")
Species_RDMC_VUOC<-Trait_rstanarm_indv(Species_splits$VUOC,"RDMC_100")
Species_RDMC_MACA<-Trait_rstanarm_indv(Species_splits$MACA,"RDMC_100")
Species_RDMC_HEVI<-Trait_rstanarm_indv(Species_splits$HEVI,"RDMC_100")
Species_RDMC_MUPO<-Trait_rstanarm_indv(Species_splits$MUPO,"RDMC_100")

Species_RTD_ACMI<-Trait_rstanarm_indv(Species_splits$ACMI,"RTD_100")
Species_RTD_ARTR<-Trait_rstanarm_indv(Species_splits$ARTR,"RTD_100")
Species_RTD_ELTR<-Trait_rstanarm_indv(Species_splits$ELTR,"RTD_100")
Species_RTD_HEAN<-Trait_rstanarm_indv(Species_splits$HEAN,"RTD_100")
Species_RTD_HECO<-Trait_rstanarm_indv(Species_splits$HECO,"RTD_100")
Species_RTD_PLPA<-Trait_rstanarm_indv(Species_splits$PLPA,"RTD_100")
Species_RTD_PAMU<-Trait_rstanarm_indv(Species_splits$PAMU,"RTD_100")
Species_RTD_VUOC<-Trait_rstanarm_indv(Species_splits$VUOC,"RTD_100")
Species_RTD_MACA<-Trait_rstanarm_indv(Species_splits$MACA,"RTD_100")
Species_RTD_HEVI<-Trait_rstanarm_indv(Species_splits$HEVI,"RTD_100")
Species_RTD_MUPO<-Trait_rstanarm_indv(Species_splits$MUPO,"RTD_100")

Species_SRL_ACMI<-Trait_rstanarm_indv(Species_splits$ACMI,"SRL")
Species_SRL_ARTR<-Trait_rstanarm_indv(Species_splits$ARTR,"SRL")
Species_SRL_ELTR<-Trait_rstanarm_indv(Species_splits$ELTR,"SRL")
Species_SRL_HEAN<-Trait_rstanarm_indv(Species_splits$HEAN,"SRL")
Species_SRL_HECO<-Trait_rstanarm_indv(Species_splits$HECO,"SRL")
Species_SRL_PLPA<-Trait_rstanarm_indv(Species_splits$PLPA,"SRL")
Species_SRL_PAMU<-Trait_rstanarm_indv(Species_splits$PAMU,"SRL")
Species_SRL_VUOC<-Trait_rstanarm_indv(Species_splits$VUOC,"SRL")
Species_SRL_MACA<-Trait_rstanarm_indv(Species_splits$MACA,"SRL")
Species_SRL_HEVI<-Trait_rstanarm_indv(Species_splits$HEVI,"SRL")
Species_SRL_MUPO<-Trait_rstanarm_indv(Species_splits$MUPO,"SRL")


# --------------------------------- Functions 
# To get posteriors and get them into long format
get_posts<-function(model, TRYmodel) {
  tmp1<-as.data.frame(model, regex_pars = "H_num")
  tmp2<-reshape(tmp1, 
                direction = "long",
                varying = c("H_numH1","H_numH2","H_numH3","H_numH4"),
                v.names = "Value",
                timevar = "Harvest",
                times = c("H_numH1","H_numH2","H_numH3","H_numH4"))
  tmp2<-tmp2[,!names(tmp2) %in% c("id")]
  rownames(tmp2)<-NULL
  TRYmodel_1<-as.data.frame(TRYmodel, regex_pars = "Intercept")
  TRYmodel_1$Harvest <-"TRY"
  names(TRYmodel_1)[names(TRYmodel_1) == "(Intercept)"]<-"Value"
  TRYmodel_1<-TRYmodel_1[c("Harvest","Value")]
  tmp3<-rbind(tmp2,TRYmodel_1)
  tmp3$Harvest<-factor(tmp3$Harvest, levels = c("TRY","H_numH1", "H_numH2", "H_numH3","H_numH4"))
  return(tmp3)
}


# To get posteriors and get them into long format
get_posts_noTRY_grasses<-function(model) {
  tmp1<-as.data.frame(model, regex_pars = "H_num")
  tmp2<-reshape(tmp1, 
                direction = "long",
                varying = c("H_numH1","H_numH2","H_numH3","H_numH4"),
                v.names = "Value",
                timevar = "Harvest",
                times = c("H_numH1","H_numH2","H_numH3","H_numH4"))
  tmp2<-tmp2[,!names(tmp2) %in% c("id")]
  rownames(tmp2)<-NULL
  tmp2$Harvest<-factor(tmp2$Harvest, levels = c("H_numH1", "H_numH2", "H_numH3","H_numH4"))
  return(tmp2)
}

# To get posteriors and get them into long format
get_posts_noTRY<-function(model) {
  tmp1<-as.data.frame(model, regex_pars = "H_num")
  tmp1<-tmp1[,!names(tmp1) %in% c("H_numH1")]
  tmp2<-reshape(tmp1, 
                direction = "long",
                varying = c("H_numH2","H_numH3","H_numH4"),
                v.names = "Value",
                timevar = "Harvest",
                times = c("H_numH2","H_numH3","H_numH4"))
  tmp2<-tmp2[,!names(tmp2) %in% c("id")]
  rownames(tmp2)<-NULL
  tmp2$Harvest<-factor(tmp2$Harvest, levels = c( "H_numH2", "H_numH3","H_numH4"))
  return(tmp2)
}


get_posts_forbs<-function(model, TRYmodel) {
  tmp1<-as.data.frame(model, regex_pars = "H_num")
  tmp1<-tmp1[,!names(tmp1) %in% c("H_numH1")]
  tmp2<-reshape(tmp1, 
                direction = "long",
                varying = c("H_numH2","H_numH3","H_numH4"),
                v.names = "Value",
                timevar = "Harvest",
                times = c("H_numH2","H_numH3","H_numH4"))
  tmp2<-tmp2[,!names(tmp2) %in% c("id")]
  rownames(tmp2)<-NULL
  TRYmodel_1<-as.data.frame(TRYmodel, regex_pars = "Intercept")
  TRYmodel_1$Harvest <-"TRY"
  names(TRYmodel_1)[names(TRYmodel_1) == "(Intercept)"]<-"Value"
  TRYmodel_1<-TRYmodel_1[c("Harvest","Value")]
  tmp3<-rbind(tmp2,TRYmodel_1)
  tmp3$Harvest<-factor(tmp3$Harvest, levels = c("TRY", "H_numH2", "H_numH3","H_numH4"))
  return(tmp3)
}


density_plot<-function(df, lims, adjust, yname, xname){ggplot(df, aes(x = Value, y = Harvest, fill = Harvest)) + 
    scale_x_continuous(limits = lims,
                       name = xname)+
    scale_y_discrete(breaks = c("H_numH4","H_numH3","H_numH2","H_numH1", "TRY"),
                     labels = c("H_numH1" = "10", "H_numH2" = "24", "H_numH3" = "42", "H_numH4" = "84", "TRY"), 
                     name = yname,
                     expand = expand_scale(mult = adjust))+
    geom_density_ridges(scale = 10, bandwidth = 5, rel_min_height = 0.01,
                        alpha = .2, color = c("grey40")) + 
    theme_ridges()+
    scale_fill_cyclical(values = c("grey90","gray70","grey50","grey30", "grey10"))+
    coord_cartesian(clip = "off")+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())+
    theme_bw()
}

# ------------------------- Get values

posts_ACMI<-get_posts_forbs(Species_SLA_ACMI, SLA_TRY_ACMI)
posts_ARTR<-get_posts_forbs(Species_SLA_ARTR, SLA_TRY_ARTR)
posts_ELTR<-get_posts(Species_SLA_ELTR, SLA_TRY_ELTR)
posts_HEAN<-get_posts_forbs(Species_SLA_HEAN, SLA_TRY_HEAN)
posts_HECO<-get_posts(Species_SLA_HECO, SLA_TRY_HECO)
posts_PLPA<-get_posts_forbs(Species_SLA_PLPA, SLA_TRY_PLPA)
posts_PAMU<-get_posts_forbs(Species_SLA_PAMU, SLA_TRY_PAMU)
posts_VUOC<-get_posts(Species_SLA_VUOC, SLA_TRY_VUOC)
posts_MACA<-get_posts_noTRY(Species_SLA_MACA)
posts_MUPO<-get_posts_noTRY_grasses(Species_SLA_MUPO)
posts_HEVI<-get_posts_noTRY(Species_SLA_HEVI)

posts_ACMI_RMR<-get_posts_noTRY_grasses(Species_RMR_ACMI)
posts_ARTR_RMR<-get_posts_noTRY_grasses(Species_RMR_ARTR)
posts_ELTR_RMR<-get_posts_noTRY_grasses(Species_RMR_ELTR)
posts_HEAN_RMR<-get_posts_noTRY_grasses(Species_RMR_HEAN)
posts_HECO_RMR<-get_posts_noTRY_grasses(Species_RMR_HECO)
posts_PLPA_RMR<-get_posts_noTRY_grasses(Species_RMR_PLPA)
posts_PAMU_RMR<-get_posts_noTRY_grasses(Species_RMR_PAMU)
posts_VUOC_RMR<-get_posts_noTRY_grasses(Species_RMR_VUOC)
posts_MACA_RMR<-get_posts_noTRY_grasses(Species_RMR_MACA)
posts_MUPO_RMR<-get_posts_noTRY_grasses(Species_RMR_MUPO)
posts_HEVI_RMR<-get_posts_noTRY_grasses(Species_RMR_HEVI)

posts_ACMI_LDMC<-get_posts_noTRY(Species_LDMC_ACMI)
posts_ARTR_LDMC<-get_posts_noTRY(Species_LDMC_ARTR)
posts_ELTR_LDMC<-get_posts_noTRY_grasses(Species_LDMC_ELTR)
posts_HEAN_LDMC<-get_posts_noTRY(Species_LDMC_HEAN)
posts_HECO_LDMC<-get_posts_noTRY_grasses(Species_LDMC_HECO)
posts_PLPA_LDMC<-get_posts_noTRY(Species_LDMC_PLPA)
posts_PAMU_LDMC<-get_posts_noTRY(Species_LDMC_PAMU)
posts_VUOC_LDMC<-get_posts_noTRY_grasses(Species_LDMC_VUOC)
posts_MACA_LDMC<-get_posts_noTRY(Species_LDMC_MACA)
posts_MUPO_LDMC<-get_posts_noTRY_grasses(Species_LDMC_MUPO)
posts_HEVI_LDMC<-get_posts_noTRY(Species_LDMC_HEVI)

posts_ACMI_RDMC<-get_posts_noTRY_grasses(Species_RDMC_ACMI)
posts_ARTR_RDMC<-get_posts_noTRY_grasses(Species_RDMC_ARTR)
posts_ELTR_RDMC<-get_posts_noTRY_grasses(Species_RDMC_ELTR)
posts_HEAN_RDMC<-get_posts_noTRY_grasses(Species_RDMC_HEAN)
posts_HECO_RDMC<-get_posts_noTRY_grasses(Species_RDMC_HECO)
posts_PLPA_RDMC<-get_posts_noTRY_grasses(Species_RDMC_PLPA)
posts_PAMU_RDMC<-get_posts_noTRY_grasses(Species_RDMC_PAMU)
posts_VUOC_RDMC<-get_posts_noTRY_grasses(Species_RDMC_VUOC)
posts_MACA_RDMC<-get_posts_noTRY_grasses(Species_RDMC_MACA)
posts_MUPO_RDMC<-get_posts_noTRY_grasses(Species_RDMC_MUPO)
posts_HEVI_RDMC<-get_posts_noTRY_grasses(Species_RDMC_HEVI)

posts_ACMI_RTD<-get_posts_noTRY_grasses(Species_RTD_ACMI)
posts_ARTR_RTD<-get_posts_noTRY_grasses(Species_RTD_ARTR)
posts_ELTR_RTD<-get_posts_noTRY_grasses(Species_RTD_ELTR)
posts_HEAN_RTD<-get_posts_noTRY_grasses(Species_RTD_HEAN)
posts_HECO_RTD<-get_posts_noTRY_grasses(Species_RTD_HECO)
posts_PLPA_RTD<-get_posts_noTRY_grasses(Species_RTD_PLPA)
posts_PAMU_RTD<-get_posts_noTRY_grasses(Species_RTD_PAMU)
posts_VUOC_RTD<-get_posts_noTRY_grasses(Species_RTD_VUOC)
posts_MACA_RTD<-get_posts_noTRY_grasses(Species_RTD_MACA)
posts_MUPO_RTD<-get_posts_noTRY_grasses(Species_RTD_MUPO)
posts_HEVI_RTD<-get_posts_noTRY_grasses(Species_RTD_HEVI)

posts_ACMI_SRL<-get_posts_noTRY_grasses(Species_SRL_ACMI)
posts_ARTR_SRL<-get_posts_noTRY_grasses(Species_SRL_ARTR)
posts_ELTR_SRL<-get_posts_noTRY_grasses(Species_SRL_ELTR)
posts_HEAN_SRL<-get_posts_noTRY_grasses(Species_SRL_HEAN)
posts_HECO_SRL<-get_posts_noTRY_grasses(Species_SRL_HECO)
posts_PLPA_SRL<-get_posts_noTRY_grasses(Species_SRL_PLPA)
posts_PAMU_SRL<-get_posts_noTRY_grasses(Species_SRL_PAMU)
posts_VUOC_SRL<-get_posts_noTRY_grasses(Species_SRL_VUOC)
posts_MACA_SRL<-get_posts_noTRY_grasses(Species_SRL_MACA)
posts_MUPO_SRL<-get_posts_noTRY_grasses(Species_SRL_MUPO)
posts_HEVI_SRL<-get_posts_noTRY_grasses(Species_SRL_HEVI)



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

RMR_plot_ACMI<-density_plot(posts_ACMI_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ACMI (F)", size = 4)
RMR_plot_ARTR<-density_plot(posts_ARTR_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ARTR (S)", size = 4)
RMR_plot_ELTR<-density_plot(posts_ELTR_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ELTR (G)", size = 4)
RMR_plot_HEAN<-density_plot(posts_HEAN_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEAN (AF)", size = 4)
RMR_plot_HECO<-density_plot(posts_HECO_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HECO (G)", size = 4)
RMR_plot_PAMU<-density_plot(posts_PAMU_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PAMU (F)", size = 4)
RMR_plot_PLPA<-density_plot(posts_PLPA_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PLPA (AF)", size = 4)
RMR_plot_VUOC<-density_plot(posts_VUOC_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "VUOC (AG)", size = 4)
RMR_plot_MACA<-density_plot(posts_MACA_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MACA (BF)", size = 4)
RMR_plot_HEVI<-density_plot(posts_HEVI_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEVI (F)", size = 4)
RMR_plot_MUPO<-density_plot(posts_MUPO_RMR, c(-10, 75), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MUPO (G)", size = 4)

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

RDMC_plot_ACMI<-density_plot(posts_ACMI_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ACMI (F)", size = 4)
RDMC_plot_ARTR<-density_plot(posts_ARTR_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ARTR (S)", size = 4)
RDMC_plot_ELTR<-density_plot(posts_ELTR_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "ELTR (G)", size = 4)
RDMC_plot_HEAN<-density_plot(posts_HEAN_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEAN (AF)", size = 4)
RDMC_plot_HECO<-density_plot(posts_HECO_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HECO (G)", size = 4)
RDMC_plot_PAMU<-density_plot(posts_PAMU_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PAMU (F)", size = 4)
RDMC_plot_PLPA<-density_plot(posts_PLPA_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "PLPA (AF)", size = 4)
RDMC_plot_VUOC<-density_plot(posts_VUOC_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "VUOC (AG)", size = 4)
RDMC_plot_MACA<-density_plot(posts_MACA_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MACA (BF)", size = 4)
RDMC_plot_HEVI<-density_plot(posts_HEVI_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "HEVI (F)", size = 4)
RDMC_plot_MUPO<-density_plot(posts_MUPO_RDMC, c(-10, 40), c(0,3.5), "", "") + annotate("text", x = 0, y = 12, label = "MUPO (G)", size = 4)

RTD_plot_ACMI<-density_plot(posts_ACMI_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "ACMI (F)", size = 4)
RTD_plot_ARTR<-density_plot(posts_ARTR_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "ARTR (S)", size = 4)
RTD_plot_ELTR<-density_plot(posts_ELTR_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "ELTR (G)", size = 4)
RTD_plot_HEAN<-density_plot(posts_HEAN_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "HEAN (AF)", size = 4)
RTD_plot_HECO<-density_plot(posts_HECO_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "HECO (G)", size = 4)
RTD_plot_PAMU<-density_plot(posts_PAMU_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "PAMU (F)", size = 4)
RTD_plot_PLPA<-density_plot(posts_PLPA_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "PLPA (AF)", size = 4)
RTD_plot_VUOC<-density_plot(posts_VUOC_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "VUOC (AG)", size = 4)
RTD_plot_MACA<-density_plot(posts_MACA_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "MACA (BF)", size = 4)
RTD_plot_HEVI<-density_plot(posts_HEVI_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "HEVI (F)", size = 4)
RTD_plot_MUPO<-density_plot(posts_MUPO_RTD, c(-10, 32), c(0,3.5), "", "") + annotate("text", x = -5, y = 12, label = "MUPO (G)", size = 4)

SRL_plot_ACMI<-density_plot(posts_ACMI_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "ACMI (F)", size = 4)
SRL_plot_ARTR<-density_plot(posts_ARTR_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "ARTR (S)", size = 4)
SRL_plot_ELTR<-density_plot(posts_ELTR_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "ELTR (G)", size = 4)
SRL_plot_HEAN<-density_plot(posts_HEAN_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "HEAN (AF)", size = 4)
SRL_plot_HECO<-density_plot(posts_HECO_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "HECO (G)", size = 4)
SRL_plot_PAMU<-density_plot(posts_PAMU_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "PAMU (F)", size = 4)
SRL_plot_PLPA<-density_plot(posts_PLPA_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "PLPA (AF)", size = 4)
SRL_plot_VUOC<-density_plot(posts_VUOC_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "VUOC (AG)", size = 4)
SRL_plot_MACA<-density_plot(posts_MACA_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "MACA (BF)", size = 4)
SRL_plot_HEVI<-density_plot(posts_HEVI_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "HEVI (F)", size = 4)
SRL_plot_MUPO<-density_plot(posts_MUPO_SRL, c(-25, 110), c(0,3.5), "", "") + annotate("text", x = -10, y=13, label = "MUPO (G)", size = 4)




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
