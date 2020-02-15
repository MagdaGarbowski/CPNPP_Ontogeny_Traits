## Correlation plots between traits at different timepoint
## Does covariation among traits remain consistant across time point? Are these graphs showing that
## Put r2 value onto plots at each time point? 
setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
source("Rscripts/Functions/Functions_TraitCorrelation_Plots.R")

library(psych)
library(plyr)
library(ggplot2)
library(gridExtra)
library(grid)

# Pop avg. data 
Pop_avg_data<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")
Pop_avg_data_2<-as.data.frame(Pop_avg_data[!colnames(Pop_avg_data) %in% c("X","LOCATION_CODE","POP_CODE","GrowthForm")])
Pop_avg_data_2$value[Pop_avg_data_2$value == 0] <-NA

# Pops growth rates data 
SpeciesData_GrowthRates<-read.csv("Data_Generated/TraitData_GrowthRates_2017.csv")
SpeciesData_GrowthRates$H_num<-as.factor(SpeciesData_GrowthRates$H_num)
levels(SpeciesData_GrowthRates$H_num)<-c("H1","H2","H3","H4")

# long to wide for pop_avg_data
Pop_avg_data_2<- reshape(Pop_avg_data, idvar = c("POP_ID","H_num","GrowthForm", "SPECIES"), timevar = "trait", direction = "wide")

ag2<-merge(Pop_avg_data_2, SpeciesData_GrowthRates, by=c("POP_ID","H_num"))
ag2$value.SLA<-ifelse((ag2$GrowthForm.x %in% c("FORB","SHRUB") & ag2$H_num == "H1"), NA, ag2$value.SLA) # Add NA for H1 SLA for forb species 
ag2$value.LDMC<-ifelse((ag2$GrowthForm.x %in% c("FORB","SHRUB") & ag2$H_num == "H1"), NA, ag2$value.LDMC) # Add NA for H1 SLA for forb species 

# Data transformations to improve normality 
cols<- c("value.SLA","value.LDMC","value.RMR","value.SRL","value.RDMC","value.RTD","Tot_weight_avg")
ag2[cols]<-log(ag2[cols])

H_splits = split(ag2, paste(ag2$H_num))


Harvest_pairs_panels<-lapply(H_splits[c(2,3,4)], pairs_panels_function, c("value.SLA","value.LDMC","value.RMR","value.SRL","value.RDMC","value.RTD", "RGR_Tot_ln", "RER_ln"))

pairs_panels_function(ag2,c("value.SLA","value.LDMC","value.RMR","value.SRL","value.RDMC","value.RTD", "RGR_Tot_ln", "RER_ln"))

# -----------------------------Get correlations coeffecients -----------------------
# Need to look at pearsons (panels plots output) to see significance

H2<-Cor_function(H_splits$H2)
H3<-Cor_function(H_splits$H3)
H4<-Cor_function(H_splits$H4)
rbind(H2, H3, H4)


#------------------------------- Plotting - -------------------------------

ag2_FG_splits<-split(ag2, paste(ag2$GrowthForm.x))
Grasses<-ag2_FG_splits$GRASS
Forbs<-ag2_FG_splits$FORB

Grasses_splits<-split(Grasses, paste(Grasses$H_num))
Forbs_splits<-split(Forbs, paste(Forbs$H_num))

legend_plot<-ggplot(Forbs) + 
  geom_point(aes(x = value.SLA, y = value.SRL, color = H_num), size = 5)+
  scale_color_manual("", 
                     breaks = c("H1","H2","H3","H4"),
                     values = c("#E69F00", "#CC79A7", "#56B4E9", "#009E73"))+
  theme(legend.position = "bottom", 
        legend.key = element_blank(),
        legend.key.size = unit(1,"cm"),
        legend.text=element_text(size=14))

leg<-get_legend(legend_plot)

SLA_RGR_grass<-cor_plot_all_H2H4(Grasses[!(Grasses$H_num %in% c("H1")),], "value.SLA", "RGR_Tot_ln", c(1.8, 6), c(-0.05, .25))
RER_RGR_grass<-cor_plot_all_H2H4(Grasses[!(Grasses$H_num %in% c("H1")),], "RER_ln", "RGR_Tot_ln", c(-0.05, 0.2), c(-0.15, .25))
SRL_RGR_grass<-cor_plot_all_H2H4(Grasses[!(Grasses$H_num %in% c("H1")),], "value.SRL", "RGR_Tot_ln", c(0.9,5.5), c(-0.15, 0.25))
SLA_LDMC_grass<-cor_plot(Grasses, "value.SLA", "value.LDMC", c(1,6), c(-3.2, -0.5))
SLA_SRL_grass<-cor_plot(Grasses, "value.SLA", "value.SRL", c(1.5, 6), c(0.5, 6))
SRL_RTD_grass<-cor_plot(Grasses, "value.SRL", "value.RTD", c(1, 5.5), c(-5, -0.5))

SLA_RGR_forb<-cor_plot_all_H2H4(Forbs[!(Forbs$H_num %in% c("H1")),], "value.SLA", "RGR_Tot_ln", c(1.8, 6), c(-0.05, .25))
RER_RGR_forb<-cor_plot_all_H2H4(Forbs[!(Forbs$H_num %in% c("H1")),], "RER_ln", "RGR_Tot_ln", c(-0.05, 0.2), c(-0.15, .25))
SRL_RGR_forb<-cor_plot_all_H2H4(Forbs[!(Forbs$H_num %in% c("H1")),], "value.SRL", "RGR_Tot_ln", c(0.9,5.5), c(-0.15, 0.25))
SLA_LDMC_forb<-cor_plot_all_H2H4(Forbs[!(Forbs$H_num %in% c("H1")),], "value.SLA", "value.LDMC", c(1,6), c(-3.2, -0.5))
SLA_SRL_forb<-cor_plot_all_H2H4(Forbs[!(Forbs$H_num %in% c("H1")),], "value.SLA", "value.SRL", c(1.5, 6), c(0.5, 6))
SRL_RTD_forb<-cor_plot(Forbs, "value.SRL", "value.RTD", c(1, 5.5), c(-5, -0.5))

grid.arrange(arrangeGrob(SLA_RGR_grass, SLA_RGR_forb, SRL_RGR_grass, SRL_RGR_forb,RER_RGR_grass,RER_RGR_forb,ncol = 2), leg, nrow = 2, heights = c(10,1))
grid.arrange(arrangeGrob(SLA_LDMC_grass, SLA_LDMC_forb, SLA_SRL_grass, SLA_SRL_forb, SRL_RTD_grass, SRL_RTD_forb, ncol = 2), 
             leg, nrow = 2, heights = c(12,1))

H1_G<-Cor_function_H1_Grass(Grasses_splits$H1)
H2_G<-Cor_function(Grasses_splits$H2)
H3_G<-Cor_function(Grasses_splits$H3)
H4_G<-Cor_function(Grasses_splits$H4)
rbind(H2_G, H3_G, H4_G)

H1_F<-Cor_function_Forb_H1(Forbs_splits$H1)
H2_F<-Cor_function(Forbs_splits$H2)
H3_F<-Cor_function(Forbs_splits$H3)
H4_F<-Cor_function(Forbs_splits$H4)
rbind(H2_F, H3_F, H4_F)

pdf("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Figures/Correlations.pdf", height = 11, width = 7)
grid.arrange(arrangeGrob(textGrob("Grasses"), textGrob("Forbs"), ncol = 2),
             arrangeGrob(textGrob ("Correlations with Growth Rate"), ncol = 1),
             arrangeGrob(SLA_RGR_grass, SLA_RGR_forb, SRL_RGR_grass, SRL_RGR_forb,RER_RGR_grass,RER_RGR_forb,ncol = 2),
             arrangeGrob(textGrob("Trait correlations"), ncol = 1),
             arrangeGrob(SLA_LDMC_grass, SLA_LDMC_forb, SLA_SRL_grass, SLA_SRL_forb, SRL_RTD_grass, SRL_RTD_forb, ncol = 2),
             arrangeGrob(leg, ncol = 1), heights = c(0.4,0.6,8,0.6,8,0.6))
dev.off()

