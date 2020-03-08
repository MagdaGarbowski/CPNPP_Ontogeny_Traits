## Correlation plots between traits at different timepoint
## Does covariation among traits remain consistant across time point? Are these graphs showing that
## Put r2 value onto plots at each time point? 
setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
source("Rscripts/Functions/Functions_TraitCorrelation_Plots.R")

library(psych)
library(psychTools)
library(plyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(factoextra)

# Load data 
Pop_avg_data<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")
seed_weights<-read.csv("Data_Raw_TRY_Seedweights/traits_population_seedweights.csv")

# Pop data
Pop_avg_data_2<-as.data.frame(Pop_avg_data[!colnames(Pop_avg_data) %in% c("X","LOCATION_CODE","POP_CODE")])
Pop_avg_data_2$value[Pop_avg_data_2$value == 0] <-NA
Pop_avg_data_2<-Pop_avg_data_2[!Pop_avg_data_2$H_num =="H4",]

# Long to wide format? 
Pop_avg_data_2_wide<-reshape(Pop_avg_data_2, 
                             idvar = c("H_num","POP_ID", "H_num","SPECIES","GrowthForm"),
                             timevar = "trait",
                             direction = "wide")

Pop_avg_data_2_wide_keep<-Pop_avg_data_2_wide[c("POP_ID","GrowthForm","H_num","value.ln.SLA_w_cots","value.ln.LDMC_w_cots","value.ln.HT","value.ln.RASARatio","value.ln.RDMC","value.ln.RMR",
                                                "value.ln.RTD", "value.ln.value.SumOfAvgDiam.mm.","value.RER_ln", "value.ln.SRL")]

Pop_avg_data_2_wide_keep_2<-Pop_avg_data_2_wide_keep
colnames(Pop_avg_data_2_wide_keep_2)<-c("POP_ID","GrowthForm","H_num","SLA","LDMC","HT","RASA","RDMC","RMR","RTD","Diam","RER","SRL")


# Seed data 
seed_weights<-seed_weights[,c("POP_ID","AVG_SEED_WEIGHT")]
seed_weights$Tot_weight_avg_mg <- seed_weights$AVG_SEED_WEIGHT * 1000
seed_weights$Tot_weight_avg_mg_ln<-log(seed_weights$Tot_weight_avg_mg)
seed_weights$H_num1 = "H1"
seed_weights$H_num2 = "H2"
seed_weights$H_num3 = "H3"

seed_weights_long<-reshape(seed_weights,
                           direction = "long",
                           varying = c("H_num1","H_num2","H_num3"),
                           v.names = "H_num",
                           idvar = "POP_ID",
                           timevar = "x",
                           times = c("H_num1","H_num2","H_num3"))
rownames(seed_weights_long) = NULL
seed_weights_long$SMass<-seed_weights_long$Tot_weight_avg_mg_ln
seed_weights_long[,c("x", "Tot_weight_avg_mg","AVG_SEED_WEIGHT","Tot_weight_avg_mg_ln")]<-NULL

seed_weights_long$H_num<-as.factor(seed_weights_long$H_num)

#### Combine seed and trait data 

Pop_data_w_seeds<-merge(Pop_avg_data_2_wide_keep_2, seed_weights_long, by = c("POP_ID","H_num"))
Pop_data_w_seeds<-Pop_data_w_seeds[complete.cases(Pop_data_w_seeds),]
Pop_data_w_seeds2 = Pop_data_w_seeds
Pop_data_w_seeds2_splits<-split(Pop_data_w_seeds2, paste(Pop_data_w_seeds2$H_num))


#Pop_data_w_seeds2[,c("POP_ID","GrowthForm")] = NULL

H_splits = split(Pop_data_w_seeds2, paste(Pop_data_w_seeds2$H_num))

H1_corr.out<-corr.test(H_splits$H1[,-c(1,2,3)], y = NULL, use = "pairwise", method = "pearson", adjust = "holm", alpha = 0.1, ci = TRUE, minlength = 10)
H1_corr<-print(corr.p(H1_corr.out$r, n = 51), short = FALSE, digits = 4)

H2_corr.out<-corr.test(H_splits$H2[,-1], y = NULL, use = "pairwise", method = "pearson", adjust = "holm", alpha = 0.1, ci = TRUE, minlength = 10)
H2_corr<-print(corr.p(H2_corr.out$r, n = 51), short = FALSE, digits = 4)

H3_corr.out<-corr.test(H_splits$H3[,-1], y = NULL, use = "pairwise", method = "pearson", adjust = "holm", alpha = 0.05, ci = TRUE, minlength = 10)
H3_corr<-print(corr.p(H3_corr.out$r, n = 51), short = FALSE, digits = 4)

corr_all<-cbind(H1_corr,H2_corr,H3_corr)

#Harvest_pairs_panels<-lapply(H_splits[c(1,2,3)], pairs_panels_function, c("SLA","LDMC","HT","RASA","RDMC","RMR","RTD","Diam","RER","RGR"))

# PCA - Attempt 1 
H1_dat<-H_splits$H1[complete.cases(H_splits$H1),][,c("POP_ID","GrowthForm","SLA","LDMC","HT","RASA","RDMC","RMR","RTD","Diam","RER","SRL","SMass")]
H2_dat<-H_splits$H2[complete.cases(H_splits$H2),][,c("POP_ID","GrowthForm","SLA","LDMC","HT","RASA","RDMC","RMR","RTD","Diam","RER","SRL","SMass")]
H3_dat<-H_splits$H3[complete.cases(H_splits$H3),][,c("POP_ID","GrowthForm","SLA","LDMC","HT","RASA","RDMC","RMR","RTD","Diam","RER","SRL","SMass")]


PCA_H1<- principal(H1_dat[,-c(1,2)], 3, rotate = "varimax", scores = TRUE)
PCA_H2<- principal(H2_dat[,-c(1,2)], 3, rotate = "varimax", scores = TRUE)
PCA_H3<- principal(H3_dat[,-c(1,2)], 3, rotate = "varimax", scores = TRUE)

RC_H1_1_2<-biplot(PCA_H1,choose=c(1,2), col = c( "blue", "black"), pch=c(24,21)[H1_dat[1:50,"GrowthForm"]], group = H1_dat[1:50, "GrowthForm"])
RC_H1_1_3<-biplot(PCA_H1,choose=c(1,3), col = c( "blue", "black"), pch=c(24,21)[H1_dat[1:50,"GrowthForm"]], group = H1_dat[1:50, "GrowthForm"])

RC_H2_1_2<-biplot(PCA_H2,choose=c(1,2), col = c( "blue", "black"), pch=c(24,21)[H1_dat[1:50,"GrowthForm"]], group = H1_dat[1:50, "GrowthForm"])
RC_H2_1_3<-biplot(PCA_H2,choose=c(1,3), col = c( "blue", "black"), pch=c(24,21)[H1_dat[1:50,"GrowthForm"]], group = H1_dat[1:50, "GrowthForm"])

RC_H3_1_2<-biplot(PCA_H3,choose=c(1,2), col = c( "blue", "black"), pch=c(24,21)[H1_dat[1:50,"GrowthForm"]], group = H1_dat[1:50, "GrowthForm"])
RC_H3_1_3<-biplot(PCA_H3,choose=c(1,3), col = c( "blue", "black"), pch=c(24,21)[H1_dat[1:50,"GrowthForm"]], group = H1_dat[1:50, "GrowthForm"])



pdf("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Figures/PCAs.pdf", height = 6, width = 10)
grid.arrange(PCA_h1, PCA_h2, PCA_h3, ncol = 3)
dev.off()

# H1: RC1: Neg loadings of SLA pos. loadings LDMC, RASA, RTD # cheap leaf, expensive/extensive root? 
#     RC2: Positive loadings HT, Diameter, RER. Neg : RGR # Aquisitive tissues, low mass 
#     RC3: Positive loadings of RDMC, RMR # Root investment 
# H2: RC1: Positive loadings of Hieght, diamter, RGR # Growth axis 
#     RC2: Neg SLA positive LDMC, RASA #Leaf tissue construction
#     RC3: positive RMR, RER # Root investment 
# H3: RC1: Neg SLA, Positive RASA, RMR, RTD #Similar to H1
#     RC2: Positive: LDMC, HT, RDMC, Diam # Construction 
#     RC3: Positive RER, RGR #Root investment 

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
#
#
#


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




