### Traits through time plotting
### Janurary 24, 2019 
### Plots to be made: 



# (1) Traits through time: SLA, LDMC, RMR, RDMC, RTD, SRL 
# (a) Organize species along x-axis by functional group? 
# (2) Relative change of traits through time for each species 
# (3) Plasticity across time-points 

setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/") 
source("Rscripts/Functions/Functions_Figures.R")

library(ggplot2)

SpeciesxTimeMedians<- read.csv("Data_Generated/SpeciesxTimeMedians.csv")
GrowthrateMedians<- read.csv("Data_Generated/RelativegrowthrateMedians.csv")
RelativeDifferencesMedians<- read.csv("Data_Generated/RelativeDifferencesMedians.csv")
PlasticityMedians<- read.csv("Data_Generated/PlasticityMedians.csv")

# Specify species levels for plotting 
# VUOC missing from plasticity? 

SpeciesxTimeMedians$Species<-factor(SpeciesxTimeMedians$Species, levels = c("HEAN","MACA","PLPA","VUOC","ELTR","HECO","MUPO","ARTR", "ACMI","HEVI","PAMU"))
GrowthrateMedians$Species<-factor(GrowthrateMedians$Species, levels = c("HEAN","MACA","PLPA","VUOC","ELTR","HECO","MUPO","ARTR", "ACMI","HEVI","PAMU"))
RelativeDifferencesMedians$Species<-factor(RelativeDifferencesMedians$Species, levels = c("HEAN","MACA","PLPA","VUOC","ELTR","HECO","MUPO","ARTR", "ACMI","HEVI","PAMU"))
PlasticityMedians$trait<-factor(PlasticityMedians$trait, levels = c("SLA","LDMC","RMR","SRL","RDMC","RTD", "RER","RGR_Below_ln", "RGR_Above_ln", "RGR_Tot_ln"))
PlasticityMedians$Species<-factor(PlasticityMedians$Species, levels = c("HEAN","MACA","PLPA","ELTR","HECO","MUPO","ACMI","ARTR", "ACMI","HEVI","PAMU"))


SLA<-plot_function(SpeciesxTimeMedians, "SLA", 180)
LDMC<-plot_function(SpeciesxTimeMedians, "LDMC", 0.38)
RMR<-plot_function(SpeciesxTimeMedians, "RMR", 0.63)
RDMC<-plot_function(SpeciesxTimeMedians, "RDMC", 0.32)
RTD<-plot_function(SpeciesxTimeMedians, "RTD", 0.20)
SRL<-plot_function(SpeciesxTimeMedians, "SRL", 75)

rel_SLA<-plot_function(RelativeDifferencesMedians, "rel_SLA2", 480)
rel_LDMC<-plot_function(RelativeDifferencesMedians, "rel_LDMC2", 495)
rel_RMR<-plot_function(RelativeDifferencesMedians, "rel_RMR2", 150)
rel_RDMC<-plot_function(RelativeDifferencesMedians, "rel_RDMC2", 200)
rel_RTD<-plot_function(RelativeDifferencesMedians, "rel_RTD2", 180)
rel_SRL<-plot_function(RelativeDifferencesMedians, "rel_SRL2", 200)

RER<-plot_function(GrowthrateMedians, "RER", 0.12)
Tot_GR<-plot_function(GrowthrateMedians, "RGR_Tot", 0.14)
Above_GR<-plot_function(GrowthrateMedians, "RGR_Above", 0.15)
Below_GR<-plot_function(GrowthrateMedians, "RGR_Below", 0.14)


PI<-plot_function_PI(PlasticityMedians)

pdf("~/CommitteeMeeting/Figures_Reports/Traits_relative.pdf", height = 8, width = 12)
grid.arrange(rel_SLA, rel_SRL, rel_LDMC, rel_RDMC,rel_RMR, rel_RTD,  ncol = 2)
dev.off()

pdf("~/CommitteeMeeting/Figures_Reports/Traits_CI.pdf", height = 8, width = 12)
grid.arrange(SLA, SRL,LDMC, RDMC,RMR, RTD,  ncol = 2)
dev.off()

pdf("~/CommitteeMeeting/Figures_Reports/Traits_RGR.pdf", height = 6, width = 12)
grid.arrange(Tot_GR, RER, Above_GR,Below_GR, ncol = 2)
dev.off()

pdf("~/CommitteeMeeting/Figures_Reports/Traits_PlasticityALL.pdf", height = 6, width = 12)
PI
dev.off()

# ----------------------------------------  FIGURES TO IDENTIFY OUTLIERS -----------------------------------
SpeciesData_all_traits<-SpeciesData_all[colnames(SpeciesData_all) %in% c("SAMPLE_ID", "POP_ID","SPECIES", "LDMC","SLA","SRL","RDMC","RTD","RMR", "H_num")]
SpeciesData_all_traits_long<-reshape(data =SpeciesData_all_traits, 
                                     idvar = c("SAMPLE_ID", "SPECIES", "POP_ID", "H_num"), 
                                     varying = c("LDMC","SLA","SRL","RDMC","RTD","RMR"),
                                     v.names = c("value"),
                                     times = c("LDMC","SLA","SRL","RDMC","RTD","RMR"),
                                     direction = "long")
rownames(SpeciesData_all_traits_long) <- NULL
names(SpeciesData_all_traits_long)[names(SpeciesData_all_traits_long) %in% c("time")]<-"trait"
SpeciesData_all_traits_split<-split(SpeciesData_all_traits_long, paste(SpeciesData_all_traits_long$SPECIES))

ggplot(SpeciesData_all_traits_split$ACMI, aes (y=value, x = H_num)) + 
  geom_point(aes(color = POP_ID), position = position_jitter(width = 0.3, height = 0.1))+
  facet_wrap(~trait, ncol = 2, scales = "free")


# ---------------------------------------- FIGURE FOR SCB ABSTRACT -------------------------
library("RColorBrewer")

head(SpeciesxTimeMedians)
SpeciesxTimeMedians$GrowthForm<- ifelse(SpeciesxTimeMedians$Species %in% c("ACMI","ARTR","HEAN","HEVI","MACA","PAMU","PLPA"), "FORB", "GRASS")
Forbs_SLA<-SpeciesxTimeMedians[(SpeciesxTimeMedians$GrowthForm == "FORB" & SpeciesxTimeMedians$trait == "SLA"),]
Forbs_LDMC<-SpeciesxTimeMedians[(SpeciesxTimeMedians$GrowthForm == "FORB" & SpeciesxTimeMedians$trait == "LDMC"),]
Grass_SLA<-SpeciesxTimeMedians[(SpeciesxTimeMedians$GrowthForm == "GRASS" & SpeciesxTimeMedians$trait == "SLA"),]
Grass_LDMC<-SpeciesxTimeMedians[(SpeciesxTimeMedians$GrowthForm == "GRASS" & SpeciesxTimeMedians$trait == "LDMC"),]

SCB_Fig_Forb<-function(df, pd, y_text, lims, Title_text){ggplot(df, aes (y = CI_50, x = Harvest, group = Species)) + 
  geom_errorbar(aes(min = CI_10, max = CI_90), position = pd, size = 0.4, width = 0.2) +
  geom_errorbar(aes(color = Species, min = CI_20, max = CI_80), position = pd, size = 1.2, width = 0) +
  geom_point(aes(fill = Species, shape = Species), position = pd, size = 5) + 
  scale_color_manual(values =  c("#666666","#666666","#666666","#666666","#666666","#666666","#666666"))+
  scale_fill_manual(values =  c("white","white","white","#666666","white","white","white")) +
  scale_shape_manual(values = c(21,8,22,10,24,4,7,9))+
  scale_x_discrete(labels = c("10 days", "24 days", "42 days", "84 days"))+
  ylim(lims)+
  ylab(y_text) +
  ggtitle(Title_text)+
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size = 14), 
        plot.title = element_text(size = 18, hjust = 0.5))
}

Forbs_SLA_plot<-SCB_Fig_Forb(Forbs_SLA, position_dodge(0.70), expression(paste("SLA [m" ^{2}," g" ^{-1},"]")), c(-50,190), "Forbs") +
  annotate("text", x = 0.6, y = 190, label = "a)", size = 5)
Forbs_LDMC_plot<-SCB_Fig_Forb(Forbs_LDMC, position_dodge(0.70), expression(paste("LDMC [g" ," g" ^{-1},"]")), c(-0.05,0.45), "Forbs")+
  annotate("text", x = 0.6, y = .44, label = "c)", size = 5)
Grasses_SLA_plot<-SCB_Fig_Forb(Grass_SLA, position_dodge(0.45), expression(paste("SLA [m" ^{2}," g" ^{-1},"]")), c(-50,190), "Grasses")+
  annotate("text", x = 0.6, y = 190, label = "b)", size = 5)
Grasses_LDMC_plot<-SCB_Fig_Forb(Grass_LDMC, position_dodge(0.45),expression(paste("LDMC [g" ," g" ^{-1},"]")), c(-0.05,0.45), "Grasses")+
  annotate("text", x = 0.6, y = .44, label = "d)", size = 5)

pdf("~/Downloads/NACCB.pdf", height=8, width=14)
grid.arrange(Forbs_SLA_plot, Grasses_SLA_plot, 
             Forbs_LDMC_plot, Grasses_LDMC_plot, ncol = 2)
dev.off()




#  -----------------------ORIGINAL FIGURES BROKEN UP BY FUNCTIONAL GROUP-----------------------------------
# (1) Need in import correct datasets for plotting 

Grass_colors<-c("blue3", "violetred4","navajowhite2","grey65")
Grass_shapes<-c(21,21,21)

Forb_colors<-c("coral2","chartreuse4","goldenrod3", "plum4")
Forb_shapes<-c(22,22,22,22,22)

Annual_colors<-c("grey25","grey45", "grey65", "grey85")
Annual_shapes<-c(23,23,23,23)

# Turn Medians_intervals_SLA into a list and lapply plot function. How to name? 

SLA_GRASS_plot<-plot_function(Medians_intervals_SLA, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("SLA [mm" ^{2}," g" ^{-1},"]")))
SRL_GRASS_plot<-plot_function(Medians_intervals_SRL, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("SRL [m" ^{2}," g" ^{-1},"]")))
LDMC_GRASS_plot<-plot_function(Medians_intervals_LDMC, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("LDMC [g" ," g" ^{-1},"]")))
RDMC_GRASS_plot<-plot_function(Medians_intervals_RDMC, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("RDMC [g" ," g" ^{-1},"]")))
RTD_GRASS_plot<-plot_function(Medians_intervals_RTD, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("RTD [g" ," cm" ^{-3},"]")))
RMR_GRASS_plot<-plot_function(Medians_intervals_RMR, "PER_GRASS", Grass_colors, Grass_shapes, "RMR")
#RGR_GRASS_plot<-plot_function_rates(Medians_intervals_RGR, "PER_GRASS", Grass_colors, Grass_shapes, "RGR")
AB_RGR_GRASS_plot<-plot_function_rates(Medians_intervals_AB_RGR, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("Shoot GR [g" ," g" ^{-1}, " days" ^{-1},"]")), c(-0.075, 0.2), "none")
BG_RGR_GRASS_plot<-plot_function_rates(Medians_intervals_BG_RGR, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("Root GR [g" ," g" ^{-1}, " days" ^{-1},"]")), c(-0.075, 0.2), "none")
RER_GRASS_plot<-plot_function_rates(Medians_intervals_RER, "PER_GRASS", Grass_colors, Grass_shapes, expression(paste("RER [cm" ," cm" ^{-1}, " days" ^{-1},"]")), c(-0.01, 0.15), "none")
RER_GRASS_plot_legend<-plot_function_rates(Medians_intervals_RER, "PER_GRASS", Grass_colors, Grass_shapes, "RER", c(-0.01, 0.15), "right")


SLA_FORB_plot<-plot_function(Medians_intervals_SLA, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("SLA [mm" ^{2}," g" ^{-1},"]")))
SRL_FORB_plot<-plot_function(Medians_intervals_SRL, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("SRL [m" ^{2}," g" ^{-1},"]")))
LDMC_FORB_plot<-plot_function(Medians_intervals_LDMC, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("LDMC [g" ," g" ^{-1},"]")))
RDMC_FORB_plot<-plot_function(Medians_intervals_RDMC, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("RDMC [g" ," g" ^{-1},"]")))
RTD_FORB_plot<-plot_function(Medians_intervals_RTD, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("RTD [g" ," cm" ^{-3},"]")))
RMR_FORB_plot<-plot_function(Medians_intervals_RMR, "PER_FORB", Forb_colors, Forb_shapes, "RMR")
#RGR_FORB_plot<-plot_function_rates(Medians_intervals_RGR, "PER_FORB", Forb_colors, Forb_shapes, "RGR")
AB_RGR_FORB_plot<-plot_function_rates(Medians_intervals_AB_RGR, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("Shoot GR [g" ," g" ^{-1}, " days" ^{-1},"]")), c(-0.075, 0.2), "none")
BG_RGR_FORB_plot<-plot_function_rates(Medians_intervals_BG_RGR, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("Root GR [g" ," g" ^{-1}, " days" ^{-1},"]")), c(-0.075, 0.2), "none")
RER_FORB_plot<-plot_function_rates(Medians_intervals_RER, "PER_FORB", Forb_colors, Forb_shapes, expression(paste("RER [cm" ," cm" ^{-1}, " days" ^{-1},"]")), c(-0.01, 0.15), "none") 
RER_FORB_plot_legend<-plot_function_rates(Medians_intervals_RER, "PER_FORB", Forb_colors, Forb_shapes, "RER", c(-0.01, 0.15), "right")


SLA_ANNUAL_plot<-plot_function(Medians_intervals_SLA, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("SLA [mm" ^{2}," g" ^{-1},"]")))
SRL_ANNUAL_plot<-plot_function(Medians_intervals_SRL, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("SRL [m" ^{2}," g" ^{-1},"]")))
LDMC_ANNUAL_plot<-plot_function(Medians_intervals_LDMC, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("LDMC [g" ," g" ^{-1},"]")))
RDMC_ANNUAL_plot<-plot_function(Medians_intervals_RDMC, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("RDMC [g" ," g" ^{-1},"]")))
RTD_ANNUAL_plot<-plot_function(Medians_intervals_RTD, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("RTD [g" ," cm" ^{-3},"]")))
RMR_ANNUAL_plot<-plot_function(Medians_intervals_RMR, "ANNUAL", Annual_colors, Annual_shapes, "RMR")
#RGR_ANNUAL_plot<-plot_function_rates(Medians_intervals_RGR, "ANNUAL", Annual_colors, Annual_shapes, "RGR")
AB_RGR_ANNUAL_plot<-plot_function_rates(Medians_intervals_AB_RGR, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("Shoot GR [g" ," g" ^{-1}, " days" ^{-1},"]")), c(-0.01, 0.25), "none")
BG_RGR_ANNUAL_plot<-plot_function_rates(Medians_intervals_BG_RGR, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("Root GR [g" ," g" ^{-1}, " days" ^{-1},"]")), c(-0.01, 0.25), "none")
RER_ANNUAL_plot<-plot_function_rates(Medians_intervals_RER, "ANNUAL", Annual_colors, Annual_shapes, expression(paste("RER [cm" ," cm" ^{-1}, " days" ^{-1},"]")), c(-0.01, 0.2), "none") 
RER_ANNUAL_plot_legend<-plot_function_rates(Medians_intervals_RER, "ANNUAL", Annual_colors, Annual_shapes, "RER", c(-0.01, 0.2), "right")

# Get legends 
GRASS_legend<-get_legend(RER_GRASS_plot_legend)
FORB_legend<-get_legend(RER_FORB_plot_legend)
ANNUAL_legend<-get_legend(RER_ANNUAL_plot_legend)

# Plasticity plots 
Plasticity_annuals<-plot_function_PI(plasticity_df, 
                                     "ANNUAL", c("grey95","gray80","grey60","grey45","grey30","black"), 
                                     c(21,22,23,24,25,16), c(-0.1,0.64), "right")
Plasticity_PER_GRASS<-plot_function_PI(plasticity_df, 
                                       "PER_GRASS", c("grey95","gray80","grey60","grey45","grey30","black"), 
                                       c(21,22,23,24,25,16), c(-0.1,0.64), "right")
Plasticity_PER_FORB<-plot_function_PI(plasticity_df, 
                                      "PER_FORB", c("grey95","gray80","grey60","grey45","grey30","black"), 
                                      c(21,22,23,24,25,16), c(-0.1,0.64), "right")

# Plasticity plots - without H 1 harvest 
Plasticity_annuals_10day<-plot_function_PI(plasticity_df_10day, 
                                           "ANNUAL", c("grey95","gray80","grey60","grey45","grey30","black"), 
                                           c(21,22,23,24,25,16), c(-0.1,0.64), "right")
Plasticity_PER_GRASS_10day<-plot_function_PI(plasticity_df_10day, 
                                             "PER_GRASS", c("grey95","gray80","grey60","grey45","grey30","black"), 
                                             c(21,22,23,24,25,16), c(-0.1,0.64), "right")
Plasticity_PER_FORB_10day<-plot_function_PI(plasticity_df_10day, 
                                            "PER_FORB", c("grey95","gray80","grey60","grey45","grey30","black"), 
                                            c(21,22,23,24,25,16), c(-0.1,0.64), "right")




pdf("~/CommitteeMeeting/Figures_Reports/Grass_figs_1.pdf", height=8, width=8)
grid.arrange(SLA_GRASS_plot, LDMC_GRASS_plot, 
             SRL_GRASS_plot, RMR_GRASS_plot,
             RDMC_GRASS_plot, RTD_GRASS_plot, ncol = 2, 
             top = textGrob ("Grasses", gp = gpar(fontsize=14)))

dev.off()

pdf("~/CommitteeMeeting/Figures_Reports/Forb_figs_1.pdf", height=8, width=8)
grid.arrange(SLA_FORB_plot, LDMC_FORB_plot, 
             SRL_FORB_plot, RMR_FORB_plot,
             RDMC_FORB_plot, RTD_FORB_plot, ncol = 2, 
             top = textGrob ("Forbs", gp = gpar(fontsize=14)))
dev.off()

pdf("~/CommitteeMeeting/Figures_Reports/Annual_figs_1.pdf", height=8, width=8)
grid.arrange(SLA_ANNUAL_plot, LDMC_ANNUAL_plot, 
             SRL_ANNUAL_plot, RMR_ANNUAL_plot,
             RDMC_ANNUAL_plot, RTD_ANNUAL_plot, ncol = 2, 
             top = textGrob ("Annuals", gp = gpar(fontsize=14)))
dev.off()

pdf("~/CommitteeMeeting/Figures_Reports/Rates_figs_1.pdf", height=7, width=9)
grid.arrange(arrangeGrob(AB_RGR_GRASS_plot, BG_RGR_GRASS_plot, RER_GRASS_plot, GRASS_legend, ncol = 4, nrow  = 1, layout_matrix = rbind(c(1,2,3,4)),widths = c(3, 3, 3, 1.5), heights = c(2), 
                         top = textGrob ("Perennial Grasses", gp = gpar (fontsize = 12))),
             arrangeGrob(AB_RGR_FORB_plot, BG_RGR_FORB_plot, RER_FORB_plot, FORB_legend, ncol = 4, nrow  = 1, layout_matrix = rbind(c(1,2,3,4)),widths = c(3, 3, 3, 1.5), heights = c(2), 
                         top = textGrob ("Perennial Forbs", gp = gpar (fontsize = 12))),
             arrangeGrob(AB_RGR_ANNUAL_plot, BG_RGR_ANNUAL_plot, RER_ANNUAL_plot, ANNUAL_legend, ncol = 4, nrow  = 1, layout_matrix = rbind(c(1,2,3,4)),widths = c(3, 3, 3, 1.5), heights = c(2), 
                         top = textGrob ("Annuals", gp = gpar (fontsize = 12))),
             top = textGrob ("Growth Rates", gp = gpar(fontsize=14, face = "bold")))
dev.off()


pdf("~/CommitteeMeeting/Figures_Reports/PlasticityIndex.pdf", height = 8, width = 12)
grid.arrange(arrangeGrob(Plasticity_annuals, Plasticity_PER_GRASS, Plasticity_PER_FORB, ncol = 3, top = textGrob("All Harvest Plasticity Index")),
             arrangeGrob(Plasticity_annuals_10day, Plasticity_PER_GRASS_10day, Plasticity_PER_FORB_10day, ncol = 3, top = textGrob("Without 10-day harvest Plasticity Index")))
dev.off()



