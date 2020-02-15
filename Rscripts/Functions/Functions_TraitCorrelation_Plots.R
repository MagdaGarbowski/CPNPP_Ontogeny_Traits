## CPNPP Ontogeny Traits - Functions for correlations between traits at given timepoints
## Associated code: 

# ---------------------------------------Pearson correlation panels--------------------------

pairs_panels_function<-function(df, vars){
  pairs.panels(df[,vars],
               bg=c("blue3", "violetred4","navajowhite2","grey65","coral2","chartreuse4","goldenrod3", "plum4",
                    "grey25","grey45", "grey65", "grey85")[df$SPECIES],  stars = TRUE, main = paste("Harvest", df$H_num))} 

# ---------------------------- Get correlation coeffecients at given timepoints --------------------

Cor_function<-function(df){
  SLA_LDMC<-cor(df$value.SLA, df$value.LDMC, use = "complete.obs")
  SLA_SRL<-cor(df$value.SLA, df$value.SRL, use = "complete.obs")
  SLA_RDMC<-cor(df$value.SLA, df$value.RDMC, use = "complete.obs")
  SLA_RTD<-cor(df$value.SLA, df$value.RTD, use = "complete.obs")
  SLA_RGR<-cor(df$value.SLA, df$RGR_Tot_ln, use = "complete.obs")
  SLA_RER<-cor(df$value.SLA, df$RER_ln, use = "complete.obs")
  LDMC_SRL<-cor(df$value.LDMC, df$value.SRL, use = "complete.obs")
  LDMC_RDMC<-cor(df$value.LDMC, df$value.RDMC, use = "complete.obs")
  LDMC_RTD<-cor(df$value.LDMC, df$value.RTD, use = "complete.obs")
  SRL_RDMC<-cor(df$value.SRL, df$value.RDMC, use = "complete.obs")
  SRL_RTD<-cor(df$value.SRL, df$value.RTD, use = "complete.obs")
  SRL_RGR<-cor(df$value.SRL, df$RGR_Tot_ln, use = "complete.obs")
  SRL_RER<-cor(df$value.SRL, df$RER_ln, use = "complete.obs")
  RDMC_RTD<-cor(df$value.RDMC, df$value.RTD, use = "complete.obs")
  RDMC_RGR<-cor(df$value.RDMC, df$RGR_Tot_ln, use = "complete.obs")
  RDMC_RER<-cor(df$value.RDMC, df$RER_ln, use = "complete.obs")
  RTD_RGR<-cor(df$value.RTD, df$RGR_Tot_ln, use = "complete.obs")
  RTD_RER<-cor(df$value.RTD, df$RER_ln, use = "complete.obs")
  RER_RGR<-cor(df$RER_ln, df$RGR_Tot_ln, use = "complete.obs")
  out<-cbind(SLA_LDMC,SLA_SRL,SLA_RDMC,SLA_RTD,SLA_RGR,SLA_RER,LDMC_SRL,LDMC_RDMC,LDMC_RTD,SRL_RDMC,SRL_RTD,SRL_RGR,SRL_RER,
             RDMC_RTD, RDMC_RGR, RDMC_RER, RTD_RGR,RTD_RER,RER_RGR)
  out
}

Cor_function_H1_Grass<-function(df){
  SLA_LDMC<-cor(df$value.SLA, df$value.LDMC, use = "complete.obs")
  SLA_SRL<-cor(df$value.SLA, df$value.SRL, use = "complete.obs")
  SLA_RDMC<-cor(df$value.SLA, df$value.RDMC, use = "complete.obs")
  SLA_RTD<-cor(df$value.SLA, df$value.RTD, use = "complete.obs")
  LDMC_SRL<-cor(df$value.LDMC, df$value.SRL, use = "complete.obs")
  LDMC_RDMC<-cor(df$value.LDMC, df$value.RDMC, use = "complete.obs")
  LDMC_RTD<-cor(df$value.LDMC, df$value.RTD, use = "complete.obs")
  SRL_RDMC<-cor(df$value.SRL, df$value.RDMC, use = "complete.obs")
  SRL_RTD<-cor(df$value.SRL, df$value.RTD, use = "complete.obs")
  RDMC_RTD<-cor(df$value.RDMC, df$value.RTD, use = "complete.obs")
  out<-cbind(SLA_LDMC,SLA_SRL,SLA_RDMC,SLA_RTD,LDMC_SRL,LDMC_RDMC,LDMC_RTD,SRL_RDMC,SRL_RTD,
             RDMC_RTD)
  out
}

Cor_function_Forb_H1<-function(df){
  SRL_RDMC<-cor(df$value.SRL, df$value.RDMC, use = "complete.obs")
  SRL_RTD<-cor(df$value.SRL, df$value.RTD, use = "complete.obs")
  RDMC_RTD<-cor(df$value.RDMC, df$value.RTD, use = "complete.obs")
  out<-cbind(SRL_RDMC,SRL_RTD,RDMC_RTD)
  out
}

#---------------------------------------- Plotting functions ------------------------------------------

cor_plot<-function(df, var1, var2,xlim, ylim){
  ggplot(df) +
    labs(x = var1, y = var2) +
    geom_point(aes(x = df[[var1]], y = df[[var2]], alpha = 0.2, color = H_num, fill = H_num))+
    stat_ellipse(geom = "polygon", alpha = 0.2, aes(x = df[[var1]], y = df[[var2]], color = H_num, fill = H_num))+
    stat_ellipse(geom ="polygon", alpha = 0.2, aes (x = df[[var1]], y = df[[var2]], linetype = "All"), color = "black") + 
    ylim (ylim) +
    xlim (xlim) + 
    scale_fill_manual("", 
                      guide = FALSE, 
                      breaks = c("H1","H2","H3","H4"),
                      values = c("#E69F00", "#CC79A7", "#56B4E9", "#009E73"))+
    scale_color_manual("", 
                       guide = FALSE,
                       breaks = c("H1","H2","H3","H4"),
                       values = c("#E69F00", "#CC79A7", "#56B4E9", "#009E73"))+
    scale_linetype_manual("", 
                          guide = FALSE,
                          values = c("dashed"))+
    scale_alpha(guide = 'none')+
    theme(legend.position = "none") +
    theme_bw() 
}


cor_plot_all_H2H4<-function(df, var1, var2, xlim, ylim){
  ggplot(df) +
    labs(x = var1, y = var2)+
    geom_point(aes(x = df[[var1]], y = df[[var2]], alpha = 0.2, color = H_num, fill = H_num))+
    stat_ellipse(geom = "polygon", alpha = 0.2, aes(x = df[[var1]], y = df[[var2]], color = H_num, fill = H_num))+
    stat_ellipse(geom ="polygon", alpha = 0.2, aes (x = df[[var1]], y = df[[var2]], linetype = "All"), color = "black") + 
    ylim (ylim) +
    xlim (xlim) + 
    scale_fill_manual("", 
                      guide = FALSE, 
                      breaks = c("H2","H3","H4"),
                      values = c("#CC79A7", "#56B4E9", "#009E73"))+
    scale_color_manual("", 
                       guide = FALSE,
                       breaks = c("H2","H3","H4"),
                       values = c("#CC79A7", "#56B4E9", "#009E73"))+
    scale_linetype_manual("", 
                          guide = FALSE,
                          values = c("dashed"))+
    scale_alpha(guide = 'none')+
    theme(legend.position = "none") +
    theme_bw() 
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
