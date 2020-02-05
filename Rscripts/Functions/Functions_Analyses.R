# CPNPP Ontogenty Traits
# ---------------------------------- Functions (Analyses) --------------------------------------

# Posteriors to dataframe functions 
# Traits through time
df_function<-function(df){
  df<-as.data.frame(summary(df, digits = 4, probs = c(0.1,0.2, 0.5, 0.8, 0.9)))
  df$names<-as.factor(rownames(df))
  tmp<-df[!grepl("Intercept|sigma|mean_PPD|log-posterior", df$names),]
  levels(tmp$names)[levels(tmp$names) %in% c("SPECIESACMI","SPECIESARTR","SPECIESELTR","SPECIESHEAN",
                                             "SPECIESHECO","SPECIESHEVI","SPECIESMACA","SPECIESMUPO",
                                             "SPECIESPAMU","SPECIESPLPA","SPECIESVUOC")] <- c("SPECIESACMI:H_numH1","SPECIESARTR:H_numH1","SPECIESELTR:H_numH1",
                                                                                              "SPECIESHEAN:H_numH1","SPECIESHECO:H_numH1","SPECIESHEVI:H_numH1","SPECIESMACA:H_numH1",
                                                                                              "SPECIESMUPO:H_numH1","SPECIESPAMU:H_numH1","SPECIESPLPA:H_numH1","SPECIESVUOC:H_numH1")
  tmp2<-tmp
  levels(tmp2$names)[levels(tmp2$names) %in% c("H_numH2","H_numH3","H_numH4")]<-c("SPECIESACMI:H_numH2","SPECIESACMI:H_numH3","SPECIESACMI:H_numH4")
  tmp3<-strsplit(as.character(tmp2$names), ":")
  tmp4<-do.call(rbind, tmp3)
  colnames(tmp4)<- c("Species","Harvest")
  tmp5<-cbind(tmp4,tmp2)
  tmp5$Species<-gsub("SPECIES","",tmp5$Species)
  names(tmp5)[names(tmp5) %in% c( "10%","20%","50%","80%","90%")]<-c("CI_10","CI_20","CI_50","CI_80","CI_90")
  tmp5<-tmp5[,!names(tmp5) %in% c("n_eff","Rhat","names")]
  return(tmp5)
}

# Plasticity data frame function 
dat_frame_function_PI<-function(df){
  df<-as.data.frame(summary(df, digits = 4, probs = c(0.1,0.2, 0.5, 0.8, 0.9)))
  df$Species<-as.factor(rownames(df))
  tmp<-df[!grepl(c("sigma|log-fit_ratio|R2|mean_PPD|log-posterior"), df$Species),]
  tmp$Species<-gsub("Species","",tmp$Species)
  names(tmp)[names(tmp) %in% c( "10%","20%","50%","80%","90%")]<-c("CI_10","CI_20","CI_50","CI_80","CI_90")
  return(tmp)
}

Median_CI_calc_function<-function(df){
  d<-reshape(df, idvar = "Species", timevar = "Harvest", direction = "wide")
  d$CI_10_cal_h1<-((d$CI_10.H_numH1))
  d$CI_20_cal_h1<-((d$CI_20.H_numH1))
  d$CI_50_cal_h1<-((d$CI_50.H_numH1))
  d$CI_80_cal_h1<-((d$CI_80.H_numH1))
  d$CI_90_cal_h1<-((d$CI_90.H_numH1))
  d$CI_10_cal_h2<-((d$mean.H_numH1)+(d$CI_10.H_numH2))
  d$CI_20_cal_h2<-((d$mean.H_numH1)+(d$CI_20.H_numH2))
  d$CI_50_cal_h2<-((d$mean.H_numH1)+(d$CI_50.H_numH2))
  d$CI_80_cal_h2<-((d$mean.H_numH1)+(d$CI_80.H_numH2))
  d$CI_90_cal_h2<-((d$mean.H_numH1)+(d$CI_90.H_numH2))
  d$CI_10_cal_h3<-((d$mean.H_numH1)+(d$CI_10.H_numH3))
  d$CI_20_cal_h3<-((d$mean.H_numH1)+(d$CI_20.H_numH3))
  d$CI_50_cal_h3<-((d$mean.H_numH1)+(d$CI_50.H_numH3))
  d$CI_80_cal_h3<-((d$mean.H_numH1)+(d$CI_80.H_numH3))
  d$CI_90_cal_h3<-((d$mean.H_numH1)+(d$CI_90.H_numH3))
  d$CI_10_cal_h4<-((d$mean.H_numH1)+(d$CI_10.H_numH4))
  d$CI_20_cal_h4<-((d$mean.H_numH1)+(d$CI_20.H_numH4))
  d$CI_50_cal_h4<-((d$mean.H_numH1)+(d$CI_50.H_numH4))
  d$CI_80_cal_h4<-((d$mean.H_numH1)+(d$CI_80.H_numH4))
  d$CI_90_cal_h4<-((d$mean.H_numH1)+(d$CI_90.H_numH4))
  d2<-d[names(d) %in% c("Species",
                        "CI_10_cal_h1","CI_20_cal_h1","CI_50_cal_h1","CI_80_cal_h1","CI_90_cal_h1",
                        "CI_10_cal_h2","CI_20_cal_h2","CI_50_cal_h2","CI_80_cal_h2","CI_90_cal_h2",
                        "CI_10_cal_h3","CI_20_cal_h3","CI_50_cal_h3","CI_80_cal_h3","CI_90_cal_h3",
                        "CI_10_cal_h4","CI_20_cal_h4","CI_50_cal_h4","CI_80_cal_h4","CI_90_cal_h4")]
  d2_long<- gather(d2, ConfInt, value, CI_10_cal_h1:CI_90_cal_h4, factor_key=TRUE)
  d3_long<-strsplit(as.character(d2_long$ConfInt), "_")
  d4_long<-do.call(rbind, d3_long)
  colnames(d4_long)<- c("CI","Prob","Cal","Harvest")
  d5_long<-cbind(d4_long,d2_long)
  d6_long<-d5_long[!names(d5_long) %in% c("CI","Cal","ConfInt")]
  data_wide <- spread(d6_long, Prob, value)
  names(data_wide)[names(data_wide) %in% c( "10","20","50","80","90")]<-c("CI_10","CI_20","CI_50","CI_80","CI_90")
  data_wide<-data_wide[,!names(data_wide) %in% c("n_eff","Rhat","names")]
  data_wide
}

