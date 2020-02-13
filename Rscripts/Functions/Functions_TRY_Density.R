# Functions for CPNPP_Ontogeny_Density_TRY_Figures
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
