# CPNPP Ontogenty Traits
# ---------------------------------- Functions (Figures) --------------------------------------


# Plot function - traits through time 
plot_function<-function(df, var, y) {
  ggplot(data = df[df$trait %in% c(var),], aes (y = CI_50, x = Species, group = Harvest)) + 
    geom_errorbar(aes(min = CI_10, max = CI_90), position = position_dodge(0.45), size = 0.4, width = 0.2) +
    geom_errorbar(aes(color = Harvest, min = CI_20, max = CI_80), position = position_dodge(0.45), size = 1.2, width = 0) +
    geom_point(aes(fill = Harvest, shape = Harvest), position = position_dodge(0.45), size = 4) + 
    scale_fill_manual(values =  c("#99CCFF", "#6699CC","#336699","#003399"))+
    scale_color_manual(values =  c("#99CCFF", "#6699CC","#336699","#003399")) +
    scale_shape_manual(values = c(21,21,21,21))+
    theme_bw() + 
    theme(axis.text=element_text(size=9))+
    annotate("text", x=1.5, y = y, label = paste(var), size = 5)
}


# Plot function - Plasticity
plot_function_PI<-function(df) {
  ggplot(data = df, aes (y = CI_50, x = Species, group = trait)) + 
    geom_errorbar(aes(min = CI_10, max = CI_90), position = position_dodge(0.65), size = 0.4, width = 0.2) +
    geom_errorbar(aes(color = trait, min = CI_20, max = CI_80), position = position_dodge(0.65), size = 1.2, width = 0) +
    geom_point(aes(fill = trait, shape = trait), position = position_dodge(0.65), size = 4) + 
    scale_fill_manual(values =  c("#66CC66","#336633","#99CCFF", "#6699CC","#336699","#003399","#FFCC66","#CC9933","#996600","#FFCC99"))+
    scale_color_manual(values =  c("#66CC66","#336633","#99CCFF", "#6699CC","#336699","#003399","#FFCC66","#CC9933","#996600","#FFCC99")) +
    scale_shape_manual(values = c(21,21,21,22,22,22,22,23,23,23))+
    theme_bw() + 
    theme(axis.text=element_text(size=9))
}