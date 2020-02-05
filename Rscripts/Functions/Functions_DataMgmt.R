### Functions for data cleaning and calcultating trait values for 2017 CPNPP trait data
# Scripts that use these functions: 

# -------------------------- Functions -------------------------------
# Get mean for appropriate group 
agg_mean_fun<-function(df, var, time, group, colnames){
  tmp<-aggregate(df[[var]], by = list (df[[time]], df[[group]]), 
                 FUN = mean, na.rm = TRUE)
  tmp$trait<-paste(var)
  colnames(tmp)<-colnames
  print(tmp)
}

# Growth rate function for RER and RGR (above and below)
growth_function<- function(trait, days){
  c(NA, diff(log(trait))/diff(days))
}

# Relative change function 
relative_change_function<- function(trait){
  c(NA,diff(trait))
}

# Distance function for Plasticity Index 
dist_fun<- function(df){
  combs <- as.data.frame(as.table(combn(df$value, 2)))
  combs_1 <- reshape(combs, idvar = "Var2", timevar = "Var1", direction = "wide") 
  combs_1$dist <- (abs(combs_1$Freq.A - combs_1$Freq.B))/(combs_1$Freq.A + combs_1$Freq.B)
  combs_1$Trait <- paste((df[1,4]))
  combs_1$Species<- paste((df[1,2]))
  print(combs_1)
}

pop_id_function<-function(df,ID, character, columnnames){
  tmp<-strsplit(as.character(df[[ID]]), character)
  tmp2<-do.call(rbind, tmp)
  colnames(tmp2)<- columnnames
  tmp3<-cbind(df,tmp2)
  tmp3$POP_ID<-paste(tmp3$SPECIES,tmp3$LOCATION_CODE,tmp3$POP_CODE, sep = "_")
  tmp3$H_num<-as.factor(tmp3$H_num)
  tmp3$GrowthForm<- ifelse(tmp3$SPECIES %in% c("ACMI","ARTR","HEAN","HEVI","MACA","PAMU","PLPA"), "FORB", "GRASS")
  tmp3
}
