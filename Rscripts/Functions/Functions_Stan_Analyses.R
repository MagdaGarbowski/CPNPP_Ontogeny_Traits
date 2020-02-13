# Functions for CPNPP_Ontogeny_Traits_Stan_Analyses

na.omit.fun = function (df, cols) {
  na.omit(df[,cols])
}

mk_data_function = function(df, y_var){
  mk_dat = list (N = nrow(df), 
                 K = nlevels(as.factor(df$H_num)), 
                 x = as.integer(as.factor(df$H_num)), 
                 y = df[[y_var]], 
                 J = nlevels(as.factor(as.character(df$POP_ID))), 
                 p = as.integer(as.factor(as.character(df$POP_ID))))
  mk_dat
}
