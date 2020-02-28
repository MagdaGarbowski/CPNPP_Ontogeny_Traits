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

mk_data_traitvar_function = function(df){
  mk_dat = list (mu1 = 0,
                 sigma1 = 1, 
                 mu2 = 0, 
                 sigma2 = 1, 
                 N = nrow(df), 
                 K = nlevels(as.factor(df$trait)), 
                 J = nlevels(as.factor(as.character(df$SPECIES))),
                 y = df$value, 
                 P = as.integer(as.factor(as.character(df$trait))),
                 s = as.integer(as.factor(as.character(df$SPECIES))))
  mk_dat
}

mk_data_all_function = function(df){
  mk_dat = list(N = nrow(df), 
                K = nlevels(as.factor(df$H_num)), 
                x = as.integer(as.factor(as.character(df$H_num))), 
                S = nlevels(as.factor(df$SPECIES)),
                s = as.integer(as.factor(as.character(df$SPECIES))),
                y = df$value, 
                J = nlevels(as.factor(as.character(df$POP_ID))), 
                p = as.integer(as.factor(as.character(df$POP_ID))))
}

make_matrix_function<-function(df){
  mat_dat = list(
    M = model.matrix(~SPECIES * H_num, df),
    y = df$value,
    n = nrow (model.matrix(~SPECIES * H_num, df)),
    s = ncol (model.matrix(~SPECIES * H_num, df)))
}

################################################################################
## Matt's functions for prepping data

to_index = function(x) as.integer(as.factor(x))

prep_data = function(df){
    Hnum = to_index(df$H_num)
    sp = to_index(df$SPECIES)
    pop = to_index(df$POP_ID)
    sp_Hnum = to_index(paste(df$SPECIES, df$H_num))
    list(N = nrow(df),
         n_Hnum = max(Hnum),
         n_sp = max(sp),
         n_pop = max(pop),
         n_sp_Hnum = max(sp_Hnum),
         Hnum = Hnum,
         sp = sp,
         pop = pop,
         sp_Hnum = sp_Hnum,
         y = df$value)
}

# -------------------------- Functions to estimate differences among groups ----------------------
Hnum_difference_function<-function(df){
  df_Hnum<-as.data.frame(rstan::extract(df, pars = "beta_Hnum"))
  df_Hnum$H1_neg<-df_Hnum[,1]*-1
  df_Hnum$H2_neg<-df_Hnum[,2]*-1
  df_Hnum$H3_neg<-df_Hnum[,3]*-1
  df_Hnum$H4_neg<-df_Hnum[,4]*-1
  
  df_Hnum$diff_h1_h2<-rowSums(df_Hnum[,c(2,5)])
  df_Hnum$diff_h1_h3<-rowSums(df_Hnum[,c(3,5)])
  df_Hnum$diff_h1_h4<-rowSums(df_Hnum[,c(4,5)])
  df_Hnum$diff_h2_h3<-rowSums(df_Hnum[,c(3,6)])
  df_Hnum$diff_h2_h4<-rowSums(df_Hnum[,c(4,6)])
  df_Hnum$diff_h3_h4<-rowSums(df_Hnum[,c(4,7)])
  
  dff<-data.frame(matrix(NA, nrow = 2, ncol = 0))
  
  dff$h1h2_90<-hdi(df_Hnum$diff_h1_h2, credMass = 0.9)
  dff$h1h3_90<-hdi(df_Hnum$diff_h1_h3, credMass = 0.9)
  dff$h1h4_90<-hdi(df_Hnum$diff_h1_h4, credMass = 0.9)
  dff$h2h3_90<-hdi(df_Hnum$diff_h2_h3, credMass = 0.9)
  dff$h2h4_90<-hdi(df_Hnum$diff_h2_h4, credMass = 0.9)
  dff$h3h4_90<-hdi(df_Hnum$diff_h3_h4, credMass = 0.9)
  
  dff$h1h2_80<-hdi(df_Hnum$diff_h1_h2, credMass = 0.8)
  dff$h1h3_80<-hdi(df_Hnum$diff_h1_h3, credMass = 0.8)
  dff$h1h4_80<-hdi(df_Hnum$diff_h1_h4, credMass = 0.8)
  dff$h2h3_80<-hdi(df_Hnum$diff_h2_h3, credMass = 0.8)
  dff$h2h4_80<-hdi(df_Hnum$diff_h2_h4, credMass = 0.8)
  dff$h3h4_80<-hdi(df_Hnum$diff_h3_h4, credMass = 0.8)
  return(dff)
  
}

Hnum_difference_function_noH1<-function(df){
  df_Hnum<-as.data.frame(rstan::extract(df, pars = "beta_Hnum"))
  df_Hnum$H1_neg<-df_Hnum[,1]*-1
  df_Hnum$H2_neg<-df_Hnum[,2]*-1
  df_Hnum$H3_neg<-df_Hnum[,3]*-1
  
  df_Hnum$diff_h2_h3<-rowSums(df_Hnum[,c(4,2)])
  df_Hnum$diff_h2_h4<-rowSums(df_Hnum[,c(4,3)])
  df_Hnum$diff_h3_h4<-rowSums(df_Hnum[,c(5,3)])
  
  dff<-data.frame(matrix(NA, nrow = 2, ncol = 0))
  
  dff$h2h3_90<-hdi(df_Hnum$diff_h2_h3, credMass = 0.9)
  dff$h2h4_90<-hdi(df_Hnum$diff_h2_h4, credMass = 0.9)
  dff$h3h4_90<-hdi(df_Hnum$diff_h3_h4, credMass = 0.9)
  
  dff$h2h3_80<-hdi(df_Hnum$diff_h2_h3, credMass = 0.8)
  dff$h2h4_80<-hdi(df_Hnum$diff_h2_h4, credMass = 0.8)
  dff$h3h4_80<-hdi(df_Hnum$diff_h3_h4, credMass = 0.8)
  return(dff)
}

#-------------------Functions for backtransforming parameter values ---------------------

back_trans_function_H_num<-function(df){
  df_Hnum<-as.data.frame(summary(df)[["summary"]])
  df_Hnum<-df_Hnum[c("alpha", "beta_Hnum[1]","beta_Hnum[2]","beta_Hnum[3]","beta_Hnum[4]"),]
  df_bt<-data.frame(matrix(NA, nrow = 4, ncol = 0))
  df_bt$H_num<-c("H1","H2","H3","H4")
  df_bt$bt_CI2.5<-c(exp(df_Hnum[1,4] + df_Hnum[2,4]), 
                    exp(df_Hnum[1,4] + df_Hnum[3,4]),
                    exp(df_Hnum[1,4] + df_Hnum[4,4]),
                    exp(df_Hnum[1,4] + df_Hnum[5,4]))
  df_bt$bt_CI50<-c(exp(df_Hnum[1,6] + df_Hnum[2,6]), 
                   exp(df_Hnum[1,6] + df_Hnum[3,6]),
                   exp(df_Hnum[1,6] + df_Hnum[4,6]),
                   exp(df_Hnum[1,6] + df_Hnum[5,6]))
  df_bt$bt_CI97.5<-c(exp(df_Hnum[1,8] + df_Hnum[2,8]), 
                     exp(df_Hnum[1,8] + df_Hnum[3,8]),
                     exp(df_Hnum[1,8] + df_Hnum[4,8]),
                     exp(df_Hnum[1,8] + df_Hnum[5,8]))
  return(df_bt)
}
