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

# Difference function
diff_by_group = function(samples, diff_var = "_H[0-9]")
{
    ## browser()
    grps = factor(gsub(diff_var, "", colnames(samples)))
    ans = lapply(seq_along(levels(grps)), function(i){
        pairwise_diff(samples[,as.integer(grps) == i])
    })
    names(ans) = levels(grps)
    ans
}

pairwise_diff = function(samples)
{
    i = combn(seq(ncol(samples)), 2, simplify = FALSE)
    ans = lapply(i, function(j) samples[,j[1]] - samples[,j[2]])
    names(ans) = sapply(i, paste, collapse="-")
    ans
}

includes_zero = function (x){
 all(x<0) | all(x>0) 
}


# -------------------------- Functions to estimate differences among groups ----------------------
Hnum_difference_function_noH4<-function(df){
  df_Hnum<-as.data.frame(rstan::extract(df, pars = "beta_Hnum"))
  df_Hnum$H1_neg<-df_Hnum[,1]*-1
  df_Hnum$H2_neg<-df_Hnum[,2]*-1
  df_Hnum$H3_neg<-df_Hnum[,3]*-1
  
  df_Hnum$diff_h1_h2<-rowSums(df_Hnum[,c(4,2)])
  df_Hnum$diff_h1_h3<-rowSums(df_Hnum[,c(4,3)])
  df_Hnum$diff_h2_h3<-rowSums(df_Hnum[,c(5,3)])
  
  dff<-data.frame(matrix(NA, nrow = 2, ncol = 0))
  
  dff$h1h2_90<-hdi(df_Hnum$diff_h1_h2, credMass = 0.9)
  dff$h1h3_90<-hdi(df_Hnum$diff_h1_h3, credMass = 0.9)
  dff$h2h3_90<-hdi(df_Hnum$diff_h2_h3, credMass = 0.9)
  
  dff$h1h2_80<-hdi(df_Hnum$diff_h1_h2, credMass = 0.8)
  dff$h1h3_80<-hdi(df_Hnum$diff_h1_h3, credMass = 0.8)
  dff$h2h3_80<-hdi(df_Hnum$diff_h2_h3, credMass = 0.8)
  return(dff)
}

Hnum_difference_function_wTRY<-function(df){
  df_Hnum<-as.data.frame(rstan::extract(df, pars = "beta_Hnum"))
  df_Hnum$H1_neg<-df_Hnum[,1]*-1
  df_Hnum$H2_neg<-df_Hnum[,2]*-1
  df_Hnum$H3_neg<-df_Hnum[,3]*-1
  
  df_Hnum$diff_h1_h2<-rowSums(df_Hnum[,c(4,2)])
  df_Hnum$diff_h1_h3<-rowSums(df_Hnum[,c(4,3)])
  df_Hnum$diff_h2_h3<-rowSums(df_Hnum[,c(5,3)])
  
  dff<-data.frame(matrix(NA, nrow = 2, ncol = 0))
  
  dff$h1h2_90<-hdi(df_Hnum$diff_h1_h2, credMass = 0.9)
  dff$h1h3_90<-hdi(df_Hnum$diff_h1_h3, credMass = 0.9)
  dff$h2h3_90<-hdi(df_Hnum$diff_h2_h3, credMass = 0.9)
  
  dff$h1h2_80<-hdi(df_Hnum$diff_h1_h2, credMass = 0.8)
  dff$h1h3_80<-hdi(df_Hnum$diff_h1_h3, credMass = 0.8)
  dff$h2h3_80<-hdi(df_Hnum$diff_h2_h3, credMass = 0.8)
  return(dff)
}


sp_Hnum_diff<-function(datframe){
  df<-datframe
  df$H1_neg<-df[,1]*-1
  df$H2_neg<-df[,2]*-1
  df$H3_neg<-df[,3]*-1
  df$diff_h1_h2<-rowSums(df[,c(4,2)])
  df$diff_h1_h3<-rowSums(df[,c(4,3)])
  df$diff_h2_h3<-rowSums(df[,c(5,3)])
  
  dff<-data.frame(matrix(NA, nrow = 2, ncol = 0))
  
  dff$h1h2_90<-hdi(df$diff_h1_h2, credMass = 0.9)
  dff$h1h3_90<-hdi(df$diff_h1_h3, credMass = 0.9)
  dff$h2h3_90<-hdi(df$diff_h2_h3, credMass = 0.9)
  return(dff)
}

# -------------------- Functions for pulling sp_Hnum samples from mod output------------------

sample_species_function<-function(df){
  dat<-reshape(df, 
               direction = "long",
               varying = names(df)[1:33],
               v.names = "sample",
               idvar = "x",
               timevar = "sp_Hnum",
               times = names(df)[1:33])
  dat$species<-substr(dat$sp_Hnum, start = 1, stop = 4)
  return(dat)
}

samples_wide_fun<-function(df){
  sp_df<-as.data.frame(df[c("sp_Hnum","sample")])
  rownames(sp_df)<-NULL
  sp_df$ID <- rep(1:2000, times = 3)
  sp_df_wide<-reshape(sp_df, 
                      timevar = "sp_Hnum",
                      idvar = "ID",
                      direction = "wide")
  sp_df_wide<-sp_df_wide[-c(1)]
  sp_df_wide
}

#-------------------Functions for backtransforming parameter values ---------------------
# Back transform from log 
back_trans_function_H_num_noH4<-function(df, var){
  df_Hnum<-as.data.frame(summary(df, pars = c("alpha","beta_Hnum"), probs = c(.05,.5,.95))[["summary"]])
  df_bt<-data.frame(matrix(NA, nrow = 3, ncol = 0))
  df_bt$H_num<-c("H1","H2","H3")
  df_bt$bt_CI05<-c(exp(df_Hnum[1,4] + df_Hnum[2,4]), 
                    exp(df_Hnum[1,4] + df_Hnum[3,4]),
                    exp(df_Hnum[1,4] + df_Hnum[4,4]))
  df_bt$bt_CI50<-c(exp(df_Hnum[1,5] + df_Hnum[2,5]), 
                   exp(df_Hnum[1,5] + df_Hnum[3,5]),
                   exp(df_Hnum[1,5] + df_Hnum[4,5]))
  df_bt$bt_CI95<-c(exp(df_Hnum[1,6] + df_Hnum[2,6]), 
                     exp(df_Hnum[1,6] + df_Hnum[3,6]),
                     exp(df_Hnum[1,6] + df_Hnum[4,6]))
  df_bt$trait <- var
  return(df_bt)
}

# Not log transformed traits 
back_function_H_num_noH4<-function(df, var){
  df_Hnum<-as.data.frame(summary(df, pars = c("alpha","beta_Hnum"), probs = c(.05,.5,.95))[["summary"]])
  df_bt<-data.frame(matrix(NA, nrow = 3, ncol = 0))
  df_bt$H_num<-c("H1","H2","H3")
  df_bt$bt_CI05<-c(df_Hnum[1,4] + df_Hnum[2,4], 
                   df_Hnum[1,4] + df_Hnum[3,4],
                   df_Hnum[1,4] + df_Hnum[4,4])
  df_bt$bt_CI50<-c(df_Hnum[1,5] + df_Hnum[2,5], 
                   df_Hnum[1,5] + df_Hnum[3,5],
                   df_Hnum[1,5] + df_Hnum[4,5])
  df_bt$bt_CI95<-c(df_Hnum[1,6] + df_Hnum[2,6], 
                   df_Hnum[1,6] + df_Hnum[3,6],
                   df_Hnum[1,6] + df_Hnum[4,6])
  df_bt$trait <- var
  return(df_bt)
}

