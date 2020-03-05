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

########################  Difference function - Species x H_num Interactions #######################
diff_by_group = function(samples,
                         name_var_num = 1, diff_var_num = 2, split_char = "_",
                         use_names = TRUE, names = NULL)
{
    nms = do.call(rbind, strsplit(colnames(samples), split_char))
    grps = factor(nms[,diff_var_num])

    if(use_names)
        names = levels(factor(nms[,name_var_num]))
    
    ans = lapply(seq_along(levels(grps)), function(i){
        pairwise_diff(samples[,as.integer(grps) == i], names)
    })
        
    names(ans) = levels(grps)
    ans
}

pairwise_diff = function(samples, var_names = NULL)
{
    i = combn(seq(ncol(samples)), 2, simplify = FALSE)
    ans = lapply(i, function(j) samples[,j[1]] - samples[,j[2]])
    if(!is.null(var_names)){
        names(ans) = sapply(i, function(j) paste(var_names[j], collapse="-"))
    } else {
        names(ans) = sapply(i, paste, collapse="-")
    }
    ans
}

summarize_diffs = function(x, FUN = quantile, only_sig_zero = FALSE, ...){
    ans = lapply(x, function(y) t(sapply(y, FUN, ...)))
}

includes_zero = function(x)
    all(x < 0) | all(x > 0)

################### Difference function - species at each H_num  ####################
diff_within_H_nums<-function(mod_data){
  dat<-as.data.frame(extract(mod_data, pars = c("alpha", "beta_Hnum","beta_sp","beta_sp_Hnum")))
  colnames(dat)<-c("alpha","H1","H2","H3","ACMI","ARTR","ELTR","HEAN","HECO","HEVI","MACA","MUPO","PAMU", "PLPA", "VUOC", inter_names$ln.LDMC_w_cots)
  out<-data.frame(ACMI_H1 = apply(dat[,c("alpha","H1","ACMI","ACMI_H1")], 1, sum), 
                  ACMI_H2 = apply(dat[,c("alpha","H2","ACMI","ACMI_H2")], 1, sum), 
                  ACMI_H3 = apply(dat[,c("alpha","H3","ACMI","ACMI_H3")], 1, sum),
                  ARTR_H1 = apply(dat[,c("alpha","H1","ARTR","ARTR_H1")], 1, sum), 
                  ARTR_H2 = apply(dat[,c("alpha","H2","ARTR","ARTR_H2")], 1, sum), 
                  ARTR_H3 = apply(dat[,c("alpha","H3","ARTR","ARTR_H3")], 1, sum),
                  ELTR_H1 = apply(dat[,c("alpha","H1","ELTR","ELTR_H1")], 1, sum), 
                  ELTR_H2 = apply(dat[,c("alpha","H2","ELTR","ELTR_H2")], 1, sum), 
                  ELTR_H3 = apply(dat[,c("alpha","H3","ELTR","ELTR_H3")], 1, sum),
                  HEAN_H1 = apply(dat[,c("alpha","H1","HEAN","HEAN_H1")], 1, sum), 
                  HEAN_H2 = apply(dat[,c("alpha","H2","HEAN","HEAN_H2")], 1, sum), 
                  HEAN_H3 = apply(dat[,c("alpha","H3","HEAN","HEAN_H3")], 1, sum),
                  HECO_H1 = apply(dat[,c("alpha","H1","HECO","HECO_H1")], 1, sum), 
                  HECO_H2 = apply(dat[,c("alpha","H2","HECO","HECO_H2")], 1, sum), 
                  HECO_H3 = apply(dat[,c("alpha","H3","HECO","HECO_H3")], 1, sum),
                  HEVI_H1 = apply(dat[,c("alpha","H1","HEVI","HEVI_H1")], 1, sum), 
                  HEVI_H2 = apply(dat[,c("alpha","H2","HEVI","HEVI_H2")], 1, sum), 
                  HEVI_H3 = apply(dat[,c("alpha","H3","HEVI","HEVI_H3")], 1, sum),
                  MACA_H1 = apply(dat[,c("alpha","H1","MACA","MACA_H1")], 1, sum), 
                  MACA_H2 = apply(dat[,c("alpha","H2","MACA","MACA_H2")], 1, sum), 
                  MACA_H3 = apply(dat[,c("alpha","H3","MACA","MACA_H3")], 1, sum),
                  MUPO_H1 = apply(dat[,c("alpha","H1","MUPO","MUPO_H1")], 1, sum), 
                  MUPO_H2 = apply(dat[,c("alpha","H2","MUPO","MUPO_H2")], 1, sum), 
                  MUPO_H3 = apply(dat[,c("alpha","H3","MUPO","MUPO_H3")], 1, sum),
                  PAMU_H1 = apply(dat[,c("alpha","H1","PAMU","PAMU_H1")], 1, sum), 
                  PAMU_H2 = apply(dat[,c("alpha","H2","PAMU","PAMU_H2")], 1, sum), 
                  PAMU_H3 = apply(dat[,c("alpha","H3","PAMU","PAMU_H3")], 1, sum),
                  PLPA_H1 = apply(dat[,c("alpha","H1","PLPA","PLPA_H1")], 1, sum), 
                  PLPA_H2 = apply(dat[,c("alpha","H2","PLPA","PLPA_H2")], 1, sum), 
                  PLPA_H3 = apply(dat[,c("alpha","H3","PLPA","PLPA_H3")], 1, sum),
                  VUOC_H1 = apply(dat[,c("alpha","H1","VUOC","VUOC_H1")], 1, sum), 
                  VUOC_H2 = apply(dat[,c("alpha","H2","VUOC","VUOC_H2")], 1, sum), 
                  VUOC_H3 = apply(dat[,c("alpha","H3","VUOC","VUOC_H3")], 1, sum))
  return(out)
}

######################  Differences function - H_num ############################ 
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
################ Back transform - species at each H_num #####################
# df come from diff_within_H_nums function 

back_function_sp_atHnum = function(df){
  dat_back<-exp(df)
  dat_back_quant <- apply(dat_back,2, quantile, p = c(0.025, 0.5, 0.975))
  dat_back_quant<-as.data.frame(t(dat_back_quant))
  dat_back_quant<-round(dat_back_quant, digits =3)
  dat_back_quant$CI = paste(dat_back_quant[,1], dat_back_quant[,3], sep = ", ")
  dat_back_quant$names<-rownames(dat_back_quant)
  rownames(dat_back_quant)<-NULL
  dat_back_quant<-dat_back_quant[,-c(1,3)]
  dat_back_quant_long<-reshape(dat_back_quant,
          direction = "long",
          varying = names(dat_back_quant[1:2]),
          v.names = "val",
          timevar = "CI",
          times = names(dat_back_quant[1:2]))
  rownames(dat_back_quant_long)<-NULL
  dat_back_quant_long <- dat_back_quant_long[order(dat_back_quant_long$names),]
  y<-strsplit(dat_back_quant_long$names, "_" )
  y <- do.call(rbind, y)
  colnames(y)<-c("SPECIES","H_num")
  z <- cbind(dat_back_quant_long,y)
  z<-reshape (z, 
              idvar = c("SPECIES", "CI"),
              timevar = c("H_num"),
              direction = "wide")
  zz<-z[,c(1,2,4,7,10)]
  return(zz)
}






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

bt_function_sp_Hnum_log<-function(df){
  df_Hnum_sp<-as.data.frame(summary(df, pars = c("alpha","beta_sp_Hnum"), 
                                    probs = c(.05,.5,.95))[["summary"]])
  df_Hnum_sp$names<-c("alpha",inter_names$ln.LDMC_w_cots)
  
  df_bt<-data.frame(matrix(NA, nrow = 34, ncol = 0))
  df_bt$sp_H_num<-c("alpha",inter_names$ln.LDMC_w_cots)
  
  exp_fun_CI05<-function(x){ exp(df_Hnum_sp[1,4] + x)}
  exp_fun_CI50<-function(x){ exp(df_Hnum_sp[1,5] + x)}
  exp_fun_CI95<-function(x){ exp(df_Hnum_sp[1,6] + x)}
  
  df_bt$bt_CI05<-apply(df_Hnum_sp[4],2, exp_fun_CI05)
  df_bt$bt_CI50<-apply(df_Hnum_sp[5],2, exp_fun_CI50)
  df_bt$bt_CI95<-apply(df_Hnum_sp[6],2, exp_fun_CI95)
  
  df_bt_Hnum_Sp<-strsplit(df_bt$sp_H_num,  "_") 
  df_bt_Hnum_Sp<-do.call(rbind, df_bt_Hnum_Sp)
  colnames(df_bt_Hnum_Sp)<-c("Sp","Hnum")
  df_bt_Hnum_Sp<-cbind(df_bt,df_bt_Hnum_Sp )
  df_bt_Hnum_Sp<-df_bt_Hnum_Sp[-1,]
  return(df_bt_Hnum_Sp)
}

bt_function_sp_Hnum<-function(df){
  df_Hnum_sp<-as.data.frame(summary(df, pars = c("alpha","beta_sp_Hnum"), 
                                    probs = c(.05,.5,.95))[["summary"]])
  df_Hnum_sp$names<-c("alpha",inter_names$ln.LDMC_w_cots)
  
  df_bt<-data.frame(matrix(NA, nrow = 34, ncol = 0))
  df_bt$sp_H_num<-c("alpha",inter_names$ln.LDMC_w_cots)
  
  exp_fun_CI05<-function(x){ (df_Hnum_sp[1,4] + x)}
  exp_fun_CI50<-function(x){ (df_Hnum_sp[1,5] + x)}
  exp_fun_CI95<-function(x){ (df_Hnum_sp[1,6] + x)}
  
  df_bt$bt_CI05<-apply(df_Hnum_sp[4],2, exp_fun_CI05)
  df_bt$bt_CI50<-apply(df_Hnum_sp[5],2, exp_fun_CI50)
  df_bt$bt_CI95<-apply(df_Hnum_sp[6],2, exp_fun_CI95)
  
  df_bt_Hnum_Sp<-strsplit(df_bt$sp_H_num,  "_") 
  df_bt_Hnum_Sp<-do.call(rbind, df_bt_Hnum_Sp)
  colnames(df_bt_Hnum_Sp)<-c("Sp","Hnum")
  df_bt_Hnum_Sp<-cbind(df_bt,df_bt_Hnum_Sp )
  df_bt_Hnum_Sp<-df_bt_Hnum_Sp[-1,]
  return(df_bt_Hnum_Sp)
}

