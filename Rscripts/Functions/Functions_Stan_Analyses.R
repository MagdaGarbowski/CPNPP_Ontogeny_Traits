# Functions for CPNPP_Ontogeny_Traits_Stan_Analyses

na.omit.fun = function (df, cols) {
  na.omit(df[,cols])
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


##################  Differences - H_num main effects (means) #############################

diff_H_num<-function(mod_data){
  dat<-as.data.frame(rstan::extract(mod_data, pars = c("alpha", "beta_Hnum")))
  colnames(dat)<-c("alpha","H1","H2","H3")  
  out<-data.frame(H1 = apply(dat[,c("alpha","H1")], 1, sum), 
                  H2 = apply(dat[,c("alpha","H2")], 1, sum), 
                  H3 = apply(dat[,c("alpha","H3")], 1, sum))
  return(out)
}


#### Fix this##### 

Hnum_difference_function_noH4<-function(df){
  df_Hnum<-df
  df_Hnum$H1_neg<-df_Hnum[,1]*-1
  df_Hnum$H2_neg<-df_Hnum[,2]*-1
  df_Hnum$H3_neg<-df_Hnum[,3]*-1
  
  df_Hnum$diff_h1_h2<-rowSums(df_Hnum[,c(4,2)])
  df_Hnum$diff_h1_h3<-rowSums(df_Hnum[,c(4,3)])
  df_Hnum$diff_h2_h3<-rowSums(df_Hnum[,c(5,3)])
  
  dff<-data.frame(matrix(NA, nrow = 2, ncol = 0))
  
  dff$h1h2_95<-hdi(df_Hnum$diff_h1_h2, credMass = 0.90)
  dff$h1h3_95<-hdi(df_Hnum$diff_h1_h3, credMass = 0.90)
  dff$h2h3_95<-hdi(df_Hnum$diff_h2_h3, credMass = 0.90)

  return(dff)
}

###################### Differences - Species main effects (means)  ####################

diff_among_species <- function(mod_data,
                                spNames = c("ACMI","ARTR","ELTR",
                                            "HEAN","HECO","HEVI",
                                            "MACA","MUPO","PAMU",
                                            "PLPA", "VUOC"))
{
  dat<-as.data.frame(rstan::extract(mod_data,
                                    pars = c("alpha", "beta_sp")))
  
  colnames(dat)<-c("alpha", spNames)  
  combos = expand.grid(list(alpha = "alpha", Sp = spNames),
                       stringsAsFactors = FALSE)
  combos$inter = paste(combos$Sp, sep="_")
  
  out = as.data.frame(lapply(seq(nrow(combos)), function(i) {
    rowSums(dat[,as.character(combos[i,])])
  }))
  
  colnames(out) = combos$inter
  out
}

diff_by_sps = function(samples)
{
  nms = colnames(samples)
  grps = factor(nms)
  ans = lapply(seq_along(levels(grps)), function(i){
    pairwise_diff(samples[,as.integer(grps) == i], names)
  })
  names(ans) = levels(grps)
  ans
}


###################### Differences - Species at each H_num  ####################

diff_within_H_nums2 <- function(mod_data,
                                spNames = c("ACMI","ARTR","ELTR",
                                            "HEAN","HECO","HEVI",
                                            "MACA","MUPO","PAMU",
                                            "PLPA", "VUOC"))
{
    dat<-as.data.frame(rstan::extract(mod_data,
                                      pars = c("alpha", "beta_Hnum","beta_sp",
                                               "beta_sp_Hnum")))
    
    colnames(dat)<-c("alpha","H1","H2","H3", spNames, inter_names$ln.LDMC_w_cots)  
    combos = expand.grid(list(alpha = "alpha",
                              Hnum = c("H1","H2", "H3"), Sp = spNames),
                         stringsAsFactors = FALSE)
    combos$inter = paste(combos$Sp, combos$Hnum, sep="_")

    out = as.data.frame(lapply(seq(nrow(combos)), function(i) {
        rowSums(dat[,as.character(combos[i,])])
    }))

    colnames(out) = combos$inter
    out
}
 
############### Backtransforming parameter values - Species at each H_num ########
# df comes from diff_within_H_nums function 

back_function_sp_atHnum = function(df){
  dat_back<-exp(df)
  dat_back_quant = apply(dat_back,2, quantile, p = c(0.025, 0.5, 0.975))
  dat_back_quant = as.data.frame(t(dat_back_quant))
  dat_back_quant = round(dat_back_quant, digits = 2)
  dat_back_quant$CI = paste(dat_back_quant[,1], dat_back_quant[,3], sep = ", ")
  dat_back_quant$names = rownames(dat_back_quant)
  rownames(dat_back_quant) = NULL
  dat_back_quant = dat_back_quant[,-c(1,3)]
  dat_back_quant_long = reshape(dat_back_quant,
                               direction = "long",
                               varying = names(dat_back_quant[1:2]),
                               v.names = "val",
                               timevar = "CI",
                               times = names(dat_back_quant[1:2]))
  rownames(dat_back_quant_long) = NULL
  dat_back_quant_long = dat_back_quant_long[order(dat_back_quant_long$names),]
  y = strsplit(dat_back_quant_long$names, "_" )
  yy = do.call(rbind, y)
  colnames(yy) = c("SPECIES","H_num")
  z = cbind(dat_back_quant_long,yy)
  z = reshape (z, 
              idvar = c("SPECIES", "CI"),
              timevar = c("H_num"),
              direction = "wide")
  zz = z[,c(1,2,4,7,10)]
  return(zz)
}

back_function_sp_atHnum_noexp = function(df ){
  dat_back<-(df)
  dat_back_quant = apply(dat_back,2, quantile, p = c(0.025, 0.5, 0.975))
  dat_back_quant = as.data.frame(t(dat_back_quant))
  dat_back_quant = round(dat_back_quant, digits = 2)
  dat_back_quant$CI = paste(dat_back_quant[,1], dat_back_quant[,3], sep = ", ")
  dat_back_quant$names = rownames(dat_back_quant)
  rownames(dat_back_quant) = NULL
  dat_back_quant = dat_back_quant[,-c(1,3)]
  dat_back_quant_long = reshape(dat_back_quant,
                                direction = "long",
                                varying = names(dat_back_quant[1:2]),
                                v.names = "val",
                                timevar = "CI",
                                times = names(dat_back_quant[1:2]))
  rownames(dat_back_quant_long) = NULL
  dat_back_quant_long = dat_back_quant_long[order(dat_back_quant_long$names),]
  y = strsplit(dat_back_quant_long$names, "_" )
  yy = do.call(rbind, y)
  colnames(yy) = c("SPECIES","H_num")
  z = cbind(dat_back_quant_long,yy)
  z = reshape (z, 
               idvar = c("SPECIES", "CI"),
               timevar = c("H_num"),
               direction = "wide")
  zz = z[,c(1,2,4,7,10)]
  return(zz)
}

back_function_sp_atHnum_plotting = function(df){
  dat_back<-exp(df)
  dat_back_quant = apply(dat_back,2, quantile, p = c(0.025, 0.5, 0.975))
  dat_back_quant = as.data.frame(t(dat_back_quant))
  dat_back_quant$names<-rownames(dat_back_quant)
  y = strsplit(dat_back_quant$names, "_" )
  yy = do.call(rbind, y)
  colnames(yy) = c("SPECIES","H_num")
  z = cbind(dat_back_quant,yy)
  colnames(z) = c("ci_lower","median", "ci_upper", "names","SPECIES","H_num")
  z
}

#



