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
