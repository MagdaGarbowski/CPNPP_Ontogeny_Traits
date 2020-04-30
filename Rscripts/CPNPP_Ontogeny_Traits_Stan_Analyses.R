if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

Pop_avg_data_wrates<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(rstanarm)
library(bayesplot)
library(HDInterval)
library(gridExtra)
library(ggplot2)
library(cowplot)

# Data for full model - mean estimates of traits by species and time 
Pop_avg_data_wrates$H_num<-as.character(as.factor(Pop_avg_data_wrates$H_num))
na.omit_pop<-na.omit.fun(Pop_avg_data_wrates, c("POP_ID", "value","trait","SPECIES", "H_num"))
Traits_all<-na.omit_pop[na.omit_pop$trait %in% c ("ln.RDMC","ln.RTD","ln.RMR","ln.SRL", "ln.value.SumOfAvgDiam.mm.", "ln.SLA_w_cots", "ln.LDMC_w_cots","ln.RASARatio","RRER_ln","RGR_Tot", "RGR_Tot_ln"),] 

Traits_all_noH4<-Traits_all[!Traits_all$H_num %in% c ("H4"),] 
Trait_splits_noH4 = split(Traits_all_noH4, paste(Traits_all_noH4$trait))
data_noH4<-lapply(Trait_splits_noH4, prep_data)

################################################################################
# Full model (Species and H_num together) 
# Will still need to add random effect of population 
mod = stan_model("stan_models/All_Traits_Random_Effects_HS.stan")
  
mods_function_all <- function(df,
                          mod_file = "stan_models/All_Traits_Random_Effects.stan",
                          iter = 1000,
                          cores = 2,
                          mod = stan_model(mod_file), ...){
  
  sampling(mod, df, iter = iter, cores = cores, ...) 
}

all_mods_noH4 = lapply(data_noH4, mods_function_all, mod = mod,
                       pars = c("beta_Hnum_raw", "beta_sp_raw",
                                "beta_pop_raw", "beta_sp_Hnum_raw"),
                       include = FALSE,
                       warmup = 1000, iter = 1500,
                       control = list(adapt_delta = 0.95)) # Divergent transitions and Bulk Effective Sample Size too low 


################################################################################
# Get names 
# For one: beta_names = levels(factor(paste(Traits_all_noH4$SPECIES, Traits_all_noH4$H_num)))

# So, to get the correct names for each model, run this
inter_names = lapply(Trait_splits_noH4, function(df)
    levels(factor(paste(df$SPECIES, df$H_num, sep = "_"))))

# And to put those names on the model, do this
tmp = lapply(seq_along(inter_names), function(i) {
    x = all_mods_noH4[[i]]
    names(x)[grep("beta_sp_Hnum", names(x))] = inter_names[[i]]
    x
    })

interactions = lapply(all_mods_noH4, function(x) {
    ans = as.data.frame(x)
    ans[,grep("beta_sp_Hnum", colnames(ans))]
    })

# apply the names again
interactions2 = mapply(function(x, nm) {colnames(x) = nm; x},
                      x = interactions, nm = inter_names, SIMPLIFY = FALSE)

# big list of lists with all differences for species x h_num effects 
tmp_sp_at_Hnum = lapply(interactions2, diff_by_group)
## Summarize it all - species at H_num 
ans_sp_at_Hnum = lapply(tmp_sp_at_Hnum, function(x) lapply(x, function(y) t(sapply(y, quantile, p=c(0.05, 0.5, 0.95)))))
lapply(ans_sp_at_Hnum, function(x) lapply(x, function(y) y [apply(y,1, includes_zero),]))

############################## Differences among species within H_num ##################################################
### Question for Matt: Is there a way to display letters or numbers with species for which CI do not contain zero (aka are differnt)
### This would be similar to the cld function in emmeans.

diff_within_H_nums_est<-lapply(all_mods_noH4,diff_within_H_nums2)

# Now we want to differences for all column combinations... 
SLA_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.SLA_w_cots)
SLA_sp_at_Hnum_vals_CI = lapply(SLA_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(SLA_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

SRL_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.SRL)
SRL_sp_at_Hnum_vals_CI = lapply(SRL_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(SRL_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

LDMC_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.LDMC_w_cots)
LDMC_sp_at_Hnum_vals_CI = lapply(LDMC_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(LDMC_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RDMC_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.RDMC)
RDMC_sp_at_Hnum_vals_CI = lapply(RDMC_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RDMC_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RTD_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.RTD)
RTD_sp_at_Hnum_vals_CI = lapply(RTD_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RTD_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RDiam_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.value.SumOfAvgDiam.mm.)
RDiam_sp_at_Hnum_vals_CI = lapply(RDiam_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RDiam_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RMR_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.RMR)
RMR_sp_at_Hnum_vals_CI = lapply(RMR_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RMR_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RASA_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.RASA)
RASA_sp_at_Hnum_vals_CI = lapply(RASA_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RASA_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RGR_Tot_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$RGR_Tot_ln)
RGR_Tot_sp_at_Hnum_vals_CI = lapply(RGR_Tot_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RGR_Tot_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

RER_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$RRER_ln)
RER_sp_at_Hnum_vals_CI = lapply(RER_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RER_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

############################## Differences among H_num within Species ##################################################
diff_within_H_nums_est<-lapply(all_mods_noH4,diff_within_H_nums2)

# Now we want to differences for all column combinations... 
SLA_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.SLA_w_cots, name_var_num = 2, diff_var_num = 1)
SLA_Hnum_at_Sps_vals_CI = lapply(SLA_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(SLA_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

LDMC_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.LDMC_w_cots, name_var_num = 2, diff_var_num = 1)
LDMC_Hnum_at_Sps_vals_CI = lapply(LDMC_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(LDMC_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

SRL_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.SRL, name_var_num = 2, diff_var_num = 1)
SRL_Hnum_at_Sps_vals_CI = lapply(SRL_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(SRL_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RDMC_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.RDMC, name_var_num = 2, diff_var_num = 1)
RDMC_Hnum_at_Sps_vals_CI = lapply(RDMC_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RDMC_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RTD_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.RTD, name_var_num = 2, diff_var_num = 1)
RTD_Hnum_at_Sps_vals_CI = lapply(RTD_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RTD_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RDiam_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.value.SumOfAvgDiam.mm., name_var_num = 2, diff_var_num = 1)
RDiam_Hnum_at_Sps_vals_CI = lapply(RDiam_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RDiam_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RMR_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.RMR, name_var_num = 2, diff_var_num = 1)
RMR_Hnum_at_Sps_vals_CI = lapply(RMR_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RMR_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RASA_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$ln.RASA, name_var_num = 2, diff_var_num = 1)
RASA_Hnum_at_Sps_vals_CI = lapply(RASA_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RASA_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RGR_Tot_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$RGR_Tot_ln, name_var_num = 2, diff_var_num = 1)
RGR_Tot_Hnum_at_Sps_vals_CI = lapply(RGR_Tot_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RGR_Tot_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])

RER_Hnum_at_Sps_vals = diff_by_group(diff_within_H_nums_est$RRER_ln, name_var_num = 2, diff_var_num = 1)
RER_Hnum_at_Sps_vals_CI = lapply(RER_Hnum_at_Sps_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RER_Hnum_at_Sps_vals_CI, function(y) y [apply(y,1, includes_zero),])


##################################  Plotting #######################################################################
sp_level_fun<-function(df){
  dat=df
  dat$SPECIES = factor(df$SPECIES, levels = c("ACMI","ARTR","HEAN","HEVI","MACA","PAMU","PLPA","ELTR","HECO","MUPO","VUOC"))
  dat
}

GrowthForm_fun<-function(df){
  dat = df
  dat$GrowthForm = ifelse(df$SPECIES %in% c("ACMI","ARTR","HEAN","HEVI","MACA","PAMU","PLPA"), "Forb", "Grass")
  dat
}

scaleFUN <- function(x) sprintf("%.2f", x)

bt_w_CI_all_plotting<-lapply(diff_within_H_nums_est, back_function_sp_atHnum_plotting)
bt_w_CI_all_plotting2<-lapply(bt_w_CI_all_plotting,sp_level_fun )
bt_w_CI_all_plotting2<-lapply(bt_w_CI_all_plotting2, GrowthForm_fun)

############################ Table for species comparisons within timepoint ########################################
TrueFalse_CI<-function(df){
  xx = as.data.frame(df)
  xx$zero <- apply(xx, 1, function (x) all(x < 0) | all(x > 0))
  return(xx)
}

lapply(SLA_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(LDMC_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RASA_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RDMC_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RMR_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RTD_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(SRL_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RDiam_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RGR_Tot_sp_at_Hnum_vals_CI, TrueFalse_CI)
lapply(RER_sp_at_Hnum_vals_CI, TrueFalse_CI)


# Long to wide based on H_num to make table of difference at each H_num 
tab_function<-function(df){
  dr<-df[,c("median","SPECIES","H_num")]
  drr<-reshape(dr, idvar = "SPECIES", timevar = "H_num", direction = "wide")
}

list_for_comps<-lapply(bt_w_CI_all_plotting2,tab_function)
table_for_comps<-do.call(cbind, list_for_comps)

write.csv(table_for_comps, "/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/MediansTable.csv")

################################## Back to Plotting #######################################################################

# How to save a list so can open in a new "graphing" script? 
# save(bt_w_CI_all_plotting, file="Data_Generated/Med_CI_trait.RData")

plot_fun_med_ci_stan<-function(df, y_lab,  ...){
ggplot(df, aes (x = SPECIES, y = median, group = H_num, colour = H_num, shape = GrowthForm)) + 
    geom_errorbar(aes (x = SPECIES, ymin = ci_lower, ymax = ci_upper),position = position_dodge(.5), size = .3, width = .2, color = "black" )+
    geom_point(aes(x = SPECIES), position = position_dodge(.5), size = 4.3)+
    scale_shape_manual(name = "Growth Form:", 
                       values = c(20,18))+
    labs(y = y_lab) + 
    scale_fill_manual(name = "Harvest:", 
                      values = c("grey70","grey40","grey25"),
                      labels = c("10-day", "24-day","42-day"))+
    scale_colour_manual(name = "Harvest:", 
                        values = c("grey70","grey40","grey25"), 
                        labels = c("10-day", "24-day","42-day"))+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12), 
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = "hide")
  }

plots_SLA<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.SLA_w_cots, expression(paste("SLA [mm"^{2}, "g"^{-1}, "]"))) + scale_y_continuous(trans='log2') + 
  annotate("text", 
           x = c(1,2,3,4,5,6,7,8,9,10,11),
           y = 165,
           label = c("BCD","BCD","BCD","BC","BCD","BC","B","B","CD","D","A")) + 
  annotate("text", 
           x = c(2,4,7),
           y = 135,
           label = c("a/ab/b","a/a/b","a/ab/b")) + 
  annotate("text", 
           x = 0.75,
           y = 15,
           label = "(a)", 
           size = 5) 

plots_LDMC<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.LDMC, expression(paste("LDMC [mg g"^{-1}, "]"))) + 
  annotate("text", 
           x = c(1,2,3,4,5,6,7,8,9,10,11),
           y = .3,
           label = c("B","B","AB","AB","AB","AB","B","AB","A","AB","AB")) +
  annotate("text", 
           x = c(2,4,7),
           y = .28,
           label = c("b/a/a","b/a/a","b/ab/a"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.05,
           label = "(b)", 
           size = 5) 


plots_SRL<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.SRL, expression(paste("SRL [m g"^{-1}, "]"))) + scale_y_continuous(trans='log2') + 
  annotate("text", 
           x = c(1,2,3,4,5,6,7,8,9,10,11),
           y = 1300,
           label = c("B","BC","CD","BC","BC","B","B","B","D","C","A")) + 
  annotate("text", 
           x = c(1,2,3,5,6,7,9,10),
           y = 1000,
           label = c("a/b/ab","a/b/a","a/b/ab","a/b/ab","a/b/ab","a/b/ab","a/b/ab","a/b/a"))+ 
  annotate("text", 
           x = 0.75,
           y = 115,
           label = "(d)", 
           size = 5) 

plots_RTD<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.RTD, expression(paste("RTD [mg mm"^{-3},"]"))) + 
  annotate("text", 
             x = c(1,2,3,4,5,6,7,8,9,10,11), 
             y = 0.145,
             label = c("AB","AB","AB","AB","AB","C","D","CD","BC","A","D"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.026,
           label = "(f)", 
           size = 5)

plots_RDMC<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.RDMC, expression(paste("RDMC [mg g"^{-1}, "]"))) + 
  annotate("text", 
           x = c(1,2,3,4,5,6,7,8,9,10,11), 
           y = 0.2,
           label = c("AB","AB","B","AB","AB","B","AB","AB","AB","AB","A"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.03,
           label = "(e)", 
           size = 5)

plots_Diam<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.value.SumOfAvgDiam.mm., expression(paste("Diameter [mm]"))) + 
  annotate("text", 
           x = c(1,2,3,4,5,6,7,8,9,10,11), 
           y = 0.45,
           label = c("CD","C","BC","C","C","BC","B","BC","A","D", "CD"))+ 
  annotate("text", 
           x = c(6),
           y = 0.43,
           label = c("a/b/ab"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.14,
           label = "(g)", 
           size = 5)

plots_RMR<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.RMR, "RMR") + 
  annotate("text", 
           x = c(7,11),
           y = 0.67,
           label = c("a/b/b", "a/ab/b"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.11,
           label = "(c)", 
           size = 5)

plots_RASA<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$ln.RASA, "R:S Area") + scale_y_continuous(trans='log2')

plots_RGR<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$RGR_Tot_ln, expression(paste("RGR [mg day"^{-1}, "]")))  + 
  scale_y_continuous(trans='log2', breaks = c(0.004, 0.031,0.250,2.000)) + 
  annotate("text", 
           x = c(2,3,4,5,6,7,8,9,10,11),
           y = 2,
           label = c("a/b/b","a/b/b","a/b/b","a/b/b","a/a/b","a/b/b","a/b/b","a/b/b","a/ab/b", "a/b/b"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.001,
           label = "(h)", 
           size = 5)

plots_RER<-plot_fun_med_ci_stan (bt_w_CI_all_plotting2$RRER_ln, expression(paste("RRER [cm day"^{-1}, "]")))  +
  scale_y_continuous(trans='log2') + 
  annotate("text", 
           x = c(1,2,3,4,5,6,7,8,9,10,11),
           y = 22,
           label = c("a/b/b","a/b/b","a/b/b","a/b/c","a/b/b","a/b/b","a/c/b","a/b/c","a/b/b","a/b/b", "a/b/b"))+ 
  annotate("text", 
           x = 0.75,
           y = 0.02,
           label = "(i)", 
           size = 5)

#legend<-get_legend(plots_SLA)
p_grid<-plot_grid(plots_SLA, plots_LDMC, plots_RMR, plots_SRL,plots_RDMC, plots_RTD,plots_Diam, plots_RGR, plots_RER, nrow=3, align='v', axis = "l")
p_grd_2<-plot_grid(p_grid, legend, ncol =1, rel_heights = c(8,1))

pdf("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Figures/Traits_xSPxHnum.pdf", height = 12, width = 18)
p_grd_2
dev.off()

pdf("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Figures/SLA_LDMC.pdf", height = 10, width = 10)
grid.arrange(plots_SLA, plots_LDMC, ncol = 1)
dev.off()

#
#
############################## Species main effects ##################################################
diff_among_sp_est<-lapply(all_mods_noH4,diff_among_species )

# Now we want to differences for all column combinations... 
quant_sp = function (df){
  quantile(df, p = c(0.025, 0.5, 0.975))
}
SLA_Sps_diff = pairwise_diff(diff_among_sp_est$ln.SLA_w_cots)
SLA_SpsCI = lapply(SLA_Sps_diff, quant_sp)

LDMC_Sps_diff = pairwise_diff(diff_among_sp_est$ln.LDMC_w_cots)
LDMC_SpsCI = lapply(LDMC_Sps_diff, quant_sp)

SRL_Sps_diff = pairwise_diff(diff_among_sp_est$ln.SRL)
SRL_SpsCI = lapply(SRL_Sps_diff, quant_sp)

RTD_Sps_diff = pairwise_diff(diff_among_sp_est$ln.RTD)
RTD_SpsCI = lapply(RTD_Sps_diff, quant_sp)

RDMC_Sps_diff = pairwise_diff(diff_among_sp_est$ln.RDMC)
RDMC_SpsCI = lapply(RDMC_Sps_diff, quant_sp)

Diam_Sps_diff = pairwise_diff(diff_among_sp_est$ln.value.SumOfAvgDiam.mm.)
Diam_SpsCI = lapply(Diam_Sps_diff, quant_sp)

RMR_Sps_diff = pairwise_diff(diff_among_sp_est$ln.RMR)
RMR_SpsCI = lapply(RMR_Sps_diff, quant_sp)

RASA_Sps_diff = pairwise_diff(diff_among_sp_est$ln.RASARatio)
RASA_SpsCI = lapply(RASA_Sps_diff, quant_sp)

RGR_Sps_diff = pairwise_diff(diff_among_sp_est$RGR_Tot_ln)
RGR_SpsCI = lapply(RGR_Sps_diff, quant_sp)

RER_Sps_diff = pairwise_diff(diff_among_sp_est$RRER_ln)
RER_SpsCI = lapply(RER_Sps_diff, quant_sp)

#
########################### Differences at H_num for all traits together ###########################################
diff_Hnum_est<-lapply(all_mods_noH4, diff_H_num)

H_num_diff_SLA<-Hnum_difference_function_noH4(diff_Hnum_est$ln.SLA_w_cots) # no effect 
H_num_diff_LDMC<-Hnum_difference_function_noH4(diff_Hnum_est$ln.LDMC_w_cots) # no effect 
H_num_diff_RER<-Hnum_difference_function_noH4(diff_Hnum_est$RRER_ln) # yes effect 
H_num_diff_RGR<-Hnum_difference_function_noH4(diff_Hnum_est$RGR_Tot_ln) # no effect 
H_num_diff_RASA<-Hnum_difference_function_noH4(diff_Hnum_est$ln.RASARatio) # no effect 
H_num_diff_RMR<-Hnum_difference_function_noH4(diff_Hnum_est$ln.RMR) # no effect 
H_num_diff_SRL<-Hnum_difference_function_noH4(diff_Hnum_est$ln.SRL) # yes effect 
H_num_diff_RTD<-Hnum_difference_function_noH4(diff_Hnum_est$ln.RTD) # no effect 
H_num_diff_RDMC<-Hnum_difference_function_noH4(diff_Hnum_est$ln.RDMC) # no effect 
H_num_diff_RDIam<-Hnum_difference_function_noH4(diff_Hnum_est$ln.value.SumOfAvgDiam.mm.) # no effect 

#################### Back transformed values for all species together at each H_num ###############################
SLA_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.SLA_w_cots, "SLA")
LDMC_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.LDMC_w_cots, "LDMC")
RASA_bt_Hnum<-back_trans_function_H_num_noH1(all_mods_noH4$ln.RASARatio, "RASA")
RMR_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.RMR, "RMF")

Rdiam_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.value.SumOfAvgDiam.mm., "R.Diam")
RDMC_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.RDMC, "RDMC")
RTD_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.RTD, "RTD")
SRL_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.SRL, "SRL")

RGR_bt_Hnum<-back_function_H_num_noH4(all_mods_noH4$RGR_Tot, "RGR")
RER_bt_Hnum<-back_function_H_num_noH4(all_mods_noH4$RER, "RER")

#

# Plotting H_num main effects ------------------------------------------------------------------------

plot_Hnum <- function(df, yscale, var, traitunits){
  ggplot(df, aes(x = H_num, y = bt_CI50))+
  geom_errorbar(aes(ymin = bt_CI05, ymax = bt_CI95), size =.4, width = .1)+
  geom_point() + 
  ylim(yscale)+ 
  labs(title = var, y = traitunits, x = " ")+
  scale_x_discrete(labels = c("10","24","42"))+
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        title = element_text(size = 9),
        panel.grid = element_blank()
  )
}

SLA_plot<-plot_Hnum(SLA_bt_Hnum, c(25,65), "Specific leaf area", expression(paste("mm "^{2}, "mg"^{-1}))) + annotate("text", x = c(1,2,3), y = c(60,58, 50), label = c("a","ab","b"), size = 4.5)
LDMC_plot<-plot_Hnum(LDMC_bt_Hnum, c(100, 175), "Leaf dry matter content", expression(paste("mg g"^{-1})))  
SRL_plot<-plot_Hnum(SRL_bt_Hnum, c(200,510), "Specific root length", expression(paste("m g"^{-1}))) + annotate("text", x = c(1,2,3), y = c(370,490,385), label = c("b","a","b"), size = 4.5)
RDMC_plot<-plot_Hnum(RDMC_bt_Hnum, c(7,15), "Root dry matter content", expression(paste("mg g"^{-1}))) + annotate("text", x = c(1,2,3), y = c(14.5,12.5,12), label = c("a","ab","b"), size = 4.5)
RTD_plot<-plot_Hnum(RTD_bt_Hnum, c(0.05,0.09), "Root tissue density", expression(paste("mg mm"^{-3})))
RootDiam_plot<-plot_Hnum(Rdiam_bt_Hnum, c(0.225, 0.3), "Root Diameter", "mm")

RMR_plot<-plot_Hnum(RMR_bt_Hnum, c(0.19, 0.39), "Root:Shoot Mass Ratio", " ")
RASA_plot<-plot_Hnum(RASA_bt_Hnum, c(0.5, 1.5), "Root:Shoot Area Ratio", " ")

RER_plot<-plot_Hnum(RER_bt_Hnum, c(0, 8), "Root elongation rate", expression(paste("cm day"^{-1}))) + annotate("text", x = c(1,2,3), y = c(4.3,4.65,7.65), label = c("a","a","b"), size = 4.5)
RGR_plot<-plot_Hnum(RGR_bt_Hnum, c(0,2), "Relative growth rate", expression(paste("mg day"^{-1}))) + annotate("text", x = c(1,2,3), y = c(0.3,0.7,1.7), label = c("a","a","b"), size = 4.5)

# Plotting RER for species ----------------------------------------------------------------------------
df_bt_Hnum_Sp$Sp<-factor(df_bt_Hnum_Sp$Sp, levels = c("ACMI","ARTR","HEVI","PAMU",
                                                      "ELTR","HECO","MUPO","VUOC",
                                                      "HEAN","MACA","PLPA"))

RER_fig<-ggplot(df_bt_Hnum_Sp, aes(x = Hnum, y = bt_CI50)) + geom_point() + 
geom_errorbar(aes(ymin =  bt_CI05, ymax = bt_CI95), size = 0.5, width = 0.3)+ 
scale_x_discrete(labels = c("10","24","42"))+
labs(x = "Days", y = "Root elongation rate (cm/day)")+
facet_wrap(.~Sp, ncol = 4) + 
ggtitle("Root elongation rates")+
theme_bw()



# Exporting plots ------------------------------------------------------------------
pdf("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Figures/Trait_H_num_Figs.pdf", height = 10, width = 7)
grid.arrange(SLA_plot,LDMC_plot, 
             SRL_plot,RDMC_plot,
             RTD_plot, RootDiam_plot,
             RMR_plot,RASA_plot,
             RGR_plot,RER_plot, ncol = 2)
dev.off()

