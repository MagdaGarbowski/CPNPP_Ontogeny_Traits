if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

Pop_avg_data_wrates<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(rstanarm)
library(bayesplot)
library(HDInterval)
library(gridExtra)

bayesplot::color_scheme_set("brightblue")

# Data for full model - mean estimates of traits by species and time 
Pop_avg_data_wrates$H_num<-as.character(as.factor(Pop_avg_data_wrates$H_num))
na.omit_pop<-na.omit.fun(Pop_avg_data_wrates, c("POP_ID", "value","trait","SPECIES", "H_num"))
Traits_all<-na.omit_pop[na.omit_pop$trait %in% c ("ln.RDMC","ln.RTD","ln.RMR","ln.SRL", "ln.value.SumOfAvgDiam.mm.", "ln.SLA_w_cots", "ln.LDMC_w_cots","ln.RASARatio","RER_ln","RGR_Tot_ln"),] 

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

#How to adjust function to get "main effect" differences of H_num? 
#ans<-lapply(interactions, diff_by_group, diff_var = "[A-Z]{4}_")
#ans_out<-lapply(ans, function(x) lapply(x, function(y) t(sapply(y, quantile, p=c(0.05, 0.5, 0.95)))))
#lapply(ans_out, function(x) lapply(x, function(y) y [apply(y,1, includes_zero),]))

################################################################################
diff_within_H_nums_est<-lapply(all_mods_noH4,diff_within_H_nums)

# Now we want to differences for all column combinations... 
SLA_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$ln.SLA_w_cots)
SLA_sp_at_Hnum_vals_CI = lapply(SLA_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(SLA_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])

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

RER_sp_at_Hnum_vals = diff_by_group(diff_within_H_nums_est$RER_ln)
RER_sp_at_Hnum_vals_CI = lapply(RER_sp_at_Hnum_vals, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975))))
lapply(RER_sp_at_Hnum_vals_CI, function(y) y [apply(y,1, includes_zero),])


############################ Backtransform values for each species at each H_num #####################################
bt_w_CI_all<-lapply(diff_within_H_nums_est, back_function_sp_atHnum)
bt_w_CI_all_ord <- bt_w_CI_all[c("ln.SLA_w_cots","ln.LDMC_w_cots","ln.RASARatio","ln.RMR","ln.SRL","ln.RDMC","ln.RTD", "ln.value.SumOfAvgDiam.mm.", "RER_ln","RGR_Tot_ln")]
bt_w_CI_all_out<-do.call(cbind,bt_w_CI_all_ord)

write.csv(bt_w_CI_all_out, "/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Tables/Stan_Medians_CI.csv")
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
########################### Differences at H_num for all traits together ###########################################
H_num_diff_SLA<-Hnum_difference_function_noH4(all_mods_noH4$ln.SLA_w_cots)
H_num_diff_LDMC<-Hnum_difference_function_noH4(all_mods_noH4$ln.LDMC_w_cots)
H_num_diff_RER<-Hnum_difference_function_noH4(all_mods_noH4$RER)
H_num_diff_RGR<-Hnum_difference_function_noH4(all_mods_noH4$RGR_Tot)
H_num_diff_RASA<-Hnum_difference_function_noH4(all_mods_noH4$ln.RASARatio)
H_num_diff_RMR<-Hnum_difference_function_noH4(all_mods_noH4$ln.RMR)
H_num_diff_SRL<-Hnum_difference_function_noH4(all_mods_noH4$ln.SRL)
H_num_diff_RTD<-Hnum_difference_function_noH4(all_mods_noH4$ln.RTD)
H_num_diff_RDMC<-Hnum_difference_function_noH4(all_mods_noH4$ln.RDMC)
H_num_diff_RDIam<-Hnum_difference_function_noH4(all_mods_noH4$ln.value.SumOfAvgDiam.mm.)

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
#
#
##
#
#
#
#
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

pdf("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/Output/Figures/RER_Fig.pdf", height = 10, width = 7)

dev.off()


#
#
#
#
#
#
#
#
#
#
#
#


