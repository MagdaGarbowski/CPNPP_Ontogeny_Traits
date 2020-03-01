if(Sys.info()["login"] == "MagdaGarbowski")
    setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")

TraitData_Pop_Avg_2017<-read.csv("Data_Generated/TraitData_PopAvg_2017.csv")
Pop_avg_data_wrates<-read.csv("Data_Generated/TraitData_PopAvg_wRates_2017.csv")

source("Rscripts/Functions/Functions_Stan_Analyses.R")
library(rstan)
library(rstanarm)
library(bayesplot)
library(HDInterval)
library(gridExtra)
library(tidyr)

bayesplot::color_scheme_set("brightblue")

Pop_avg_data_wrates$H_num<-as.character(as.factor(Pop_avg_data_wrates$H_num))
# Pop_avg_data_wrates$value<-ifelse((Pop_avg_data_wrates$GrowthForm %in% c("FORB","SHRUB", "GRASS") & Pop_avg_data_wrates$H_num == "H1" & Pop_avg_data_wrates$trait %in% c("ln.SLA", "ln.LDMC", "value.SLA","value.LDMC")), NA, Pop_avg_data_wrates$value)

# Data for variance estimates of traits
# na.omit_traits<-na.omit.fun(TraitData_Pop_Avg_2017, c("POP_ID", "value","trait","SPECIES"))
# traits_stan_df<-mk_data_traitvar_function(na.omit_traits)

# Data for full model - mean estimates of traits by species and time 
na.omit_species<-na.omit.fun(Pop_avg_data_wrates, c("POP_ID", "value","trait","SPECIES", "H_num"))
Traits_all<-na.omit_species[na.omit_species$trait %in% c ("ln.RDMC","ln.RTD","ln.RMR","ln.SRL", "ln.value.SumOfAvgDiam.mm.", "ln.SLA_w_cots", "ln.LDMC_w_cots","ln.RASARatio","RER","RGR_Tot"),] 
Traits_all[Traits_all$trait == "RGR_Tot",]$value<-Traits_all[Traits_all$trait == "RGR_Tot",]$value*1000


Traits_all_noH4<-Traits_all[!Traits_all$H_num %in% c ("H4"),] 
Trait_splits_noH4 = split(Traits_all_noH4, paste(Traits_all_noH4$trait))

################################################################################

data_noH4<-lapply(Trait_splits_noH4, prep_data)


# Full model (Species and H_num together) --------------------------------------#
# Will still need to add random effect of population 
mod = stan_model("stan_models/All_Traits_Random_Effects.stan")
  
mods_function_all <- function(df,
                          mod_file = "stan_models/All_Traits_Random_Effects.stan",
                          iter = 1000,
                          cores = 2,
                          mod = stan_model(mod_file), ...){
  
  sampling(mod, df, iter = iter, cores = cores, ...) 
}

# The new model does not sample as efficiently as the matrix
# representation, so higher adapt delta and iterations are needed

# all_mods_full = lapply(data_full, mods_function_all, mod = mod,
#                       pars = c("beta_Hnum_raw", "beta_sp_raw",
#                                "beta_pop_raw", "beta_sp_Hnum_raw"),
#                       include = FALSE,
#                       warmup = 1500, iter = 2000,
#                       control = list(adapt_delta = 0.95))  # Bulk Effective Sample Size too low 

all_mods_noH4 = lapply(data_noH4, mods_function_all, mod = mod,
                       pars = c("beta_Hnum_raw", "beta_sp_raw",
                                "beta_pop_raw", "beta_sp_Hnum_raw"),
                       include = FALSE,
                       warmup = 1000, iter = 1500,
                       control = list(adapt_delta = 0.95)) # Divergent transitions and Bulk Effective Sample Size too low 

# Estimate differences between H_nums for all traits ----------------------------

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

# Back transformed values for plotting (H_num main effects) --------------------------------------------
SLA_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.SLA_w_cots, "SLA")
LDMC_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.LDMC_w_cots, "LDMC")
RASA_bt_Hnum<-back_trans_function_H_num_noH1(all_mods_noH4$ln.RASARatio, "RASA")
RMR_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.RMR, "RMF")

Rdiam_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.value.SumOfAvgDiam.mm., "R.Diam")
RDMC_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.RDMC, "R.Diam")
RTD_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.RTD, "RTD")
SRL_bt_Hnum<-back_trans_function_H_num_noH4(all_mods_noH4$ln.SRL, "SRL")

RGR_bt_Hnum<-back_function_H_num_noH4(all_mods_noH4$RGR_Tot, "RGR")
RER_bt_Hnum<-back_function_H_num_noH4(all_mods_noH4$RER, "RER")


# Get names  ----------------------------------------------------------# 

## beta_names_noh4<-unique(Traits_all_noH4[,c("SPECIES", "H_num")])
## beta_names_noh4$names_noh4<-paste(beta_names_noh4$SPECIES, beta_names_noh4$H_num, sep = "_")
## beta_names_noh4<- beta_names_noh4[order(beta_names_noh4$names_noh4),]
## beta_names_noh4<-beta_names_noh4$names_noh4

# Not quite right - the order is based on creating a factor but it for
# each split, so if a species or H_num is not present in one, all the
# levels are shifted. This is how you would do for one
# beta_names = levels(factor(paste(Traits_all_noH4$SPECIES, Traits_all_noH4$H_num)))

# So, to get the correct names for each model, run this
inter_names = lapply(Trait_splits_noH4, function(df)
    levels(factor(paste(df$SPECIES, df$H_num, sep = "_"))))

# And to put those names on the model, do this
tmp = lapply(seq_along(inter_names), function(i) {
    x = all_mods_noH4[[i]]
    names(x)[grep("beta_sp_Hnum", names(x))] = inter_names[[i]]
    x
    })

# Estimate differences between H_num for given species x trait  --------------------------
# Come back and figure out how to make this more effecient

SLA_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.SLA_w_cots, pars = "beta_sp_Hnum"))
## names(SLA_sp_Hnum)[grep("beta_sp_Hnum",names(SLA_sp_Hnum))]<-beta_names_noh4
SLA_sp_samples<-sample_species_function(SLA_sp_Hnum)
SLA_sps_Hnum_splits<-split(SLA_sp_samples, paste(SLA_sp_samples$species))

SLA_dfs<-lapply(SLA_sps_Hnum_splits, samples_wide_fun)
SLA_sp_Hnum_diffs<-lapply(SLA_dfs, sp_Hnum_diff)
SLA_sp_Hnum_diffs_df<-do.call(rbind, SLA_sp_Hnum_diffs)
################################################################################
interactions = lapply(all_mods_noH4, function(x) {
    ans = as.data.frame(x)
    ans[,grep("beta_sp_Hnum", colnames(ans))]
    })

# apply the names again
interactions = mapply(function(x, nm) {colnames(x) = nm; x},
                      x = interactions, nm = inter_names, SIMPLIFY = FALSE)

# big list of lists with all differences for each species by every h_num
tmp = lapply(interactions, diff_by_group)

## Summarize it all
lapply(tmp, function(x) lapply(x, function(y) t(sapply(y, quantile, p=c(0.025, 0.5, 0.975)))))

################################################################################



LDMC_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.LDMC_w_cots, pars = "beta_sp_Hnum"))
names(LDMC_sp_Hnum)[grep("beta_sp_Hnum",names(LDMC_sp_Hnum))]<-beta_names_noh4
LDMC_sp_samples<-sample_species_function(LDMC_sp_Hnum)
LDMC_sps_Hnum_splits<-split(LDMC_sp_samples, paste(LDMC_sp_samples$species))

LDMC_dfs<-lapply(LDMC_sps_Hnum_splits, samples_wide_fun)
LDMC_sp_Hnum_diffs<-lapply(LDMC_dfs, sp_Hnum_diff)
LDMC_sp_Hnum_diffs_df<-do.call(rbind, LDMC_sp_Hnum_diffs)

SRL_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.SRL, pars = "beta_sp_Hnum"))
names(SRL_sp_Hnum)[grep("beta_sp_Hnum",names(SRL_sp_Hnum))]<-beta_names_noh4
SRL_sp_samples<-sample_species_function(SRL_sp_Hnum)
SRL_sps_Hnum_splits<-split(SRL_sp_samples, paste(SRL_sp_samples$species))

SRL_dfs<-lapply(SRL_sps_Hnum_splits, samples_wide_fun)
SRL_sp_Hnum_diffs<-lapply(SRL_dfs, sp_Hnum_diff)
SRL_sp_Hnum_diffs_df<-do.call(rbind, SRL_sp_Hnum_diffs)

RDMC_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.RDMC, pars = "beta_sp_Hnum"))
names(RDMC_sp_Hnum)[grep("beta_sp_Hnum",names(RDMC_sp_Hnum))]<-beta_names_noh4
RDMC_sp_samples<-sample_species_function(RDMC_sp_Hnum)
RDMC_sps_Hnum_splits<-split(RDMC_sp_samples, paste(RDMC_sp_samples$species))

RDMC_dfs<-lapply(RDMC_sps_Hnum_splits, samples_wide_fun)
RDMC_sp_Hnum_diffs<-lapply(RDMC_dfs, sp_Hnum_diff)
RDMC_sp_Hnum_diffs_df<-do.call(rbind, RDMC_sp_Hnum_diffs)

RTD_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.RTD, pars = "beta_sp_Hnum"))
names(RTD_sp_Hnum)[grep("beta_sp_Hnum",names(RTD_sp_Hnum))]<-beta_names_noh4
RTD_sp_samples<-sample_species_function(RTD_sp_Hnum)
RTD_sps_Hnum_splits<-split(RTD_sp_samples, paste(RTD_sp_samples$species))

RTD_dfs<-lapply(RTD_sps_Hnum_splits, samples_wide_fun)
RTD_sp_Hnum_diffs<-lapply(RTD_dfs, sp_Hnum_diff)
RTD_sp_Hnum_diffs_df<-do.call(rbind, RTD_sp_Hnum_diffs)

RDiam_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.value.SumOfAvgDiam.mm., pars = "beta_sp_Hnum"))
names(RDiam_sp_Hnum)[grep("beta_sp_Hnum",names(RDiam_sp_Hnum))]<-beta_names_noh4
RDiam_sp_samples<-sample_species_function(RDiam_sp_Hnum)
RDiam_sps_Hnum_splits<-split(RDiam_sp_samples, paste(RDiam_sp_samples$species))

RDiam_dfs<-lapply(RDiam_sps_Hnum_splits, samples_wide_fun)
RDiam_sp_Hnum_diffs<-lapply(RDiam_dfs, sp_Hnum_diff)
RDiam_sp_Hnum_diffs_df<-do.call(rbind, RDiam_sp_Hnum_diffs)

RMR_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.RMR, pars = "beta_sp_Hnum"))
names(RMR_sp_Hnum)[grep("beta_sp_Hnum",names(RMR_sp_Hnum))]<-beta_names_noh4
RMR_sp_samples<-sample_species_function(RMR_sp_Hnum)
RMR_sps_Hnum_splits<-split(RMR_sp_samples, paste(RMR_sp_samples$species))

RMR_dfs<-lapply(RMR_sps_Hnum_splits, samples_wide_fun)
RMR_sp_Hnum_diffs<-lapply(RMR_dfs, sp_Hnum_diff)
RMR_sp_Hnum_diffs_df<-do.call(rbind, RMR_sp_Hnum_diffs)

RASA_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$ln.RASARatio, pars = "beta_sp_Hnum"))
names(RASA_sp_Hnum)[grep("beta_sp_Hnum",names(RASA_sp_Hnum))]<-beta_names_noh4
RASA_sp_samples<-sample_species_function(RASA_sp_Hnum)
RASA_sps_Hnum_splits<-split(RASA_sp_samples, paste(RASA_sp_samples$species))

RASA_dfs<-lapply(RASA_sps_Hnum_splits, samples_wide_fun)
RASA_sp_Hnum_diffs<-lapply(RASA_dfs, sp_Hnum_diff)
RASA_sp_Hnum_diffs_df<-do.call(rbind, RASA_sp_Hnum_diffs)

RGR_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$RGR_Tot, pars = "beta_sp_Hnum"))
names(RGR_sp_Hnum)[grep("beta_sp_Hnum",names(RGR_sp_Hnum))]<-beta_names_noh4
RGR_sp_samples<-sample_species_function(RGR_sp_Hnum)
RGR_sps_Hnum_splits<-split(RGR_sp_samples, paste(RGR_sp_samples$species))

RGR_dfs<-lapply(RGR_sps_Hnum_splits, samples_wide_fun)
RGR_sp_Hnum_diffs<-lapply(RGR_dfs, sp_Hnum_diff)
RGR_sp_Hnum_diffs_df<-do.call(rbind, RGR_sp_Hnum_diffs)

RER_sp_Hnum<-as.data.frame(rstan::extract(all_mods_noH4$RER, pars = "beta_sp_Hnum"))
names(RER_sp_Hnum)[grep("beta_sp_Hnum",names(RER_sp_Hnum))]<-beta_names_noh4
RER_sp_samples<-sample_species_function(RER_sp_Hnum)
RER_sps_Hnum_splits<-split(RER_sp_samples, paste(RER_sp_samples$species))

RER_dfs<-lapply(RER_sps_Hnum_splits, samples_wide_fun)
RER_sp_Hnum_diffs<-lapply(RER_dfs, sp_Hnum_diff)
RER_sp_Hnum_diffs_df<-do.call(rbind, RER_sp_Hnum_diffs)


# Back transformed values for plotting H_num x species (RER only) -------------------------------------------
df_Hnum_sp<-as.data.frame(summary(all_mods_noH4$RER_ln, 
                                  pars = c("alpha",
                                           "ACMI_H1","ACMI_H2","ACMI_H3",
                                           "ARTR_H1","ARTR_H2","ARTR_H3",
                                           "ELTR_H1","ELTR_H2","ELTR_H3",
                                           "HEAN_H1","HEAN_H2","HEAN_H3",
                                           "HECO_H1","HECO_H2","HECO_H3",
                                           "HEVI_H1","HEVI_H2","HEVI_H3",
                                           "MACA_H1","MACA_H2","MACA_H3",
                                           "MUPO_H1","MUPO_H2","MUPO_H3",
                                           "PAMU_H1","PAMU_H2","PAMU_H3",
                                           "PLPA_H1","PLPA_H2","PLPA_H3",
                                           "VUOC_H1","VUOC_H2","VUOC_H3"), 
                                  probs = c(.05,.5,.95))[["summary"]])
df_bt<-data.frame(matrix(NA, nrow = 34, ncol = 0))
df_bt$sp_H_num<-c("alpha",
  "ACMI_H1","ACMI_H2","ACMI_H3",
  "ARTR_H1","ARTR_H2","ARTR_H3",
  "ELTR_H1","ELTR_H2","ELTR_H3",
  "HEAN_H1","HEAN_H2","HEAN_H3",
  "HECO_H1","HECO_H2","HECO_H3",
  "HEVI_H1","HEVI_H2","HEVI_H3",
  "MACA_H1","MACA_H2","MACA_H3",
  "MUPO_H1","MUPO_H2","MUPO_H3",
  "PAMU_H1","PAMU_H2","PAMU_H3",
  "PLPA_H1","PLPA_H2","PLPA_H3",
  "VUOC_H1","VUOC_H2","VUOC_H3")
  
exp_fun_CI05<-function(x){ exp(1.18511 + x)}
exp_fun_CI50<-function(x){ exp(1.19133 + x)}
exp_fun_CI95<-function(x){ exp(1.198069 + x)}

df_bt$bt_CI05<-apply(df_Hnum_sp[4],2, exp_fun_CI05)
df_bt$bt_CI50<-apply(df_Hnum_sp[5],2, exp_fun_CI50)
df_bt$bt_CI95<-apply(df_Hnum_sp[6],2, exp_fun_CI95)
  
df_bt_Hnum_Sp<-strsplit(df_bt$sp_H_num,  "_") 
df_bt_Hnum_Sp<-do.call(rbind, df_bt_Hnum_Sp)
colnames(df_bt_Hnum_Sp)<-c("Sp","Hnum")
df_bt_Hnum_Sp<-cbind(df_bt,df_bt_Hnum_Sp )
df_bt_Hnum_Sp<-df_bt_Hnum_Sp[-1,]


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

# Look at results ---------------------------------------------------# 
names(all_mods_noH4$ln.SLA_w_cots)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.SLA_w_cots))]<-beta_names_noh4
plot(all_mods_noH4$ln.SLA_w_cots, pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.LDMC_w_cots)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.LDMC_w_cots))]<-beta_names_noh4
plot(all_mods_noH4$ln.LDMC_w_cots, pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.SRL)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.SRL))]<-beta_names_noh4
plot(all_mods_noH4$ln.SRL, pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.RDMC)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.RDMC))]<-beta_names_noh4
plot(all_mods_noH4$ln.RDMC, pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.RTD)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.RTD))]<-beta_names_noh4
plot(all_mods_noH4$ln.RTD, pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.value.SumOfAvgDiam.mm.)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.value.SumOfAvgDiam.mm.))]<-beta_names_noh4
plot(all_mods_noH4$ln.value.SumOfAvgDiam.mm., pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.RMR)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.RMR))]<-beta_names_noh4
plot(all_mods_noH4$ln.RMR, pars = "beta_sp_Hnum")

names(all_mods_noH4$ln.RASARatio)[grep("beta_sp_Hnum",names(all_mods_noH4$ln.RASARatio))]<-beta_names_noh4
plot(all_mods_noH4$ln.RASARatio, pars = "beta_sp_Hnum")

names(all_mods_noH4$RGR_Tot)[grep("beta_sp_Hnum",names(all_mods_noH4$RGR_Tot))]<-beta_names_noh4
plot(all_mods_noH4$RGR_Tot, pars = "beta_sp_Hnum")

names(all_mods_noH4$RER)[grep("beta_sp_Hnum",names(all_mods_noH4$RER))]<-beta_names_noh4
plot(all_mods_noH4$RER, pars = "beta_sp_Hnum")

