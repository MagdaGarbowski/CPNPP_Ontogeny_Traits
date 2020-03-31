### Data Mgmt Code for Trait Data 2017 
### Janurary 29, 2020
### 2017 Trait data - Traits through time 

# (1) Create dataset with species calculations for SLA, LDMC, RMR, RDMC, RTD, SRL (for each sample)
# (2) Create dataset with growth rates and root elongation (for each population (POP_ID))
# (3) Create dataset with relative change for each trait (for each population (POP_ID))
# (4) Create dataset with Plasticity Index through time (for each species) based on distances between avg. values at each timepoint x species 

setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
library(tidyr)
library(data.table)
source("Rscripts/Functions/Functions_DataMgmt.R")

seed_weights<-read.csv("Data_Raw_TRY_Seedweights/traits_population_seedweights.csv")
seed_weights<-seed_weights[,c("POP_ID","AVG_SEED_WEIGHT")]
seed_weights$Tot_weight_avg_mg <- seed_weights$AVG_SEED_WEIGHT * 1000
seed_weights$Tot_weight_avg_mg_ln<-log(seed_weights$Tot_weight_avg_mg)
seed_weights$H_num = 0.01
seed_weights$RL_avg_cm = 0.3
seed_weights$RL_avg_cm_ln = log(0.2)
seed_weights <- seed_weights[c("H_num","POP_ID","Tot_weight_avg_mg","Tot_weight_avg_mg_ln","RL_avg_cm", "RL_avg_cm_ln")]

# Load data 
filedirectory = "Data_Raw/"
file_list = list.files(path=filedirectory, pattern="\\.csv", full.names = TRUE)
ll = lapply(file_list, read.csv, stringsAsFactor = FALSE)
ll = do.call(rbind, ll)
Trait_data<-(ll)

# Remove "_IN_" from ELTR_NMNW 
Trait_data$SAMPLE_ID<-gsub("_IN_","_",Trait_data$SAMPLE_ID)
Trait_data$SAMPLE_ID<-gsub("ARTRWY","ARTR",Trait_data$SAMPLE_ID)
names(Trait_data)[names(Trait_data) == "H_num"]<-"H_num_og"

Trait_data<-pop_id_function(Trait_data, "SAMPLE_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE","SAMPLE_NUM","H_num","LETTER"))

# -------------------------- Values entered wrong from data sheets (usually extra zeros) ---------------------- 
# Zero values entered wrong 
Trait_data [Trait_data$SAMPLE_ID == "MACA_UTNE_BON_9_H1_O", "LWF_A"] <-NA  # Entered for wrong sample 
Trait_data [Trait_data$SAMPLE_ID == "MACA_UTNE_BON_9_H1_O", "LWS_A"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "MACA_UTNE_BON_9_H1_O", "LWD_A"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "MACA_UTNE_BON_8_H4_O", "LWF_A"] <-0.1269
Trait_data [Trait_data$SAMPLE_ID == "MACA_UTNE_BON_8_H4_O", "LWS_A"] <-0.1270
Trait_data [Trait_data$SAMPLE_ID == "MACA_UTNE_BON_8_H4_O", "LWD_A"] <-0.0179
Trait_data [Trait_data$SAMPLE_ID == "HEAN_AZSE_HWD_1_H1_O", "CWD"] <-0.00048
Trait_data [Trait_data$SAMPLE_ID == "ACMI_UTNW_CCC_1_H1_R", "LWD"] <-0.0001273
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTSC_B_1_H1_O", "LWD"] <-0.0001
Trait_data [Trait_data$SAMPLE_ID == "VUOC_UTSE_IM_4_H2_B", "LWD"] <-0.0024
Trait_data [Trait_data$SAMPLE_ID == "ACMI_UTNW_VE_2_H3_G", "LWD"] <-0.00075
Trait_data [Trait_data$SAMPLE_ID == "ELTR_NMNW_C_18_H1_O", "LWD"] <-0.0021
Trait_data [Trait_data$SAMPLE_ID == "MUPO_AZSE_SIT_8_H3_O", "LWD"] <-0.0005
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTSW_CP_10_H3_G", "LWD"] <-0.001  
Trait_data [Trait_data$SAMPLE_ID == "PAMU_AZNC_WMR_2_H2_R", "CWD"] <-0.0005 
Trait_data [Trait_data$SAMPLE_ID == "PAMU_AZNC_WMR_2_H2_R", "CWS"] <-0.005  
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTEC_MC_10_H2_O", "CWD"] <-0.0004
Trait_data [Trait_data$SAMPLE_ID == "ACMI_UTNW_VE_3_H3_O", "LWS"] <-0.0518 
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTC_MC_6_H2_G", "LWS"] <-0.0120 
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTSW_BJ_11_H2_G", "LWS"] <-0.0120 
Trait_data [Trait_data$SAMPLE_ID == "ELTR_UTSW_DNF_18_H2_R", "RWD"] <-0.014
Trait_data [Trait_data$SAMPLE_ID == "ELTR_UTC_TMC_8_H1_R", "RWD"] <-0.00096
Trait_data [Trait_data$SAMPLE_ID == "HEVI_AZNC_KV_9_H1_G", "RWD"] <-0.000155
Trait_data [Trait_data$SAMPLE_ID == "HEVI_COSW_DE_15_H3_G", "RWD"] <-0.093
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTSE_MO_15_H2_G", "RWD"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTSW_MV_5_H1_G", "RWD"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTSW_OCR_1_H1_O", "RWD"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "ARTR_UTSC_B_13_H1_G", "RWD"] <-.00015
Trait_data [Trait_data$SAMPLE_ID == "ELTR_UTC_TMC_13_H3_B", "RWF"] <-0.013
Trait_data [Trait_data$SAMPLE_ID == "HEVI_COSW_DE_18_H2_O", "RWF"] <-0.00252
Trait_data [Trait_data$SAMPLE_ID == "PLPA_UTSC_BTR_6_H2_R", "RWD"] <-0.0001
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTEC_MC_12_H2_O", "TRANSP_DATE"] <-as.character("8/3")
Trait_data [Trait_data$SAMPLE_ID == "ELTR_UTC_TMC_3_H1_R", "TRANSP_DATE"] <-as.character("10/5")
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTC_WR_1_H1_O", "TRANSP_DATE"] <-as.character("7/31")
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTC_WR_5_H1_O", "TRANSP_DATE"] <-as.character("7/31")
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTSW_CP_1_H1_O", "TRANSP_DATE"] <-as.character("7/27")
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTSW_CP_4_H2_O", "TRANSP_DATE"] <-as.character("7/27")
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTSW_CP_4_H2_G", "HD"] <-as.character("8/21")
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTSW_CP_4_H2_O", "TRANSP_DATE"] <-as.character("7/27")
Trait_data [Trait_data$SAMPLE_ID == "VUOC_AZSW_WB_1_H1_R", "TRANSP_DATE"] <-as.character("8/16")
Trait_data [Trait_data$SAMPLE_ID == "MACA_NMNW_FC_19_H4_G", "TRANSP_DATE"] <-as.character("7/31")
Trait_data [Trait_data$SAMPLE_ID == "HEVI_NMNW_NC_7_H2_O", "HD"] <-as.character("7/27")

# --------Check leaf weight dry (LWD) values and impute where missing and appropriate -----------# 
# Maybe don't do this... often gives negative numbers? Just accept NA? 

is.na(Trait_data$LWD)<-!(Trait_data$LWD) # Make zeros NA - this is okay because will assign based on saturated weight next (i.e. there was actually a leaf there)
table(is.na(Trait_data$LWD[(Trait_data$GrowthForm == "FORB" & Trait_data$H_num == "H1")]))  # out of 120 values 61 are true 
table(is.na(Trait_data$CWD[(Trait_data$GrowthForm == "FORB" & Trait_data$H_num == "H1")]))  # out of 120 values 26 are true 
table(is.na(Trait_data %in% c("CWD", "LWD")[(Trait_data$GrowthForm == "FORB" & Trait_data$H_num == "H1")]))  # out of 120 values 61 are true 
table(is.na(Trait_data$LWD[(Trait_data$GrowthForm == "FORB" & Trait_data$H_num %in% c("H2","H3","H4"))]))  # out of 327 values, 12 have NA, these will be imputed
table(is.na(Trait_data$LWD[(Trait_data$GrowthForm == "GRASS" & Trait_data$H_num == "H1")]))  # out of 65 values, 1 are NA, will be imputed
table(is.na(Trait_data$LWD[(Trait_data$GrowthForm == "GRASS" & Trait_data$H_num %in% c("H2","H3","H4"))]))  # All have values? really??? 

# If SAMPLE_ID has a LWS but not a LWD for GRASS GrowthForm add imputed value from species at the harvest 
# If SAMPLE_ID has a LWS but not a LWD for FORB GrowthForm at H_num = 2,3,4 add imputed value from species at the harvest 

mod_LWD<-lm(LWD ~ LWS + SPECIES + H_num, data = Trait_data)
X_LWD=predict(mod_LWD,Trait_data)
Trait_data$LWD_tmp = ifelse(Trait_data$GrowthForm == "GRASS" & Trait_data$LWS > 0 & (is.na(Trait_data$LWD) | Trait_data$LWD == 0), X_LWD , Trait_data$LWD)
Trait_data$LWD_tmp = ifelse(Trait_data$GrowthForm == "FORB" & Trait_data$H_num %in% c("H2","H3","H4") & Trait_data$LWS >0 & (is.na(Trait_data$LWD) | Trait_data$LWD == 0), X_LWD , Trait_data$LWD)

# --------Check root weight dry (RWD) values and impute where missing and appropriate -----------# 

# If SAMPLE_ID has a RWF (root weight fresh) add imputed value from species at the harvest 

is.na(Trait_data$RWD)<-!Trait_data[(Trait_data$H_num %in% c("H1", "H2", "H3")),]$RWD # Make zeros NA for H1-H3 plants. H4 sometimes did not have root weights taken 
table(is.na(Trait_data$RWD)) # 6 NA for H1-H3, 11 NA for H4 
mod_RWD<-lm(RWD ~ RWF + SPECIES + H_num, data = Trait_data)
X_RWD=predict(mod_RWD, Trait_data)
Trait_data$RWD_tmp = ifelse (Trait_data$RWF>0 & (is.na(Trait_data$RWD) | Trait_data$RWD == 0), X_RWD , Trait_data$RWD) # 11 NA remain from H4 values 

# Predictions resulted in some negative values - remove these
Trait_data$RWD_tmp[Trait_data$RWD_tmp < 0]<-NA
Trait_data$LWD_tmp[Trait_data$LWD_tmp < 0]<-NA

# Remove "late" measurements of PLPA - these obviously didn't grow as the first set and I think they should be removed from analyses
Trait_data<-Trait_data [!(Trait_data$NOTES%in%c("Late")), ]

# ---------------------- Calculate SLA, LDMC, RMR, RDMC, RTD, SRL for "Traits" dataset ----------------------------

# LDMC & SLA  - For H1, H2, H3 all values being used, for H4 leaf values being used unless all leaves scanned and weighed 


# ------------------------------------------------ SLA -------------------------------------------------------------
Trait_data$SLA <-ifelse((Trait_data$H_num  == "H4" & Trait_data$LEAF_TOTAL > 0), 
                        (Trait_data$LEAF_TOTAL/rowSums(Trait_data[,c("LWD_A","LWD_B", "BWD_A")], na.rm=TRUE)), 
                        (Trait_data$LEAVES_TOTAL/Trait_data$LWD_tmp))

# Species cases of SLA in which the wrong leaf scans were divided by the wrong weights (e.g. LEAF_TOTAL/LWD_A should be LEAF_TOTAL/LWD)
# Likely resulted from the wrong identifier in the WinRhizo Scan

# LEAF_TOTAL/LWD
Trait_data$SLA<-ifelse(Trait_data$SAMPLE_ID %in% c("ACMI_UTNW_CCC_2_H2_B", "HECO_UTC_KRRS_10_H4_G", "HECO_UTC_KRRS_14_H4_G",
                                                   "HECO_UTC_TMR_12_H4_G", "HECO_UTC_TMR_13_H4_G", "HECO_UTC_TMR_7_H4_O",
                                                   "HECO_UTEC_GRC_14_H4_G", "HECO_UTEC_GRC_4_H4_O", "HECO_UTEC_GRC_5_H4_O"), 
                       (Trait_data$LEAF_TOTAL/Trait_data$LWD_tmp),
                       (Trait_data$SLA))

# Total.Of.PROJ_AREA_AG/LWD_A + BWD_A
Trait_data$SLA<-ifelse(Trait_data$SAMPLE_ID %in% c("ARTR_UTC_MC_4_H4_O"), 
                       ((Trait_data$Total.Of.PROJ_AREA_AG)/(Trait_data$LWD_A + Trait_data$BWD_A)),
                       (Trait_data$SLA))

# Total.Of.PROJ_AREA_AG/LWD_A
Trait_data$SLA<-ifelse(Trait_data$SAMPLE_ID %in% c("ARTR_UTSW_BJ_9_H4_G","ARTR_UTSW_MV_20_H4_G","MACA_UTSW_BJ_13_H3_G",
                                                   "MUPO_AZSE_BRR_17_H4_O","MUPO_AZSE_BRR_19_H4_G"), 
                       ((Trait_data$Total.Of.PROJ_AREA_AG)/(Trait_data$LWD_A)),
                       (Trait_data$SLA))

# Total.Of.PROJ_AREA_AG/LWD
Trait_data$SLA<-ifelse(Trait_data$SAMPLE_ID %in% c("HEVI_NMNW_NC_5_H3_G","HEVI_NMNW_ND_1_H3_G"), 
                       ((Trait_data$Total.Of.PROJ_AREA_AG)/(Trait_data$LWD)),
                       (Trait_data$SLA))

# LEAVES_TOTAL/LWD
Trait_data$SLA<-ifelse(Trait_data$SAMPLE_ID %in% c("HEVI_NMNW_NC_16_H4_O"), 
                       ((Trait_data$LEAVES_TOTAL)/(Trait_data$LWD)),
                       (Trait_data$SLA))

Trait_data$SLA <- Trait_data$SLA/10

Leaf_vars<-(Trait_data[names(Trait_data) %in% c("SAMPLE_ID", "LEAF_TOTAL", "Total.Of.PROJ_AREA_AG", "LEAF", "LEAVES_TOTAL","LWS","LWD", "LWD_tmp", "LWD_A", "LWD_B","LW_Max_sum","SLA")])
subset(Trait_data, Trait_data[,"SLA"] >100)
Trait_data$SLA[Trait_data$SLA > 400] <- NA 
Trait_data$SLA[Trait_data$SLA == 0] <- NA 

# ---------------------------------------- SLA based on cotyledons and leaves for H1 forbs ------------------------

Trait_data$CWD[Trait_data$CWD == 0] <- NA
Trait_data$LWD[Trait_data$LWD == 0] <- NA
Trait_data$CWD[Trait_data$SWD == 0] <- NA

Trait_data$SLA_w_cots<-ifelse((Trait_data$GrowthForm == "FORB" & Trait_data$H_num %in% c("H1","H2")),
                              ((Trait_data$Total.Of.PROJ_AREA_AG/rowSums(Trait_data[,c("CWD","LWD")], na.rm=TRUE))/10),
                              Trait_data$SLA)

Trait_data$SLA_w_cots[Trait_data$SLA_w_cots > 400] <- NA 
Trait_data$SLA_w_cots[Trait_data$SLA_w_cots == 0] <- NA

# View(Trait_data[,c("SAMPLE_ID","SLA","SLA_w_cots")])
# View(Trait_data[,c("SAMPLE_ID","Total.Of.PROJ_AREA_AG","COTS_TOTAL", "LEAF_TOTAL","LEAF","LEAVES_TOTAL","CWD","LWD","SAMPLE_ID","SLA","SLA_w_cots")])


# ---------------------------------------------- LDMC -------------------------------------------------------------
# Select "max" weight from LWF_A or LWS_A, LWF_B or LWS_B), BWS or BWF... if NA on those three then do LWS and LWD 

Trait_data$LW_Max_A<-pmax(Trait_data$LWS_A,Trait_data$LWF_A, na.rm = TRUE)
Trait_data$LW_Max_B<-pmax(Trait_data$LWS_B,Trait_data$LWF_B, na.rm = TRUE)
Trait_data$BW_Max<-pmax(Trait_data$BWS_A,Trait_data$BWF_A, na.rm = TRUE)

Trait_data$LW_Max_Sum<-rowSums(Trait_data[,c("LW_Max_A","LW_Max_B", "BW_Max")], na.rm = TRUE)
Trait_data$LW_Max_Sum[Trait_data$LW_Max_Sum == 0] <- NA

Trait_data$Grass_LW_Max<-pmax(Trait_data$L.SFW, Trait_data$LWS, na.rm = TRUE)

Trait_data$LDMC <- ifelse((Trait_data$GrowthForm == "GRASS" & is.na(Trait_data[,c("LW_Max_Sum")])), 
                          (Trait_data$LWD/Trait_data$Grass_LW_Max), NA)


Trait_data$LDMC <- ifelse((is.na(Trait_data[,c("LW_Max_Sum")]) & is.na(Trait_data[,c("LDMC")])),
                          ((Trait_data$LWD)/(Trait_data$LWS)), (Trait_data$LDMC))
                          
Trait_data$LDMC <- ifelse((Trait_data[,c("LW_Max_Sum")]>0 & is.na(Trait_data[,c("LDMC")])),
                          ((rowSums(Trait_data[,c("LWD_A","LWD_B","BWD_A")], na.rm = TRUE))/(Trait_data$LW_Max_Sum)), 
                          (Trait_data$LDMC))

Trait_data$LDMC<-ifelse(Trait_data$SAMPLE_ID %in% c("ARTR_UTSW_MV_20_H4_G"), 
                       (Trait_data$LWD_A/Trait_data$LWS_A),
                       (Trait_data$LDMC))

Trait_data$LDMC<-ifelse(Trait_data$SAMPLE_ID %in% c("ACMI_UTNW_CCC_4_H4_R"), 
                        (Trait_data$LWD/Trait_data$L.SFW),
                        (Trait_data$LDMC))
                        
Leaf_vars<-(Trait_data[names(Trait_data) %in% c("SAMPLE_ID","LWS","L.SFW","LWD","LWS_A","LWS_B","BWS_A","LWD_tmp", "LWD_A", "LWD_B","BWD_A","LW_Max_Sum","LDMC")])
Trait_data$LDMC[Trait_data$LDMC > 0.8] <- NA 
Trait_data$LDMC<-Trait_data$LDMC

# ------------------------------------------ LDMC with cotelydons for H1 and H2 Forbs ---------------------
Trait_data$LDMC_w_cots<-ifelse((Trait_data$GrowthForm == "FORB" & Trait_data$H_num %in% c("H1","H2")),
                              ((rowSums(Trait_data[,c("CWD","LWD")], na.rm = TRUE))/(rowSums(Trait_data[,c("CWS","LWS")], na.rm = TRUE))),
                              (Trait_data$LDMC))


# View(Trait_data[,c("SAMPLE_ID","SLA","SLA_w_cots")])
# View(Trait_data[,c("SAMPLE_ID","Total.Of.PROJ_AREA_AG","CWS","LWS","CWD","LWD","SAMPLE_ID","LDMC","LDMC_w_cots")])
Trait_data$LDMC_w_cots[Trait_data$LDMC_w_cots > .5] <- NA 

# ------------------------------------------ SRL (m g^-1) -------------------------------------------------

Trait_data$SRL <- ((Trait_data$SumOfLength.cm./Trait_data$RWD))/100 # SRL - Specific root length
# Trait_data$SRL[Trait_data$SRL > 1300] <- NA 

#------------------------------------------ RDMC  -------------------------------------------------

Trait_data$RDMC <- ((Trait_data$RWD)/(Trait_data$RWF))  # RDMC (Root Dry Matter Content)

Trait_data [Trait_data$SAMPLE_ID == "MUPO_AZSE_BRR_5_H1_G", "RDMC"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "MUPO_AZSE_SIT_3_H1_O", "RDMC"] <-NA
Trait_data$RDMC[Trait_data$RDMC > .4] <- NA 

#---------------------------------------------RTD-------------------------------------------

Trait_data$RTD <- (Trait_data$RWD)/(Trait_data$SumOfRootVolume.cm3.) # RTD (Root tissue density - mg / mm ^3) 
Trait_data$RTD[Trait_data$RTD > 0.6] <- NA 

# ----------------------------------- RMR ------------------------------------------------------------------

Trait_data$RMR <- (Trait_data$RWD)/(rowSums(Trait_data[,c( "CWD", "SWD","LWD", "RWD", "LWD_A", "LWD_B", "BWD_A")], na.rm = TRUE)) # RMR (Root Mass Ratio)
Trait_data$RMR[Trait_data$RMR> 0.9] <- NA 

#Root_vars<-(Trait_data[names(Trait_data) %in% c("SAMPLE_ID","SumOfLength.cm.","RWF","RWD","SRL", "RDMC", "RTD", "RMR")])


# ----------------------------------- HEIGHT ------------------------------------------------------------------

Trait_data$HT[Trait_data$HT<0.05] <- 0.1

#Root_vars<-(Trait_data[names(Trait_data) %in% c("SAMPLE_ID","SumOfLength.cm.","RWF","RWD","SRL", "RDMC", "RTD", "RMR")])

# ------------------------------- One sided root area/one sided shoot area -----------------------------------

Trait_data$RASARatio<-(Trait_data$SumOfProjArea.cm2./Trait_data$Total.Of.PROJ_AREA_AG)
Trait_data$RASARatio<-ifelse((Trait_data$H_num == "H4"), NA, Trait_data$RASARatio)
Trait_data$RASARatio[Trait_data$RASARatio > 40] <- NA

# --------------- log transform traits of interest at the individual level --------------------------------------

# Data into long format 
Trait_data_3<-Trait_data[c("SAMPLE_ID","POP_ID","SPECIES","GrowthForm", "H_num",
                                     "SumOfAvgDiam.mm.","HT","SLA_w_cots","SRL","RTD","RASARatio", "RMR", "LDMC_w_cots","RDMC")] 

Trait_data_3_long<-reshape(Trait_data_3, 
                               direction = "long", 
                               varying = names(Trait_data_3)[6:14],
                               v.names = "value",
                               idvar = c("SAMPLE_ID"),
                               timevar = "trait",
                               times = names(Trait_data_3)[6:14])


# ------------------------------ Calculate Population Avgs ------------------------

Trait_avg_pop_df<-lapply(c("SumOfLength.cm.","SumOfAvgDiam.mm.","HT","SLA","LDMC", "RMR","SRL","RDMC","RTD","RASARatio", "SLA_w_cots", "LDMC_w_cots"), function (y) 
  agg_mean_fun(Trait_data, y, "H_num","POP_ID", c("H_num","POP_ID","value","trait")))

Trait_avg_pop_df<-do.call(rbind, Trait_avg_pop_df)

Trait_avg_pop_df_3<-pop_id_function(Trait_avg_pop_df, "POP_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE"))

Trait_avg_pop_df_3<-as.data.frame(Trait_avg_pop_df_3[!colnames(Trait_avg_pop_df_3) %in% c("LOCATION_CODE","POP_CODE")])

# ------------------------------ Calculate Species Avgs ------------------------

Trait_avg_sps_df<-lapply(c("SLA","LDMC", "RMR","SRL","RDMC","RTD"), function (y) 
  agg_mean_fun(Trait_data, y, "H_num","SPECIES", c("H_num","SPECIES","value","trait")))

Trait_avg_sps_df<-do.call(rbind, Trait_avg_sps_df)


# ------------------------Calculate Relative Growth Rates ----------------------

# Sum of aboveground for each sample 
Trait_data$Above_weight_sum <- rowSums(Trait_data[,c("CWD","LWD", "SWD","LWD_A","LWD_B", "BWD_A", "SWD_A")], na.rm = TRUE)
Trait_data$Tot_weight_sum <- rowSums(Trait_data[,c("Above_weight_sum","RWD")], na.rm = TRUE)
Trait_data$Tot_weight_sum_mg<-Trait_data$Tot_weight_sum *1000
Trait_data$Tot_weight_sum_mg[Trait_data$Tot_weight_sum_mg == 0] <- NA

Trait_data$Tot_weight_sum_mg_ln<-log(Trait_data$Tot_weight_sum_mg)
Trait_data$SumOfLength.cm._ln<-log(Trait_data$SumOfLength.cm.)


# Add a "Days" column to dataset to calculate growth rates 
Trait_data$Days <- as.factor(Trait_data$H_num)
levels(Trait_data$Days) <- c("10","24","42","84")

agg_ls_RGR<-lapply(c("Tot_weight_sum_mg","Tot_weight_sum_mg_ln","SumOfLength.cm.", "SumOfLength.cm._ln"), function (y) 
  agg_mean_fun(Trait_data, y, "Days","POP_ID", c("H_num","POP_ID","value","trait")))

# Merge into one dataset
agg_df<-do.call(cbind, agg_ls_RGR)
agg_df<-agg_df[ c(1,2,3,7,11, 15) ] # column names are similar so indexing 

colnames(agg_df)<-c("H_num","POP_ID","Tot_weight_avg_mg","Tot_weight_avg_mg_ln","RL_avg_cm", "RL_avg_cm_ln")
GrowthRates_Traits2017_3<-agg_df
GrowthRates_Traits2017_3$H_num<-as.numeric(as.character(GrowthRates_Traits2017_3$H_num))

GrowthRates_Traits2017_3<-rbind(GrowthRates_Traits2017_3,seed_weights)
GrowthRates_Traits2017_3<-GrowthRates_Traits2017_3[order( GrowthRates_Traits2017_3[,2], GrowthRates_Traits2017_3[,1] ),]
rownames(GrowthRates_Traits2017_3)<-NULL

Growth_calc<- by(GrowthRates_Traits2017_3, GrowthRates_Traits2017_3$POP_ID, function(df){
  df$RER_rel = relative_change_function(df$RL_avg_cm, df$H_num)
  df$GR_rel = relative_change_function(df$Tot_weight_avg_mg, df$H_num)
  #df$RER = growth_function_not_log(df$RL_avg_cm, df$H_num)
  #df$RGR_Tot = growth_function_not_log(df$Tot_weight_avg_mg, df$H_num)
  #df$RGR_Tot_lnvals<-growth_function_not_log(df$Tot_weight_avg_mg_ln, df$H_num)
  df
})

Growth_calc2<-lapply(Growth_calc, shift_down_function_RER)
Growth_calc3<-lapply(Growth_calc2, shift_down_function_GR)

Growth_calc_df<-do.call(rbind, Growth_calc3)
Growth_calc_df$RER_rel2<-ifelse((Growth_calc_df$RER_rel2 <0), 0.001, Growth_calc_df$RER_rel2)
Growth_calc_df$GR_rel2<-ifelse((Growth_calc_df$GR_rel2 <0), 0.001, Growth_calc_df$GR_rel2)

Growth_calc_df$RRER_ln<-log(Growth_calc_df$RER_rel2)
Growth_calc_df$RGR_Tot_ln<-log(Growth_calc_df$GR_rel2)

GrowthRates_Traits2017_3<-pop_id_function(Growth_calc_df, "POP_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE"))
GrowthRates_Traits2017_3$Days<-GrowthRates_Traits2017_3$H_num
levels(GrowthRates_Traits2017_3$H_num)[(GrowthRates_Traits2017_3$H_num) %in% c("10","24","42","84")]<-c("H1","H2","H3","H4")
rownames(GrowthRates_Traits2017_3)<-NULL
GrowthRates_Traits2017_3[!GrowthRates_Traits2017_3$H_num  == 0.01, ]


GrowthRates_Traits2017_3<-GrowthRates_Traits2017_3[c("H_num","POP_ID","Tot_weight_avg_mg","Tot_weight_avg_mg_ln","RL_avg_cm","RL_avg_cm_ln","RER_rel2", "GR_rel2", "RRER_ln", "RGR_Tot_ln")]


# -------------------- Combine growth rate data and trait population data --------------------
Pop_avg_data_2<- reshape(Trait_avg_pop_df_3, idvar = c("POP_ID","H_num","GrowthForm", "SPECIES"), timevar = "trait", direction = "wide")
ag2<-Pop_avg_data_2
# Data transformations to improve normality 
cols<- c("value.SumOfLength.cm.","value.SumOfAvgDiam.mm.","value.HT","value.SLA","value.LDMC","value.RMR","value.SRL","value.RDMC","value.RTD","value.RASARatio", "value.SLA_w_cots","value.LDMC_w_cots")
ag2[cols]<-log(ag2[cols])
colnames(ag2)[5:16]<-c("ln.SumOfLength.cm.","ln.value.SumOfAvgDiam.mm.","ln.HT", "ln.SLA", "ln.LDMC", "ln.RMR", "ln.SRL", "ln.RDMC", "ln.RTD","ln.RASARatio", "ln.SLA_w_cots","ln.LDMC_w_cots")

Pop_avg_data_3<-merge(Pop_avg_data_2,ag2, by = c("H_num","POP_ID", "SPECIES", "GrowthForm"))
Pop_avg_data_wrates<-merge(Pop_avg_data_3,GrowthRates_Traits2017_3, by = c("H_num","POP_ID"))

#Back into long format
Pop_avg_data_wrates_long<-reshape(Pop_avg_data_wrates, 
  direction = "long",
  varying = names(Pop_avg_data_wrates)[c(5:28,33:36)],
  v.names = "value",
  idvar = c("H_num","POP_ID"),
  timevar = "trait",
  times = names(Pop_avg_data_wrates)[c(5:28,33:36)])

Pop_avg_data_wrates_long<-Pop_avg_data_wrates_long[,!names(Pop_avg_data_wrates_long) %in% c("Tot_weight_avg_mg", "Tot_weight_avg_mg_ln", "RL_avg_cm","RL_avg_cm_ln")]

rownames(Pop_avg_data_wrates_long) <- NULL

### Write complete dataset to CSV 
write.csv(Trait_data_w_ln3_long, file = "Data_Generated/TraitData_2017.csv")
write.csv(GrowthRates_Traits2017_3, "Data_Generated/TraitData_GrowthRates_2017.csv")
write.csv(distances_all, "Data_Generated/TraitData_Plasticity_2017.csv")
write.csv(distances_all_10day, "Data_Generated/TraitData_Plasticity_2017_10day.csv")
write.csv(Trait_avg_sps_df, "Data_Generated/TraitData_SpeciesAvg_2017.csv")
write.csv(Trait_avg_pop_df_3, "Data_Generated/TraitData_PopAvg_2017.csv")
write.csv(Rel_cal_df_3, "Data_Generated/TraitData_RelativeDiff_2017.csv")
write.csv(Pop_avg_data_wrates_long, "Data_Generated/TraitData_PopAvg_wRates_2017.csv")

######## Growth plots 
agg_gr<-aggregate(x = GrowthRates_Traits2017_3$Tot_weight_avg_mg,
                  by = list(GrowthRates_Traits2017_3$SPECIES, GrowthRates_Traits2017_3$H_num), 
                  FUN = mean)
colnames(agg_gr) <- c("Species","Harvest","mg")

GR_plot_ALL<-ggplot(agg_gr, aes(x = Harvest, y = mg))+geom_point() + facet_wrap(.~Species, ncol = 3, scales = "free") # Look at the data 
GR_plot_noH4<-ggplot(agg_gr[!agg_gr$Harvest %in% c("H4"),], aes(x = Harvest, y = mg))+geom_point() + facet_wrap(.~Species, ncol = 3, scales = "free") # Look at the data 


