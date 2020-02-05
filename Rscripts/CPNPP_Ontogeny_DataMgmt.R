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
Trait_data [Trait_data$SAMPLE_ID == "VUOC_UTSE_IM_4_H2_B", "LWD"] <-0.0024
Trait_data [Trait_data$SAMPLE_ID == "ACMI_UTNW_VE_2_H3_G", "LWD"] <-0.00075
Trait_data [Trait_data$SAMPLE_ID == "ELTR_NMNW_C_18_H1_O", "LWD"] <-0.0021
Trait_data [Trait_data$SAMPLE_ID == "MUPO_AZSE_SIT_8_H3_O", "LWD"] <-0.0005
Trait_data [Trait_data$SAMPLE_ID == "PAMU_UTSW_CP_10_H3_G", "LWD"] <-0.001   
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


# ------------------------------------------ SRL -------------------------------------------------

Trait_data$SRL <- ((Trait_data$SumOfLength.cm./Trait_data$RWD)/1000) # SRL - Specific root length
Trait_data$SRL[Trait_data$SRL > 300] <- NA 

#------------------------------------------ RDMC & RTD  -------------------------------------------------

Trait_data$RDMC <- (Trait_data$RWD)/(Trait_data$RWF) # RDMC (Root Dry Matter Content)

Trait_data [Trait_data$SAMPLE_ID == "MUPO_AZSE_BRR_5_H1_G", "RDMC"] <-NA
Trait_data [Trait_data$SAMPLE_ID == "MUPO_AZSE_SIT_3_H1_O", "RDMC"] <-NA
Trait_data$RDMC[Trait_data$RDMC > 0.6] <- NA 

#---------------------------------------------RTD-------------------------------------------

Trait_data$RTD <- (Trait_data$RWD)/(Trait_data$SumOfRootVolume.cm3.) # RTD (Root tissue density - g / cm ^3) 
Trait_data$RTD[Trait_data$RTD > 0.6] <- NA 

# ----------------------------------- RMR ------------------------------------------------------------------

Trait_data$RMR <- (Trait_data$RWD)/(rowSums(Trait_data[,c( "CWD", "SWD","LWD", "RWD", "LWD_A", "LWD_B", "BWD_A")], na.rm = TRUE)) # RMR (Root Mass Ratio)
Trait_data$RMR[Trait_data$RMR> 0.9] <- NA 

Root_vars<-(Trait_data[names(Trait_data) %in% c("SAMPLE_ID","SumOfLength.cm.","RWF","RWD","SRL", "RDMC", "RTD", "RMR")])
# View(Root_vars) 

# ------------------------------ Calculate Population Avgs ------------------------

Trait_avg_pop_df<-lapply(c("SLA","LDMC", "RMR","SRL","RDMC","RTD"), function (y) 
  agg_mean_fun(Trait_data, y, "H_num","POP_ID", c("H_num","POP_ID","value","trait")))

Trait_avg_pop_df<-do.call(rbind, Trait_avg_pop_df)

Trait_avg_pop_df_3<-pop_id_function(Trait_avg_pop_df, "POP_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE"))


# ------------------------------ Calculate Species Avgs ------------------------

Trait_avg_sps_df<-lapply(c("SLA","LDMC", "RMR","SRL","RDMC","RTD"), function (y) 
  agg_mean_fun(Trait_data, y, "H_num","SPECIES", c("H_num","SPECIES","value","trait")))

Trait_avg_sps_df<-do.call(rbind, Trait_avg_sps_df)

Trait_avg_pop_df_3<-pop_id_function(Trait_avg_pop_df, "POP_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE"))


# ------------------------Calculate Relative Growth Rates ----------------------
# Tried to actually calculate based on days but messing up in the "function" section. 

#View((Trait_data[names(Trait_data) %in% c("SAMPLE_ID","TRANSP_DATE","TLD","HD","TRANSP_DATE_X", "TRANSP_DATE_Y", "TRANSP_DATE_Z", "HD_Z", "Days")]))
# Convert dates in 28-Jul format to 7/29 format 

# Trait_data$TRANSP_DATE_X <- as.Date(Trait_data$TRANSP_DATE,format="%d-%B")
# Trait_data$TRANSP_DATE_Y <- as.Date(Trait_data$TRANSP_DATE,format="%m/%d")
# Trait_data$TRANSP_DATE_Z<-paste(Trait_data$TRANSP_DATE_X, Trait_data$TRANSP_DATE_Y)
# Trait_data$TRANSP_DATE_Z<-gsub("NA", "", Trait_data$TRANSP_DATE_Z)
# Trait_data$TRANSP_DATE_Z<-as.Date(Trait_data$TRANSP_DATE_Z)
# 
# Trait_data$HD_X <- as.Date(Trait_data$HD,format="%d-%B")
# Trait_data$HD_Y <- as.Date(Trait_data$HD,format="%m/%d")
# Trait_data$HD_Z<-paste(Trait_data$HD_X, Trait_data$HD_Y)
# Trait_data$HD_Z<-gsub("NA", "", Trait_data$HD_Z)
# Trait_data$HD_Z<-as.Date(Trait_data$HD_Z)
# 
# Trait_data$TRANSP_DATE_Z<-gsub("2020", "2017",Trait_data$TRANSP_DATE_Z)
# Trait_data$HD_Z<-gsub("2020", "2017",Trait_data$HD_Z)
# 
# Trait_data [Trait_data$SAMPLE_ID == "ELTR_UTC_TMC_4_H4_B", "HD_Z"] <-as.character("2018-01-03")
# Trait_data [Trait_data$SAMPLE_ID == "ELTR_UTC_TMC_7_H4_B", "HD_Z"] <-as.character("2018-01-03")
# Trait_data [Trait_data$SAMPLE_ID == "HEAN_UTEC_CI_18_H4_B", "HD_Z"] <-as.character("2018-01-02")
# 
# Trait_data$Days<-round(as.numeric(difftime(Trait_data$HD_Z, Trait_data$TRANSP_DATE_Z, units = "days"), digits = 0))
# Trait_data<-Trait_data[!names(Trait_data) %in% c("TRANSP_DATE_X","TRANSP_DATE_Y", "HD_X", "HD_Y") ]
# 
# Traitdata_days_avg<-agg_mean_fun(Trait_data, "Days", "H_num", "POP_ID", c("H_num","POP_ID","value", "Days_avg"))

# Sum of aboveground for each sample 
Trait_data$Above_weight_sum <- rowSums(Trait_data[,c("CWD","LWD", "SWD","LWD_A","LWD_B", "BWD_A", "SWD_A")], na.rm = TRUE)
Trait_data$Tot_weight_sum <- rowSums(Trait_data[,c("Above_weight_sum","RWD")], na.rm = TRUE)

# Add a "Days" column to dataset to calculate growth rates 
Trait_data$Days <- as.factor(Trait_data$H_num)
levels(Trait_data$Days) <- c("10","24","42","84")

agg_ls_RGR<-lapply(c("Above_weight_sum","Tot_weight_sum", "RWD","SumOfLength.cm."), function (y) 
  agg_mean_fun(Trait_data, y, "Days","POP_ID", c("H_num","POP_ID","value","trait")))

# Merge into one dataset
agg_df<-do.call(cbind, agg_ls_RGR)
agg_df<-agg_df[ c(1,2,3,7,11,15) ] # column names are similar so indexing 
colnames(agg_df)<-c("H_num","POP_ID","Above_weight_avg","Tot_weight_avg","Below_weight_avg","RE_avg")
GrowthRates_Traits2017_3<-agg_df
GrowthRates_Traits2017_3$H_num<-as.numeric(as.character(GrowthRates_Traits2017_3$H_num))

Growth_calc<- by(GrowthRates_Traits2017_3, GrowthRates_Traits2017_3$POP_ID, function(df){
  df$RGR_Above_ln = growth_function(df$Above_weight_avg, df$H_num )
  df$RGR_Below_ln = growth_function(df$Below_weight_avg, df$H_num)
  df$RGR_Tot_ln = growth_function(df$Tot_weight_avg, df$H_num)
  df$RER_ln = growth_function(df$RE_avg, df$H_num)
  df
})

Growth_calc_df<-do.call(rbind, Growth_calc)

GrowthRates_Traits2017_3<-pop_id_function(Growth_calc_df, "POP_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE"))
GrowthRates_Traits2017_3$Days<-GrowthRates_Traits2017_3$H_num
levels(GrowthRates_Traits2017_3$H_num)[(GrowthRates_Traits2017_3$H_num) %in% c("10","24","42","84")]<-c("H1","H2","H3","H4")

# ------------------------Calculate Plasticity Index ----------------------

agg_ls_PI_traits<-lapply(c("SLA","LDMC", "RMR","SRL","RDMC","RTD"), function (y) 
  agg_mean_fun(Trait_data, y, "Days","SPECIES", c("Days","SPECIES","value","trait")))

agg_ls_PI_growthrates<-lapply(c("RER_ln","RGR_Tot_ln","RGR_Above_ln","RGR_Below_ln"), function (y) 
  agg_mean_fun(GrowthRates_Traits2017_3, y, "Days","SPECIES", c("Days","SPECIES","value","trait")))

## Full dataset 
Trait_avg<-do.call(rbind,agg_ls_PI_traits)
RGR_avg<-do.call(rbind, agg_ls_PI_growthrates)

All_Species_avg<-rbind(Trait_avg,RGR_avg)
Trait_splits_PI = split(All_Species_avg, paste(All_Species_avg$SPECIES, All_Species_avg$trait))

distances<-lapply(Trait_splits_PI, dist_fun)
distances_all<-do.call(rbind.data.frame, distances)
distances_all<-distances_all[complete.cases(distances_all), ]

## Without H1 datast for SLA and LDMC for FORB species (most didn't really have true leaves)
Trait_avg_10day_ForbsDrops<-Trait_avg[!((Trait_avg$Days) %in% c("10") & 
                               (Trait_avg$trait) %in% c("SLA","LDMC") &
                               (Trait_avg$SPECIES) %in% c("ACMI","ARTR","HEVI","MACA","PLPA","PAMU","HEAN")),]
Species_splits_PI_10day = split(Trait_avg_10day_ForbsDrops, paste(Trait_avg_10day_ForbsDrops$SPECIES, Trait_avg_10day_ForbsDrops$Trait))

distances_10day<-lapply(Species_splits_PI_10day, dist_fun)
distances_all_10day<-do.call(rbind.data.frame, distances_10day)


### ----------------------- Calculate relative change (base traits ) ----------------------------

agg_ls_relative_traits<-lapply(c("SLA","LDMC", "RMR","SRL","RDMC","RTD"), function (y) 
  agg_mean_fun(Trait_data, y, "Days","POP_ID", c("Days","POP_ID","value","trait")))

agg_relative_traits <-do.call(rbind, agg_ls_relative_traits)
agg_relative_traits[agg_relative_traits == 0] <- NA 

agg_relative_traits_wide<-spread(agg_relative_traits, trait, value)

RC_calc<- by(agg_relative_traits_wide, agg_relative_traits_wide$POP_ID, function(df){
  df$rel_LDMC = relative_change_function(df$LDMC)
  df$rel_LDMC2 = df$rel_LDMC/shift(df$LDMC, 1L, type = "lag")
  df$rel_RDMC = relative_change_function(df$RDMC)
  df$rel_RDMC2 = df$rel_RDMC/shift(df$RDMC, 1L, type = "lag")
  df$rel_RMR = relative_change_function(df$RMR)
  df$rel_RMR2 = df$rel_RMR/shift(df$RMR, 1L, type = "lag")
  df$rel_RTD = relative_change_function(df$RTD)
  df$rel_RTD2 = df$rel_RTD/shift(df$RTD, 1L, type = "lag")
  df$rel_SLA = relative_change_function(df$SLA)
  df$rel_SLA2 = df$rel_SLA/shift(df$SLA, 1L, type = "lag")
  df$rel_SRL = relative_change_function(df$SRL)
  df$rel_SRL2 = df$rel_SRL/shift(df$SRL, 1L, type = "lag")
  df
}) 

Rel_cal_df<-do.call(rbind, RC_calc)
Rel_cal_df<-Rel_cal_df[,!names(Rel_cal_df) %in% c("rel_LDMC", "rel_RDMC","rel_RMR", "rel_RTD","rel_SLA", "rel_SRL")]
Rel_cal_df$H_num<-as.factor(Rel_cal_df$Days)
Rel_cal_df_3<-pop_id_function(Rel_cal_df, "POP_ID", "_", c("SPECIES","LOCATION_CODE","POP_CODE"))

levels(Rel_cal_df_3$H_num)[(Rel_cal_df_3$H_num) %in% c("10","24","42","84")]<-c("H1","H2","H3","H4")

### Write complete dataset to CSV 
write.csv(Trait_data, file = "Data_Generated/TraitData_2017.csv")
write.csv(GrowthRates_Traits2017_3, "Data_Generated/TraitData_GrowthRates_2017.csv")
write.csv(distances_all, "Data_Generated/TraitData_Plasticity_2017.csv")
write.csv(distances_all_10day, "Data_Generated/TraitData_Plasticity_2017_10day.csv")
write.csv(Trait_avg_pop_df_3, "Data_Generated/TraitData_SpeciesAvg_2017.csv")
write.csv(Trait_avg_pop_df_3, "Data_Generated/TraitData_PopAvg_2017.csv")
write.csv(Rel_cal_df_3, "Data_Generated/TraitData_RelativeDiff_2017.csv")




