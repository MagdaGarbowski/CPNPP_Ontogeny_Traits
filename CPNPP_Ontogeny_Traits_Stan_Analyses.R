setwd("/Users/MagdaGarbowski/CPNPP_Ontogeny_Traits/")
SpeciesData<-read.csv("Data_Generated/TraitData_2017.csv")
source("Rscripts/Functions/Functions_Stan_Analyses.R")

SpeciesData$RMR_100<-SpeciesData$RMR*100
SpeciesData$RDMC_100<-SpeciesData$RDMC*100
SpeciesData$LDMC_100<-SpeciesData$LDMC*100
SpeciesData$RTD_100<-SpeciesData$RTD*100
SpeciesData$H_num<-as.character(as.factor(SpeciesData$H_num))

# Add NA for H1 SLA for forb species 
SpeciesData$SLA<-ifelse((SpeciesData$GrowthForm %in% c("FORB","SHRUB") & SpeciesData$H_num == "H1"), NA, SpeciesData$SLA)
SpeciesData$SLA<-ifelse((SpeciesData$SLA == 0), NA, SpeciesData$SLA)

Species_splits = split(SpeciesData, paste(SpeciesData$SPECIES))

dat_na.omit_SLA<-lapply(Species_splits, na.omit.fun, c("SAMPLE_ID","POP_ID","H_num","SLA"))
mk_data_SLA<-lapply(dat_na.omit_SLA, mk_data_function, "SLA")

dat_na.omit_RMR<-lapply(Species_splits, na.omit.fun, c("SAMPLE_ID","POP_ID","H_num","RMR"))
mk_data_RMR<-lapply(dat_na.omit_RMR, mk_data_function, "RMR")

# ----------------------------- Run SLA models one by one to see if they are working. Divergent transitions for most -------------------------------------------------------

fit_SLA_ACMI <- stan(file = "stan_models/example.stan", data = mk_data_SLA$ACMI, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) #  Divergent transitions
summary(fit_SLA_ACMI) # MATT - why is it printing summaries for all chains? 
aggregate(.~H_num, dat_na.omit_SLA$ACMI, FUN = mean)

fit_SLA_ARTR <- stan(file = "stan_models/example.stan", data = mk_data_SLA$ARTR, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) #  Divergent transtitions
summary(fit_SLA_ARTR) 
aggregate(.~H_num, dat_na.omit_SLA$ARTR, FUN = mean)

fit_SLA_ELTR <- stan(file = "stan_models/example.stan", data = mk_data_SLA$ELTR, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # OK 
summary(fit_SLA_ELTR) 
aggregate(.~H_num, dat_na.omit_SLA$ELTR, FUN = mean)

fit_SLA_HECO <- stan(file = "stan_models/example.stan", data = mk_data_SLA$HECO, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) #  Divergent transtitions
summary(fit_SLA_HECO) 
aggregate(.~H_num, dat_na.omit_SLA$HECO, FUN = mean)

fit_SLA_HEAN <- stan(file = "stan_models/example.stan", data = mk_data_SLA$HEAN, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) #  Divergent transtitions
summary(fit_SLA_HEAN) 
aggregate(.~H_num, dat_na.omit_SLA$HEAN, FUN = mean)

fit_SLA_HEVI <- stan(file = "stan_models/example.stan", data = mk_data_SLA$HEVI, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # Divergent transtitions
summary(fit_SLA_HEVI) 
aggregate(.~H_num, dat_na.omit_SLA$HEVI, FUN = mean)

fit_SLA_MUPO <- stan(file = "stan_models/example.stan", data = mk_data_SLA$MUPO, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # Divergent transtitions
summary(fit_SLA_MUPO) 
aggregate(.~H_num, dat_na.omit_SLA$MUPO, FUN = mean)

fit_SLA_MACA <- stan(file = "stan_models/example.stan", data = mk_data_SLA$MACA, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # Divergent transitions
summary(fit_SLA_MACA) 
aggregate(.~H_num, dat_na.omit_SLA$MACA, FUN = mean)

fit_SLA_PAMU <- stan(file = "stan_models/example.stan", data = mk_data_SLA$PAMU, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # Divergent transitions
summary(fit_SLA_PAMU) 
aggregate(.~H_num, dat_na.omit_SLA$PAMU, FUN = mean)

fit_SLA_PLPA <- stan(file = "stan_models/example.stan", data = mk_data_SLA$PLPA, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # Divergent transitions
summary(fit_SLA_PLPA) 
aggregate(.~H_num, dat_na.omit_SLA$PLPA, FUN = mean)

fit_SLA_VUOC <- stan(file = "stan_models/example.stan", data = mk_data_SLA$VUOC, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) # Divergent transitions
summary(fit_SLA_VUOC) 
aggregate(.~H_num, dat_na.omit_SLA$VUOC, FUN = mean)


# ----------------------------- Run RMR models one by one to see if they are working. Divergent transitions for most -------------------------------------------------------

fit_RMR_ACMI <- stan(file = "stan_models/example.stan", data = mk_data_RMR$ACMI, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99))  
summary(fit_RMR_ACMI) 
aggregate(.~H_num, dat_na.omit_RMR$ACMI, FUN = mean)

fit_RMR_ARTR <- stan(file = "stan_models/example.stan", data = mk_data_RMR$ARTR, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_ARTR) 
aggregate(.~H_num, dat_na.omit_RMR$ARTR, FUN = mean)

fit_RMR_ELTR <- stan(file = "stan_models/example.stan", data = mk_data_RMR$ELTR, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99))  
summary(fit_RMR_ELTR) 
aggregate(.~H_num, dat_na.omit_RMR$ELTR, FUN = mean)

fit_RMR_HECO <- stan(file = "stan_models/example.stan", data = mk_data_RMR$HECO, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_HECO) 
aggregate(.~H_num, dat_na.omit_RMR$HECO, FUN = mean)

fit_RMR_HEAN <- stan(file = "stan_models/example.stan", data = mk_data_RMR$HEAN, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_HEAN) 
aggregate(.~H_num, dat_na.omit_RMR$HEAN, FUN = mean)

fit_RMR_HEVI <- stan(file = "stan_models/example.stan", data = mk_data_RMR$HEVI, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_HEVI) 
aggregate(.~H_num, dat_na.omit_RMR$HEVI, FUN = mean)

fit_RMR_MUPO <- stan(file = "stan_models/example.stan", data = mk_data_RMR$MUPO, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_MUPO) 
aggregate(.~H_num, dat_na.omit_RMR$MUPO, FUN = mean)

fit_RMR_MACA <- stan(file = "stan_models/example.stan", data = mk_data_RMR$MACA, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_MACA) 
aggregate(.~H_num, dat_na.omit_RMR$MACA, FUN = mean)

fit_RMR_PAMU <- stan(file = "stan_models/example.stan", data = mk_data_RMR$PAMU, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_PAMU) 
aggregate(.~H_num, dat_na.omit_RMR$PAMU, FUN = mean)

fit_RMR_PLPA <- stan(file = "stan_models/example.stan", data = mk_data_RMR$PLPA, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_PLPA) 
aggregate(.~H_num, dat_na.omit_RMR$PLPA, FUN = mean)

fit_RMR_VUOC <- stan(file = "stan_models/example.stan", data = mk_data_RMR$VUOC, warmup = 500, iter = 3000, chains = 4, cores = 2, thin = 1, control = list(adapt_delta = 0.99)) 
summary(fit_RMR_VUOC) 
aggregate(.~H_num, dat_na.omit_RMR$VUOC, FUN = mean)