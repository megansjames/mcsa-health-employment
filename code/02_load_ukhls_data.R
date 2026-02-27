
## LOAD UKHLS DATA ##

# Data: Understanding Society, the UK Household Longitudinal Study (UKHLS)

# Citation: 
# University of Essex, Institute for Social and Economic Research. (2025). 
# Understanding Society: Waves 1-14, 2009-2023 and Harmonised BHPS: Waves 1-18, 1991-2009. [data collection]. 19th Edition. UK Data Service. 
# SN: 6614, DOI: http://doi.org/10.5255/UKDA-SN-6614-20

# Data is accessible from the UK Data Service (UKDS) and downloaded for use in Stata format (.dta)

# Clear R environment
rm(list=ls())


#### 1. READ IN ALL WAVES OF DATA ####

# 14 waves (a-n)
waves <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n") 
length(waves)

# Read in the individual respondent data files (_indresp.dta)

for (i in 1:length(waves)){
  
  # Wave letter from wave number
  wave <- waves[i]
  
  # Load data
  indresp <- haven::read_dta(file=paste0("data/raw/", wave, "_indresp.dta"))
  indresp <- select(indresp, 
                    pidp, contains("hidp"), # identifiers
                    ends_with(c("_jbstat", "_jbft_dv",  # economic activity
                                "_sf1", "_scsf1",  # self-rated health
                                "_age_dv", "_sex_dv", # demographics
                                "_hiqual_dv", "_hhtype_dv", "_marstat", "_gor_dv", # socioeconomic information
                                "_scghq1_dv", "_scghq2_dv", # GHQ-12 
                                "_sf12pcs_dv", "_sf12mcs_dv", # SF-12 summary scores
                                "_intdaty_dv", "_intdatm_dv", "_intdatd_dv" # interview date
                                ))) 
  
  # Remove wave prefixes
  indresp <- indresp %>%
    rename_with(~stringr::str_replace(.x, paste0(waves[i], "_"), ""),
                starts_with(paste0(waves[i], "_")))
  
  # Remove labels
  indresp <- indresp %>% mutate(across(everything(), labelled::remove_labels))
  
  # Add in wave information 
  indresp$wave_letter <- waves[i]
  indresp$wave_number <- as.numeric(paste(i)) # need numbers for the lags!
  
  # Export to the global environment 
  assign(paste0(wave, "_indresp"), indresp)
  rm(indresp)
}

# Read in the household respondent data files (_hhresp.dta)

for (i in 1:length(waves)){
  
  # Wave letter from wave number
  wave <- waves[i]
  
  # Load data
  hhresp <- haven::read_dta(file=paste0("data/raw/", wave, "_hhresp.dta"))
  hhresp <- select(hhresp,
                    contains(c("hidp", # identifiers
                               "fihhmngrs_dv", "ieqmoecd_dv", # household income and equivalisation
                               "hhdenu"))) # household cross-sectional weight

  # Remove wave prefixes
  hhresp <- hhresp %>%
     rename_with(~stringr::str_replace(.x, paste0(wave, "_"), ""),
                 starts_with(paste0(wave, "_")))
  
   # Remove labels
  hhresp <- hhresp %>% mutate(across(everything(), labelled::remove_labels))
  
   # Add in wave information
   hhresp$wave_letter <- waves[i]
   hhresp$wave_number <- as.numeric(paste(i))
  
   # Equivalise gross monthly household income
   hhresp$equivincome <- hhresp$fihhmngrs_dv / hhresp$ieqmoecd_dv
  
  # # Remove negative values
   hhresp$equivincome[hhresp$equivincome < 0] <- NA
  
  # # Split into weighted quintiles of equivalised gross monthly household income
   hhresp$equivquint <- cut(
    hhresp$equivincome,
     breaks = Hmisc::wtd.quantile(hhresp$equivincome, probs=seq(0,1,0.2), na.rm=TRUE),
     labels = c("Q1", "Q2", "Q3", "Q4", "Q5"),
     include.lowest = TRUE)

  # Export to the global environment
  assign(paste0(wave, "_hhresp"), hhresp)
  
  rm(hhresp)
}


#### 2. JOIN ALL WAVES TO CREATE A LONGITUDINAL DATASET ####

# For individual respondents
for (i in 1:length(waves)){
  
  indresp <- get(paste0(waves[i], "_indresp"))
  
  # Create a base dataframe and add on additional follow-up waves
  if(i == 1){ 
    ind_data <- indresp
  } else {
    ind_data <- bind_rows(ind_data, indresp) 
    }
  
  rm(indresp)
}

# For households
for (i in 1:length(waves)){ 
  
  hhresp <- get(paste0(waves[i], "_hhresp"))
  
  # Create a base dataframe and add on additional follow-up waves
  if(i == 1){ 
    hh_data <- hhresp
  } else {
    hh_data <- bind_rows(hh_data, hhresp) }

  rm(hhresp)
}

# Join individual respondents with relevant household information
data_combined <- left_join(ind_data, hh_data) # Join by hidp, wave_letter, wave_number

# Remove original data files 
rm(list = ls()[grepl("indresp", ls())])
rm(list = ls()[grepl("hhresp", ls())])
rm(ind_data, hh_data)

# Save the longitudinal dataset
saveRDS(data_combined, file = "data/processed/raw_longitudinal_data.rds")

