
## PREPROCESS LONGITUDINAL DATA ##

# Clear R environment
rm(list=ls())


#### 1. CLEAN VARIABLES IN RAW LONGITUDINAL DATASET ####

data <- readRDS("data/processed/raw_longitudinal_data.rds")

# Economic activity
data <- data %>% 
  mutate(jb = case_when(jbstat %in% c(1,2) ~ 1,    # employed
                        jbstat %in% c(3,13) ~ 2,   # unemployed
                        jbstat == 8 ~ 3,           # long-term sickness
                        jbstat == 4 ~ 4,           # retired
                        jbstat == 6 ~ 5,           # care
                        jbstat %in% c(7,11,5,14,15,12,10,9,97) ~ 6, # other
                        TRUE ~ NA))  # missing
data$jb <- factor(data$jb,
                  levels=c(1,2,3,4,5,6),
                  labels=c("employed", "unemployed", "long-term sick", "retired", "care", "other"))

# Self-rated health 
data <- data %>% 
  mutate(srh = ifelse(!is.na(scsf1) & scsf1>0, scsf1, #  Use self-completed self-rated health (scsf1) where available (i.e., Wave 2 onwards for all non-proxy respondents)
                      ifelse(sf1>0, sf1, NA))) # If not, use face-to-face self-rated health if given (i.e., all in Wave 1 and proxy respondents)                      
data$srh <- factor(data$srh, 
                   levels=c(1,2,3,4,5), 
                   labels=c("excellent", "very good", "good", "fair", "poor"))

# Sex
data <- data %>% 
  mutate(sex = case_when(sex_dv == 1 ~ "male", 
                         sex_dv == 2 ~ "female", 
                         TRUE ~ NA)) # missing

# Age
data <- data %>% 
  mutate(age = case_when(age_dv >= 0 ~ age_dv, 
                         TRUE ~ NA))

# Age group - 10 year for working-age population (16-64)
data <- data %>% 
  mutate(agegrp_10yr = case_when(age %in% c(16:24) ~ "16-24", 
                                 age %in% c(25:34) ~ "25-34",
                                 age %in% c(35:44) ~ "35-44",
                                 age %in% c(45:54) ~ "45-54",
                                 age %in% c(55:64) ~ "55-64", 
                                 TRUE ~ NA))
data %>% group_by(agegrp_10yr) %>% summarise(min = min(age), max = max(age))

# Highest educational qualification (3 category)
data <- data %>% 
  mutate(hiqual_3cat = case_when(hiqual_dv %in% c(1,2) ~ "high", 
                                 hiqual_dv %in% c(3,4) ~ "medium", 
                                 hiqual_dv %in% c(5,9) ~ "low",
                                 TRUE ~ NA))

# Highest educational qualification (4 category)
data <- data %>% 
  mutate(hiqual_4cat = case_when(hiqual_dv %in% c(1,2) ~ "high", 
                                 hiqual_dv %in% c(3,4) ~ "medium", 
                                 hiqual_dv %in% c(5) ~ "other",
                                 hiqual_dv %in% c(9) ~ "no",
                                 TRUE ~ NA))

# Household type (4 category)
data <- data %>% 
  mutate(hhtype_4cat = case_when(hhtype_dv %in% c(1,2,3,16,17,22) ~ 1, 
                                 hhtype_dv %in% c(4,5,18,23) ~ 2, 
                                 hhtype_dv %in% c(6,8,19) ~ 3, 
                                 hhtype_dv %in% c(10,11,12,20,21) ~ 4))
data$hhtype_4cat <- factor(data$hhtype_4cat, 
                           levels=c(1,2,3,4), 
                           labels=c("Single without children", "Single with children", "Couple without children", "Couple with children"))

# Marital/Partnership status (3 category)
data <- data %>% 
  mutate(marstat_3cat = case_when(marstat == 1 ~ 1, 
                                  marstat %in% c(2,3) ~ 2, 
                                  marstat %in% c(4:9) ~ 3, 
                                  TRUE ~ NA))
data$marstat_3cat <- factor(data$marstat_3cat, 
                            levels=c(1,2,3), 
                            labels=c("Single", "Partnered", "Separated"))

# Marital/Partnership status (2 category)
data <- data %>% 
  mutate(marstat_2cat = case_when(marstat %in% c(1,4:9) ~ 1, 
                                  marstat %in% c(2,3) ~ 2, 
                                  TRUE ~ NA))
data$marstat_2cat <- factor(data$marstat_2cat, 
                            levels=c(1,2), 
                            labels=c("Single", "Partnered"))

# Government Office Region (GOR)
data <- data %>% 
  mutate(gor = ifelse(gor_dv < 0, NA, gor_dv))
data$gor <- factor(data$gor, 
                   level=c(1,2,3,4,5,6,7,8,9,10,11,12), 
                   labels=c("North East", "North West", "Yorkshire and the Humber", "East Midlands", 
                            "West Midlands", "East of England", "London", "South East", "South West", 
                            "Wales", "Scotland", "Northern Ireland"))

# GHQ-12 Likert Scale (0-36)
data <- data %>% 
  mutate(ghqscale = ifelse(scghq1_dv < 0, NA, scghq1_dv))

# GHQ-12 Caseness (4+ on 12-point scale)
data <- data %>% 
  mutate(scghq2 = ifelse(scghq2_dv < 0, NA, scghq2_dv), 
         ghqcase = ifelse(scghq2 >= 4, 1, 0))
# Make GHQ-12 caseness a factor variable 
data$ghqcase <- factor(data$ghqcase, 
                       levels = c(0,1), 
                       labels = c("No", "Yes"))

# SF-12 MCS
data <- data %>% 
  mutate(sf12mcs = ifelse(sf12mcs_dv >= 0, sf12mcs_dv, NA))

# SF-12 PCS
data <- data %>% 
  mutate(sf12pcs = ifelse(sf12pcs_dv >= 0, sf12pcs_dv, NA))


#### 2. IDENTIFY TRANSITIONS IN ECONOMIC ACTIVITY STATUS ####

# Create a 2-digit summary variable of transitions in economic activity status between waves 
# This code is heavily based on the UKHLS Introduction to Understanding Society R training course materials

data <- data %>%
  arrange(pidp, wave_number) %>% # sorting by pidp and wave to allow for the lag
  mutate(jb_num = as.numeric(jb)) %>% # creating a numerical version of the jb (makes it 1,2,3,4,5,6)
  group_by(pidp) %>% # grouping by the individual identifier pidp
  mutate(jb_lag = lag(jb_num, n = 1, default = NA)) %>% # creating a one period lagged variable by one row (n=1)
  mutate(jbch = factor(jb_lag * 10 + as.numeric(jb_num), # creating change variable
                       level=c(11,12,13,14,15,16,
                               21,22,23,24,25,26,
                               31,32,33,34,35,36,
                               41,42,43,44,45,46,
                               51,52,53,54,55,56,
                               61,62,63,64,65,66),
                       label=c("Remains employed","Employed to unemployed","Employed to long-term sick","Employed to retired","Employed to care","Employed to other",
                               "Unemployed to employed","Remains unemployed","Unemployed to long-term sick","Unemployed to retired","Unemployed to care","Unemployed to other",
                               "Long-term sick to employed","Long-term sick to unemployed","Remains long-term sick","Long-term sick to retired","Long-term sick to care","Long-term sick to other",
                               "Retired to employed","Retired to unemployed","Retired to long-term sick","Remains retired","Retired to care","Retired to other",
                               "Care to employed","Care to unemployed","Care to long-term sick","Care to retired","Remains care","Care to other",
                               "Other to employed","Other to unemployed","Other to long-term sick","Other to retired","Other to care","Remains other"))) %>%
  ungroup() %>% select(-jb_num, -jb_lag) # remove the lag variables


# Identify those with single year difference in wave numbers as this was not specified in the lagging process
data <- data %>%
  arrange(pidp, wave_number) %>% # sorting by pidp and wave to allow for the lag 
  group_by(pidp) %>% # grouping by the individual identifier pidp
  mutate(wave_number_lag = lag(wave_number, n = 1, default = NA)) %>% # creating a one period lagged variable by one row (n=1)
  mutate(wave_diff = wave_number - wave_number_lag) %>%   # creating wave difference variable
  ungroup() %>% select(-wave_number_lag)

# Set as NA for all non-consecutive wave transitions to remove non-one-year lags
data$jbch[data$wave_diff > 1] <- NA

# Reduce data file to only processed variables (i.e. those with NAs)
longitudinal_transitions_data <- data %>% select(pidp, wave_number,  # IDs
                                                 jb, jbch, # Economic activity and transitions in economic activity
                                                 srh, # Self-rated health
                                                 age, agegrp_10yr, sex, gor, # Demographics
                                                 hiqual_3cat, hiqual_4cat, hhtype_4cat, marstat_2cat, marstat_3cat, equivquint, # Socioeconomics
                                                 ghqcase, ghqscale, # GHQ-12
                                                 sf12mcs, sf12pcs, # SF-12
                                                 intdaty_dv, intdatm_dv, intdatd_dv) # Interview date

# Save as the processed longitudinal data
saveRDS(longitudinal_transitions_data, file = "data/processed/longitudinal_transitions_data.rds")

                             