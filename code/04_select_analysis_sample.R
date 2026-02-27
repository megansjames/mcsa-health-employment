
## SELECT ANALYSIS SAMPLE ##

# Clear R environment
rm(list=ls()) 

# Read in the processed longitudinal transitions data
data <- readRDS("data/processed/longitudinal_transitions_data.rds")


#### 1. SAMPLE SELECTION PROCESS ####

# Number of individual respondents of all ages in UKHLS
length(unique(data$pidp))

# Filter to those aged 16 to 64 years old (the working-age population)
data_16_64 <- filter(data, age > 15 & age < 65)
length(unique(data_16_64$pidp))

# Filter to those who experience transitions to long-term sickness in the working-age population
transitions_to_lts_names <- c("Employed to long-term sick", 
                              "Unemployed to long-term sick", 
                              "Retired to long-term sick", 
                              "Care to long-term sick", 
                              "Other to long-term sick")

lts_transitions <- filter(data_16_64, jbch %in% transitions_to_lts_names)
length(unique(lts_transitions$pidp))

# Filter to those transitions to long-term sickness prior to Wave 11 so that 5 years of follow-up are possible (regardless if they are observed)
lts_transitions_5year <- filter(lts_transitions, wave_number < 11)
length(unique(lts_transitions_5year$pidp))

# Label the first wave of the transition to long-term sickness (regardless of other transitions)
lts_transitions_first <- lts_transitions_5year %>% 
  group_by(pidp) %>% 
  summarise(lts_first_wave = first(wave_number)) %>% 
  ungroup()

# Get all the observed data for all waves for these individuals who experience a transition for whom 5 years of follow-up are possible
# Create variable 'transition_time' which is the time since the first wave with the transition
# Filter to the wave prior to transition and all waves following the transition to long-term sickness (transition_time >= -1)

lts_transitions_1yr_prior_5yr_followup_possible <- data %>%
  filter(pidp %in% unique(lts_transitions_5year$pidp)) %>%
  left_join(lts_transitions_first) %>%
  mutate(transition_time = wave_number-lts_first_wave) %>% 
  filter(transition_time >= -1) 

length(unique(lts_transitions_1yr_prior_5yr_followup_possible$pidp))
  
# Save for missing data analysis
saveRDS(lts_transitions_1yr_prior_5yr_followup_possible, 
        file = "data/processed/lts_transitions_1yr_prior_5yr_followup_possible.rds")

# Tidy
rm(data, data_16_64, lts_transitions_first, lts_transitions_5year, lts_transitions)

# Remove those with missing covariates

# Baseline/year of onset of long-term sickness 
lts_transitions_baseline <- filter(lts_transitions_1yr_prior_5yr_followup_possible, transition_time == 0)
lts_transitions_baseline_noNAs <- lts_transitions_baseline %>% 
  select(pidp, jb, srh, age, sex, hiqual_4cat, marstat_2cat, gor, ghqcase) %>% # Variables of concern for prior prediction
  tidyr::drop_na()

# one year prior to long-term sickness 
lts_transitions_1yr_prior <- filter(lts_transitions_1yr_prior_5yr_followup_possible, transition_time == -1) 
lts_transitions_1yr_prior_noNAs <- lts_transitions_1yr_prior %>% 
  select(pidp, jb, srh, ghqcase, equivquint) %>% # Variables of concern for prior prediction
  tidyr::drop_na()

# Filter out those with missing data for baseline or prediction covariates
lts_transitions_1yr_prior_noNAs_5yr_followup_possible <- lts_transitions_1yr_prior_5yr_followup_possible %>% 
  filter(pidp %in% unique(lts_transitions_baseline_noNAs$pidp)) %>% 
  filter(pidp %in% unique(lts_transitions_1yr_prior_noNAs$pidp))

# Number of individuals after removing NAs from covariates
length(unique(lts_transitions_1yr_prior_noNAs_5yr_followup_possible$pidp))

# Save dataset
saveRDS(lts_transitions_1yr_prior_noNAs_5yr_followup_possible, 
        file = "data/processed/lts_transitions_1yr_prior_noNAs_5yr_followup_possible.rds")

# Filter to those with 5 years of complete consecutive follow-up
lts_transitions_5yr_max_followup <- filter(lts_transitions_1yr_prior_noNAs_5yr_followup_possible, transition_time < 5)
tally <- lts_transitions_5yr_max_followup %>% group_by(pidp) %>% tally()
table(tally$n)
lts_transitions_1yr_prior_5yr_followup_complete_incNAs <- left_join(lts_transitions_5yr_max_followup, tally) %>% 
  filter(n == 6) %>% select(-n) # includes the wave prior to transition

# Number of individuals with 5 years of consequent follow-up
length(unique(lts_transitions_1yr_prior_5yr_followup_complete_incNAs$pidp))

# Save for missing data analysis due to missing states in sequences
saveRDS(lts_transitions_1yr_prior_5yr_followup_complete_incNAs, 
        file = "data/processed/lts_transitions_1yr_prior_5yr_followup_complete_incNAs.rds")

# Tidy
rm(tally, lts_transitions_baseline, lts_transitions_baseline_noNAs, 
   lts_transitions_1yr_prior, lts_transitions_1yr_prior_noNAs, 
   lts_transitions_5yr_max_followup, lts_transitions_1yr_prior_5yr_followup_possible, 
   lts_transitions_1yr_prior_noNAs_5yr_followup_possible)

# Mark out who has missing economic activity data for the sequence follow-up
lts_transitions_jb_NA <- lts_transitions_1yr_prior_5yr_followup_complete_incNAs %>% 
  filter(transition_time >= 0) %>% 
  group_by(pidp) %>% 
  summarise(jb_NA = sum(as.numeric(jb))) %>% # convert to numeric and summarise as sum will be NA if one value is NA
  mutate(jb_NA = ifelse(is.na(jb_NA), 1, 0)) # indicator on whether any value of JB is missing

# Table of how many have missing economic activity data
table(lts_transitions_jb_NA$jb_NA)

# Mark out who has missing self-rated health data for the sequence follow-up
lts_transitions_srh_NA <- lts_transitions_1yr_prior_5yr_followup_complete_incNAs %>% 
  filter(transition_time >= 0) %>% 
  group_by(pidp) %>% 
  summarise(srh_NA = sum(as.numeric(srh))) %>% # convert to numeric and summarise as sum will be NA if one value is NA
  mutate(srh_NA = ifelse(is.na(srh_NA), 1, 0)) # indicator on whether any value of SRH is missing

# Table of how many have missing self-rated health data
table(lts_transitions_srh_NA$srh_NA)

# Join missing indicators to the main data to allow for a filter to remove those with missing states
lts_transitions_1yr_prior_5yr_followup_complete_noNAs <- lts_transitions_1yr_prior_5yr_followup_complete_incNAs %>% 
  left_join(lts_transitions_jb_NA) %>% filter(jb_NA == 0) %>% 
  left_join(lts_transitions_srh_NA) %>% filter(srh_NA == 0) %>% 
  select(-jb_NA, -srh_NA)

# Number of individuals after removing the NAs
length(unique(lts_transitions_1yr_prior_5yr_followup_complete_noNAs$pidp))

# Save dataset
saveRDS(lts_transitions_1yr_prior_5yr_followup_complete_noNAs, 
        file = "data/processed/lts_transitions_1yr_prior_5yr_followup_complete_noNAs.rds")

# Tidy
rm(lts_transitions_jb_NA, lts_transitions_srh_NA, lts_transitions_1yr_prior_5yr_followup_complete_incNAs, lts_transitions_1yr_prior_5yr_followup_complete_noNAs)


#### 2. DESCRIPTIVE ANALYSIS OF SAMPLE ####

## DESCRIPTION OF STATES FOR SEQUENCE ANALYSIS ##

sample <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_complete_noNAs.rds")

# Economic activity states in the sample
table(sample$jb)
sum(table(sample$jb))
round(prop.table(table(sample$jb))*100,1)
sum(round(prop.table(table(sample$jb))*100,1))

# Self-rated health states in the sample
table(sample$srh)
sum(table(sample$srh))
round(prop.table(table(sample$srh))*100,1)
sum(round(prop.table(table(sample$srh))*100,1))

## DESCRIPTION OF COVARIATES ## 

# Baseline = initial transition into long-term sickness 
sample_at_baseline <- sample %>% filter(transition_time == 0)
# One-year prior = one-year prior to transition to long-term sickness
sample_1yr_prior <- sample %>% filter(transition_time == -1) %>%
  rename_with(~ paste0(.x, "_prior"), .cols = !pidp)
sample_for_table <- left_join(sample_at_baseline, sample_1yr_prior)
rm(sample_at_baseline, sample_1yr_prior)

##### Table: Analytical sample by covariate #####

table <- sample_for_table %>%
  gtsummary::tbl_summary(
    include = c(sex, agegrp_10yr, marstat_2cat, hiqual_4cat, gor, srh, ghqcase, 
                jb_prior, srh_prior, ghqcase_prior, equivquint_prior), 
    statistic = list(all_categorical() ~ c("{n} ({p}%)") , 
                     missing = "no"), 
    digits = list(all_categorical() ~ c(0,1)), 
    label = list(sex = "Sex", 
                 agegrp_10yr = "Age group", 
                 marstat_2cat = "Partnership status", 
                 hiqual_4cat = "Highest educational qualification", 
                 gor = "Region", 
                 srh = "Self-rated health", 
                 ghqcase = "GHQ-12 caseness", 
                 jb_prior = "Prior economic activity", 
                 srh_prior = "Prior self-rated health", 
                 ghqcase_prior = "Prior GHQ-12 caseness", 
                 equivquint_prior = "Prior equvalised household income quintile"
                 )) %>% 
  modify_header(label = "**Variable**") %>%
  bold_labels() %>% 
  as_gt() %>% gt::tab_source_note(gt::md("*Percentages may not sum to 100% due to rounding*"))

# Print as HTML in the Viewer
table

# Export to a Word document for future use
gt::gtsave(table, filename = "output/main/table_overall_sample.docx")


## TIME BETWEEN INTERVIEW DATES ## 
sample <- sample %>% 
  mutate(interview_date = lubridate::make_date(year = intdaty_dv, 
                                               month = intdatm_dv, 
                                               day = intdatd_dv)) %>% 
  arrange(pidp, transition_time) %>% 
  group_by(pidp) %>% 
  mutate(gap_days = as.numeric(interview_date - lag(interview_date))) %>% 
  ungroup()

# Plot of the variation in the time between waves
ggplot(data = sample, aes(x = gap_days)) + 
  geom_histogram(color = "darkgrey", fill="white", binwidth = 1) +
  theme_minimal() +  
  geom_vline(aes(xintercept = 365), color = "darkblue", linewidth=1, linetype = "dashed") + 
  labs(x="Number of days between wave interviews", 
       y="Count") + 
  scale_x_continuous(breaks=c(seq(150,550,50)))