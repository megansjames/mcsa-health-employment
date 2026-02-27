
## RUN MISSING DATA ANALYSIS ##

rm(list=ls())


#### 1. SUMMARY OF MISSING DATA ####

# Those with missing covariate information prior to possible 5 years of follow-up
sample_1yr_NAs_5yr_possible <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_possible.rds") %>% 
  mutate(missing_1yr = "1yr_NAs")
length(unique(sample_1yr_NAs_5yr_possible$pidp))

# Those with possible 5 year follow-up following transition to long-term sickness (but no missing covariates)
sample_5yr_possible <- readRDS("data/processed/lts_transitions_1yr_prior_noNAs_5yr_followup_possible.rds") %>% 
  mutate(missing_5yr = "5yr_possible")
length(unique(sample_5yr_possible$pidp))

# Those with 5 year follow-up occurring (but with missing information for sequences)
sample_5yr_complete_incNAs <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_complete_incNAs.rds") %>% 
  mutate(missing_NAs = "5yr_incNAs")
length(unique(sample_5yr_complete_incNAs$pidp))

# Those with 5 year follow-up with no missing information for sequences (i.e., the final sample)
sample_5yr_complete_noNAs <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_complete_noNAs.rds") %>% 
  mutate(missing_final = "final")
length(unique(sample_5yr_complete_noNAs$pidp))

# Combine all together and create indicator of when they exit from the sample
missing_data <- list(sample_1yr_NAs_5yr_possible, sample_5yr_possible, sample_5yr_complete_incNAs, sample_5yr_complete_noNAs) %>% reduce(left_join)
missing_indicator <- missing_data %>% 
  mutate(missing = case_when(missing_final == "final" ~ "final", 
                             is.na(missing_final) & missing_NAs == "5yr_incNAs" ~ "5yr_incNAs", 
                             is.na(missing_final) & is.na(missing_NAs) & missing_5yr == "5yr_possible" ~ "5yr_possible", 
                             TRUE ~ "1yr_NAs")) %>% 
  select(-missing_final, -missing_NAs, -missing_5yr, -missing_1yr) 
         
# Separate into baseline and one-year prior and rejoin together 
missing_baseline <- filter(missing_indicator, transition_time == 0)
missing_prior <- filter(missing_indicator, transition_time == -1) %>% 
  rename_with(~ paste0("prior_", .x), .cols = !pidp)
missing_combined <- left_join(missing_baseline, missing_prior)
table(missing_combined$missing)

##### Appendix Table: Sample by reason for missingness #####

table <- missing_combined %>%
  gtsummary::tbl_summary(
    by = missing, 
    include = c(sex, agegrp_10yr, marstat_2cat, hiqual_4cat, gor, srh, ghqcase, 
                prior_jb, prior_srh, prior_ghqcase, prior_equivquint, wave_number), 
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
                 prior_jb = "Prior economic activity", 
                 prior_srh = "Prior self-rated health", 
                 prior_ghqcase = "Prior GHQ-12 caseness", 
                 prior_equivquint = "Prior equivalised household income quintile"
    )) %>% 
  modify_header(label = "**Variable**") %>%
  bold_labels() %>% 
  as_gt() %>% gt::tab_source_note(gt::md("*Percentages may not sum to 100% due to rounding*"))
# Print as HTML in the Viewer
table
# Export to a Word document for future use
gt::gtsave(table, filename = "output/missing_sample.docx")


#### 2. MISSING DUE TO ITEM NON-RESPONSE ####

# Conduct MSCA with an additional state for missing states
# Use a 'missing' state for the NAs for economic activity (jb) and self-rated health (srh)
SA <- "missing_state"

# Create a new state for missing data for both sequences
sample_5yr_complete_incNAs <- sample_5yr_complete_incNAs %>% 
  mutate(jb = ifelse(is.na(jb), "missing", jb), 
         srh = ifelse(is.na(srh), "missing", srh))

# Redefine the factor variables and palettes including the missing states         
jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other", "missing")         
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854", "#B3B3B3") # select and reorder colours to match jb_levels
sample_5yr_complete_incNAs$jb <- factor(sample_5yr_complete_incNAs$jb, levels = c(1:6, "missing"), labels = jb_levels)

srh_levels <- c("excellent", "very good", "good", "fair", "poor", "missing")
srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C", "#B3B3B3") # reorder colours to match srh_levels
sample_5yr_complete_incNAs$srh <- factor(sample_5yr_complete_incNAs$srh, levels = c(1:5, "missing"), labels = srh_levels)

# Create economic activity (jb) sequence dataframe
jb_seq_df <- sample_5yr_complete_incNAs %>% select(pidp, jb, transition_time) %>% filter(transition_time >= 0) %>% arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "jb")
jb_seq <- seqdef(jb_seq_df, 2:6, alphabet = jb_levels, cpal = jb_palette)

# Create self-rated health (srh) sequence dataframe
srh_seq_df <- sample_5yr_complete_incNAs %>% select(pidp, srh, transition_time) %>% filter(transition_time >= 0) %>% arrange(pidp, transition_time) |>
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "srh")
srh_seq <- seqdef(srh_seq_df, 2:6, alphabet = srh_levels, cpal = srh_palette)

# Create a multi-channel sequence combining economic activity and self-rated health sequences
multi_seq <- seqHMM::mc_to_sc_data(list(jb_seq, srh_seq))
alphabet(multi_seq)

# Create substitution-cost matrices using transition rates for each of the channels
sub_mats <- list(employment = TraMineR::seqsubm(jb_seq, method = "TRATE"),
                 health = TraMineR::seqsubm(srh_seq, method = "TRATE")) 

# Build multi-domain sequence 
channels <- list(employment = jb_seq, health = srh_seq)
mc_dist <- TraMineR::seqMD(channels, 
                           sm = sub_mats,        # subs = substitution matrix based on transition rates for both channels 
                           indel = 1,            # indels = 1
                           what = "diss", 
                           method = "OM",        # optimal matching algorithm
                           cweight = c(1, 1))    # equal weighting of both channels 

# Conduct clustering using Ward's agglomerative hierarchical method of clustering
Clustered <- hclust(as.dist(mc_dist), method = "ward.D2")

# Estimate cluster quality statistics using those included in the WeightedCluster package
Clustered_range <- WeightedCluster::as.clustrange(Clustered, diss = mc_dist, ncluster = 10)

# Table of cluster quality statistics
cluster_stats <- Clustered_range[["stats"]]
write.csv(x = cluster_stats, file = paste0("output/sensitivity/cluster_stats_", SA, ".csv"))

# Plot of cluster quality statistics
plot(Clustered_range, stat = "all", norm = "zscoremed", lwd = 2)
png(paste0("output/sensitivity/cluster_stats_", SA, ".png"), width = 2400, height = 1600, res = 300)
plot(Clustered_range, stat = c("PBC", "HG", "ASW", "HC", "CH"), norm = "zscoremed", lwd = 2)
dev.off()

# 4 clusters
Cuts <- cutree(Clustered, k = 4)
Groups <- factor(Cuts, labels = paste("Cluster", 1:4))

# Plot using seqHMM and gridplot
ssp_k <- list()
for(k in 1:4) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k], 
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}

# Export the seqHMM plot
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, 
                 ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()

# 3 clusters
Cuts <- cutree(Clustered, k = 3)
Groups <- factor(Cuts, labels = paste("Cluster", 1:3))

# Plot using seqHMM and gridplot
ssp_k <- list()
for(k in 1:3) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k], 
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}


# Export the seqHMM plot
png(paste0("output/sensitivity/plot_3clusters_", SA, ".png"), width = 900, height = 1000)
seqHMM::gridplot(ssp_k, 
                 ncol = 2, byrow = TRUE, 
                 row.prop = c(0.3, 0.3, 0.3,0.1), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()


#### 3. MISSING DUE TO LOSS-TO-FOLLOW-UP ####

# Inverse probability of attrition weighting (IAPW) #

# Create a binary indicator variable of inclusion in the sample
attrition_sample <- missing_combined %>% 
  filter(missing != "1yr_NAs") %>% # remove those who are excluded due to missing covariates
  mutate(followed_up = as.factor(case_when(missing %in% c("5yr_incNAs", "final") ~ 1,  # have follow-up in the study even if it has some NAs within it
                                 TRUE ~ 0))) # not followed up 
table(attrition_sample$followed_up)
class(attrition_sample$followed_up)

# Relevel education level 

attrition_sample$hiqual_4cat <-  factor(attrition_sample$hiqual_4cat,  levels = c("high", "medium", "low", "no"))

# Using WeightIt to conduct a weighting of the sample using ATO to balance as overlap of the two samples
weightit_object <- WeightIt::weightit(followed_up ~ sex + agegrp_10yr + fct_relevel(hiqual_4cat,  "high", "medium", "low", "no") + 
                                        prior_jb + prior_srh + prior_equivquint + prior_ghqcase,
                   data = attrition_sample, 
                   estimand = "ATO")  # ATO = weights to produces the same as the other method (this is the overlap sample method)

# Add these weights to the data
attrition_sample$weights <- weightit_object$weights

##### Appendix Figure: Love plot of covariate balance #####

# Love plot to show the balance of covariates
new_names <- c(prop.score = "Propensity score", 
               sex_male = "Sex: Male", 
               `agegrp_10yr_16-24` = "Age group: 16-24", 
               `agegrp_10yr_25-34` = "Age group: 25-34", 
               `agegrp_10yr_35-44` = "Age group: 35-44", 
               `agegrp_10yr_45-54` = "Age group: 45-54", 
               `agegrp_10yr_55-64` = "Age group: 55-64", 
               hiqual_4cat_high = "Educational level: High", 
               hiqual_4cat_medium = "Educational level: Medium", 
               hiqual_4cat_low = "Educational level: Low", 
               hiqual_4cat_no = "Educational level: No", 
               prior_jb_employed = "Prior economic activity: Employed", 
               prior_jb_unemployed = "Prior economic activity: Unemployed", 
               prior_jb_retired = "Prior economic activity: Retired", 
               prior_jb_care = "Prior economic activity: Care", 
               prior_jb_other = "Prior economic activity: Other", 
               prior_srh_excellent = "Prior self-rated health: Excellent", 
               `prior_srh_very good` = "Prior self-rated health: Very good", 
               prior_srh_good = "Prior self-rated health: Good", 
               prior_srh_fair = "Prior self-rated health: Fair", 
               prior_srh_poor = "Prior self-rated health: Poor", 
               prior_equivquint_Q1 = "Prior household income quintile: Q1", 
               prior_equivquint_Q2 = "Prior household income quintile: Q2", 
               prior_equivquint_Q3 = "Prior household income quintile: Q3", 
               prior_equivquint_Q4 = "Prior household income quintile: Q4", 
               prior_equivquint_Q5 = "Prior household income quintile: Q5", 
               prior_ghqcase_Yes = "Prior GHQ-12 caseness: Yes"
               )

cobalt::set.cobalt.options(binary = "std") # standardised mean differences
love_plot <- cobalt::love.plot(weightit_object, var.names = new_names)
love_plot

png("output/sensitivity/love_plot.png", width = 2400, height = 1600, res = 300)
love_plot
dev.off()

# Additional plot of propensity score distributions
bal.plot(weightit_object, un = TRUE)

# Add weights to the overall data
attrition_weights <- attrition_sample %>% select(pidp, weights)

# Add these weights to the main data
sample_weighted <- filter(missing_data, pidp %in% unique(attrition_weights$pidp)) %>% left_join(attrition_weights) %>% 
  filter(transition_time > -1 & transition_time < 5) %>% filter(!is.na(missing_NAs))
length(unique(sample_weighted$pidp))                            

## Conduct MCSA with IAPW weighting ###
SA <- "IAPW"

# Create a new state for missing data for both sequences
sample_weighted <- sample_weighted %>% 
  mutate(jb = ifelse(is.na(jb), "missing", jb), 
         srh = ifelse(is.na(srh), "missing", srh))

# Redefine the factor variables and palettes including the missing states         
jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other", "missing") 
jb_labels <- c("Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other", "Missing")
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854", "#B3B3B3") # select and reorder colours to match jb_levels
sample_weighted$jb <- factor(sample_weighted$jb, levels = c(1:6, "missing"), labels = jb_levels)

srh_levels <- c("excellent", "very good", "good", "fair", "poor", "missing")
srh_labels <- c("Excellent", "Very good", "Good", "Fair", "Poor", "Missing")
srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C", "#B3B3B3") # reorder colours to match srh_levels
sample_weighted$srh <- factor(sample_weighted$srh, levels = c(1:5, "missing"), labels = srh_levels)

# Create economic activity (jb) sequence dataframe
jb_seq_df <- sample_weighted %>% select(pidp, weights, jb, transition_time) %>% filter(transition_time >= 0) %>% arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = c("pidp", "weights"), names_from = "transition_time", values_from = "jb")
jb_seq <- seqdef(jb_seq_df, 3:7, alphabet = jb_levels, cpal = jb_palette, labels = jb_labels, weights = jb_seq_df$weights)

# Create self-rated health (srh) sequence dataframe
srh_seq_df <- sample_weighted %>% select(pidp, weights, srh, transition_time) %>% filter(transition_time >= 0) %>% arrange(pidp, transition_time) |>
  pivot_wider(id_cols = c("pidp", "weights"), names_from = "transition_time", values_from = "srh")
srh_seq <- seqdef(srh_seq_df, 3:7, alphabet = srh_levels, cpal = srh_palette, labels = srh_labels, weights = srh_seq_df$weights)

# Create a multi-channel sequence combining economic activity and self-rated health sequences
multi_seq <- seqHMM::mc_to_sc_data(list(jb_seq, srh_seq))
alphabet(multi_seq)

# Create substitution-cost matrices using transition rates for each of the channels
sub_mats <- list(employment = TraMineR::seqsubm(jb_seq, method = "TRATE", weighted = TRUE),
                 health = TraMineR::seqsubm(srh_seq, method = "TRATE", weighted = TRUE)) 

# Build multi-domain sequence 
channels <- list(employment = jb_seq, health = srh_seq)
mc_dist <- TraMineR::seqMD(channels, 
                           sm = sub_mats,        # subs = substitution matrix based on transition rates for both channels - having an error!
                           indel = 1,            # indels = 1
                           what = "diss", 
                           method = "OM",        # optimal matching algorithm
                           cweight = c(1, 1))    # equal weighting of both channels 

# Conduct clustering using Ward's agglomerative hierarchical method of clustering
Clustered <- hclust(as.dist(mc_dist), method = "ward.D2", members = jb_seq_df$weights)

# Estimate cluster quality statistics using those included in the WeightedCluster package
Clustered_range <- WeightedCluster::as.clustrange(Clustered, diss = mc_dist, ncluster = 10, weights = jb_seq_df$weights)

# Table of cluster quality statistics
cluster_stats <- Clustered_range[["stats"]]
write.csv(x = cluster_stats, file = paste0("output/sensitivity/cluster_stats_", SA, ".csv"))

# Plot of cluster quality statistics
plot(Clustered_range, stat = "all", norm = "zscoremed", lwd = 2)
png(paste0("output/sensitivity/cluster_stats_", SA, ".png"), width = 2400, height = 1600, res = 300)
plot(Clustered_range, stat = c("PBC", "HG", "ASW", "HC", "CH"), norm = "zscoremed", lwd = 2)
dev.off()

# 4 clusters 
Cuts <- cutree(Clustered, k = 4)
Groups <- factor(Cuts, labels = paste("Cluster", 1:4))

# Plot using seqHMM
ssp_k <- list()
for(k in 1:4) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k])
}

# Export the seqHMM plot
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, 
                 ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()

##### Appendix Figure A3: State distribution plots for each cluster for IPAW and missing states #####

# Reorder the clusters to be more informative and logical with cluster names
new_labels <- c("Cluster 4", "Cluster 2", "Cluster 3", "Cluster 1") 
new_Groups <- factor(Groups, levels = paste("Cluster", 1:4), labels = new_labels)
new_Cuts <- as.numeric(gsub("Cluster", "", new_Groups))

# Cluster 1: Entering employment
ssp_k[[1]] <- ssp(list(jb_seq[new_Cuts == 1,], 
                       srh_seq[new_Cuts == 1,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 1: Entering employment"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Cluster 2: Healthy but inactive
ssp_k[[2]] <- ssp(list(jb_seq[new_Cuts == 2,], 
                       srh_seq[new_Cuts == 2,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 2: Healthy but inactive"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Cluster 3: Continued long-term sickness
ssp_k[[3]] <- ssp(list(jb_seq[new_Cuts == 3,], 
                       srh_seq[new_Cuts == 3,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 3: Continued long-term sickness"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Cluster 4: Transition to retirement
ssp_k[[4]] <- ssp(list(jb_seq[new_Cuts == 4,], 
                       srh_seq[new_Cuts == 4,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 4: Transition to retirement"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Save publication-quality plot
png(paste0("output/sensitivity/plot_4clusters_", SA, "_publication.png"), width = 4500, height = 5500, res=500)
gridplot(ssp_k, 
         ncol = 2, byrow = TRUE, 
         row.prop = c(0.41, 0.41, 0.18), 
         cex.legend = 1, 
         title.legend = "auto", 
)
dev.off()





  









