
## ADDITIONAL ANALYSIS ##

# These are additional analyses conducted for the interest for my PhD rather than the paper

#### 1. SF-12 ####

# Read in the regression sample for the cluster numbers
regression_sample <- readRDS("data/processed/regression_sample.rds")
cluster_numbers <- select(regression_sample, pidp, cluster, group)

# Read in the longitudinal data
sample <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_complete_noNAs.rds")

# Add on cluster numbers
sequence_sf12 <- left_join(sample, cluster_numbers) # PIDP only 

# Reorder the cluster group 
sequence_sf12$group <- factor(sequence_sf12$group, 
                              c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))

# Calculate mean values for SF-12 PCS and MCS summary scores 
mean_sequence_sf12 <- sequence_sf12 %>% 
  group_by(transition_time, group) %>% 
  dplyr::summarize(mean_sf12pcs = mean(sf12pcs, na.rm=TRUE), 
                   mean_sf12mcs = mean(sf12mcs, na.rm=TRUE))

# Basic plot with facets for each cluster
mean_sequence_sf12 %>% 
  pivot_longer(cols = mean_sf12pcs:mean_sf12mcs, 
               names_to = "type", 
               values_to = "value") %>% 
  ggplot(aes(x=transition_time, y=value, color = type)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~group)

# Basic plot with facets for PCS and MCS seperately
mean_sequence_sf12 %>% 
  pivot_longer(cols = mean_sf12pcs:mean_sf12mcs, 
               names_to = "type", 
               values_to = "value") %>% 
  ggplot(aes(x=transition_time, y=value, color = group)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~type)


## Sample size ##

# Total in main analysis
sequence_sf12 %>% distinct(pidp) %>% tally() # 919

# Those with some observed SF-12 information
sequence_sf12 %>% 
  filter(!is.na(sf12pcs) | !is.na(sf12mcs)) %>% 
  distinct(pidp) %>%  
  tally() # 918

# Those with complete SF-12 throughout the follow-up period
sf12_complete <- sequence_sf12 %>% 
  select(pidp, wave_number, transition_time, group, sf12mcs, sf12pcs) %>% 
  drop_na() %>% 
  dplyr::group_by(pidp) %>% dplyr::summarize(n = n()) %>% 
  filter(n == 6) 
sf12_complete %>% distinct(pidp) %>% tally() 

# Recalculate mean values with only those with complete information for SF-12
mean_sequence_sf12_complete <- sequence_sf12 %>% 
  select(pidp, wave_number, transition_time, group, sf12mcs, sf12pcs) %>% 
  filter(pidp %in% sf12_complete$pidp) %>% 
  group_by(transition_time, group) %>% 
  dplyr::summarize(mean_sf12pcs = mean(sf12pcs, na.rm=TRUE), 
                   mean_sf12mcs = mean(sf12mcs, na.rm=TRUE))


## Mean SF-12 over follow-up period, by cluster ##

mean_sequence_sf12 %>% 
  pivot_longer(cols = mean_sf12pcs:mean_sf12mcs, 
               names_to = "type", 
               values_to = "value") %>%
  mutate(type = recode(type, "mean_sf12mcs" = "SF-12 MCS", "mean_sf12pcs" = "SF-12 PCS")) %>%
  ggplot(aes(x=transition_time, y=value, color = group)) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~type) + 
  theme_bw()+ 
  geom_vline(xintercept = -1, linetype = 2, color = "darkgrey") + 
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") + 
  geom_vline(xintercept = 1, linetype = 2, color = "darkgrey") + 
  labs(y = "SF-12 Score", color = "Cluster", x = "Year since onset of long-term sickness", 
       caption = "N = 811\nOnly those with complete SF-12 for all time-points") + 
  scale_color_discrete(breaks = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
                       labels = c("Cluster 1: Entering employment", 
                                  "Cluster 2: Precarious non-employment", 
                                  "Cluster 3: Continued long-term sickness", 
                                  "Cluster 4: Rapid retirement")) 

## Adding confidence intervals ##

# MCS
summary_mcs <- sequence_sf12 %>% as.data.frame() %>% 
  select(pidp, wave_number, transition_time, group, sf12mcs, sf12pcs) %>% 
  filter(pidp %in% sf12_complete$pidp) %>% 
  Rmisc::summarySE(measurevar= "sf12mcs", groupvars=c("group", "transition_time"), na.rm=TRUE) %>% 
  mutate(sf12 = sf12mcs, type = "SF-12 MCS") %>% select(-sf12mcs)
# PCS
summary_pcs <- sequence_sf12 %>% as.data.frame() %>% 
  select(pidp, wave_number, transition_time, group, sf12mcs, sf12pcs) %>% 
  filter(pidp %in% sf12_complete$pidp) %>% 
  Rmisc::summarySE(measurevar= "sf12pcs", groupvars=c("group", "transition_time"), na.rm=TRUE) %>% 
  mutate(sf12 = sf12pcs, type = "SF-12 PCS") %>% select(-sf12pcs)

summary <- rbind(summary_mcs, summary_pcs)

##### Additional Figure: Mean SF-12 over follow-up period with 95% CIs, by cluster #####

sf12_plot <- ggplot(summary, aes(x=transition_time, y=sf12, colour=group)) + 
  geom_errorbar(aes(ymin=sf12-se, ymax=sf12+se), width=0.3) + 
  geom_line() + 
  geom_point() + 
  facet_wrap(~type) + 
  theme_bw()+ 
  geom_vline(xintercept = 0, linetype = 2, color = "darkgrey") + 
  geom_vline(xintercept = 1, linetype = 2, color = "darkgrey") + 
  labs(y = "SF-12 Score", color = "Cluster", x = "Year since onset of long-term sickness", 
       caption = "N = 811\nOnly those with complete SF-12 for all time-points") + 
  scale_color_discrete(breaks = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
                       labels = c("Cluster 1: Entering employment", 
                                  "Cluster 2: Healthy but inactive", 
                                  "Cluster 3: Continued long-term sickness", 
                                  "Cluster 4: Transition to retirement")) 

png("output/extra/sf12_plot.png", width = 1211, height = 1000, res=200)
sf12_plot
dev.off()



#### 2. FOLLOW-UP TIME VARIATION ####

# Those with possible follow-up following transition to long-term sickness (but no missing covariates)
sample_baseline_complete <- readRDS("data/processed/lts_transitions_1yr_prior_noNAs_5yr_followup_possible.rds")

# Filter to those with 4 years of complete consecutive follow-up
followup_times <- c(1:4)
followup_list <- list()

for (time in followup_times) { 
  
  # Filter to those observations at transition times less than the follow-up
  lts_transitions_max_followup <- sample_baseline_complete %>% filter(transition_time <= time)
  
  # Tally of how many timepoints is an individual observed at
  tally <- lts_transitions_max_followup %>% group_by(pidp) %>% tally()
  
  # Filter to those with complete observation history throughout this follow-up period
  lts_transitions_1yr_prior_followup_complete_incNAs <- left_join(lts_transitions_max_followup, tally, by = "pidp") %>% 
    filter(n == (time + 2)) %>% # include times 0 and time - 1 as 2 additional timepoints
    select(-n)
  
  # Count the number of individuals
  followup_with_nas <- lts_transitions_1yr_prior_followup_complete_incNAs %>% distinct(pidp) %>% count() %>% pull()
  
  # Export out the complete data
  assign(paste0("followup_with_NAs_", time), lts_transitions_1yr_prior_followup_complete_incNAs)
  
  # Mark out who has missing economic activity data for the sequence follow-up
  lts_transitions_jb_NA <- lts_transitions_1yr_prior_followup_complete_incNAs %>% 
    filter(transition_time >= 0) %>% 
    group_by(pidp) %>% 
    summarise(jb_NA = sum(as.numeric(jb))) %>% # convert to numeric and summarise as sum will be NA if one value is NA
    mutate(jb_NA = ifelse(is.na(jb_NA), 1, 0)) # indicator on whether any value of JB is missing
  
  # Table of how many have missing economic activity data
  table(lts_transitions_jb_NA$jb_NA)
  
  # Mark out who has missing self-rated health data for the sequence follow-up
  lts_transitions_srh_NA <- lts_transitions_1yr_prior_followup_complete_incNAs %>% 
    filter(transition_time >= 0) %>% 
    group_by(pidp) %>% 
    summarise(srh_NA = sum(as.numeric(srh))) %>% # convert to numeric and summarise as sum will be NA if one value is NA
    mutate(srh_NA = ifelse(is.na(srh_NA), 1, 0)) # indicator on whether any value of SRH is missing
  
  # Table of how many have missing self-rated health data
  table(lts_transitions_srh_NA$srh_NA)
  
  # Join missing indicators to the main data to allow for a filter to remove those with missing states
  lts_transitions_1yr_prior_followup_complete_noNAs <- lts_transitions_1yr_prior_followup_complete_incNAs %>% 
    left_join(lts_transitions_jb_NA, by = "pidp") %>% filter(jb_NA == 0) %>% 
    left_join(lts_transitions_srh_NA, by = "pidp") %>% filter(srh_NA == 0) %>% 
    select(-jb_NA, -srh_NA)
  
  # Count the number of individuals
  followup_without_nas <- lts_transitions_1yr_prior_followup_complete_noNAs %>% distinct(pidp) %>% count() %>% pull()
  
  # Export out the complete data
  assign(paste0("followup_complete_noNAs_", time), lts_transitions_1yr_prior_followup_complete_noNAs)
  
  # Store in a list
  followup_list[[time]] <- tibble(
    followup_time = time, 
    followup_count_with_nas = followup_with_nas, 
    followup_count_without_nas = followup_without_nas
  )
  
} 

# Make a little table
followup_table <- bind_rows(followup_list)
followup_table <- followup_table %>% 
  mutate(percent_remaining_with_nas = followup_count_with_nas/max(followup_count_with_nas), 
         percent_remaining_without_nas = followup_count_without_nas/max(followup_count_without_nas))
# Export
write.csv(x = followup_table, file = "output/eztra/followup_table_9.csv")

# Make a little plot (because why not but not sure which is better sometimes?)
followup_plot <- followup_table %>% 
  select(followup_time, followup_count_with_nas, followup_count_without_nas) %>% 
  pivot_longer(!followup_time, names_to = "missingness", values_to = "count") %>% 
  ggplot(aes(x=as.factor(followup_time), y=count, colour=missingness, group = missingness)) + 
  geom_point(size = 2) + 
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  labs(x="Total number of years of follow-up since year of onset of long-term sickness", 
       y = "Number of individuals", 
       colour = "Missingness") + 
  scale_color_discrete(palette=scales::pal_brewer(palette="Set1"),
                       label=c("NAs due to item non-response", 
                               "Complete responses with no NAs"))
followup_plot
# Export
png("output/extra/followup_plot.png", width = 1211, height = 1000, res=200)
followup_plot
dev.off()

# Plot of sample remaining
sample_remaining_plot <- followup_table %>% 
  select(followup_time, percent_remaining_with_nas, percent_remaining_without_nas) %>% 
  pivot_longer(!followup_time, names_to = "missingness", values_to = "count") %>% 
  ggplot(aes(x=as.factor(followup_time), y=count, colour=missingness, group = missingness)) + 
  geom_point(size = 2) + 
  geom_line(linewidth = 1) + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  labs(x="Total number of years of follow-up since year of onset of long-term sickness", 
       y = "Number of individuals", 
       colour = "Missingness") + 
  scale_color_discrete(palette=scales::pal_brewer(palette="Set1"),
                       label=c("NAs due to item non-response", 
                               "Complete responses with no NAs")) + 
  scale_y_continuous(labels=scales::percent)
sample_remaining_plot
# Export
png("output/extra/sample_remaining_plot.png", width = 1211, height = 1000, res=200)
sample_remaining_plot
dev.off()

# Run the sequence analysis with each of the follow-up times
for (time in followup_times) { 
  
  sample <- get(paste0("followup_complete_noNAs_", time))
  print(length(unique(sample$pidp))) 
  
  # Create economic activity (jb) sequence dataframe
  jb_seq_df <- sample %>% select(pidp, jb, transition_time) %>% filter(transition_time >= 0) %>% arrange(pidp, transition_time) %>% 
    pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "jb")
  
  jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other")
  jb_labels <- c("Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other")
  brewer.pal(n=6, "Set2")
  jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854") # select and reorder colours to match jb_levels
  
  jb_seq <- seqdef(jb_seq_df, 2:(time+2), alphabet = jb_levels, labels = jb_labels, cpal = jb_palette)
  
  # Create self-rated health (srh) sequence dataframe
  srh_seq_df <- sample %>% select(pidp, srh, transition_time) %>% filter(transition_time >= 0) %>% arrange(pidp, transition_time) |>
    pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "srh")
  
  srh_levels <- c("excellent", "very good", "good", "fair", "poor")
  srh_labels <- c("Excellent", "Very good", "Good", "Fair", "Poor")
  brewer.pal(n=5, "RdYlBu")
  srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C") # reorder colours to match srh_levels
  
  srh_seq <- seqdef(srh_seq_df, 2:(time+2), alphabet = srh_levels, labels = srh_labels, cpal = srh_palette) 
  
  # Create a multi-channel sequence combining economic activity and self-rated health sequences
  multi_seq <- seqHMM::mc_to_sc_data(list(jb_seq, srh_seq))
  alphabet(multi_seq)
  
  # Subs = Matrices based on transition rates
  # Indels = 1
  # Method = Optimal Matching
  
  # Create a multi-channel sequence combining economic activity and self-rated health sequences
  multi_seq <- seqHMM::mc_to_sc_data(list(jb_seq, srh_seq))
  alphabet(multi_seq)
  
  # Create substitution-cost matrices using transition rates for each of the channels
  sub_mats_basecase <- list(employment = TraMineR::seqsubm(jb_seq, method = "TRATE"),
                            health = TraMineR::seqsubm(srh_seq, method = "TRATE")) 
  
  # Build multi-domain sequence 
  channels <- list(employment = jb_seq, health = srh_seq)
  mc_dist <- TraMineR::seqMD(channels, 
                             sm = sub_mats_basecase, # subs = substitution matrix based on transition rates for both channels 
                             indel = 1,              # indels = 1
                             what = "diss", 
                             method = "OM",        # optimal matching algorithm
                             cweight = c(1, 1))    # equal weighting of both channels 
  
  # Conduct clustering using Ward's agglomerative hierarchical method of clustering
  Clustered <- hclust(as.dist(mc_dist), method = "ward.D2")
  
  # Assign into the global environment
  assign(paste0("Clustered_", time), Clustered)
  
  # Estimate cluster quality statistics using those included in the WeightedCluster package
  Clustered_range <- WeightedCluster::as.clustrange(Clustered, diss = mc_dist, ncluster = 10)
  
  # Table of cluster quality statistics
  cluster_stats <- Clustered_range[["stats"]]
  write.csv(x = cluster_stats, file = paste0("output/eztra/cluster_stats_followup_", time, ".csv"))
  assign(paste0("cluster_stats_", time), cluster_stats)
  
  # Plot of cluster quality statistics
  plot(Clustered_range, stat = "all", norm = "zscoremed", lwd = 2)
  png(paste0("output/extra/cluster_stats_followup_", time, ".png"), width = 2400, height = 1600, res = 300)
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
                              title = levels(Groups)[k], 
                              ylab = c("Economic activity", "Self-rated health"),
                              xlab = "Year since onset of long-term sickness")
  }
  # Export the seqHMM plot
  png(paste0("output/extra/plot_4clusters_followup_", time, ".png"), width = 900, height = 700)
  seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
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
  png(paste0("output/extra/plot_3clusters_followup_", time, ".png"), width = 900, height = 1000)
  seqHMM::gridplot(ssp_k, 
                   ncol = 2, byrow = TRUE, 
                   row.prop = c(0.3, 0.3, 0.3,0.1), cex.legend = 0.8, 
                   title.legend = "auto")
  dev.off()
  
}

# Little plot comparing the ASW of the cluster stats
asw_1 <- cluster_stats_1 %>% mutate(followup = 1, cluster = row_number()+1) %>% select(cluster, followup, ASW)
asw_2 <- cluster_stats_2 %>% mutate(followup = 2, cluster = row_number()+1) %>% select(cluster, followup, ASW)
asw_3 <- cluster_stats_3 %>% mutate(followup = 3, cluster = row_number()+1) %>% select(cluster, followup, ASW)
asw_4 <- cluster_stats_4 %>% mutate(followup = 4, cluster = row_number()+1) %>% select(cluster, followup, ASW)
asw_to_plot <- rbind(asw_1, asw_2, asw_3, asw_4)

asw_followup_plot <- ggplot(asw_to_plot, aes(x=cluster, y=ASW, color=as.factor(followup))) + 
  geom_point(size = 2) + 
  geom_line(linewidth = 1) + 
  scale_color_discrete(palette = scales::brewer_pal(palette = "PuOr")) + 
  theme_minimal() + 
  theme(legend.position="top") + 
  labs(x="Number of clusters",
       y="Average Silhouette Width (ASW)", 
       color="Years of follow-up after the onset of long-term sickness")
asw_followup_plot
# Export
png("output/extra/asw_followup_plot.png", width = 1211, height = 1000, res=200)
asw_followup_plot
dev.off()


