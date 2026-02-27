
## RUN MULTI-CHANNEL SEQUENCE ANALYSIS ##

rm(list=ls())


#### 1. SEQUENCE DATA ####

# Read in analysis sample data
sample <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_complete_noNAs.rds")

# Create economic activity (jb) sequence dataframe
jb_seq_df <- sample %>% select(pidp, jb, transition_time) %>% 
  filter(transition_time >= 0) %>% arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "jb")

jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other")
jb_labels <- c("Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other")
brewer.pal(n=6, "Set2")
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854") # select and reorder colours to match jb_levels

jb_seq <- TraMineR::seqdef(jb_seq_df, 2:6, alphabet = jb_levels, labels = jb_labels, cpal = jb_palette)

# Plot economic activity sequences
TraMineR::seqdplot(jb_seq, border = NA,  ncol = 3, legend.prop = 0.2) # distribution plot
TraMineR::seqIplot(jb_seq, border = NA,  ncol = 3, legend.prop = 0.2, sortv = "from.start") # index plot

# Create self-rated health (srh) sequence dataframe
srh_seq_df <- sample %>% select(pidp, srh, transition_time) %>% 
  filter(transition_time >= 0) %>% arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "srh")

srh_levels <- c("excellent", "very good", "good", "fair", "poor")
srh_labels <- c("Excellent", "Very good", "Good", "Fair", "Poor")
brewer.pal(n=5, "RdYlBu")
srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C") # reorder colours to match srh_levels

srh_seq <- TraMineR::seqdef(srh_seq_df, 2:6, alphabet = srh_levels, labels = srh_labels, cpal = srh_palette)

# Plot self-rated health sequences
TraMineR::seqdplot(srh_seq, border = NA,  ncol = 3, legend.prop = 0.2) # distribution plot
TraMineR::seqIplot(srh_seq, border = NA,  ncol = 3, legend.prop = 0.2, sortv = "from.start") # index plot

# Create a multi-channel sequence combining economic activity and self-rated health sequences
multi_seq <- seqHMM::mc_to_sc_data(list(jb_seq, srh_seq))
alphabet(multi_seq)

# Plot for the extended alphabet approach (i.e., categories for each combination)
TraMineR::seqdplot(multi_seq, border = NA, legend.prop = 0.2) # distribution plot
TraMineR::seqIplot(multi_seq, sortv = "from.start") # index plot

##### Figures: Distribution and index plots for the overall sample #####

# Plots showing the channels separately using seqHMM package functions

# State distribution plot for overall sample
png("output/main/overall_distribution_plot.png", width = 2200, height = 1500, res = 300)
plot(seqHMM::ssp(list("Economic activity" = jb_seq, 
                      "Self-rated health" = srh_seq),
                 type = "d", # state distribution plot
                 xlab = "Year since onset of long-term sickness",
                 ylab = c("Economic activity", "Self-rated health"), , 
                 with.legend = TRUE, 
                 border = NA, 
                 xtlab = 0:4))
dev.off()

# Sequence index plot for overall sample
png("output/main/overall_index_plot.png", width = 2200, height = 1500, res = 300)
plot(seqHMM::ssp(list("Economic activity" = jb_seq, 
                      "Self-rated health" = srh_seq),
                 type = "I", # sequence index plot
                 sortv = "from.start", 
                 sort.channel = 1, # sort using the first channel (economic activity)
                 xlab = "Year since onset of long-term sickness",
                 ylab = c("Economic activity", "Self-rated health"), 
                 with.legend = TRUE, 
                 border = NA, 
                 xtlab = 0:4))
dev.off()

##### Appendix Tables: Matrices of substitution costs #####

# Create substitution cost matrices using transition rates for each of the channels
sub_mats <- list(employment = TraMineR::seqsubm(jb_seq, method = "TRATE"),
                 health = TraMineR::seqsubm(srh_seq, method = "TRATE")) 

# Export these matrices
write.csv(round(sub_mats$employment, 2), file = "output/main/jb_transition_matrix.csv")
write.csv(round(sub_mats$health, 2), file = "output/main/srh_transition_matrix.csv")


# Build multi-domain sequence using this matrices for the substition costs
channels <- list(employment = jb_seq, health = srh_seq)
mc_dist <- TraMineR::seqMD(channels, 
                           sm = sub_mats,        # subs = substitution matrix based on transition rates for both channels 
                           indel = 1,            # indels = 1
                           what = "diss", 
                           method = "OM",        # optimal matching algorithm
                           cweight = c(1, 1))    # equal weighting of both channels 


#### 2. CLUSTERING ####

# Conduct clustering using Ward's agglomerative hierarchical method of clustering
Clustered <- hclust(as.dist(mc_dist), method = "ward.D2")

# Estimate cluster quality statistics using those included in the WeightedCluster package
Clustered_range <- WeightedCluster::as.clustrange(Clustered, diss = mc_dist, ncluster = 10)

# Table: Cluster quality statistics
cluster_stats <- Clustered_range[["stats"]]
write.csv(x = cluster_stats, file = "output/main/cluster_stats_basecase.csv")

##### Appendix Figure: Cluster quality statistics #####

# Plot of cluster quality statistics
plot(Clustered_range, stat = "all", norm = "zscoremed", lwd = 2) # all measures on a plot
png("output/main/cluster_stats_basecase.png", width = 2400, height = 1600, res = 300)
plot(Clustered_range, stat = c("PBC", "HG", "ASW", "HC", "CH"), norm = "zscoremed", lwd = 2) # selected measures
dev.off()


#### 3. PLOTS ####

# 4 clusters
Cuts <- cutree(Clustered, k = 4)
Groups <- factor(Cuts, labels = paste("Cluster", 1:4))

# Individual plots for each cluster
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

# Combine plots together
png("output/main/plot_4clusters_basecase.png", width = 900, height = 700)
gridplot(ssp_k, 
         ncol = 2, byrow = TRUE, 
         row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
         title.legend = "auto" 
)
dev.off()


##### Figure: State distribution plots for each cluster ####

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
png("output/main/plot_4clusters_basecase_publication.png", width = 4300, height = 4816, res=500)
gridplot(ssp_k, 
         ncol = 2, byrow = TRUE, 
         row.prop = c(0.41, 0.41, 0.18), 
         cex.legend = 1, 
         title.legend = "auto")
dev.off()

##  Plots of each individual cluster for presentations ##

# Cluster 1: Entering employment
png("output/main/plot_distribution_cluster1.png", width = 900, height = 600, res = 125)
plot(ssp(list(jb_seq[new_Cuts == 1,], 
         srh_seq[new_Cuts == 1,]),
    type = "d", border = NA, 
    with.legend = TRUE, 
    title = paste0("Cluster 1: Entering employment"), 
    ylab = c("Economic activity", "Self-rated health"),
    xlab = "Year since onset of long-term sickness", 
    xtlab = 0:4))
dev.off()

# Cluster 2: Healthy but inactive
png("output/main/plot_distribution_cluster2.png", width = 900, height = 600, res = 125)
plot(ssp(list(jb_seq[new_Cuts == 2,], 
              srh_seq[new_Cuts == 2,]),
         type = "d", border = NA, 
         with.legend = TRUE, 
         title = paste0("Cluster 2: Healthy but inactive"), 
         ylab = c("Economic activity", "Self-rated health"),
         xlab = "Year since onset of long-term sickness", 
         xtlab = 0:4))
dev.off()

# Cluster 3: Continued long-term sickness
png("output/main/plot_distribution_cluster3.png", width = 900, height = 600, res = 125)
plot(ssp(list(jb_seq[new_Cuts == 3,], 
              srh_seq[new_Cuts == 3,]),
         type = "d", border = NA, 
         with.legend = TRUE, 
         title = paste0("Cluster 3: Continued long-term sickness"), 
         ylab = c("Economic activity", "Self-rated health"),
         xlab = "Year since onset of long-term sickness", 
         xtlab = 0:4))
dev.off()

# Cluster 4: Transition to retirement
png("output/main/plot_distribution_cluster4.png", width = 900, height = 600, res = 125)
plot(ssp(list(jb_seq[new_Cuts == 4,], 
              srh_seq[new_Cuts == 4,]),
         type = "d", border = NA, 
         with.legend = TRUE, 
         title = paste0("Cluster 4: Transition to retirement"), 
         ylab = c("Economic activity", "Self-rated health"),
         xlab = "Year since onset of long-term sickness", 
         xtlab = 0:4))
dev.off()

##### Appendix Figure: Sequence index plot by cluster #####

# Cluster 1: Entering employment
ssp_k[[1]] <- ssp(list(jb_seq[new_Cuts == 1,], 
                       srh_seq[new_Cuts == 1,]),
                  type = "I", sortv = "from.start", sort.channel = 1, 
                  border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 1: Entering employment"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Cluster 2: Healthy but inactive
ssp_k[[2]] <- ssp(list(jb_seq[new_Cuts == 2,], 
                       srh_seq[new_Cuts == 2,]),
                  type = "I", sortv = "from.start", sort.channel = 1,
                  border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 2: Healthy but inactive"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Cluster 3: Continued long-term sickness
ssp_k[[3]] <- ssp(list(jb_seq[new_Cuts == 3,], 
                       srh_seq[new_Cuts == 3,]),
                  type = "I", sortv = "from.start", sort.channel = 1, 
                  border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 3: Continued long-term sickness"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Cluster 4: Transition to retirement
ssp_k[[4]] <- ssp(list(jb_seq[new_Cuts == 4,], 
                       srh_seq[new_Cuts == 4,]),
                  type = "I", sortv = "from.start", sort.channel = 1, 
                  border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 4: Transition to retirement"), 
                  ylab = c("Economic activity", "Self-rated health"),
                  xlab = "Year since onset of long-term sickness", 
                  xtlab = 0:4)

# Save publication-quality plot
png("output/main/plot_index_4clusters_publication.png", width = 4300, height = 4816, res=500)
gridplot(ssp_k, 
         ncol = 2, byrow = TRUE, 
         row.prop = c(0.41, 0.41, 0.18), 
         cex.legend = 1, 
         title.legend = "auto")
dev.off()


#### 4. PREPROCESS FOR MULTINOMIAL LOGISTIC REGRESSION ####

# Filter to baseline sample 
sample_baseline <- filter(sample, transition_time == 0)

# Add on the cluster numbers 
sample_baseline$cluster <- new_Cuts
sample_baseline$group <- new_Groups

# Filter to characteristics prior to the sample and add a prefix
sample_prior <- filter(sample, transition_time == -1) %>% 
  select(pidp, jb, srh, ghqcase, equivquint) %>% 
  dplyr::rename_with(.fn = ~paste0("prior_", .x), .cols = !pidp)
  
# Join with baseline 
regression_sample <- left_join(sample_baseline, sample_prior)
  
# Number and percentages in each of the clusters
table(regression_sample$cluster)
round(prop.table(table(regression_sample$cluster))*100,0)

# Save as an RDS format (.rds)
saveRDS(regression_sample, file = "data/processed/regression_sample.rds")


#### 5. DESCRIPTION OF COVARIATES BY CLUSTER ####

regression_sample <- readRDS("data/processed/regression_sample.rds")

##### Table: Analytical sample by covariate by cluster #####

table <- regression_sample %>%
  gtsummary::tbl_summary(
    by = cluster, 
    include = c(sex, agegrp_10yr, marstat_2cat, hiqual_4cat, gor, srh, ghqcase, 
                prior_jb, prior_srh, prior_ghqcase, prior_equivquint), 
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
                 prior_equivquint = "Prior equvalised household income quintile"
    )) %>% 
  modify_header(label = "**Variable**") %>%
  bold_labels() %>% 
  as_gt() %>% gt::tab_source_note(gt::md("*Percentages may not sum to 100% due to rounding*"))
# Print as HTML in the Viewer
table
# Export to a Word document for future use
gt::gtsave(table, filename = "output/main/table_sample_by_cluster.docx")

