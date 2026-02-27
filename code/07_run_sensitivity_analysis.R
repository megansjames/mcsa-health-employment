
## RUN SENSITIVITY ANALYSES FOR MULTI-CHANNEL SEQUENCE ANALYSIS ##

rm(list=ls())

# Read in analysis sample data used in the main analysis
sample <- readRDS("data/processed/lts_transitions_1yr_prior_5yr_followup_complete_noNAs.rds")

#### 1. DIFFERENT MCSA SPECIFICATIONS ####

# Set-up of sequence objects for different MCSA specifications (costs for subs and indels and dissimilarity methods)

# Create economic activity (jb) sequence dataframe
jb_seq_df <- sample %>% select(pidp, jb, transition_time) %>% 
  filter(transition_time >= 0) %>% 
  arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "jb")

jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other")
jb_labels <- c("Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other")
brewer.pal(n=6, "Set2")
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854") # select and reorder colours to match jb_levels

jb_seq <- seqdef(jb_seq_df, 2:6, alphabet = jb_levels, labels = jb_labels, cpal = jb_palette)


# Create self-rated health (srh) sequence dataframe
srh_seq_df <- sample %>% select(pidp, srh, transition_time) %>% 
  filter(transition_time >= 0) %>% 
  arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "srh")

srh_levels <- c("excellent", "very good", "good", "fair", "poor")
srh_labels <- c("Excellent", "Very good", "Good", "Fair", "Poor")
brewer.pal(n=5, "RdYlBu")
srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C") # reorder colours to match srh_levels

srh_seq <- seqdef(srh_seq_df, 2:6, alphabet = srh_levels, labels = srh_labels, cpal = srh_palette)

# Create a multi-channel sequence combining economic activity and self-rated health sequences
multi_seq <- seqHMM::mc_to_sc_data(list(jb_seq, srh_seq))
alphabet(multi_seq)


##### 0: BASECASE #####

SA <- "basecase"

# Subs = Matrices based on transition rates
# Indels = 1
# Method = Optimal Matching

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
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()

# 6 clusters 
Cuts <- cutree(Clustered, k = 6)
Groups <- factor(Cuts, labels = paste("Cluster", 1:6))
# Plot using seqHMM
ssp_k <- list()
for(k in 1:6) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k], 
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}
png(paste0("output/sensitivity/plot_6clusters_", SA, ".png"), width = 900, height = 900)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.3, 0.3, 0.3, 0.1), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()


##### 1: SUBSTITUTION COSTS = 2 #####

SA <- "subs_constant"

# Subs = 2
# Indels = 1
# Method = Optimal Matching

# Create substitution-cost matrices with constant values of "2"
sub_mats_constant <- list(employment = TraMineR::seqsubm(jb_seq, method = "CONSTANT"),
                           health = TraMineR::seqsubm(srh_seq, method = "CONSTANT")) 

# Build multi-domain sequence 
channels <- list(employment = jb_seq, health = srh_seq)
mc_dist <- TraMineR::seqMD(channels, 
                           sm = sub_mats_constant,        # subs = substitution matrix based on transition rates for both channels 
                           indel = 1,                     # indels = 1
                           what = "diss", 
                           method = "OM",                 # optimal matching algorithm
                           cweight = c(1, 1))             # equal weighting of both channels 


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
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()

# 5 clusters 
Cuts <- cutree(Clustered, k = 5)
Groups <- factor(Cuts, labels = paste("Cluster", 1:5))
# Plot using seqHMM
ssp_k <- list()
for(k in 1:5) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k], 
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}
png(paste0("output/sensitivity/plot_5clusters_", SA, ".png"), width = 900, height = 900)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.3, 0.3, 0.3, 0.1), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()


##### 2: METHOD = DYNAMIC HAMMING DISTANCE #####

SA <- "DHD_subs_basecase"

# Subs = Matrices based on transition rates
# Indels = 1
# Method = DYNAMIC HAMMING DISTANCE (might just be Hamming Distance - need to check through this)

# Build multi-domain sequence 
mc_dist <- TraMineR::seqMD(channels, 
                           sm = sub_mats_basecase,        # subs = substitution matrix based on transition rates for both channels 
                           indel = 1,                     # indels = 1
                           what = "diss", 
                           method = "HAM",                # Dynamic Hamming Distance - I think this is just Hamming distance not Dynamic Hamming distance!!!
                           cyouweight = c(1, 1))             # equal weighting of both channels 

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
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()

# 5 clusters 
Cuts <- cutree(Clustered, k = 5)
Groups <- factor(Cuts, labels = paste("Cluster", 1:5))
# Plot using seqHMM
ssp_k <- list()
for(k in 1:5) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k], 
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}
png(paste0("output/sensitivity/plot_5clusters_", SA, ".png"), width = 900, height = 900)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.3, 0.3, 0.3, 0.1), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()


# 6 clusters 
Cuts <- cutree(Clustered, k = 6)
Groups <- factor(Cuts, labels = paste("Cluster", 1:6))
# Plot using seqHMM
ssp_k <- list()
for(k in 1:6) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(Groups)[k], 
                            ylab = c("Economic activity", "Health"),
                            xlab = "Year since onset of long-term sickness")
}
png(paste0("output/sensitivity/plot_6clusters_", SA, ".png"), width = 900, height = 900)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.3, 0.3, 0.3, 0.1), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()


##### 3: OMSTRAN: TRANSITIONS NOT STATES #####

SA <- "OMstran"

# Subs = Matrices based on transition rates
# Indels = 1
# Method = OMstran

## CAT multidomain distances based on channel specific cost methods
MDdist <- seqMD(channels, method="OM",
                sm = list("TRATE"), 
                indel = 1, 
                what="diss")

## Providing channel specific substitution costs - this is the same as the basecase
smatrix <- list()
smatrix[[1]] <- seqsubm(jb_seq, method="TRATE")
smatrix[[2]] <- seqsubm(srh_seq, method="TRATE")

# Create substitution-cost matrices using transition rates for each of the channels
sub_mats_basecase <- list(employment = TraMineR::seqsubm(jb_seq, method = "TRATE"),
                          health = TraMineR::seqsubm(srh_seq, method = "TRATE")) 

## Retrieving the MD sequences
MDseq <- seqMD(channels)
alphabet(MDseq)

## Retrieving the CAT multidomain substitution costs
CATcost <- seqMD(channels,
                 sm=smatrix, 
                 indel = 1, 
                 cweight=c(1,1), 
                 what="cost")

## OMstran distances between MD sequences
MDdist2 <- seqdist(MDseq, 
                   method="OMstran",
                   sm = CATcost, 
                   indel=attr(CATcost,"indel"), 
                   otto = 0.3)

# Conduct clustering using Ward's agglomerative hierarchical method of clustering
Clustered <- hclust(as.dist(MDdist2), method = "ward.D2")

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
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()


##### COMPARISON OF DIFFERENT MCSA SPECIFICATIONS #####

# Loop to read in the specific cluster statistics from the sensitivity analysis
SAs <- c("basecase", "subs_constant", "DHD_subs_basecase", "OMstran")
for(SA in SAs) { 
  
  cluster_stats <- read.csv(paste0("output/sensitivity/cluster_stats_", SA, ".csv")) %>% 
    select(X, ASW) %>% 
    mutate(X = parse_number(X), 
           SA = paste(SA))
  
  if (SA == "basecase") { 
    cluster_stats_ASW <- cluster_stats
  } else { 
    cluster_stats_ASW <- rbind(cluster_stats_ASW, cluster_stats)
    }
}

# Plot of ASW statistics by sensitivity analysis for different MCSA specifications
asw_plot <- cluster_stats_ASW %>% 
  ggplot(aes(x = X, y = ASW, color = SA)) + 
  geom_line() + 
  geom_point() + 
  theme_minimal() + 
  labs(x = "Number of Clusters", y = "Average Silhouette Width", color = "Sensitivity Analysis") + 
  scale_color_brewer(
    palette = "Set1",
    limits = c("basecase", "subs_constant", "DHD_subs_basecase", "OMstran"),
    labels = c("basecase" = "Main Analysis", 
               "subs_constant" = "Subs = 2", 
               "DHD_subs_basecase" = "DHD", 
               "OMstran" = "OMstran"))
asw_plot
png(paste0("output/sensitivity/asw_mcsa_spec_plot.png"), width = 2500, height = 1400, res = 300)
asw_plot
dev.off()


#### 2. EXTENDED ALPHABET APPROACH ####

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


# Original order: "Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other"
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854") 
# Alphabetical order: "Care", "Employed", "Long-term sickness", "Other", "Retired", "Unemployed"
jb_palette_reordered <- c("#FFD92F", "#8DA0CB", "#FC8D62", "#A6D854", "#E78AC3", "#66C2A5") 

# Lighten/darken the base economic activity status colour
extended_palette <- unlist(lapply(jb_palette_reordered, function(colour) {
  colorspace::lighten(colour, amount = seq(-0.3, 0.3, length.out = 5))
}))

# Apply the palette to the combined sequence
cpal(multi_seq) <- extended_palette
TraMineR::seqlegend(multi_seq, cex = 0.7, ncol = 6)

# Plot for the extended alphabet approach (i.e., categories for each combination)
TraMineR::seqdplot(multi_seq, border = NA, legend.prop = 0.2) # distribution plot
TraMineR::seqIplot(multi_seq, sortv = "from.start") # index plot

# View current labels
alphabet(multi_seq)
new_labels <- c("Care / Excellent", "Care / Very good", "Care / Good", "Care / Fair", "Care / Poor", 
                "Employed / Excellent", "Employed / Very good", "Employed / Good", "Employed / Fair", "Employed / Poor",
                "Long-term sickness / Excellent", "Long-term sickness / Very good", "Long-term sickness / Good", "Long-term sickness / Fair", "Long-term sickness / Poor",
                "Other / Excellent", "Other / Very good", "Other / Good", "Other / Fair", "Other / Poor",
                "Retired / Excellent", "Retired / Very good", "Retired / Good", "Retired / Fair", "Retired / Poor",
                "Unemployed / Excellent", "Unemployed / Very good", "Unemployed / Good", "Unemployed / Fair", "Unemployed / Poor")

# Re-assign labels
attr(multi_seq, "labels") <- new_labels
alphabet(multi_seq)


# Plotting with a condensed legend
png(paste0("output/sensitivity/overall_dist_ea.png"), width = 6000, height = 5000, res=500)
TraMineR::seqdplot(multi_seq, 
                   border = NA, 
                   with.legend = "right", 
                   legend.prop = 0.3, 
                   xlab = "Year since onset of long-term sickness") # Give the legend more room
dev.off()

# Compute transition matrix and save
tr_ea_matrix <- seqtrate(multi_seq)
write.csv(x = round(tr_ea_matrix, 3), file = paste0("output/sensitivity/transition_matrix_ea.csv"))

# Calculate distances using optimal matching
dist_om <- seqdist(multi_seq, 
                   method = "OM", 
                   sm = "TRATE") # substition matrix based on transitions


# Conduct clustering using Ward's agglomerative hierarchical method of clustering
Clustered <- hclust(as.dist(dist_om), method = "ward.D2")

# Estimate cluster quality statistics using those included in the WeightedCluster package
Clustered_range <- WeightedCluster::as.clustrange(Clustered, diss = dist_om, ncluster = 10)

# Table of cluster quality statistics
cluster_stats <- Clustered_range[["stats"]]
write.csv(x = cluster_stats, file = "output/sensitivity/cluster_stats_extended_alphabet.csv")

# Plot of cluster quality statistics
plot(Clustered_range, stat = "all", norm = "zscoremed", lwd = 2) # all measures on a plot
png("output/sensitivity/cluster_stats_extended_alphabet.png", width = 2400, height = 1600, res = 300)
plot(Clustered_range, stat = c("PBC", "HG", "ASW", "HC", "CH"), norm = "zscoremed", lwd = 2) # selected measures
dev.off()

# 4 clusters
Cuts <- cutree(Clustered, k = 4)
Groups <- factor(Cuts, labels = paste("Cluster", 1:4))


# Reorder the clusters to be more informative and logical with cluster names
new_labels <- c("Cluster 4", "Cluster 2", "Cluster 3", "Cluster 1") 
new_Groups <- factor(Groups, levels = paste("Cluster", 1:4), labels = new_labels)
new_Cuts <- as.numeric(gsub("Cluster", "", new_Groups))
new_Groups_labelled <- factor(new_Cuts, labels = c("Cluster 1: Entering employment", "Cluster 2: Healthy but inactive", "Cluster 3: Continued long-term sickness", "Cluster 4: Transition to retirement"))


# Plot using the extended alphabet approach
TraMineR::seqdplot(multi_seq, 
                   border = NA, 
                   group = new_Groups_labelled,
                   with.legend = "auto", 
                   legend.prop = 0.1, # Give the legend more room
                   cex.legend =0.8,
                   ncol = 6) 

# Save publication-quality plot
png(paste0("output/sensitivity/plot_4clusters_extended_alphabet_combined_states_publication.png"), width = 5000, height = 5000, res=500)
TraMineR::seqdplot(multi_seq, 
                   border = NA, 
                   group = new_Groups_labelled,
                   with.legend = "auto", 
                   legend.prop = 0.1, # Give the legend more room
                   cex.legend = 0.85,
                   ncol = 6, 
                   xlab = "Year since onset of long-term sickness") 
dev.off()


# Individual plots for each cluster using the same scheme as the main analysis
ssp_k <- list()
for(k in 1:4) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,], 
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA, 
                            with.legend = FALSE, 
                            title = levels(new_Groups_labelled)[k], 
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}

# Combine plots together
gridplot(ssp_k, 
         ncol = 2, byrow = TRUE, 
         row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
         title.legend = "auto" 
)


# Save publication-quality plot
png(paste0("output/sensitivity/plot_4clusters_extended_alphabet_publication.png"), width = 4500, height = 5500, res=500)
gridplot(ssp_k, 
         ncol = 2, byrow = TRUE, 
         row.prop = c(0.41, 0.41, 0.18), 
         cex.legend = 1, 
         title.legend = "auto")
dev.off()


# Comparison of cluster membership between main analysis and extended alphabet approach

# Filter to baseline sample 
sample_baseline_ea <- filter(sample, transition_time == 0)
# Add on the cluster numbers 
sample_baseline_ea$cluster_ea <- new_Cuts
sample_baseline_ea$group_ea <- new_Groups

# Read in data with the cluster numbers from the main basecase analysis
sample_baseline_main <- readRDS("data/processed/regression_sample.rds")
# Filter the main sample to those in the EA only sample
sample_baseline_main <- sample_baseline_main %>% 
  dplyr::rename(cluster_main = cluster, group_main = group)

# Join together and make a little table
cluster_comp <- full_join(sample_baseline_main, sample_baseline_ea)
cluster_comp_table <- table(cluster_comp$cluster_ea, cluster_comp$cluster_main)
write.csv(x = cluster_comp_table, file = paste0("output/sensitivity/cluster_comp_extended_alphabet_main_table.csv"))



#### 3. SUBSAMPLE: ECONOMICALLY ACTIVE PRIOR ONLY ####

# Filter to a subsample of only those who are economically active in the year prior to the onset of long-term sickness 
# Economically active = Employed or Unemployed

SA <- "ea_only"

sample_prior_ea_only <- sample %>% 
  filter(transition_time == -1, jb %in% c("employed", "unemployed"))
sample_ea_only <- sample %>% filter(pidp %in% unique(sample_prior_ea_only$pidp))
length(unique(sample_prior_ea_only$pidp))

# Create economic activity (jb) sequence dataframe
jb_seq_df_ea_only <- sample_ea_only %>% select(pidp, jb, transition_time) %>% 
  filter(transition_time >= 0) %>% 
  arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "jb")

jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other")
jb_labels <- c("Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other")
brewer.pal(n=6, "Set2")
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854") # select and reorder colours to match jb_levels

jb_seq <- seqdef(jb_seq_df_ea_only, 2:6, alphabet = jb_levels, labels = jb_labels, cpal = jb_palette)

# Create self-rated health (srh) sequence dataframe
srh_seq_df_ea_only <- sample_ea_only %>% select(pidp, srh, transition_time) %>% 
  filter(transition_time >= 0) %>% 
  arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "srh")

srh_levels <- c("excellent", "very good", "good", "fair", "poor")
srh_labels <- c("Excellent", "Very good", "Good", "Fair", "Poor")
brewer.pal(n=5, "RdYlBu")
srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C") # reorder colours to match srh_levels

srh_seq <- seqdef(srh_seq_df_ea_only, 2:6, alphabet = srh_levels, labels = srh_labels, cpal = srh_palette)

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
png(paste0("output/sensitivity/plot_4clusters_", SA, ".png"), width = 900, height = 700)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE, 
                 row.prop = c(0.4, 0.4, 0.2), cex.legend = 0.8, 
                 title.legend = "auto")
dev.off()

# 6 clusters
Cuts <- cutree(Clustered, k = 6)
Groups <- factor(Cuts, labels = paste("Cluster", 1:6))
# Plot using seqHMM
ssp_k <- list()
for(k in 1:6) {
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[Cuts == k,],
                                 srh_seq[Cuts == k,]),
                            type = "d", border = NA,
                            with.legend = FALSE,
                            title = levels(Groups)[k],
                            ylab = c("Economic activity", "Self-rated health"),
                            xlab = "Year since onset of long-term sickness")
}
png(paste0("output/sensitivity/plot_6clusters_", SA, ".png"), width = 900, height = 900)
seqHMM::gridplot(ssp_k, ncol = 2, byrow = TRUE,
                 row.prop = c(0.3, 0.3, 0.3, 0.1), cex.legend = 0.8,
                 title.legend = "auto")
dev.off()

## Compare cluster membership of EA only with the main analysis ##

# Relabel for consistency with the main analysis and the ordering of clusters already defined
# 4 clusters 
Cuts <- cutree(Clustered, k = 4)
Groups <- factor(Cuts, labels = paste("Cluster", 1:4))
new_labels <- c("Cluster 4", "Cluster 2", "Cluster 3", "Cluster 1") 
new_Groups <- factor(Groups, levels = paste("Cluster", 1:4), labels = new_labels)
new_Cuts <- as.numeric(gsub("Cluster", "", new_Groups))

# Filter to baseline sample 
sample_baseline_ea_only <- filter(sample_ea_only, transition_time == 0)
# Add on the cluster numbers 
sample_baseline_ea_only$cluster_ea_only <- new_Cuts
sample_baseline_ea_only$group_ea_only <- new_Groups

# Read in data with the cluster numbers from the main basecase analysis
sample_baseline_main <- readRDS("data/processed/regression_sample.rds")
# Filter the main sample to those in the EA only sample
sample_baseline_main <- sample_baseline_main %>% 
  filter(pidp %in% unique(sample_baseline_ea_only$pidp)) %>% 
  rename(cluster_main = cluster, group_main = group)

# Join together and make a little table
cluster_comp <- full_join(sample_baseline_main, sample_baseline_ea_only)
cluster_comp_table <- table(cluster_comp$cluster_ea_only, cluster_comp$cluster_main)
write.csv(x = cluster_comp_table, file = paste0("output/sensitivity/cluster_comp_ea_only_main_table.csv"))

### Appendix Figure ??: Sequence distribution plots for each cluster for economically active only sample ###

# Reorder the clusters to be more informative and logical with cluster names

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
         title.legend = "auto")
dev.off()
