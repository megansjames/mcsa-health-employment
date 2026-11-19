
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
Clustered_range <- WeightedCluster::as.clustrange(Clustered, diss = MDdist2, ncluster = 10) # fix this 

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

##### Appendix Figure: Plot of ASW statistics by sensitivity analysis for different MCSA specifications #####
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


##### Appendix Figure: Overall distribution plot for extended alphabet approach #####
png(paste0("output/sensitivity/overall_dist_ea.png"), width = 6000, height = 5000, res=500)
TraMineR::seqdplot(multi_seq, 
                   border = NA, 
                   with.legend = "right", 
                   legend.prop = 0.3, 
                   xlab = "Year since onset of long-term sickness") # Give the legend more room
dev.off()

##### Appendix Table: Transition matrix for extended alphabet approach #####
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
new_Groups_labelled <- factor(new_Cuts, labels = c("Cluster 1: Entering employment", "Cluster 2: Precarious non-employment", "Cluster 3: Continued long-term sickness", "Cluster 4: Rapid retirement"))


# Plot using the extended alphabet approach
TraMineR::seqdplot(multi_seq, 
                   border = NA, 
                   group = new_Groups_labelled,
                   with.legend = "auto", 
                   legend.prop = 0.1, # Give the legend more room
                   cex.legend =0.8,
                   ncol = 6) 

##### Appendix Figure: Distribution plot by cluster for extended alphabet approach #####
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
  ssp_k[[k]] <- seqHMM::ssp(list(jb_seq[new_Cuts == k,], 
                                 srh_seq[new_Cuts == k,]),
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
  dplyr::rename(cluster_main = cluster, group_main = group)

# Join together and make a little table
cluster_comp <- full_join(sample_baseline_main, sample_baseline_ea_only)
cluster_comp_table <- table(cluster_comp$cluster_ea_only, cluster_comp$cluster_main)
write.csv(x = cluster_comp_table, file = paste0("output/sensitivity/cluster_comp_ea_only_main_table.csv"))

##### Appendix Figure: Sequence distribution plots for each cluster for economically active only sample #####

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

# Cluster 2: Precarious non-employment
ssp_k[[2]] <- ssp(list(jb_seq[new_Cuts == 2,], 
                       srh_seq[new_Cuts == 2,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 2: Precarious non-employment"), 
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

# Cluster 4: Rapid retirement
ssp_k[[4]] <- ssp(list(jb_seq[new_Cuts == 4,], 
                       srh_seq[new_Cuts == 4,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 4: Rapid retirement"), 
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


#### 4. SUBSAMPLE: EARLIER YEARS ONLY (PRE-COVID) ####

# Filter to a subsample of only those who experience their transitions in the earlier waves
# This allows for all observations to be pre-COVID

SA <- "pre_covid"

sample_pre_covid <- sample %>% 
  filter(lts_first_wave <= 6)
length(unique(sample_pre_covid$pidp))

# Create economic activity (jb) sequence dataframe
jb_seq_df_pre_covid <- sample_pre_covid %>% select(pidp, jb, transition_time) %>% 
  filter(transition_time >= 0) %>% 
  arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "jb")

jb_levels <- c("employed", "unemployed", "long-term sick", "retired", "care", "other")
jb_labels <- c("Employed", "Unemployed", "Long-term sickness", "Retired", "Care", "Other")
brewer.pal(n=6, "Set2")
jb_palette <- c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3", "#FFD92F", "#A6D854") # select and reorder colours to match jb_levels

jb_seq <- seqdef(jb_seq_df_pre_covid, 2:6, alphabet = jb_levels, labels = jb_labels, cpal = jb_palette)

# Create self-rated health (srh) sequence dataframe
srh_seq_df_pre_covid <- sample_pre_covid %>% select(pidp, srh, transition_time) %>% 
  filter(transition_time >= 0) %>% 
  arrange(pidp, transition_time) %>% 
  pivot_wider(id_cols = "pidp", names_from = "transition_time", values_from = "srh")

srh_levels <- c("excellent", "very good", "good", "fair", "poor")
srh_labels <- c("Excellent", "Very good", "Good", "Fair", "Poor")
brewer.pal(n=5, "RdYlBu")
srh_palette <- c("#2C7BB6", "#ABD9E9", "#FFFFBF",  "#FDAE61", "#D7191C") # reorder colours to match srh_levels

srh_seq <- seqdef(srh_seq_df_pre_covid, 2:6, alphabet = srh_levels, labels = srh_labels, cpal = srh_palette)

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

## Compare cluster membership of EA only with the main analysis ##

# Relabel for consistency with the main analysis and the ordering of clusters already defined
# 4 clusters 
Cuts <- cutree(Clustered, k = 4)
Groups <- factor(Cuts, labels = paste("Cluster", 1:4))
new_labels <- c("Cluster 4", "Cluster 2", "Cluster 3", "Cluster 1") 
new_Groups <- factor(Groups, levels = paste("Cluster", 1:4), labels = new_labels)
new_Cuts <- as.numeric(gsub("Cluster", "", new_Groups))

# Filter to baseline sample 
sample_baseline_pre_covid <- filter(sample_pre_covid, transition_time == 0)
# Add on the cluster numbers 
sample_baseline_pre_covid$cluster_pre_covid <- new_Cuts
sample_baseline_pre_covid$group_pre_covid <- new_Groups

# Read in data with the cluster numbers from the main basecase analysis
sample_baseline_main <- readRDS("data/processed/regression_sample.rds")
# Filter the main sample to those in the pre-COVID sample
sample_baseline_main <- sample_baseline_main %>% 
  filter(pidp %in% unique(sample_baseline_pre_covid$pidp)) %>% 
  dplyr::rename(cluster_main = cluster, group_main = group)

# Join together and make a little table
cluster_comp <- full_join(sample_baseline_main, sample_baseline_pre_covid)
cluster_comp_table <- table("precovid" = cluster_comp$cluster_pre_covid, "main" = cluster_comp$cluster_main)
write.csv(x = cluster_comp_table, file = paste0("output/sensitivity/cluster_comp_", SA, "_main_table.csv"))
pre_covid_pcts <- prop.table(table(cluster_comp$cluster_pre_covid)) * 100
round(pre_covid_pcts, 1)


##### Appendix Figure: Sequence distribution plots for each cluster for pre-COVID sample #####

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

# Cluster 2: Precarious non-employment
ssp_k[[2]] <- ssp(list(jb_seq[new_Cuts == 2,], 
                       srh_seq[new_Cuts == 2,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 2: Precarious non-employment"), 
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

# Cluster 4: Rapid retirement
ssp_k[[4]] <- ssp(list(jb_seq[new_Cuts == 4,], 
                       srh_seq[new_Cuts == 4,]),
                  type = "d", border = NA, 
                  with.legend = FALSE, 
                  title = paste0("Cluster 4: Rapid retirement"), 
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


#### 5. BOOTSTRAPPED CLUSTER ALLOCATION ####

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

# Cluster stability analysis using clusterboot() 

# Select your preferred number of clusters k based on your cluster_stats (e.g., k = 4)
k_target <- 4
set.seed(123)
boot_mc <- fpc::clusterboot(
  data = mc_dist,                # Multi-channel dissimilarity matrix
  B = 100,                       # Number of bootstrap replicates (100–500)
  distances = TRUE,              # Indicates data is a dissimilarity matrix
  clustermethod = disthclustCBI, # Interface for dissimilarity-based hclust
  method = "ward.D2",            # Agglomerative method matching your main analysis
  k = k_target,                  # Target cluster count
  showplots = FALSE)

# Stability metrics
stability_results <- data.frame(
  Cluster = paste0("Cluster_", 1:k_target),
  Jaccard_Mean = boot_mc$bootmean, # Mean Jaccard similarity across bootstrap runs
  Dissolved_Count = boot_mc$bootbrd # Number of times cluster was dissolved (Jaccard < 0.5)
)
print(stability_results)
write.csv(stability_results, file = "output/sensitivity/cluster_stability_jaccard.csv", row.names = FALSE)


# Calculate average stability across all candidate k values (k = 2 to 10)
stability_by_k <- numeric(9)
names(stability_by_k) <- paste0("k_", 2:10)

for (k_i in 2:10) {
  set.seed(123)
  boot_k <- fpc::clusterboot(
    data = mc_dist, distances = TRUE, B = 100,
    clustermethod = disthclustCBI, method = "ward.D2",
    k = k_i, showplots = FALSE
  )
  stability_by_k[k_i - 1] <- mean(boot_k$bootmean)
}

# Append overall mean stability column to your cluster quality table
cluster_stats_ext <- cbind(cluster_stats, Overall_Mean_Jaccard = stability_by_k)
write.csv(cluster_stats_ext, file = "output/sensitivity/cluster_stats_with_stability.csv")



#### 6. CROSS-VALIDATION OF PREDICTION MODELS ####

# Set seed for reproducibility
set.seed(123)

# Read in the regression sample with cluster numbers
data <- readRDS("data/processed/regression_sample.rds")

# Re-level economic activity to remove long-term sickness (impossible for this to be a prior economic activity due to how the sample is defined)
levels(data$prior_jb)
data$jb <- factor(as.numeric(data$prior_jb), 
                  levels = c(1,2,4,5,6), 
                  labels = c("employed", "unemployed", "retired", "care", "other"))  
levels(data$prior_jb)

# Set reference categories for all other variables to be included in the regression
data <- data %>% 
  mutate(agegrp_10yr = fct_relevel(agegrp_10yr, "55-64"), 
         hiqual_3cat = fct_relevel(hiqual_3cat, "high", "medium", "low"), 
         hiqual_4cat = fct_relevel(hiqual_4cat, "high", "medium", "low", "no"), 
         gor = fct_relevel(gor, "London"),
         srh = fct_relevel(srh, "good"), # set good as reference category
         prior_srh = fct_relevel(prior_srh, "good"), # good as reference category
         prior_jb = fct_relevel(jb, "employed"), 
         prior_equivquint = fct_relevel(prior_equivquint, "Q3")) # set Q3 as the reference category

# Define the model formulas for each of the models

# Model 1: Baseline only
f_baseline <- group ~ agegrp_10yr + sex + hiqual_4cat + marstat_2cat + gor + 
  as.factor(intdaty_dv) + srh + ghqcase
# Model 2: Baseline and prior
f_prior <- group ~ agegrp_10yr + sex + hiqual_4cat + marstat_2cat + gor + 
  as.factor(intdaty_dv) + srh + ghqcase + 
  prior_srh + prior_ghqcase + prior_jb + prior_equivquint # prior characteristics
# Model 3: Prior only 
f_prior_only <- group ~ agegrp_10yr + sex + hiqual_4cat + marstat_2cat + gor + 
  as.factor(intdaty_dv) +
  prior_srh + prior_ghqcase + prior_jb + prior_equivquint # prior characteristics

# Create 10-fold cross-validation repeated 5 times (stratified by outcome)
cv_folds <- rsample::vfold_cv(data, v = 10, repeats = 5, strata = group)

# Start counter for fitting of models to each of the folds
progress_env <- new.env()
progress_env$count <- 0
progress_env$total <- nrow(cv_folds)

# Function to fit brmultinom and evaluate out-of-sample metrics per fold
eval_brmultinom <- function(split, formula_obj) {
  
  # Increment counter and print current status
  progress_env$count <- progress_env$count + 1
  cat(sprintf("[%d/%d] Fitting brmultinom fold...\n", progress_env$count, progress_env$total))
  
  train <- rsample::analysis(split)
  test  <- rsample::assessment(split)
  
  # Fit bias-reduced multinomial model on training set
  fit <- brglm2::brmultinom(formula_obj, 
                            data = train, 
                            type = "AS_mean", 
                            slowit = 0.5,   # Slow down steps to guarantee convergence
                            maxit = 500,    # Allow up to 500 iterations
                            epsilon = 1e-5)  # Slightly relax tolerance if necessary (default 1e-6))
  
  # Predict classes and probabilities on test set
  pred_class <- predict(fit, newdata = test, type = "class")
  pred_probs <- predict(fit, newdata = test, type = "probs")
  
  # Out-of-sample accuracy
  acc <- mean(pred_class == test$group, na.rm = TRUE)
  
  # Out-of-sample Multiclass ROC AUC (using Hand-Till)
  mc_auc <- as.numeric(
    pROC::multiclass.roc(response = test$group, 
                         predictor = pred_probs, 
                         levels = levels(test$group))$auc)
  
  tibble::tibble(accuracy = acc, roc_auc = mc_auc)
}

# Run Cross-Validation for Model 1 (Baseline)
cv_baseline_br <- cv_folds %>%
  mutate(metrics = map(splits, ~eval_brmultinom(.x, f_baseline))) %>%
  unnest(metrics) %>%
  mutate(model = "Model 1: Baseline (BR)")
saveRDS(cv_baseline_br, file = "output/sensitivity/cv_baseline_br.rds")

# Run Cross-Validation for Model 2 (Baseline + Prior)
progress_env$count <- 0 # re-set counter
cv_prior_br <- cv_folds %>%
  mutate(metrics = map(splits, ~eval_brmultinom(.x, f_prior))) %>%
  unnest(metrics) %>%
  mutate(model = "Model 2: Baseline + Prior (BR)")
saveRDS(cv_prior_br, file = "output/sensitivity/cv_prior_br.rds")

# Run Cross-Validation for Model 3 (Baseline + Prior)
progress_env$count <- 0 # re-set counter
cv_prior_only_br <- cv_folds %>%
  mutate(metrics = map(splits, ~eval_brmultinom(.x, f_prior_only))) %>%
  unnest(metrics) %>%
  mutate(model = "Model 3: Prior Only (BR)")
saveRDS(cv_prior_only_br, file = "output/sensitivity/cv_prior_only_br.rds")


##### Appendix Table: Summary of cross-validation of models #####

# Nadeau-Bengio corrected resampled t-test (Nadeau & Bengio 2003)
J     <- nrow(cv_baseline_br)    # 50
ratio <- 1/9                     # test size / training size for 10-fold
tcrit <- qt(0.975, df = J - 1)

nb_test <- function(d) {
  m  <- mean(d)
  se <- sd(d) * sqrt(1/J + ratio)
  c(mean = m, lo = m - tcrit * se, hi = m + tcrit * se,
    p = 2 * pt(-abs(m / se), df = J - 1))
}

# Calculate out-of-sample estimates with 95% CIs
oos_summary_br <- bind_rows(cv_baseline_br, cv_prior_br, cv_prior_only_br) %>%
  pivot_longer(cols = c(accuracy, roc_auc), names_to = "metric", values_to = "estimate") %>%
  dplyr::group_by(model, metric) %>%
  dplyr::summarize(
    n = n(), 
    Mean = mean(estimate, na.rm = TRUE),
    SE = sd(estimate, na.rm = TRUE) * sqrt(1/J + ratio),
    CI_Lower = Mean - tcrit * SE,
    CI_Upper = Mean + tcrit * SE,
    .groups  = "drop")

print(oos_summary_br)
write.csv(oos_summary_br, file = "output/sensitivity/oos_summary_br_nb.csv")

# Join fold-by-fold results to test difference between baseline and baseline + prior models
diff_br_2v1 <- cv_prior_br %>%
  select(id, id2, acc_prior = accuracy, auc_prior = roc_auc) %>%
  inner_join(cv_baseline_br %>% select(id, id2, acc_base = accuracy, auc_base = roc_auc),
             by = c("id", "id2")) %>%
  mutate(acc_diff = acc_prior - acc_base,
         auc_diff = auc_prior - auc_base)

# Paired t-tests across resamples
r_acc <- nb_test(diff_br_2v1$acc_diff)
r_auc <- nb_test(diff_br_2v1$auc_diff)

# Summary of differences between out-of-sample estimates by model
oos_diff_summary_2v1 <- tibble::tibble(
  Metric    = c("Accuracy", "Multiclass AUC"),
  Mean_Diff = c(r_acc["mean"], r_auc["mean"]),
  CI_Lower  = c(r_acc["lo"],   r_auc["lo"]),
  CI_Upper  = c(r_acc["hi"],   r_auc["hi"]),
  p_value   = c(r_acc["p"],    r_auc["p"]))

print(oos_diff_summary_2v1)
write.csv(oos_diff_summary_2v1, file = "output/sensitivity/oos_diff_summary_2v1_nb.csv")

# Join fold-by-fold results to test difference between prior only and baseline + prior models
diff_br_2v3 <- cv_prior_br %>%
  select(id, id2, acc_prior = accuracy, auc_prior = roc_auc) %>%
  inner_join(cv_prior_only_br %>% select(id, id2, acc_prior_only = accuracy, auc_prior_only = roc_auc),
             by = c("id", "id2")) %>%
  mutate(acc_diff = acc_prior - acc_prior_only,
         auc_diff = auc_prior - auc_prior_only)

# Paired t-tests across resamples
r_acc <- nb_test(diff_br_2v3$acc_diff)
r_auc <- nb_test(diff_br_2v3$auc_diff)

# Summary of differences between out-of-sample estimates by model
oos_diff_summary_2v3 <- tibble::tibble(
  Metric    = c("Accuracy", "Multiclass AUC"),
  Mean_Diff = c(r_acc["mean"], r_auc["mean"]),
  CI_Lower  = c(r_acc["lo"],   r_auc["lo"]),
  CI_Upper  = c(r_acc["hi"],   r_auc["hi"]),
  p_value   = c(r_acc["p"],    r_auc["p"]))

print(oos_diff_summary_2v3)
write.csv(oos_diff_summary_2v3, file = "output/sensitivity/oos_diff_summary_2v3_nb.csv")

# Join fold-by-fold results to test difference between prior only and baseline only models
diff_br_3v1 <- cv_prior_only_br %>%
  select(id, id2, acc_prior_only = accuracy, auc_prior_only = roc_auc) %>%
  inner_join(cv_baseline_br %>% select(id, id2, acc_base = accuracy, auc_base = roc_auc),
             by = c("id", "id2")) %>%
  mutate(acc_diff = acc_prior_only - acc_base,
         auc_diff = auc_prior_only - auc_base)

# Paired t-tests across resamples
r_acc <- nb_test(diff_br_3v1$acc_diff)
r_auc <- nb_test(diff_br_3v1$auc_diff)

# Summary of differences between out-of-sample estimates by model
oos_diff_summary_3v1 <- tibble::tibble(
  Metric    = c("Accuracy", "Multiclass AUC"),
  Mean_Diff = c(r_acc["mean"], r_auc["mean"]),
  CI_Lower  = c(r_acc["lo"],   r_auc["lo"]),
  CI_Upper  = c(r_acc["hi"],   r_auc["hi"]),
  p_value   = c(r_acc["p"],    r_auc["p"]))

print(oos_diff_summary_3v1)
write.csv(oos_diff_summary_3v1, file = "output/sensitivity/oos_diff_summary_3v1_nb.csv")


##### Additional Plots: Sina Plots of Cross-Validation Folds ##### 

## Accuracy ##
cv_accuracy_plot <- bind_rows(cv_baseline_br, cv_prior_br, cv_prior_only_br) %>%
  pivot_longer(cols = c(accuracy, roc_auc), names_to = "metric", values_to = "estimate") %>%  
  filter(metric == "accuracy") %>% 
  ggplot(aes(x = model, y = estimate)) + 
  geom_violin(fill = "grey95", color = "grey50", trim = FALSE, alpha = 0.6) +  # neutral violin background
  geom_jitter(aes(color = id), width = 0.04, height = 0, size = 2, alpha = 0.8) + # coloured by repeat
  labs(
    x = "", 
    y = "Cross-Validated Accuracy", 
    color = "CV Repeat") +
  scale_x_discrete(labels = c("Model 1: Baseline Only", "Model 2: Baseline + One-Year Prior", "Model 3: Prior Only")) + 
  theme_minimal()

png("output/sensitivity/cv_accuracy_plot.png", width = 2000, height = 1400, res = 300)
cv_accuracy_plot
dev.off()

## ROC AUC ##
cv_auc_plot <- bind_rows(cv_baseline_br, cv_prior_br, cv_prior_only_br) %>%
  pivot_longer(cols = c(accuracy, roc_auc), names_to = "metric", values_to = "estimate") %>%  
  filter(metric == "roc_auc") %>% 
  ggplot(aes(x = model, y = estimate)) + 
  geom_violin(fill = "grey95", color = "grey50", trim = FALSE, alpha = 0.6) +  # neutral violin background
  geom_jitter(aes(color = id), width = 0.04, height = 0, size = 2, alpha = 0.8) + # coloured by repeat
  labs(
    x = "", 
    y = "Cross-Validated Overall Area Under Curve (AUC)", 
    color = "CV Repeat") +
  scale_x_discrete(labels = c("Model 1: Baseline Only", "Model 2: Baseline + One-Year Prior", "Model 3: Prior Only")) + 
  theme_minimal()

png("output/sensitivity/cv_auc_plot.png", width = 2000, height = 1400, res = 300)
cv_accuracy_plot
dev.off()

