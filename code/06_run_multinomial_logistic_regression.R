
## RUN MULTINOMIAL LOGISTIC REGRESSION ##

rm(list=ls())


#### 1. SET UP ####

# Read in data with the cluster numbers and prior information
data <- readRDS("data/processed/regression_sample.rds")

# Re-level clusters so that Cluster 1 will be the baseline reference group
levels(data$group)
data$group <- factor(data$group, 
                     levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))

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


#### 2. MULTINOMIAL REGRESSION ####

##### CONSTANT ONLY #####

model_constant <- parsnip::multinom_reg(engine = "nnet") %>% 
  fit(group ~ 1,
      data = data)

broom::glance(model_constant)
broom::tidy(model_constant, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate_if(is.numeric, round, 4) %>%  
  select(-std.error, -statistic)

# Add predictions
mod_preds_constant <- model_constant %>% parsnip::augment(new_data = data)

# Confusion matrix 
yardstick::conf_mat(mod_preds_constant, truth = group, estimate = .pred_class)
conf_mat_constant <- as.data.frame(conf_mat(mod_preds_constant, truth = group, estimate = .pred_class)$table) %>% 
  group_by(Prediction) %>% mutate(percent = Freq/sum(Freq))

# Plot for the confusion matrix
conf_mat_constant_plot <- ggplot(data = conf_mat_constant, mapping = aes(x=Prediction, y=Truth)) + 
  geom_tile(aes(fill = percent)) +
  geom_text(aes(label = sprintf("%0.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "#458B00") +
  theme_light() + theme(legend.position = "none") + 
  labs(subtitle = "Constant only")

# Accuracy by cluster
filter(conf_mat_constant, Prediction == Truth)

# Model accuracy
yardstick::accuracy(mod_preds_constant, truth = group, estimate = .pred_class)

# Area under the ROC 
roc_auc_constant <- yardstick::roc_auc(mod_preds_constant, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4')$.estimate
yardstick::roc_curve(mod_preds_constant, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4') |> 
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = .level)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_abline(slope = 1, linetype = "dotted") +
  coord_fixed() +
  labs(color = NULL) +
  theme_light()

##### BASELINE ONLY ######

model_baseline <- parsnip::multinom_reg(engine = "nnet") %>% 
  fit(group ~ agegrp_10yr + sex + hiqual_4cat + marstat_2cat + gor + as.factor(intdaty_dv) + # demographics
        srh + ghqcase, # current health status
      data = data)

broom::glance(model_baseline)
broom::tidy(model_baseline, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate_if(is.numeric, round, 4) %>%  
  select(-std.error, -statistic)

# Add predictions
mod_preds_baseline <- model_baseline %>% parsnip::augment(new_data = data)

# Confusion matrix 
yardstick::conf_mat(mod_preds_baseline, truth = group, estimate = .pred_class)
conf_mat_baseline <- as.data.frame(conf_mat(mod_preds_baseline, truth = group, estimate = .pred_class)$table) %>% 
  group_by(Prediction) %>% mutate(percent = (Freq/sum(Freq)*100))

# Plot for the confusion matrix
conf_mat_baseline_plot <- ggplot(data = conf_mat_baseline, mapping = aes(x=Truth, y=Prediction)) + 
  geom_tile(aes(fill = percent)) +
  geom_text(aes(label = sprintf("%0.0f", Freq)), vjust = 1, size = 6) +
  scale_fill_gradient(low = "white", high = "#458B00") +
  theme_light() + theme(legend.position = "none", 
                        text = element_text(size = 20),
                        plot.title = element_text(size = 16)) + 
  labs(subtitle = "Baseline only")

# Accuracy by cluster
filter(conf_mat_baseline, Prediction == Truth)

# Model accuracy
yardstick::accuracy(mod_preds_baseline, truth = group, estimate = .pred_class)

# Area under the ROC 
roc_auc_baseline <- yardstick::roc_auc(mod_preds_baseline, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4')$.estimate
yardstick::roc_curve(mod_preds_baseline, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4') |> 
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = .level)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_abline(slope = 1, linetype = "dotted") +
  coord_fixed() +
  labs(color = NULL) +
  theme_light() + 
  scale_color_discrete(palette = c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3"))

##### BASELINE + ONE YEAR PRIOR #####

model_prior <- parsnip::multinom_reg(engine = "nnet") %>% 
  fit(group ~ agegrp_10yr + sex + hiqual_4cat + marstat_2cat + gor + as.factor(intdaty_dv) + # demographics
        srh + ghqcase + # current health status
        prior_srh + prior_ghqcase +  # prior health 
        prior_jb + prior_equivquint, # prior economic activity and income 
      data = data)

broom::glance(model_prior)
broom::tidy(model_prior, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate_if(is.numeric, round, 4) %>%  
  select(-std.error, -statistic)

# Add predictions
mod_preds_prior <- model_prior %>% parsnip::augment(new_data = data)

# Confusion matrix 
yardstick::conf_mat(mod_preds_prior, truth = group, estimate = .pred_class)
conf_mat_prior <- as.data.frame(conf_mat(mod_preds_prior, truth = group, estimate = .pred_class)$table) %>% 
  group_by(Prediction) %>% mutate(percent = Freq/sum(Freq))

# Plot for the confusion matrix
conf_mat_prior_plot <- ggplot(data = conf_mat_prior, mapping = aes(x=Truth, y=Prediction)) + 
  geom_tile(aes(fill = percent)) +
  geom_text(aes(label = sprintf("%0.0f", Freq)), vjust = 1, size = 6) +
  scale_fill_gradient(low = "white", high = "#458B00") +
  theme_light() + theme(legend.position = "none", 
                        text = element_text(size = 20), 
                        plot.title = element_text(size = 16)) + 
  labs(subtitle = "Baseline and one-year prior")

# Accuracy by cluster
filter(conf_mat_prior, Prediction == Truth)

# Model accuracy
yardstick::accuracy(mod_preds_prior, truth = group, estimate = .pred_class)

# Area under the ROC 
roc_auc_prior <- yardstick::roc_auc(mod_preds_prior, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4')$.estimate
yardstick::roc_curve(mod_preds_prior, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4') %>%  
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = .level)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_abline(slope = 1, linetype = "dotted") +
  coord_fixed() +
  labs(color = NULL) +
  theme_light() + 
  scale_color_discrete(palette = c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3"))

#### 3. COMPARISION OF REGRESSION MODELS ####

# Save Confusion Matrix plots 
png("output/main/conf_mat_baseline.png", width = 1950, height = 1950, res = 300)
conf_mat_baseline_plot
dev.off()

png("output/main/conf_mat_prior.png", width = 1950, height = 1950, res = 300)
conf_mat_prior_plot
dev.off()

png("output/main/conf_mat_both.png", width = 4000, height = 1900, res = 300)
conf_mat_baseline_plot + conf_mat_prior_plot & patchwork::plot_annotation(tag_levels = '1')
dev.off()

##### Figure: Confusion matrices for baseline only and baseline and one-year prior characteristics ######

# Use the cvms package for much more detailed confusion matrices #

## Baseline ##

# Change labels by modifying factor levels to include cluster names 
mod_preds_baseline$group <- factor(mod_preds_baseline$group, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), 
                                   labels = c("1: Entering \nemployment", "2: Precarious \nnon-employment", "3: Continued \nlong-term sickness", "4: Rapid \nretirement"))
mod_preds_baseline$`.pred_class` <- factor(mod_preds_baseline$`.pred_class`, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), 
                                           labels = c("1: Entering \nemployment", "2: Precarious \nnon-employment", "3: Continued \nlong-term sickness", "4: Rapid \nretirement"))

# Create confusion matrix
conf_mat_baseline <- cvms::confusion_matrix(targets = mod_preds_baseline$group, 
                                            predictions = mod_preds_baseline$`.pred_class`)
plot_baseline <- cvms::plot_confusion_matrix(
  conf_mat_baseline$`Confusion Matrix`[[1]],
  add_sums = TRUE,
  add_normalized = FALSE, 
  palette = list(low = "#ffffff", high = "#000000"), 
  font_counts = font(size = 4.5),
  font_row_percentages = font(size = 3.5), 
  font_col_percentages = font(size = 3.5),
  place_x_axis_above = FALSE,
  dynamic_font_colors = dynamic_font_color_settings(
    threshold = 130,
    by = "counts",
    all = c("#000000", "#ffffff")),
  sums_settings = sum_tile_settings(
    label = "Total")) +
  labs(x = "Actual Cluster", y = "Predicted Cluster", title = "Model 1", subtitle = "Baseline characteristics only") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 1))


## Baseline and one-year prior ##

# Change labels by modifying factor levels to include cluster names
mod_preds_prior$group <- factor(mod_preds_prior$group, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), 
                                labels = c("1: Entering \nemployment", "2: Precarious \nnon-employment", "3: Continued \nlong-term sickness", "4: Rapid \nretirement"))
mod_preds_prior$`.pred_class` <- factor(mod_preds_prior$`.pred_class`, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"), 
                                        labels = c("1: Entering \nemployment", "2: Precarious \nnon-employment", "3: Continued \nlong-term sickness", "4: Rapid \nretirement"))

# Create confusion matrix
confusion_matrix_prior <- cvms::confusion_matrix(targets = mod_preds_prior$group, 
                                                 predictions = mod_preds_prior$`.pred_class`)
plot_prior <- cvms::plot_confusion_matrix(
  confusion_matrix_prior$`Confusion Matrix`[[1]],
  add_sums = TRUE,
  add_normalized = FALSE, 
  palette = list(low = "#ffffff", high = "#000000"), 
  font_counts = font(size = 4.5),
  font_row_percentages = font(size = 3.5), 
  font_col_percentages = font(size = 3.5),
  place_x_axis_above = FALSE,
  dynamic_font_colors = dynamic_font_color_settings(
    threshold = 130,
    by = "counts",
    all = c("#000000", "#ffffff")),
  sums_settings = sum_tile_settings(
    label = "Total")) +
  labs(x = "Actual Cluster", y = "Predicted Cluster", title = "Model 2", subtitle = "Baseline and one-year prior characteristics") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_blank(),       # Removes Y-axis text labels
        axis.title.y = element_blank(),      # Ensures Y-axis title space is completely removed
        axis.ticks.y = element_blank())      # Removes small tick marks on the Y-axis


# Combine into a single figure to output using patchwork package
png("output/main/conf_mat_both_detailed.png", width = 3500, height = 1750, res = 300)
plot_baseline + plot_prior
dev.off()


##### Appendix Figure: ROC curves by cluster and model #####

# Join the ROCs together into a single dataframe for comparison
roc_combined <- rbind(roc_curve(mod_preds_baseline, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4') %>% mutate(regression = "baseline"), 
                roc_curve(mod_preds_prior, truth = group, '.pred_Cluster 1', '.pred_Cluster 2', '.pred_Cluster 3', '.pred_Cluster 4') %>% mutate(regression = "prior"))

# Plot of ROC AUCs for both models together
roc_auc_plot <- ggplot(data = roc_combined, aes(x = 1 - specificity, y = sensitivity, color = .level)) +
  geom_line(linewidth = 1.1, alpha = 0.9, aes(linetype = regression)) +
  geom_abline(slope = 1, linetype = "dashed") +
  coord_fixed() +
  labs(color = "",
       linetype = "", 
       caption = paste0("Hand-Till AUC for Model 1 = ", sprintf("%.3f", round(roc_auc_baseline, 3)), 
                        "\nHand-Till AUC for Model 2 = ", sprintf("%.3f", round(roc_auc_prior, 3)))) +
  scale_linetype_manual(values=c("dotted", "solid"), 
                        labels = c("Model 1: Baseline only", 
                                   "Model 2: Baseline and one-year prior to long-term sickness")) +
  ggh4x::facet_wrap2(~.level, strip = strip_themed(background_x = elem_list_rect(fill = c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3")), 
                                                   text_x = elem_list_text(face = "bold")))
  scale_color_discrete(palette = c("#8DA0CB", "#66C2A5", "#FC8D62", "#E78AC3")) + 
  theme_bw() + 
  guides(color = "none") + 
  theme(legend.position = "top", 
        strip.text = element_text(face = "bold"))  # Bold for facet strip titles

# Save ROC-AUC plot 
png("output/main/roc_auc_plot.png", width = 1950, height = 1950, res = 300)
roc_auc_plot
dev.off()


#### 4. REGRESSION TABLE ####

##### Table: Associations between covariates and cluster memberships #####

# Baseline
model_baseline_to_save <- broom::tidy(model_baseline, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate_at(vars(matches("estimate")), round, 2) %>%  
  select(-std.error, -statistic, -conf.low, -conf.high) %>% 
  mutate(stars = ifelse(p.value < 0.001, "***", 
                        ifelse(p.value < 0.01, "**", 
                               ifelse(p.value < 0.05, "*", "")))) %>% 
  pivot_wider(names_from = y.level, values_from = c(estimate, p.value, stars))
write.csv(model_baseline_to_save, "output/main/model_baseline.csv")

# Baseline and one-year prior
model_prior_to_save <- tidy(model_prior, exponentiate = TRUE, conf.int = TRUE) %>% 
  mutate_at(vars(matches("estimate")), round, 2) %>%  
  select(-std.error, -statistic, -conf.low, -conf.high) %>% 
  mutate(stars = ifelse(p.value < 0.001, "***", 
                ifelse(p.value < 0.01, "**", 
                       ifelse(p.value < 0.05, "*", "")))) %>% 
  pivot_wider(names_from = y.level, values_from = c(estimate, p.value, stars))
write.csv(model_prior_to_save, "output/main/model_prior.csv")






         