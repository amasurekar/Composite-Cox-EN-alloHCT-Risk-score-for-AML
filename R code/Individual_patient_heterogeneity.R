# Individual Outcome Heterogeneity Analysis
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(survival)
library(devEMF)

# Set up output directory
output_dir <- "C:/Users/drash/OneDrive/Desktop/AML_Model_Comparison/Individual_Patient_Analysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load workspace
workspace_path <- "C:/Users/drash/OneDrive/Desktop/AML_Model_Comparison/model_comparison/complete_analysis.RData"
if(!file.exists(workspace_path)) {
  stop("Workspace file not found. Please run the complete analysis first.")
}
load(workspace_path)

if(!exists("scores_data")) {
  stop("scores_data object not found in workspace")
}

cat("Data loaded successfully. Rows:", nrow(scores_data), "\n")

# Check for recalibrated probability columns - adjust column names as needed
recal_columns <- c("recal_prob_Cox_TVC_24", "recal_prob_Cox_EN_24", "recal_prob_RSF_24")
alt_columns <- c("calibrated_Cox_TVC_24", "calibrated_Cox_EN_24", "calibrated_RSF_24")

# Check which recalibrated columns exist
available_recal <- recal_columns[recal_columns %in% names(scores_data)]
available_alt <- alt_columns[alt_columns %in% names(scores_data)]

if(length(available_recal) == 3) {
  # Use recal_prob_ columns
  scores_data$prob_Cox_TVC_24 <- scores_data$recal_prob_Cox_TVC_24
  scores_data$prob_Cox_EN_24 <- scores_data$recal_prob_Cox_EN_24
  scores_data$prob_RSF_24 <- scores_data$recal_prob_RSF_24
  cat("Using recalibrated probabilities (recal_prob_ columns)\n")
} else if(length(available_alt) == 3) {
  # Use calibrated_ columns
  scores_data$prob_Cox_TVC_24 <- scores_data$calibrated_Cox_TVC_24
  scores_data$prob_Cox_EN_24 <- scores_data$calibrated_Cox_EN_24
  scores_data$prob_RSF_24 <- scores_data$calibrated_RSF_24
  cat("Using recalibrated probabilities (calibrated_ columns)\n")
} else {
  # Fall back to raw predictions with warning
  scores_data$prob_Cox_TVC_24 <- scores_data$pred_prob_Cox_TVC_24
  scores_data$prob_Cox_EN_24 <- scores_data$pred_prob_Cox_EN_24
  scores_data$prob_RSF_24 <- scores_data$pred_prob_RSF_24
  cat("WARNING: Recalibrated columns not found. Using raw predictions.\n")
  cat("Available columns:", paste(names(scores_data)[grepl("prob|calibrated", names(scores_data))], collapse = ", "), "\n")
}

# Verify required columns exist
required_cols <- c("prob_Cox_TVC_24", "prob_Cox_EN_24", "prob_RSF_24", "Time", "Status")
missing_cols <- required_cols[!required_cols %in% names(scores_data)]
if(length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

cat("Required columns present:", all(required_cols %in% names(scores_data)), "\n")

# Create risk groups for each model separately using recalibrated probabilities
create_risk_groups <- function(data, prob_col, group_suffix) {
  tertiles <- quantile(data[[prob_col]], probs = c(0.333, 0.667), na.rm = TRUE)
  risk_col <- paste0("risk_", group_suffix)
  
  cat("Tertiles for", prob_col, ":", round(tertiles, 3), "\n")
  
  data[[risk_col]] <- cut(data[[prob_col]], 
                         breaks = c(0, tertiles, 1),
                         labels = c("Low", "Intermediate", "High"),
                         include.lowest = TRUE)
  return(data)
}

# Apply risk stratification for each model
scores_data <- scores_data %>%
  create_risk_groups("prob_Cox_TVC_24", "Cox_TVC") %>%
  create_risk_groups("prob_Cox_EN_24", "Cox_EN") %>%
  create_risk_groups("prob_RSF_24", "RSF")

# Filter for complete data
analysis_data <- scores_data %>%
  filter(!is.na(prob_Cox_TVC_24) & !is.na(prob_Cox_EN_24) & !is.na(prob_RSF_24) &
         !is.na(risk_Cox_TVC) & !is.na(risk_Cox_EN) & !is.na(risk_RSF) &
         !is.na(Time) & !is.na(Status))

cat("Total patients analyzed:", nrow(analysis_data), "\n")
cat("Overall event rate:", round(mean(analysis_data$Status) * 100, 1), "%\n")

# Create 24-month outcome variable
analysis_data$actual_24mo <- factor(
  ifelse(analysis_data$Time <= 24 & analysis_data$Status == 1, 
         "Died by 24mo", "Survived 24mo"),
  levels = c("Survived 24mo", "Died by 24mo")
)

# Create separate datasets for each model to preserve risk group information
cox_tvc_data <- analysis_data %>%
  select(Time, Status, actual_24mo, predicted_24mo = prob_Cox_TVC_24, risk_group = risk_Cox_TVC) %>%
  mutate(model = "Cox-TVC")

cox_en_data <- analysis_data %>%
  select(Time, Status, actual_24mo, predicted_24mo = prob_Cox_EN_24, risk_group = risk_Cox_EN) %>%
  mutate(model = "Cox-EN")

rsf_data <- analysis_data %>%
  select(Time, Status, actual_24mo, predicted_24mo = prob_RSF_24, risk_group = risk_RSF) %>%
  mutate(model = "RSF")

# Combine for density plots
density_data <- bind_rows(cox_tvc_data, cox_en_data, rsf_data) %>%
  filter(!is.na(risk_group) & !is.na(predicted_24mo))

cat("Density plot data prepared. Total observations:", nrow(density_data), "\n")

# Figure 2A: Density plots by model and risk group
p_calibration <- ggplot(density_data, aes(x = predicted_24mo, fill = actual_24mo)) +
  geom_density(alpha = 0.6, adjust = 1.2) +
  facet_grid(risk_group ~ model) +
  scale_fill_manual(values = c("Survived 24mo" = "#2E8B57", 
                              "Died by 24mo" = "#DC143C"),
                   name = "24-Month Outcome") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "A. Individual Outcome Heterogeneity Within Risk Groups by Model",
    x = "Recalibrated 24-month mortality probability",
    y = "Density"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13)
  )

# Prepare scatter plot data for each model
scatter_data_list <- list(
  Cox_TVC = analysis_data %>% 
    mutate(model = "Cox-TVC", predicted = prob_Cox_TVC_24, risk_group = risk_Cox_TVC),
  Cox_EN = analysis_data %>% 
    mutate(model = "Cox-EN", predicted = prob_Cox_EN_24, risk_group = risk_Cox_EN),
  RSF = analysis_data %>% 
    mutate(model = "RSF", predicted = prob_RSF_24, risk_group = risk_RSF)
)

scatter_data <- bind_rows(scatter_data_list)

# Figure 2B: Scatter plots by model
p_scatter <- ggplot(scatter_data, aes(x = predicted, y = Time)) +
  geom_point(aes(color = risk_group, shape = factor(Status)), 
            size = 1.5, alpha = 0.7) +
  facet_wrap(~ model, ncol = 3) +
  scale_shape_manual(values = c("0" = 16, "1" = 4),
                    labels = c("0" = "Censored", "1" = "Died"),
                    name = "Outcome") +
  scale_color_manual(values = c("Low" = "#2E8B57", 
                               "Intermediate" = "#FF8C00",
                               "High" = "#DC143C"),
                    name = "Risk Group") +
  geom_smooth(method = "loess", se = TRUE, color = "black", 
             linewidth = 0.8, alpha = 0.8) +
  labs(
    title = "B. Individual Outcome Heterogeneity Across Prediction Spectrum",
    subtitle = paste("Total patients:", nrow(analysis_data)),
    x = "Recalibrated 24-month mortality probability",
    y = "Survival time (months)"
  ) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(trans = "log10",
                    breaks = c(1, 3, 6, 12, 24, 36, 60),
                    limits = c(0.5, 80)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold")
  )

# Calculate KM mortality rates for each model's risk groups
calculate_km_by_model <- function(data, risk_col, model_name) {
  km_fit <- survfit(as.formula(paste("Surv(Time, Status) ~", risk_col)), data = data)
  km_24 <- summary(km_fit, times = 24)
  
  if(is.null(km_24$surv)) {
    warning(paste("No 24-month data available for", model_name))
    return(data.frame())
  }
  
  mortality_rates <- 1 - km_24$surv
  strata_names <- gsub(paste0(risk_col, "="), "", names(km_fit$strata))
  
  result <- data.frame(
    model = model_name,
    risk_group = strata_names,
    km_mortality_24 = mortality_rates,
    n_risk = km_24$n.risk,
    n_event = km_24$n.event,
    stringsAsFactors = FALSE
  )
  return(result)
}

# Calculate KM rates for all models
km_results <- bind_rows(
  calculate_km_by_model(analysis_data, "risk_Cox_TVC", "Cox-TVC"),
  calculate_km_by_model(analysis_data, "risk_Cox_EN", "Cox-EN"),
  calculate_km_by_model(analysis_data, "risk_RSF", "RSF")
)

cat("KM mortality rates calculated for", nrow(km_results), "model-risk group combinations\n")

# Create summary statistics
create_model_summary <- function(data, model_col, risk_col, model_name) {
  data %>%
    filter(!is.na(.data[[model_col]]) & !is.na(.data[[risk_col]])) %>%
    group_by(risk_group = .data[[risk_col]]) %>%
    summarise(
      model = model_name,
      n = n(),
      deaths = sum(Status),
      death_rate = mean(Status),
      mean_pred_24 = mean(.data[[model_col]]),
      sd_pred_24 = sd(.data[[model_col]]),
      median_pred_24 = median(.data[[model_col]]),
      q25_pred_24 = quantile(.data[[model_col]], 0.25),
      q75_pred_24 = quantile(.data[[model_col]], 0.75),
      .groups = 'drop'
    ) %>%
    left_join(km_results %>% 
              filter(model == model_name) %>% 
              select(risk_group, km_mortality_24),
              by = "risk_group") %>%
    mutate(
      calibration_diff = mean_pred_24 - km_mortality_24,
      risk_group = as.character(risk_group)
    )
}

model_summaries <- bind_rows(
  create_model_summary(analysis_data, "prob_Cox_TVC_24", "risk_Cox_TVC", "Cox-TVC"),
  create_model_summary(analysis_data, "prob_Cox_EN_24", "risk_Cox_EN", "Cox-EN"),
  create_model_summary(analysis_data, "prob_RSF_24", "risk_RSF", "RSF")
)

# Save figures
cat("Saving figures...\n")

# Save Figure 2A
tryCatch({
  emf(file.path(output_dir, "Figure_2A_Individual_Outcome_Heterogeneity.emf"), width = 12, height = 8)
  print(p_calibration)
  dev.off()
  cat("Figure 2A saved successfully\n")
}, error = function(e) {
  cat("Error saving Figure 2A:", e$message, "\n")
  if(names(dev.cur()) != "null device") dev.off()
})

# Save Figure 2B
tryCatch({
  emf(file.path(output_dir, "Figure_2B_Individual_Outcome_Heterogeneity.emf"), width = 12, height = 6)
  print(p_scatter)
  dev.off()
  cat("Figure 2B saved successfully\n")
}, error = function(e) {
  cat("Error saving Figure 2B:", e$message, "\n")
  if(names(dev.cur()) != "null device") dev.off()
})

# Save combined figure
tryCatch({
  emf(file.path(output_dir, "Figure_2_Combined_Individual_Outcome_Heterogeneity.emf"), width = 12, height = 14)
  grid.arrange(p_calibration, p_scatter, ncol = 1, heights = c(1.2, 1))
  dev.off()
  cat("Combined figure saved successfully\n")
}, error = function(e) {
  cat("Error saving combined figure:", e$message, "\n")
  if(names(dev.cur()) != "null device") dev.off()
})

# Save summary data
write.csv(model_summaries, file.path(output_dir, "Individual_Outcome_Heterogeneity_Summary.csv"), row.names = FALSE)
write.csv(km_results, file.path(output_dir, "Individual_Outcome_Heterogeneity_KM_Results.csv"), row.names = FALSE)

# Create comprehensive log file
log_file <- file.path(output_dir, "Individual_Outcome_Heterogeneity_Log.txt")
sink(log_file)

cat("INDIVIDUAL OUTCOME HETEROGENEITY\n")
cat("================================\n")
cat("Analysis Date:", format(Sys.Date(), "%B %d, %Y"), "\n")
cat("Analysis Time:", format(Sys.time(), "%H:%M:%S"), "\n\n")

cat("DATA SUMMARY\n")
cat("------------\n")
cat("Total patients analyzed:", nrow(analysis_data), "\n")
cat("Total events (deaths):", sum(analysis_data$Status), "\n")
cat("Event rate:", round(mean(analysis_data$Status) * 100, 1), "%\n")
cat("Median follow-up:", round(median(analysis_data$Time), 1), "months\n\n")

cat("RISK GROUP DISTRIBUTIONS BY MODEL\n")
cat("---------------------------------\n")
for(model in c("Cox-TVC", "Cox-EN", "RSF")) {
  risk_col <- paste0("risk_", gsub("-", "_", model))
  counts <- table(analysis_data[[risk_col]])
  cat(model, "risk groups:\n")
  print(counts)
  cat("Proportions:\n")
  print(round(prop.table(counts) * 100, 1))
  cat("\n")
}

cat("KAPLAN-MEIER MORTALITY RATES BY MODEL AND RISK GROUP\n")
cat("----------------------------------------------------\n")
print(km_results)

cat("\n\nMODEL SUMMARY STATISTICS\n")
cat("------------------------\n")
print(model_summaries)

sink()

# Print summary to console
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total patients:", nrow(analysis_data), "\n")
cat("Overall event rate:", round(mean(analysis_data$Status) * 100, 1), "%\n\n")

cat("Model Summary (Calibration within Risk Groups):\n")
print(model_summaries %>% select(model, risk_group, n, km_mortality_24, mean_pred_24, calibration_diff))

cat("\nFiles generated in:", output_dir, "\n")
files_created <- list.files(output_dir, pattern = "Individual_Outcome_Heterogeneity", full.names = FALSE)
if(length(files_created) > 0) {
  print(files_created)
} else {
  print(list.files(output_dir, full.names = FALSE))
}

cat("\nAnalysis demonstrates individual-level outcome heterogeneity within risk strata\n")
cat("using recalibrated probabilities for optimal model performance.\n")

