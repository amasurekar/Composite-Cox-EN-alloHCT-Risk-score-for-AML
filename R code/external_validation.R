################################################################################
# EXTERNAL VALIDATION 
################################################################################

# Clear workspace
rm(list = ls())
gc()

# Set seed for reproducibility
set.seed(147201)

# Load packages (minimal set)
library(readxl)
library(survival)
library(survcomp)
library(boot)
library(ggplot2)
library(openxlsx)
library(randomForestSRC)

# CHANGE THIS: Set to YOUR actual desktop path
desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
output_folder <- file.path("C:/Users/drash/OneDrive/Desktop/AML_Model_Comparison/03_Results/External_Validation_of_Models")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Initialize log file
log_file <- file.path(output_folder, paste0("validation_log_", format(Sys.Date(), "%Y%m%d"), ".txt"))
sink(log_file, split = TRUE)

cat("\n=== EXTERNAL VALIDATION ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Saving all outputs to:", output_folder, "\n\n")

# Set working directory for input files
setwd("C:/Users/drash/OneDrive/Desktop/AML_Model_Comparison/01_Raw_Data")

################################################################################
# LOAD DATA
################################################################################

cat("Loading data...\n")
external_cohort <- read_excel("External cohort.xlsx")
external_cohort$dummyid <- NULL

# Convert to numeric
vars <- c("age", "sex", "race", "kps", "hctci", "amlgp", "WC", "elngroup", 
          "time_to_CR", "mrd", "donor", "cond", "graftype", "lymphodepletion", 
          "Time", "Status")

for (v in vars) {
  if (v %in% names(external_cohort)) {
    external_cohort[[v]] <- as.numeric(external_cohort[[v]])
  }
}

# Create graft_donor
external_cohort$graft_donor <- NA
external_cohort$graft_donor[external_cohort$graftype == 1 & external_cohort$donor == 1] <- 1
external_cohort$graft_donor[external_cohort$graftype == 1 & external_cohort$donor == 2] <- 2
external_cohort$graft_donor[external_cohort$graftype == 1 & external_cohort$donor == 3] <- 3
external_cohort$graft_donor[external_cohort$graftype == 2 & external_cohort$donor == 1] <- 4
external_cohort$graft_donor[external_cohort$graftype == 2 & external_cohort$donor == 2] <- 5
external_cohort$graft_donor[external_cohort$graftype == 2 & external_cohort$donor == 3] <- 6
external_cohort$graft_donor[external_cohort$graftype == 3 & external_cohort$donor == 4] <- 7

cat("Loaded", nrow(external_cohort), "patients\n")

################################################################################
# CALCULATE SCORES
################################################################################

# Helper function
extract_var <- function(x) {
  v <- gsub("\\d.*", "", x)
  if (v == "timetoachieveCR") v <- "time_to_CR"
  if (v == "Lymphodepletion") v <- "lymphodepletion"
  return(v)
}

# Cox TVC
cat("\nCalculating Cox TVC scores...\n")
cox_tvc_data <- read_excel("cox_tvc_coefficients.xlsx")
external_cohort$Cox_TVC <- 0

for (i in 1:nrow(cox_tvc_data)) {
  if (is.na(cox_tvc_data$hr[i]) || grepl("^TVC:", cox_tvc_data$variable[i])) next
  
  parts <- strsplit(cox_tvc_data$variable[i], "\\.")[[1]]
  if (length(parts) < 2) next
  
  category <- as.numeric(parts[1])
  var_name <- extract_var(paste(parts[-1], collapse = "."))
  
  if (var_name %in% names(external_cohort)) {
    mask <- external_cohort[[var_name]] == category & !is.na(external_cohort[[var_name]])
    if (sum(mask) > 0) {
      external_cohort$Cox_TVC[mask] <- external_cohort$Cox_TVC[mask] + log(cox_tvc_data$hr[i])
    }
  }
}

# Cox EN
cat("Calculating Cox EN scores...\n")
cox_en_data <- read_excel("cox_en_coefficients.xlsx")
external_cohort$Cox_EN <- 0

for (i in 1:nrow(cox_en_data)) {
  if (is.na(cox_en_data$pooled_coef[i]) || cox_en_data$pooled_coef[i] == 0) next
  
  parts <- strsplit(cox_en_data$variable[i], "\\.")[[1]]
  if (length(parts) < 2) next
  
  category <- as.numeric(parts[1])
  var_name <- extract_var(paste(parts[-1], collapse = "."))
  
  if (var_name %in% names(external_cohort)) {
    mask <- external_cohort[[var_name]] == category & !is.na(external_cohort[[var_name]])
    if (sum(mask) > 0) {
      external_cohort$Cox_EN[mask] <- external_cohort$Cox_EN[mask] + cox_en_data$pooled_coef[i]
    }
  }
}

# RSF
cat("Calculating RSF scores...\n")
load("rsf_workspace.RData")
rsf_data <- external_cohort[, c("age", "sex", "race", "kps", "hctci", "amlgp", 
                                 "WC", "elngroup", "time_to_CR", "mrd", "donor", 
                                 "cond", "graftype", "lymphodepletion", "Time", "Status")]

predictions <- matrix(NA, nrow = nrow(rsf_data), ncol = length(final_models))
for (i in 1:length(final_models)) {
  tryCatch({
    pred <- predict(final_models[[i]], newdata = rsf_data, na.action = "na.impute")
    predictions[, i] <- pred$predicted
  }, error = function(e) {})
}
external_cohort$RSF <- rowMeans(predictions, na.rm = TRUE)

################################################################################
# CALCULATE C-INDEX
################################################################################

cat("\nCalculating C-indices...\n")

calc_cindex <- function(time, status, score, model_name) {
  valid <- complete.cases(time, status, score)
  
  c_res <- concordance.index(score[valid], time[valid], status[valid], method = "noether")
  
  # Bootstrap CI with seed
  set.seed(147201)
  boot_fun <- function(data, idx) {
    concordance.index(data$score[idx], data$time[idx], data$status[idx], method = "noether")$c.index
  }
  
  boot_data <- data.frame(time = time[valid], status = status[valid], score = score[valid])
  boot_res <- boot(boot_data, boot_fun, R = 1000)
  ci <- quantile(boot_res$t, c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    model = model_name,
    n = sum(valid),
    events = sum(status[valid]),
    c_index = c_res$c.index,
    se = c_res$se,
    ci_lower = ci[1],
    ci_upper = ci[2]
  )
}

# Calculate for each model
results <- rbind(
  calc_cindex(external_cohort$Time, external_cohort$Status, external_cohort$Cox_TVC, "Cox TVC"),
  calc_cindex(external_cohort$Time, external_cohort$Status, external_cohort$Cox_EN, "Cox EN"),
  calc_cindex(external_cohort$Time, external_cohort$Status, external_cohort$RSF, "RSF")
)

################################################################################
# SAVE OUTPUTS DIRECTLY
################################################################################

cat("\n=== RESULTS ===\n")
print(results)

# 1. CSV
csv_file <- file.path(output_folder, "cindex_results.csv")
write.csv(results, csv_file, row.names = FALSE)
cat("\nSaved:", csv_file, "\n")

# 2. Plot
p <- ggplot(results, aes(x = model, y = c_index, fill = model)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", c_index)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Cox TVC" = "#E41A1C", "Cox EN" = "#377EB8", "RSF" = "#4DAF4A")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  labs(title = "External Validation C-index", x = "Model", y = "C-index") +
  theme_minimal() + theme(legend.position = "none") + ylim(0, 0.8)

png_file <- file.path(output_folder, "cindex_plot.png")
ggsave(png_file, p, width = 6, height = 5, dpi = 300)
cat("Saved:", png_file, "\n")

# 3. Text summary
txt_file <- file.path(output_folder, "results_summary.txt")
sink(txt_file)
cat("EXTERNAL VALIDATION RESULTS\n")
cat("Date:", format(Sys.Date()), "\n\n")
for(i in 1:nrow(results)) {
  cat(sprintf("%s: %.3f [%.3f-%.3f]\n", 
              results$model[i], results$c_index[i], 
              results$ci_lower[i], results$ci_upper[i]))
}
sink()
cat("Saved:", txt_file, "\n")

# 4. Excel
wb <- createWorkbook()
addWorksheet(wb, "Results")
writeData(wb, "Results", results)
xlsx_file <- file.path(output_folder, "results.xlsx")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
cat("Saved:", xlsx_file, "\n")

# 5. Workspace
rdata_file <- file.path(output_folder, "validation_workspace.RData")
save.image(rdata_file)
cat("Saved:", rdata_file, "\n")

cat("\n=== COMPLETE ===\n")
cat("All files saved to:", output_folder, "\n")

# Close log file
sink()
