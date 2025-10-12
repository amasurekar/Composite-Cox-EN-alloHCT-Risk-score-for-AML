# ====================================================================
# AML MODEL COMPARISON - COMPLETE CORRECTED VERSION.
# ====================================================================

rm(list = ls())
gc()
set.seed(147201)

# Set working directory and create output structure
setwd("C:/Users/drash/OneDrive/Desktop/AML_Model_Comparison")
dir.create("model_comparison", recursive = TRUE, showWarnings = FALSE)
dir.create("model_comparison/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("model_comparison/Tables", recursive = TRUE, showWarnings = FALSE)

# Initialize logging
log_file <- file("model_comparison/Complete_Analysis_Log.txt", open = "wt")
sink(log_file, type = "output", split = TRUE)
sink(log_file, type = "message", append = TRUE)

cat("===============================================================\n")
cat("COMPLETE AML MODEL COMPARISON ANALYSIS\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("===============================================================\n\n")

# Load required packages
required_packages <- c("readxl", "survival", "survcomp", "timeROC", "pec",
                      "riskRegression", "ggplot2", "gridExtra", "dplyr",
                      "tidyr", "boot", "reshape2", "openxlsx", "cowplot")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ====================================================================
# PART 1: DATA LOADING AND PREPARATION
# ====================================================================

# Load scores data
scores_data <- read_excel("01_Raw_Data/internal_scores.xlsx")
scores_data$HCTCI <- round(scores_data$HCTCI)

# Load RSF workspace
rsf_workspace_path <- "01_Raw_Data/rsf_workspace.RData"
rsf_available <- FALSE

if(file.exists(rsf_workspace_path)) {
  load(rsf_workspace_path)
  if(exists("final_models")) {
    rsf_oob_predictions <- rowMeans(sapply(final_models, function(m) m$predicted.oob))
    scores_data$RSF <- rsf_oob_predictions
    rsf_available <- TRUE
  }
}

# Pre-calculated optimism-corrected C-index values
cox_optimism_corrected <- data.frame(
  model = c("Cox_TVC", "Cox_EN"),
  corrected_c = c(0.612620950, 0.613302350),
  ci_lower = c(0.609411061, 0.610091865),
  ci_upper = c(0.615830839, 0.616512835)
)

# Create survival object
surv_obj <- Surv(scores_data$Time, scores_data$Status)

# Score direction check
check_score_direction <- function(score, status, score_name) {
  mean_events <- mean(score[status == 1], na.rm = TRUE)
  mean_nonevents <- mean(score[status == 0], na.rm = TRUE)
  if(mean_events < mean_nonevents) {
    return(-score)
  } else {
    return(score)
  }
}

scores_data$Cox_TVC <- check_score_direction(scores_data$Cox_TVC, scores_data$Status, "Cox_TVC")
scores_data$Cox_EN <- check_score_direction(scores_data$Cox_EN, scores_data$Status, "Cox_EN")

if(rsf_available) {
  scores_data$RSF <- check_score_direction(scores_data$RSF, scores_data$Status, "RSF")
}

# Create binary outcome variables
scores_data$event_12 <- ifelse(scores_data$Time <= 12 & scores_data$Status == 1, 1, 0)
scores_data$event_24 <- ifelse(scores_data$Time <= 24 & scores_data$Status == 1, 1, 0)
scores_data$event_36 <- ifelse(scores_data$Time <= 36 & scores_data$Status == 1, 1, 0)

# Convert categorical variables to factors
scores_data$HCTCI <- as.factor(scores_data$HCTCI)
scores_data$ELN <- as.factor(scores_data$ELN)

# Model names
model_names <- c("Cox_TVC", "Cox_EN", "HCTCI", "ELN")
if(rsf_available) model_names <- c(model_names, "RSF")

# Define colors
model_colors <- c(
  "Cox_TVC" = "#E41A1C",
  "Cox_EN" = "#377EB8",
  "RSF" = "#4DAF4A",
  "HCTCI" = "#984EA3",
  "ELN" = "#FF7F00"
)

# ====================================================================
# PART 2: PRE-RECALIBRATION PROBABILITY CALCULATION
# ====================================================================

calculate_predicted_probabilities <- function(data, t) {
  # Cox models - convert to probabilities using baseline hazard
  cox_tvc_fit <- coxph(surv_obj ~ Cox_TVC, data = data)
  cox_en_fit <- coxph(surv_obj ~ Cox_EN, data = data)
  
  basehaz_tvc <- basehaz(cox_tvc_fit, centered = FALSE)
  basehaz_en <- basehaz(cox_en_fit, centered = FALSE)
  
  idx_tvc <- findInterval(t, basehaz_tvc$time)
  idx_en <- findInterval(t, basehaz_en$time)
  
  H0_tvc <- ifelse(idx_tvc == 0, 0, basehaz_tvc$hazard[idx_tvc])
  H0_en <- ifelse(idx_en == 0, 0, basehaz_en$hazard[idx_en])
  
  data[[paste0("pred_prob_Cox_TVC_", t)]] <- 1 - exp(-H0_tvc * exp(data$Cox_TVC))
  data[[paste0("pred_prob_Cox_EN_", t)]] <- 1 - exp(-H0_en * exp(data$Cox_EN))
  
  # RSF - use survival probabilities
  if(rsf_available) {
    surv_probs_t <- sapply(final_models, function(model) {
      time_idx <- which.min(abs(model$time.interest - t))
      return(model$survival.oob[, time_idx])
    })
    data[[paste0("pred_prob_RSF_", t)]] <- 1 - rowMeans(surv_probs_t)
  }
  
  # Categorical models (HCTCI and ELN) - use Kaplan-Meier
  for(model in c("HCTCI", "ELN")) {
    prob_col <- paste0("pred_prob_", model, "_", t)
    data[[prob_col]] <- NA
    
    for(level in levels(data[[model]])) {
      subset_idx <- which(data[[model]] == level)
      subset_data <- data[subset_idx, ]
      
      km_fit <- survfit(Surv(Time, Status) ~ 1, data = subset_data)
      km_summary <- summary(km_fit, times = t, extend = TRUE)
      
      if(length(km_summary$surv) > 0) {
        event_prob <- 1 - km_summary$surv[1]
      } else {
        idx <- max(which(km_fit$time <= t))
        event_prob <- ifelse(length(idx) == 0, 0, 1 - km_fit$surv[idx])
      }
      data[subset_idx, prob_col] <- event_prob
    }
  }
  
  return(data)
}

# Calculate predicted probabilities for all time points
for(t in c(12, 24, 36)) {
  scores_data <- calculate_predicted_probabilities(scores_data, t)
}

# ====================================================================
# PART 3: CALIBRATION ON LOGIT SCALE
# ====================================================================

calibrate_on_logit_scale <- function(pred_prob, outcome) {
  pred_logit <- qlogis(pmax(pmin(pred_prob, 0.9999), 0.0001))
  cal_fit <- glm(outcome ~ pred_logit, family = binomial())
  recal_prob <- predict(cal_fit, type = "response")
  
  return(list(
    probs = recal_prob,
    model = cal_fit,
    intercept = coef(cal_fit)[1],
    slope = coef(cal_fit)[2]
  ))
}

# Apply calibration for all models and time points
calibration_results <- list()

for(t in c(12, 24, 36)) {
  outcome <- scores_data[[paste0("event_", t)]]
  
  for(model in model_names) {
    pred_prob_col <- paste0("pred_prob_", model, "_", t)
    pred_prob <- scores_data[[pred_prob_col]]
    
    cal_result <- calibrate_on_logit_scale(pred_prob, outcome)
    scores_data[[paste0("recal_prob_", model, "_", t)]] <- cal_result$probs
    
    calibration_results[[paste0(model, "_", t)]] <- list(
      model = model,
      time = t,
      intercept = cal_result$intercept,
      slope = cal_result$slope,
      cal_fit = cal_result$model
    )
  }
}

# ====================================================================
# PART 4: CREATE MODELS FOR DISCRIMINATION ANALYSIS (CORRECTED)
# ====================================================================

cat("Creating models for discrimination analysis...\n")

cox_tvc_fit <- coxph(Surv(Time, Status) ~ Cox_TVC, data = scores_data, x = TRUE)
cox_en_fit <- coxph(Surv(Time, Status) ~ Cox_EN, data = scores_data, x = TRUE)
hctci_fit <- coxph(Surv(Time, Status) ~ as.factor(HCTCI), data = scores_data, x = TRUE)
eln_fit <- coxph(Surv(Time, Status) ~ as.factor(ELN), data = scores_data, x = TRUE)

scores_data$Cox_TVC_cal <- predict(cox_tvc_fit, type = "lp")
scores_data$Cox_EN_cal <- predict(cox_en_fit, type = "lp")
scores_data$HCTCI_cal <- predict(hctci_fit, type = "lp")
scores_data$ELN_cal <- predict(eln_fit, type = "lp")

if(rsf_available) {
  scores_data$RSF_cal <- scores_data$RSF
}

models <- list(Cox_TVC = cox_tvc_fit, Cox_EN = cox_en_fit, HCTCI = hctci_fit, ELN = eln_fit)
if(rsf_available) models$RSF <- NULL

cat("Model creation complete. _cal columns created.\n\n")

# ====================================================================
# PART 5: C-INDEX ANALYSIS
# ====================================================================

calc_apparent_cindex <- function(time, status, score, name) {
  complete_idx <- complete.cases(time, status, score)
  time <- time[complete_idx]
  status <- status[complete_idx]
  score <- score[complete_idx]
  
  c_result <- concordance.index(x = score, surv.time = time,
                               surv.event = status, method = "noether")
  c_index <- c_result$c.index
  
  set.seed(147201)
  n_boot <- 1000
  
  boot_fun <- function(data, indices) {
    d <- data[indices,]
    concordance.index(x = d$score, surv.time = d$time,
                     surv.event = d$status, method = "noether")$c.index
  }
  
  boot_data <- data.frame(time = time, status = status, score = score)
  boot_res <- boot(boot_data, boot_fun, R = n_boot)
  ci <- quantile(boot_res$t, c(0.025, 0.975), na.rm = TRUE)
  
  return(data.frame(
    model = name,
    c_index = c_index,
    ci_lower = ci[1],
    ci_upper = ci[2],
    stringsAsFactors = FALSE
  ))
}

cindex_results <- list()

# Use pre-calculated optimism-corrected values for Cox models
cindex_results[[1]] <- data.frame(
  model = "Cox_TVC",
  c_index = cox_optimism_corrected$corrected_c[1],
  ci_lower = cox_optimism_corrected$ci_lower[1],
  ci_upper = cox_optimism_corrected$ci_upper[1]
)

cindex_results[[2]] <- data.frame(
  model = "Cox_EN",
  c_index = cox_optimism_corrected$corrected_c[2],
  ci_lower = cox_optimism_corrected$ci_lower[2],
  ci_upper = cox_optimism_corrected$ci_upper[2]
)

cindex_results[[3]] <- calc_apparent_cindex(scores_data$Time, scores_data$Status,
                                          scores_data$HCTCI_cal, "HCTCI")

cindex_results[[4]] <- calc_apparent_cindex(scores_data$Time, scores_data$Status,
                                          scores_data$ELN_cal, "ELN")

if(rsf_available) {
  rsf_oob_cindex <- mean(sapply(final_models, function(m) 1 - m$err.rate[m$ntree]))
  rsf_oob_sd <- sd(sapply(final_models, function(m) 1 - m$err.rate[m$ntree]))
  
  cindex_results[[5]] <- data.frame(
    model = "RSF",
    c_index = rsf_oob_cindex,
    ci_lower = rsf_oob_cindex - 1.96 * rsf_oob_sd,
    ci_upper = rsf_oob_cindex + 1.96 * rsf_oob_sd
  )
}

cindex_df <- do.call(rbind, cindex_results)
cindex_df[, 2:4] <- round(cindex_df[, 2:4], 3)

write.csv(cindex_df, "model_comparison/Tables/cindex_results.csv", row.names = FALSE)

# ====================================================================
# PART 6: TIME-DEPENDENT AUC AND ROC CURVES (CORRECTED)
# ====================================================================

time_points <- c(12, 24, 36)
auc_results <- data.frame()
auc_results_lower <- data.frame()
auc_results_upper <- data.frame()

for(t in time_points) {
  cat("=== Processing time point:", t, "months ===\n")
  
  temp_auc <- data.frame(Time = t)
  temp_auc_lower <- data.frame(Time = t)
  temp_auc_upper <- data.frame(Time = t)
  roc_list <- list()
  
  for(model in model_names) {
    cat("Model:", model, "\n")
    
    score <- scores_data[[paste0(model, "_cal")]]
    complete_idx <- complete.cases(scores_data$Time, scores_data$Status, score)
    
    cat("  Complete cases:", sum(complete_idx), "\n")
    cat("  Events:", sum(scores_data$Status[complete_idx] == 1), "\n")
    cat("  Score range:", range(score[complete_idx], na.rm = TRUE), "\n")
    
    tryCatch({
      roc_obj <- timeROC(
        T = scores_data$Time[complete_idx],
        delta = scores_data$Status[complete_idx],
        marker = score[complete_idx],
        cause = 1,
        times = t,
        iid = FALSE
      )
      
      cat("  AUC vector:", roc_obj$AUC, "\n")
      
      # Extract first non-NA AUC value
      auc_value <- roc_obj$AUC[!is.na(roc_obj$AUC)][1]
      
      if(is.na(auc_value) || is.null(auc_value)) {
        cat("  ERROR - No valid AUC found\n")
        temp_auc[[model]] <- NA
        temp_auc_lower[[model]] <- NA
        temp_auc_upper[[model]] <- NA
        next
      }
      
      temp_auc[[model]] <- auc_value
      cat("  SUCCESS - AUC:", auc_value, "\n")
      
      # Bootstrap for confidence intervals
      set.seed(147201)
      n_boot <- 1000
      boot_auc <- numeric(n_boot)
      
      for(b in 1:n_boot) {
        boot_idx <- sample(which(complete_idx), replace = TRUE)
        
        tryCatch({
          boot_roc <- timeROC(
            T = scores_data$Time[boot_idx],
            delta = scores_data$Status[boot_idx],
            marker = scores_data[[paste0(model, "_cal")]][boot_idx],
            cause = 1,
            times = t,
            iid = FALSE
          )
          boot_auc[b] <- boot_roc$AUC[!is.na(boot_roc$AUC)][1]
        }, error = function(e) {
          boot_auc[b] <- NA
        })
      }
      
      # Calculate confidence intervals
      valid_boot <- boot_auc[!is.na(boot_auc)]
      if(length(valid_boot) > 50) {
        ci_auc <- quantile(valid_boot, c(0.025, 0.975), na.rm = TRUE)
        temp_auc_lower[[model]] <- ci_auc[1]
        temp_auc_upper[[model]] <- ci_auc[2]
        cat("  Bootstrap CI: [", round(ci_auc[1], 3), ",", round(ci_auc[2], 3), "]\n")
      } else {
        temp_auc_lower[[model]] <- NA
        temp_auc_upper[[model]] <- NA
        cat("  Bootstrap failed - insufficient valid samples\n")
      }
      
      # Create ROC data for plotting
      event_indicator <- (scores_data$Time[complete_idx] <= t) & 
                        (scores_data$Status[complete_idx] == 1)
      control_indicator <- scores_data$Time[complete_idx] > t
      
      n_thresh <- 200
      thresholds <- quantile(score[complete_idx], 
                           probs = seq(0, 1, length.out = n_thresh),
                           na.rm = TRUE)
      
      tpr <- numeric(n_thresh)
      fpr <- numeric(n_thresh)
      
      for(i in 1:n_thresh) {
        pos_pred <- score[complete_idx] >= thresholds[i]
        
        if(sum(event_indicator) > 0) {
          tpr[i] <- sum(pos_pred & event_indicator) / sum(event_indicator)
        } else {
          tpr[i] <- 0
        }
        
        if(sum(control_indicator) > 0) {
          fpr[i] <- sum(pos_pred & control_indicator) / sum(control_indicator)
        } else {
          fpr[i] <- 0
        }
      }
      
      roc_df <- data.frame(FPR = fpr, TPR = tpr)
      roc_df <- rbind(c(0, 0), roc_df, c(1, 1))
      roc_df <- roc_df[order(roc_df$FPR, roc_df$TPR), ]
      roc_df <- unique(roc_df)
      
      roc_list[[model]] <- list(
        FPR = roc_df$FPR,
        TPR = roc_df$TPR,
        AUC = auc_value
      )
      
    }, error = function(e) {
      cat("  ERROR in timeROC:", e$message, "\n")
      temp_auc[[model]] <- NA
      temp_auc_lower[[model]] <- NA
      temp_auc_upper[[model]] <- NA
    })
    
    cat("\n")
  }
  
  # Add results to main dataframes
  auc_results <- rbind(auc_results, temp_auc)
  
  if(nrow(auc_results_lower) == 0) {
    auc_results_lower <- temp_auc_lower
    auc_results_upper <- temp_auc_upper
  } else {
    auc_results_lower <- rbind(auc_results_lower, temp_auc_lower)
    auc_results_upper <- rbind(auc_results_upper, temp_auc_upper)
  }
  
  cat("Results for", t, "months:", paste(round(unlist(temp_auc[-1]), 3), collapse = ", "), "\n")
  
  # Create ROC plot if we have valid data
  if(length(roc_list) > 0) {
    plot_data <- data.frame()
    for(model in names(roc_list)) {
      temp_df <- data.frame(
        FPR = roc_list[[model]]$FPR,
        TPR = roc_list[[model]]$TPR,
        Model = model,
        AUC = roc_list[[model]]$AUC
      )
      plot_data <- rbind(plot_data, temp_df)
    }
    
    p_roc <- ggplot(plot_data, aes(x = FPR, y = TPR, color = Model)) +
      geom_line(size = 1.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      scale_color_manual(
        values = model_colors,
        labels = paste0(names(roc_list), " (AUC = ",
                       sprintf("%.3f", sapply(roc_list, function(x) x$AUC)), ")")
      ) +
      labs(
        title = sprintf("Time-Dependent ROC Curves at %d Months", t),
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()
      ) +
      coord_equal()
    
    ggsave(sprintf("model_comparison/Figures/ROC_%dmonths.png", t), p_roc,
           width = 10, height = 8, dpi = 300)
    
    cat("ROC plot saved for", t, "months\n")
  }
  
  cat("\n")
}

# Save all results
write.csv(auc_results, "model_comparison/Tables/auc_results.csv", row.names = FALSE)
write.csv(auc_results_lower, "model_comparison/Tables/auc_results_lower.csv", row.names = FALSE)
write.csv(auc_results_upper, "model_comparison/Tables/auc_results_upper.csv", row.names = FALSE)

# ====================================================================
# PART 7: CALIBRATION PLOTS (PRE- AND POST-RECALIBRATION)
# ====================================================================

assess_calibration <- function(pred_prob, obs, model_name, pre_recal = TRUE) {
  complete_idx <- complete.cases(pred_prob, obs)
  pred_prob <- pred_prob[complete_idx]
  obs <- obs[complete_idx]
  
  mean_pred <- mean(pred_prob)
  mean_obs <- mean(obs)
  
  pred_logit <- qlogis(pmax(pmin(pred_prob, 0.9999), 0.0001))
  cal_fit <- glm(obs ~ pred_logit, family = binomial)
  
  cal_intercept <- coef(cal_fit)[1]
  cal_slope <- coef(cal_fit)[2]
  E_O_ratio <- mean_pred / mean_obs
  
  n_unique <- length(unique(pred_prob))
  
  if(n_unique <= 5) {
    cal_data <- data.frame(pred = pred_prob, obs = obs) %>%
      group_by(pred) %>%
      summarise(
        n = n(),
        pred_mean = mean(pred),
        obs_mean = mean(obs),
        obs_se = sqrt(obs_mean * (1 - obs_mean) / n),
        obs_lower = pmax(0, obs_mean - 1.96 * obs_se),
        obs_upper = pmin(1, obs_mean + 1.96 * obs_se),
        .groups = 'drop'
      )
  } else {
    n_groups <- 10
    pred_groups <- cut(pred_prob,
                      breaks = quantile(pred_prob, probs = seq(0, 1, length.out = n_groups + 1)),
                      include.lowest = TRUE,
                      labels = FALSE)
    
    cal_data <- data.frame(
      pred = pred_prob,
      obs = obs,
      group = pred_groups
    ) %>%
      group_by(group) %>%
      summarise(
        n = n(),
        pred_mean = mean(pred),
        obs_mean = mean(obs),
        obs_se = sqrt(obs_mean * (1 - obs_mean) / n),
        obs_lower = pmax(0, obs_mean - 1.96 * obs_se),
        obs_upper = pmin(1, obs_mean + 1.96 * obs_se),
        .groups = 'drop'
      ) %>%
      filter(!is.na(group))
  }
  
  return(list(
    model = model_name,
    mean_predicted = mean_pred,
    mean_observed = mean_obs,
    calibration_intercept = cal_intercept,
    calibration_slope = cal_slope,
    E_O_ratio = E_O_ratio,
    cal_data = cal_data
  ))
}

standardize_cal_data <- function(cal_data, model_name) {
  data.frame(
    pred_mean = cal_data$pred_mean,
    obs_mean = cal_data$obs_mean,
    obs_lower = cal_data$obs_lower,
    obs_upper = cal_data$obs_upper,
    n = cal_data$n,
    Model = model_name
  )
}

# Pre-recalibration assessment at 24 months
cal_results_pre <- list()

for(model in model_names) {
  pred_prob_col <- paste0("pred_prob_", model, "_24")
  cal_results_pre[[model]] <- assess_calibration(
    scores_data[[pred_prob_col]],
    scores_data$event_24,
    model,
    pre_recal = TRUE
  )
}

# Create pre-recalibration plot
plot_data_pre <- data.frame()

for(model in model_names) {
  df <- standardize_cal_data(cal_results_pre[[model]]$cal_data, model)
  plot_data_pre <- rbind(plot_data_pre, df)
}

p_pre <- ggplot(plot_data_pre, aes(x = pred_mean, y = obs_mean, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 1) +
  geom_errorbar(aes(ymin = obs_lower, ymax = obs_upper), width = 0.02, alpha = 0.4) +
  geom_point(aes(size = n), alpha = 0.7) +
  scale_color_manual(values = model_colors) +
  scale_size_continuous(name = "N patients", range = c(2, 8)) +
  labs(
    title = "Pre-Recalibration Calibration at 24 Months",
    x = "Mean Predicted Probability",
    y = "Observed Event Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

for(model in c("Cox_TVC", "Cox_EN", "RSF")) {
  if(model %in% unique(plot_data_pre$Model)) {
    model_data <- plot_data_pre[plot_data_pre$Model == model,]
    if(nrow(model_data) > 3) {
      p_pre <- p_pre +
        geom_smooth(data = model_data, method = "loess", se = TRUE,
                   span = 0.8, alpha = 0.2, show.legend = FALSE)
    }
  }
}

for(model in c("HCTCI", "ELN")) {
  if(model %in% unique(plot_data_pre$Model)) {
    model_data <- plot_data_pre[plot_data_pre$Model == model,]
    p_pre <- p_pre +
      geom_line(data = model_data, linewidth = 1.2)
  }
}

ggsave("model_comparison/Figures/PreCalibration_24months.png", p_pre,
       width = 10, height = 8, dpi = 300)

# Post-recalibration assessment at 24 months
cal_results_post <- list()

for(model in model_names) {
  recal_prob_col <- paste0("recal_prob_", model, "_24")
  cal_results_post[[model]] <- assess_calibration(
    scores_data[[recal_prob_col]],
    scores_data$event_24,
    model,
    pre_recal = FALSE
  )
}

# Create post-recalibration plot
plot_data_post <- data.frame()

for(model in model_names) {
  df <- standardize_cal_data(cal_results_post[[model]]$cal_data, model)
  plot_data_post <- rbind(plot_data_post, df)
}

p_post <- ggplot(plot_data_post, aes(x = pred_mean, y = obs_mean, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 1) +
  geom_errorbar(aes(ymin = obs_lower, ymax = obs_upper), width = 0.02, alpha = 0.4) +
  geom_point(aes(size = n), alpha = 0.7) +
  scale_color_manual(values = model_colors) +
  scale_size_continuous(name = "N patients", range = c(2, 8)) +
  labs(
    title = "Post-Recalibration Calibration at 24 Months",
    x = "Mean Predicted Probability",
    y = "Observed Event Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))

for(model in c("Cox_TVC", "Cox_EN", "RSF")) {
  if(model %in% unique(plot_data_post$Model)) {
    model_data <- plot_data_post[plot_data_post$Model == model,]
    if(nrow(model_data) > 3) {
      p_post <- p_post +
        geom_smooth(data = model_data, method = "loess", se = TRUE,
                   span = 0.8, alpha = 0.2, show.legend = FALSE)
    }
  }
}

for(model in c("HCTCI", "ELN")) {
  if(model %in% unique(plot_data_post$Model)) {
    model_data <- plot_data_post[plot_data_post$Model == model,]
    p_post <- p_post +
      geom_line(data = model_data, linewidth = 1.2)
  }
}

ggsave("model_comparison/Figures/PostCalibration_24months.png", p_post,
       width = 10, height = 8, dpi = 300)

# ====================================================================
# PART 8: BRIER SCORES WITH IPCW
# ====================================================================

cens_fit <- survfit(Surv(Time, 1-Status) ~ 1, data = scores_data)

calc_ipcw_brier <- function(time, status, event_prob, t_point, cens_fit) {
  n <- length(time)
  weights <- rep(0, n)
  
  event_before_t <- time <= t_point & status == 1
  
  if(sum(event_before_t) > 0) {
    cens_prob_event <- sapply(time[event_before_t], function(ti) {
      idx <- findInterval(ti, cens_fit$time)
      if(idx == 0) return(1)
      return(cens_fit$surv[idx])
    })
    weights[event_before_t] <- 1 / pmax(cens_prob_event, 0.01)
  }
  
  after_t <- time > t_point
  
  if(sum(after_t) > 0) {
    idx_t <- findInterval(t_point, cens_fit$time)
    if(idx_t > 0) {
      cens_prob_t <- cens_fit$surv[idx_t]
      weights[after_t] <- 1 / max(cens_prob_t, 0.01)
    }
  }
  
  if(sum(weights > 0) > 0) {
    weights <- weights / mean(weights[weights > 0])
  }
  
  obs <- ifelse(time <= t_point & status == 1, 1, 0)
  valid_idx <- weights > 0
  
  if(sum(valid_idx) > 0) {
    individual_brier <- weights[valid_idx] * (obs[valid_idx] - event_prob[valid_idx])^2
    brier <- sum(individual_brier) / sum(weights[valid_idx])
    
    p0 <- sum(weights[valid_idx] * obs[valid_idx]) / sum(weights[valid_idx])
    brier_null <- p0 * (1 - p0)
    ipa <- 100 * (1 - brier / brier_null)
    
    return(list(brier = brier, ipa = ipa))
  } else {
    return(list(brier = NA, ipa = NA))
  }
}

brier_results <- data.frame()
brier_results_lower <- data.frame()
brier_results_upper <- data.frame()

for(t in c(12, 24, 36)) {
  temp_brier <- data.frame(Time = t)
  temp_brier_lower <- data.frame(Time = t)
  temp_brier_upper <- data.frame(Time = t)
  
  for(model in model_names) {
    recal_prob_col <- paste0("recal_prob_", model, "_", t)
    event_prob <- scores_data[[recal_prob_col]]
    
    brier_calc <- calc_ipcw_brier(scores_data$Time, scores_data$Status,
                                 event_prob, t, cens_fit)
    
    temp_brier[[model]] <- brier_calc$brier
    temp_brier[[paste0(model, "_IPA")]] <- brier_calc$ipa
    
    # Bootstrap for confidence intervals
    set.seed(147201)
    n_boot <- 1000
    boot_brier <- numeric(n_boot)
    boot_ipa <- numeric(n_boot)
    
    for(b in 1:n_boot) {
      boot_idx <- sample(nrow(scores_data), replace = TRUE)
      boot_prob <- event_prob[boot_idx]
      boot_time <- scores_data$Time[boot_idx]
      boot_status <- scores_data$Status[boot_idx]
      
      boot_cens_fit <- survfit(Surv(boot_time, 1-boot_status) ~ 1)
      boot_calc <- calc_ipcw_brier(boot_time, boot_status, boot_prob, t, boot_cens_fit)
      
      boot_brier[b] <- boot_calc$brier
      boot_ipa[b] <- boot_calc$ipa
    }
    
    ci_brier <- quantile(boot_brier, c(0.025, 0.975), na.rm = TRUE)
    ci_ipa <- quantile(boot_ipa, c(0.025, 0.975), na.rm = TRUE)
    
    temp_brier_lower[[model]] <- ci_brier[1]
    temp_brier_upper[[model]] <- ci_brier[2]
    temp_brier_lower[[paste0(model, "_IPA")]] <- ci_ipa[1]
    temp_brier_upper[[paste0(model, "_IPA")]] <- ci_ipa[2]
  }
  
  brier_results <- rbind(brier_results, temp_brier)
  
  if(nrow(brier_results_lower) == 0) {
    brier_results_lower <- temp_brier_lower
    brier_results_upper <- temp_brier_upper
  } else {
    brier_results_lower <- rbind(brier_results_lower, temp_brier_lower)
    brier_results_upper <- rbind(brier_results_upper, temp_brier_upper)
  }
}

write.csv(brier_results, "model_comparison/Tables/brier_results.csv", row.names = FALSE)
write.csv(brier_results_lower, "model_comparison/Tables/brier_results_lower.csv", row.names = FALSE)
write.csv(brier_results_upper, "model_comparison/Tables/brier_results_upper.csv", row.names = FALSE)

# ====================================================================
# PART 9: DECISION CURVE ANALYSIS WITH IPCW
# ====================================================================

calc_ipcw_dca <- function(time, status, pred_prob, t_point, thresh, cens_fit) {
  n <- length(time)
  weights <- rep(0, n)
  
  # Weight events before t
  event_before_t <- time <= t_point & status == 1
  
  if(sum(event_before_t) > 0) {
    cens_prob_event <- sapply(time[event_before_t], function(ti) {
      idx <- findInterval(ti, cens_fit$time)
      if(idx == 0) return(1)
      return(cens_fit$surv[idx])
    })
    weights[event_before_t] <- 1 / pmax(cens_prob_event, 0.01)
  }
  
  # Weight survivors past t
  after_t <- time > t_point
  
  if(sum(after_t) > 0) {
    idx_t <- findInterval(t_point, cens_fit$time)
    if(idx_t > 0) {
      cens_prob_t <- cens_fit$surv[idx_t]
      weights[after_t] <- 1 / max(cens_prob_t, 0.01)
    }
  }
  
  # Normalize weights
  if(sum(weights > 0) > 0) {
    weights <- weights / mean(weights[weights > 0])
  }
  
  # Calculate weighted TP and FP
  positive_pred <- pred_prob >= thresh
  obs_event <- ifelse(time <= t_point & status == 1, 1, 0)
  valid_idx <- weights > 0
  
  if(sum(valid_idx) > 0) {
    weighted_tp <- sum(weights[valid_idx] *
                      positive_pred[valid_idx] * obs_event[valid_idx]) /
                   sum(weights[valid_idx])
    
    weighted_fp <- sum(weights[valid_idx] *
                      positive_pred[valid_idx] * (1 - obs_event[valid_idx])) /
                   sum(weights[valid_idx])
    
    net_benefit <- weighted_tp - weighted_fp * (thresh / (1 - thresh))
    return(net_benefit)
  } else {
    return(NA)
  }
}

threshold_range <- seq(0.10, 0.45, by = 0.05)

for(t_dca in c(24, 36)) {
  dca_results <- data.frame()
  
  for(thresh in threshold_range) {
    for(model in model_names) {
      recal_prob_col <- paste0("recal_prob_", model, "_", t_dca)
      pred_prob <- scores_data[[recal_prob_col]]
      
      net_benefit <- calc_ipcw_dca(scores_data$Time, scores_data$Status,
                                  pred_prob, t_dca, thresh, cens_fit)
      
      dca_results <- rbind(dca_results, data.frame(
        Model = model,
        Threshold = thresh,
        NetBenefit = net_benefit,
        Time = t_dca
      ))
    }
    
    # Calculate treat all and treat none strategies
    event_weights <- rep(0, nrow(scores_data))
    event_before_t <- scores_data$Time <= t_dca & scores_data$Status == 1
    
    if(sum(event_before_t) > 0) {
      cens_prob_event <- sapply(scores_data$Time[event_before_t], function(ti) {
        idx <- findInterval(ti, cens_fit$time)
        if(idx == 0) return(1)
        return(cens_fit$surv[idx])
      })
      event_weights[event_before_t] <- 1 / pmax(cens_prob_event, 0.01)
    }
    
    after_t <- scores_data$Time > t_dca
    
    if(sum(after_t) > 0) {
      idx_t <- findInterval(t_dca, cens_fit$time)
      if(idx_t > 0) {
        cens_prob_t <- cens_fit$surv[idx_t]
        event_weights[after_t] <- 1 / max(cens_prob_t, 0.01)
      }
    }
    
    if(sum(event_weights > 0) > 0) {
      event_weights <- event_weights / mean(event_weights[event_weights > 0])
    }
    
    obs_event <- ifelse(scores_data$Time <= t_dca & scores_data$Status == 1, 1, 0)
    valid_idx <- event_weights > 0
    
    if(sum(valid_idx) > 0) {
      weighted_event_rate <- sum(event_weights[valid_idx] * obs_event[valid_idx]) /
                            sum(event_weights[valid_idx])
      net_benefit_all <- weighted_event_rate - (1 - weighted_event_rate) * (thresh / (1 - thresh))
    } else {
      net_benefit_all <- NA
    }
    
    dca_results <- rbind(dca_results,
                        data.frame(Model = "Treat All", Threshold = thresh,
                                  NetBenefit = net_benefit_all, Time = t_dca),
                        data.frame(Model = "Treat None", Threshold = thresh,
                                  NetBenefit = 0, Time = t_dca))
  }
  
  p_dca <- ggplot(dca_results[dca_results$Model != "Treat All" &
                             dca_results$Model != "Treat None",],
                  aes(x = Threshold, y = NetBenefit, color = Model)) +
    geom_line(linewidth = 1.2) +
    geom_line(data = dca_results[dca_results$Model == "Treat All",],
             color = "gray50", linetype = "dashed", linewidth = 1) +
    geom_line(data = dca_results[dca_results$Model == "Treat None",],
             color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray70") +
    labs(title = sprintf("Decision Curve Analysis at %d Months", t_dca),
         x = "Threshold Probability", y = "Net Benefit") +
    theme_minimal() +
    scale_color_manual(values = model_colors)
  
  ggsave(sprintf("model_comparison/Figures/DCA_%dmonths.png", t_dca), p_dca,
         width = 10, height = 6, dpi = 300)
  
  write.csv(dca_results,
           sprintf("model_comparison/Tables/dca_results_%dm.csv", t_dca),
           row.names = FALSE)
}

# ====================================================================
# PART 10: TIME-TO-EVENT NRI WITH IPCW
# ====================================================================

calc_ipcw_nri <- function(ref_name, new_name, time_point = 24, n_boot = 1000) {
  ref_recal_col <- paste0("recal_prob_", ref_name, "_", time_point)
  new_recal_col <- paste0("recal_prob_", new_name, "_", time_point)
  
  ref_pred <- scores_data[[ref_recal_col]]
  new_pred <- scores_data[[new_recal_col]]
  
  # IPCW weights
  weights <- rep(0, nrow(scores_data))
  event_before_t <- scores_data$Time <= time_point & scores_data$Status == 1
  
  if(sum(event_before_t) > 0) {
    cens_prob_event <- sapply(scores_data$Time[event_before_t], function(ti) {
      idx <- findInterval(ti, cens_fit$time)
      if(idx == 0) return(1)
      return(cens_fit$surv[idx])
    })
    weights[event_before_t] <- 1 / pmax(cens_prob_event, 0.01)
  }
  
  after_t <- scores_data$Time > time_point
  
  if(sum(after_t) > 0) {
    idx_t <- findInterval(time_point, cens_fit$time)
    if(idx_t > 0) {
      cens_prob_t <- cens_fit$surv[idx_t]
      weights[after_t] <- 1 / max(cens_prob_t, 0.01)
    }
  }
  
  if(sum(weights > 0) > 0) {
    weights <- weights / mean(weights[weights > 0])
  }
  
  # Define events
  event <- scores_data$Time <= time_point & scores_data$Status == 1
  valid_idx <- weights > 0
  
  # Weighted NRI calculations
  if(sum(valid_idx & event) > 0 && sum(valid_idx & !event) > 0) {
    event_up <- sum(weights[valid_idx & event] * (new_pred[valid_idx & event] > ref_pred[valid_idx & event]))
    event_down <- sum(weights[valid_idx & event] * (new_pred[valid_idx & event] < ref_pred[valid_idx & event]))
    nonevent_up <- sum(weights[valid_idx & !event] * (new_pred[valid_idx & !event] > ref_pred[valid_idx & !event]))
    nonevent_down <- sum(weights[valid_idx & !event] * (new_pred[valid_idx & !event] < ref_pred[valid_idx & !event]))
    
    total_event_weight <- sum(weights[valid_idx & event])
    total_nonevent_weight <- sum(weights[valid_idx & !event])
    
    nri_event <- (event_up - event_down) / total_event_weight
    nri_nonevent <- (nonevent_down - nonevent_up) / total_nonevent_weight
    nri_total <- nri_event + nri_nonevent
    
  } else {
    nri_total <- nri_event <- nri_nonevent <- NA
  }
  
  # Bootstrap for confidence intervals
  set.seed(147201)
  boot_nri <- replicate(n_boot, {
    idx <- sample(nrow(scores_data), replace = TRUE)
    ref_b <- ref_pred[idx]
    new_b <- new_pred[idx]
    time_b <- scores_data$Time[idx]
    status_b <- scores_data$Status[idx]
    
    # Recalculate weights for bootstrap sample
    boot_cens_fit <- survfit(Surv(time_b, 1-status_b) ~ 1)
    boot_weights <- rep(0, length(time_b))
    boot_event_before_t <- time_b <= time_point & status_b == 1
    
    if(sum(boot_event_before_t) > 0) {
      boot_cens_prob_event <- sapply(time_b[boot_event_before_t], function(ti) {
        boot_idx <- findInterval(ti, boot_cens_fit$time)
        if(boot_idx == 0) return(1)
        return(boot_cens_fit$surv[boot_idx])
      })
      boot_weights[boot_event_before_t] <- 1 / pmax(boot_cens_prob_event, 0.01)
    }
    
    boot_after_t <- time_b > time_point
    
    if(sum(boot_after_t) > 0) {
      boot_idx_t <- findInterval(time_point, boot_cens_fit$time)
      if(boot_idx_t > 0) {
        boot_cens_prob_t <- boot_cens_fit$surv[boot_idx_t]
        boot_weights[boot_after_t] <- 1 / max(boot_cens_prob_t, 0.01)
      }
    }
    
    if(sum(boot_weights > 0) > 0) {
      boot_weights <- boot_weights / mean(boot_weights[boot_weights > 0])
    }
    
    boot_event <- time_b <= time_point & status_b == 1
    boot_valid_idx <- boot_weights > 0
    
    if(sum(boot_valid_idx & boot_event) > 0 && sum(boot_valid_idx & !boot_event) > 0) {
      boot_event_up <- sum(boot_weights[boot_valid_idx & boot_event] *
                          (new_b[boot_valid_idx & boot_event] > ref_b[boot_valid_idx & boot_event]))
      boot_event_down <- sum(boot_weights[boot_valid_idx & boot_event] *
                            (new_b[boot_valid_idx & boot_event] < ref_b[boot_valid_idx & boot_event]))
      boot_nonevent_up <- sum(boot_weights[boot_valid_idx & !boot_event] *
                             (new_b[boot_valid_idx & !boot_event] > ref_b[boot_valid_idx & !boot_event]))
      boot_nonevent_down <- sum(boot_weights[boot_valid_idx & !boot_event] *
                               (new_b[boot_valid_idx & !boot_event] < ref_b[boot_valid_idx & !boot_event]))
      
      boot_total_event_weight <- sum(boot_weights[boot_valid_idx & boot_event])
      boot_total_nonevent_weight <- sum(boot_weights[boot_valid_idx & !boot_event])
      
      boot_nri_event <- (boot_event_up - boot_event_down) / boot_total_event_weight
      boot_nri_nonevent <- (boot_nonevent_down - boot_nonevent_up) / boot_total_nonevent_weight
      
      boot_nri_event + boot_nri_nonevent
    } else {
      NA
    }
  })
  
  boot_nri <- boot_nri[!is.na(boot_nri)]
  
  if(length(boot_nri) > 0) {
    nri_ci <- quantile(boot_nri, c(0.025, 0.975))
    nri_se <- sd(boot_nri)
    z_stat <- nri_total / nri_se
    p_value <- 2 * pnorm(-abs(z_stat))
  } else {
    nri_ci <- c(NA, NA)
    nri_se <- NA
    p_value <- NA
  }
  
  return(data.frame(
    comparison = paste(new_name, "vs", ref_name),
    NRI_total = nri_total,
    NRI_event = nri_event,
    NRI_nonevent = nri_nonevent,
    CI_lower = nri_ci[1],
    CI_upper = nri_ci[2],
    SE = nri_se,
    p_value = p_value,
    time_point = time_point
  ))
}

nri_results <- list()
nri_idx <- 1

for(t_nri in c(24, 36)) {
  nri_results[[nri_idx]] <- calc_ipcw_nri("HCTCI", "Cox_TVC", t_nri)
  nri_idx <- nri_idx + 1
  
  nri_results[[nri_idx]] <- calc_ipcw_nri("HCTCI", "Cox_EN", t_nri)
  nri_idx <- nri_idx + 1
  
  nri_results[[nri_idx]] <- calc_ipcw_nri("ELN", "Cox_TVC", t_nri)
  nri_idx <- nri_idx + 1
  
  nri_results[[nri_idx]] <- calc_ipcw_nri("ELN", "Cox_EN", t_nri)
  nri_idx <- nri_idx + 1
  
  if(rsf_available) {
    nri_results[[nri_idx]] <- calc_ipcw_nri("HCTCI", "RSF", t_nri)
    nri_idx <- nri_idx + 1
    
    nri_results[[nri_idx]] <- calc_ipcw_nri("ELN", "RSF", t_nri)
    nri_idx <- nri_idx + 1
  }
}

nri_df <- do.call(rbind, nri_results)
nri_df[, 2:7] <- round(nri_df[, 2:7], 3)
nri_df$p_value <- round(nri_df$p_value, 4)

write.csv(nri_df, "model_comparison/Tables/nri_results.csv", row.names = FALSE)

# ====================================================================
# PART 11: EXCEL REPORT
# ====================================================================

wb <- createWorkbook()

addWorksheet(wb, "Summary")
summary_data <- data.frame(
  Analysis = c("Total Patients", "Events", "Event Rate", "Median Follow-up",
              "Models Compared", "RSF Available", "Analysis Date"),
  Value = c(nrow(scores_data), sum(scores_data$Status),
           sprintf("%.1f%%", 100 * mean(scores_data$Status)),
           sprintf("%.1f months", median(scores_data$Time)),
           paste(model_names, collapse = ", "),
           ifelse(rsf_available, "Yes", "No"),
           format(Sys.time(), "%Y-%m-%d"))
)
writeData(wb, "Summary", summary_data)

addWorksheet(wb, "C-index")
writeData(wb, "C-index", cindex_df)

addWorksheet(wb, "AUC")
writeData(wb, "AUC", auc_results)

addWorksheet(wb, "Brier")
writeData(wb, "Brier", brier_results)

addWorksheet(wb, "NRI")
writeData(wb, "NRI", nri_df)

# Add calibration results
cal_summary <- data.frame()

for(model in model_names) {
  pre_cal <- cal_results_pre[[model]]
  post_cal <- cal_results_post[[model]]
  
  cal_summary <- rbind(cal_summary, data.frame(
    Model = model,
    Phase = "Pre-recalibration",
    Mean_Predicted = round(pre_cal$mean_predicted, 3),
    Mean_Observed = round(pre_cal$mean_observed, 3),
    Calibration_Intercept = round(pre_cal$calibration_intercept, 3),
    Calibration_Slope = round(pre_cal$calibration_slope, 3),
    E_O_Ratio = round(pre_cal$E_O_ratio, 3)
  ))
  
  cal_summary <- rbind(cal_summary, data.frame(
    Model = model,
    Phase = "Post-recalibration",
    Mean_Predicted = round(post_cal$mean_predicted, 3),
    Mean_Observed = round(post_cal$mean_observed, 3),
    Calibration_Intercept = round(post_cal$calibration_intercept, 3),
    Calibration_Slope = round(post_cal$calibration_slope, 3),
    E_O_Ratio = round(post_cal$E_O_ratio, 3)
  ))
}

addWorksheet(wb, "Calibration")
writeData(wb, "Calibration", cal_summary)

saveWorkbook(wb, "model_comparison/Tables/Complete_Model_Comparison_Results.xlsx",
            overwrite = TRUE)

# ====================================================================
# SAVE ANALYSIS
# ====================================================================

save.image(file = "model_comparison/complete_analysis.RData")

# Close logging
sink()
sink(type = "message")
close(log_file)

cat("\n===============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Results saved to model_comparison/\n")
cat("===============================================================\n")

