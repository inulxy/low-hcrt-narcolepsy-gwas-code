setwd("/Users/lixinyu/Desktop/data")

library(readxl)  
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# 1. Prepare Data
raw_data <- read_excel("compare.xlsx", skip = 2, col_names = FALSE)
colnames(raw_data) <- c(
  "Chr", "Gene", "Pos", "Marker", "OtherA", "EffectA", 
  "NatureBeta", "NatureOR", "NatureSE",    
  "MetaBeta", "MetaOR_Raw", "MetaFreq", "MetaSE", "MetaP",  
  "EUR_Beta", "EUR_OR_Raw", "EUR_SE", "EUR_P",  
  "AFR_Beta", "AFR_OR_Raw", "AFR_SE", "AFR_P",  
  "EAS_Beta", "EAS_OR_Raw", "EAS_SE", "EAS_P"   
)
base_df <- raw_data %>%
  mutate(
    NatureBeta = suppressWarnings(as.numeric(NatureBeta)),
    NatureOR   = suppressWarnings(as.numeric(NatureOR)),
    NatureSE   = suppressWarnings(as.numeric(NatureSE)), 
    
    Gene_Display = ifelse(is.na(Gene) | Gene == "", "Unknown", Gene),
    Label = paste0(Marker, " (", Gene_Display, ")")
  ) %>%
  filter(Gene != "NAB1", !grepl("proxy", Gene, ignore.case = TRUE))

fit_constrained_deming <- function(data_df) {
  mean_var_nature <- mean(data_df$Nature_SE_OR^2, na.rm = TRUE)
  mean_var_my <- mean(data_df$My_SE_OR^2, na.rm = TRUE)
  lambda_val <- mean_var_my / mean_var_nature 
  
  x_centered <- data_df$Nature_OR - 1
  y_centered <- data_df$My_OR - 1
  
  Sxx <- sum(x_centered^2)
  Syy <- sum(y_centered^2)
  Sxy <- sum(x_centered * y_centered)
  
  deming_num <- (Syy - lambda_val * Sxx) + sqrt((Syy - lambda_val * Sxx)^2 + 4 * lambda_val * Sxy^2)
  deming_den <- 2 * Sxy
  slope_val <- deming_num / deming_den
  return(list(slope = slope_val, lambda = lambda_val))
}

generate_plots <- function(data, beta_col, se_col, study_name) {
  plot_df <- data %>%
    mutate(
      My_Beta = suppressWarnings(as.numeric(!!sym(beta_col))),
      My_SE   = suppressWarnings(as.numeric(!!sym(se_col)))
    ) %>%
    filter(!is.na(My_Beta) & !is.na(NatureOR)) %>%
    mutate(
      My_OR = exp(My_Beta),
      My_SE_OR = My_OR * My_SE,            
      Nature_SE_OR = NatureOR * NatureSE, 
      
      My_Lower = exp(My_Beta - 1.96 * My_SE),
      My_Upper = exp(My_Beta + 1.96 * My_SE),
      
      Nature_OR = NatureOR,
      NatureBeta_Use = ifelse(is.na(NatureBeta), log(NatureOR), NatureBeta),
      Nature_Lower = exp(NatureBeta_Use - 1.96 * NatureSE),
      Nature_Upper = exp(NatureBeta_Use + 1.96 * NatureSE)
    )
  
  # 2. Forest Plot 
  forest_df <- plot_df %>%
    select(Label, My_OR, My_Lower, My_Upper, Nature_OR, Nature_Lower, Nature_Upper) %>%
    pivot_longer(cols = c(My_OR, Nature_OR), names_to = "Study_Type", values_to = "OR") %>%
    mutate(
      Lower = ifelse(Study_Type == "My_OR", My_Lower, Nature_Lower), 
      Upper = ifelse(Study_Type == "My_OR", My_Upper, Nature_Upper),
      Study_Label = ifelse(Study_Type == "My_OR", paste0("Refined Cohort (", study_name, ")"), "Generalized Cohort (Baseline)")
    )
  
  # 优化 Forest Plot 的学术排版
  p1 <- ggplot(forest_df, aes(x = OR, y = Label, color = Study_Label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.3, position = position_dodge(0.6), linewidth = 0.8) +
    geom_point(size = 3.5, position = position_dodge(0.6)) +
    theme_bw(base_size = 14) + # 使用更适合正式论文的主题和字号
    scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5)) + 
    labs(x = "Odds Ratio (Log Scale)", y = "") +
    scale_color_manual(values = c("#2563eb", "#10b981")) + 
    theme(
      legend.position = "bottom", 
      legend.title = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.text.y = element_text(color = "black", face = "italic") 
    )
  
  ggsave(paste0("Forest_Plot_", study_name, ".pdf"), p1, width = 10, height = 7, dpi = 300)
  
  # 3. Deming Regression & Bootstrapping for CI/P-value
  fit_data <- plot_df 
  deming_res <- fit_constrained_deming(fit_data)
  slope_val <- deming_res$slope
  intercept_val <- 1 - slope_val
  r2_val <- cor(fit_data$Nature_OR, fit_data$My_OR)^2
  
  # Use Bootstrap to calculate 95% CI & P-value
  set.seed(145810) 
  n_boot <- 1000 
  boot_slopes <- numeric(n_boot)
  
  for(i in 1:n_boot) {
    boot_sample <- fit_data[sample(nrow(fit_data), replace = TRUE), ]
    tryCatch({
      boot_res <- fit_constrained_deming(boot_sample)
      boot_slopes[i] <- boot_res$slope
    }, error = function(e) {
      boot_slopes[i] <- NA 
    })
  }
  
  boot_slopes <- na.omit(boot_slopes)
  
  # 95% CI
  ci_lower <- quantile(boot_slopes, 0.025)
  ci_upper <- quantile(boot_slopes, 0.975)

  boot_se <- sd(boot_slopes)
  z_stat <- (slope_val - 1) / boot_se
  p_val <- 2 * pnorm(-abs(z_stat))
  p_val_format <- ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val))

  cat(sprintf("\n=== %s Deming Regression Results ===\n", study_name))
  cat(sprintf("Slope: %.3f\n", slope_val))
  cat(sprintf("95%% CI: [%.3f, %.3f]\n", ci_lower, ci_upper))
  cat(sprintf("P-value (vs 1): %s\n", p_val_format))
  
  stats_text <- sprintf("Deming Slope: %.3f\n95%% CI: %.2f - %.2f\nP-value: %s\nR² = %.3f", 
                        slope_val, ci_lower, ci_upper, p_val_format, r2_val)
  
  all_vals <- c(fit_data$My_OR, fit_data$Nature_OR) 
  min_val <- min(all_vals, na.rm = TRUE) * 0.95
  max_val <- max(all_vals, na.rm = TRUE) * 1.05

  p2 <- ggplot(plot_df, aes(x = Nature_OR, y = My_OR)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "darkgray", linewidth = 1) +
    geom_abline(slope = slope_val, intercept = intercept_val, color = "black", linewidth = 1) +
    
    geom_point(size = 3.5, color = "#2563eb", alpha = 0.8) + 
    geom_text_repel(aes(label = Label), size = 4, max.overlaps = 20, fontface = "italic") +
    
    theme_bw(base_size = 14) 
    coord_fixed(ratio = 1, xlim = c(min_val, max_val), ylim = c(min_val, max_val)) +
    
    labs(x = "Generalized Cohort OR", 
         y = paste0("Refined Cohort OR (", study_name, ")")) +
    
    annotate("text", x = -Inf, y = Inf, label = stats_text, 
             hjust = -0.1, vjust = 1.2, size = 4.5, fontface = "italic") +
    
    theme(
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black")
    )
  
  ggsave(paste0("Regression_Plot_", study_name, ".pdf"), p2, width = 8, height = 8, dpi = 300) 
}

# 4. Run the loop for all analyses
analysis_list <- list(
  "Meta-All" = c("MetaBeta", "MetaSE"),
  "EUR"      = c("EUR_Beta", "EUR_SE"),
  "AFR"      = c("AFR_Beta", "AFR_SE"),
  "EAS"      = c("EAS_Beta", "EAS_SE")
)

for (study in names(analysis_list)) {
  cols <- analysis_list[[study]]
  generate_plots(base_df, cols[1], cols[2], study)
}
