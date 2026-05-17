# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 01: Stationarity Test using CADF (Pesaran, 2007)
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================

# --- 1. Preparation ---
# setwd("path/to/your/project") # Adjust this to your local directory
library(dynlm)

df <- read.csv(file.choose(), check.names = FALSE)

# Select relevant station columns
Data_Suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak", 
                    "Stasiun Meteorologi Perak I", 
                    "Stasiun Meteorologi Juanda")]

# Calculate Cross-Sectional Mean (Crucial for CADF)
CS_Mean <- rowMeans(Data_Suhu)

cat("\n=== HASIL UJI STASIONERITAS CADF PER LOKASI ===\n")

# --- 2. CADF Loop Calculation ---
# Testing for Unit Root with Cross-sectional Dependence
# Function to select optimal lag using information criteria
select_optimal_lag <- function(dy, y_lag, cs_lag, dcs, max_lag = 10, criterion = "AIC") {
  n <- length(dy)
  ic_values <- rep(NA, max_lag)
  
  for (lag in 1:max_lag) {
    # Build lagged difference matrix for dy
    lag_matrix <- matrix(NA, nrow = n, ncol = lag)
    for (j in 1:lag) {
      lag_matrix[(j+1):n, j] <- dy[1:(n-j)]
    }
    colnames(lag_matrix) <- paste0("dlag", 1:lag)
    
    df_temp <- data.frame(dy, y_lag, cs_lag, dcs, lag_matrix)
    df_temp <- df_temp[complete.cases(df_temp), ]
    
    fit_temp <- lm(dy ~ ., data = df_temp)
    
    ic_values[lag] <- if (criterion == "AIC") AIC(fit_temp) else BIC(fit_temp)
  }
  
  optimal_lag <- which.min(ic_values)
  return(list(optimal_lag = optimal_lag, ic_values = ic_values))
}

# ── Main Loop ──────────────────────────────────────────────────────────────────
max_lag    <- 10       # maximum lags to consider
criterion  <- "AIC"   # "AIC" or "BIC"

results <- data.frame(
  Stasiun    = character(),
  Optimal_Lag = integer(),
  T_Statistic = numeric(),
  P_Value     = numeric(),
  Kesimpulan  = character(),
  stringsAsFactors = FALSE
)

for (i in 1:ncol(Data_Suhu)) {
  y <- Data_Suhu[, i]
  
  # Base differences and lags
  y_lag <- c(NA, y[1:(length(y)-1)])
  dy    <- c(NA, diff(y))
  cs_lag <- c(NA, CS_Mean[1:(length(CS_Mean)-1)])
  dcs   <- c(NA, diff(CS_Mean))
  
  # ── Select optimal lag ──────────────────────────────────────────────────────
  lag_result  <- select_optimal_lag(dy, y_lag, cs_lag, dcs,
                                    max_lag = max_lag, criterion = criterion)
  optimal_lag <- lag_result$optimal_lag
  
  # ── Build augmented lag-difference matrix ───────────────────────────────────
  n           <- length(dy)
  lag_matrix  <- matrix(NA, nrow = n, ncol = optimal_lag)
  for (j in 1:optimal_lag) {
    lag_matrix[(j+1):n, j] <- dy[1:(n-j)]
  }
  colnames(lag_matrix) <- paste0("dlag", 1:optimal_lag)
  
  df_model <- data.frame(dy, y_lag, cs_lag, dcs, lag_matrix)
  df_model <- df_model[complete.cases(df_model), ]
  
  # ── Fit CADF model with optimal lag ────────────────────────────────────────
  fit_cadf <- lm(dy ~ ., data = df_model)
  
  # Extract t-statistic and p-value for y_lag
  t_stat <- summary(fit_cadf)$coefficients["y_lag", 3]
  p_val  <- summary(fit_cadf)$coefficients["y_lag", 4]
  
  kesimpulan <- if (p_val < 0.05) "Stasioner (Tolak H0)" else "Tidak Stasioner (Gagal Tolak H0)"
  
  cat("\n══════════════════════════════════════\n")
  cat("Stasiun     :", colnames(Data_Suhu)[i], "\n")
  cat("Optimal Lag :", optimal_lag, "(via", criterion, ")\n")
  cat("CADF t-stat :", round(t_stat, 4), "\n")
  cat("P-value     :", round(p_val,  4), "\n")
  cat("Kesimpulan  :", kesimpulan, "\n")
  
  results <- rbind(results, data.frame(
    Stasiun     = colnames(Data_Suhu)[i],
    Optimal_Lag = optimal_lag,
    T_Statistic = round(t_stat, 4),
    P_Value     = round(p_val,  4),
    Kesimpulan  = kesimpulan,
    stringsAsFactors = FALSE
  ))
}

cat("\n══════════════════════════════════════\n")
cat("SUMMARY RESULTS\n")
cat("══════════════════════════════════════\n")
print(results)
