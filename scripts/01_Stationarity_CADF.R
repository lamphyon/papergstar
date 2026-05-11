# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 01: Stationarity Test using CADF (Pesaran, 2007)
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================

# --- 1. Preparation ---
# setwd("path/to/your/project") # Adjust this to your local directory
df <- read.csv("data/processed/data_bmkg_combined_clean_v2.csv", check.names = FALSE)

# Select relevant station columns
Data_Suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak", 
                    "Stasiun Meteorologi Perak I", 
                    "Stasiun Meteorologi Juanda")]

# Calculate Cross-Sectional Mean (Crucial for CADF)
CS_Mean <- rowMeans(Data_Suhu)

cat("\n=== HASIL UJI STASIONERITAS CADF PER LOKASI ===\n")

# --- 2. CADF Loop Calculation ---
# Testing for Unit Root with Cross-sectional Dependence
for (i in 1:ncol(Data_Suhu)) {
  y <- Data_Suhu[, i]
  
  # Define lags and differences
  y_lag <- c(NA, y[1:(length(y)-1)])
  dy <- c(NA, diff(y))
  cs_lag <- c(NA, CS_Mean[1:(length(CS_Mean)-1)])
  dcs <- c(NA, diff(CS_Mean))
  
  # Regression model for CADF estimation
  fit_cadf <- lm(dy ~ y_lag + cs_lag + dcs)
  
  # Extract t-statistic and p-value for the y_lag coefficient
  t_stat <- summary(fit_cadf)$coefficients[2, 3]
  p_val <- summary(fit_cadf)$coefficients[2, 4]
  
  cat("\nStasiun:", colnames(Data_Suhu)[i], "\n")
  cat("CADF t-statistic:", round(t_stat, 4), "\n")
  cat("P-value:", round(p_val, 4), "\n")
  
  if(p_val < 0.05) {
    cat("Kesimpulan: Stasioner (Tolak H0)\n")
  } else {
    cat("Kesimpulan: Tidak Stasioner (Gagal Tolak H0)\n")
  }
}
