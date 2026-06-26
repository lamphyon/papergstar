# =====================================================================
# --- 1. Load Libraries ---
# =====================================================================
if (!require("xts")) install.packages("xts")
if (!require("geosphere")) install.packages("geosphere") 
if (!require("tseries")) install.packages("tseries") # Tambahan untuk ADF Test

library(xts)
library(geosphere)
library(tseries)

# =====================================================================
# --- 2. Data Preparation ---
# =====================================================================
df <- read.csv(file.choose(), check.names = FALSE) 

df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak", 
                  "Stasiun Meteorologi Perak I", 
                  "Stasiun Meteorologi Juanda")]

x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format="%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

# Bagi data Train dan Test (80:20)
s <- round(nrow(x) * 0.8)   
x_train <- x[1:s, ]
x_test  <- x[-c(1:s), ]

# Mean-Centering
mean_train <- colMeans(x_train)

x_train_centered <- x_train
for(i in 1:ncol(x_train)) {
  x_train_centered[, i] <- x_train[, i] - mean_train[i]
}

# =====================================================================
# --- 2.5. UJI STASIONERITAS (ADF TEST) ---
# =====================================================================
cat("\n============================================\n")
cat("   UJI STASIONERITAS (ADF TEST)             \n")
cat("============================================\n")

# Kita uji data training yang sudah di mean-centered
for(i in 1:ncol(x_train_centered)) {
  stasiun <- colnames(x_train_centered)[i]
  
  # Jalankan Augmented Dickey-Fuller Test
  adf_result <- adf.test(coredata(x_train_centered[, i]), alternative = "stationary")
  
  cat(sprintf("\nStasiun: %s\n", stasiun))
  cat(sprintf("Dickey-Fuller = %.4f | p-value = %.4f\n", adf_result$statistic, adf_result$p.value))
  
  # Cek nilai probabilitas (alpha = 0.05)
  if(adf_result$p.value < 0.05) {
    cat("Kesimpulan: Data STASIONER (Tolak H0)\n")
  } else {
    cat("Kesimpulan: Data TIDAK STASIONER (Terima H0)\n")
    cat("Peringatan: Data perlu didiferensiasi (d=1) sebelum lanjut ke model AR.\n")
  }
}

# =====================================================================
# --- 3. Spatial Weight Matrix (Cross-Correlation) ---
# =====================================================================
# Use training data only
corr_matrix <- cor(x_train_centered)

# Take absolute correlations
w_ccf <- abs(corr_matrix)

# Remove self-correlation
diag(w_ccf) <- 0

# Normalize rows so each row sums to 1
weight <- w_ccf / rowSums(w_ccf)

cat("\n============================================\n")
cat("--- Matriks Bobot Cross-Correlation ---\n")
print(round(weight, 4))


# =====================================================================
# --- 3.5 MACF DAN MPACF ---
# =====================================================================
cat("\n============================================\n")
cat("   IDENTIFIKASI ORDE (MACF & MPACF)         \n")
cat("============================================\n")

# Ekstrak coredata dari data train yang sudah di-centered
Z_train <- coredata(x_train_centered)

# 1. Plot MACF (Cek cut-off untuk orde MA)
par(mar=c(2, 2, 2, 2)) 
acf_model <- acf(Z_train, lag.max = 10, main = "Plot MACF (Matrix Autocorrelation Function)")

# 2. Plot MPACF (Cek cut-off untuk orde AR)
pacf_model <- pacf(Z_train, lag.max = 10, main = "Plot MPACF (Matrix Partial Autocorrelation Function)")

# Menampilkan nilai matriks MACF untuk lag 1 sampai 3 di console
cat("\n--- Nilai Matriks Autokorelasi Silang (Lag 1 hingga 3) ---\n")
for (lag_i in 1:3) {
  cat(sprintf("\nLag %d:\n", lag_i))
  print(round(acf_model$acf[lag_i + 1, , ], 4)) 
}

# =====================================================================
# --- 4. MANUAL GSTAR MODEL FITTING ---
# =====================================================================
Z_train <- coredata(x_train_centered)
N <- ncol(Z_train)
T_len <- nrow(Z_train)

# Buat variabel Target (Hari ini) dan Lag (Kemarin)
Z_t  <- Z_train[2:T_len, ]       # Observasi waktu t (Target)
Z_t1 <- Z_train[1:(T_len-1), ]   # Observasi waktu t-1 (Temporal Lag)

# Buat Spatial Lag (W * Z_t1)
V_t1 <- Z_t1 %*% t(weight)       

# List untuk menyimpan model OLS per stasiun
gstar_models <- list()
fitted_centered <- matrix(NA, nrow = T_len - 1, ncol = N)
colnames(fitted_centered) <- colnames(x)

cat("\n============================================\n")
cat("--- Ringkasan Koefisien Model GSTAR (Manual OLS) ---\n")
for (i in 1:N) {
  df_model <- data.frame(
    Y       = Z_t[, i],       
    X_Time  = Z_t1[, i],      
    X_Space = V_t1[, i]       
  )
  
  fit <- lm(Y ~ X_Time + X_Space - 1, data = df_model)
  gstar_models[[i]] <- fit
  fitted_centered[, i] <- fitted(fit)
  
  cat(sprintf("\nStasiun: %s\n", colnames(x)[i]))
  print(summary(fit)$coefficients)
}

# =====================================================================
# --- 5. ACCURACY EVALUATION ---
# =====================================================================
cat("\n============================================\n")
cat("   EVALUASI PERFORMA MODEL GSTAR            \n")
cat("============================================\n")

calc_metrics <- function(actual, predicted, label) {
  valid_idx <- complete.cases(actual, predicted)
  actual <- actual[valid_idx]
  predicted <- predicted[valid_idx]
  
  mse  <- mean((actual - predicted)^2)
  mape <- mean(abs((actual - predicted) / actual)) * 100
  ss_res <- sum((actual - predicted)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  r2 <- 1 - (ss_res / ss_tot)
  
  cat(sprintf("%-15s | MSE: %6.4f | MAPE: %5.2f%% | R²: %6.4f\n", label, mse, mape, r2))
}

cat("\n--- Performa In-Sample (Data Training) ---\n")
fitted_actual <- fitted_centered

for(i in 1:ncol(fitted_actual)) {
  fitted_actual[, i] <- fitted_actual[, i] + mean_train[i]
}

actual_train_eval <- coredata(x_train[2:T_len, ]) 

for (i in 1:N) {
  calc_metrics(actual_train_eval[, i], fitted_actual[, i], colnames(x)[i])
}

cat("\n--- Performa Out-of-Sample (Data Testing) ---\n")
Z_test <- coredata(x_test)
T_test <- nrow(Z_test)

Z_test_centered <- Z_test
for(i in 1:N) {
  Z_test_centered[, i] <- Z_test[, i] - mean_train[i]
}

last_train_val <- Z_train[T_len, ]
Z_all_centered <- rbind(last_train_val, Z_test_centered)

pred_test_centered <- matrix(NA, nrow = T_test, ncol = N)

for (t in 1:T_test) {
  Z_kemarin <- Z_all_centered[t, ]
  V_kemarin <- (Z_kemarin %*% t(weight))[1, ] 
  
  for (i in 1:N) {
    coefs <- coef(gstar_models[[i]])
    pred_test_centered[t, i] <- (coefs["X_Time"] * Z_kemarin[i]) + (coefs["X_Space"] * V_kemarin[i])
  }
}

pred_actual_test <- pred_test_centered
for(i in 1:N) {
  pred_actual_test[, i] <- pred_actual_test[, i] + mean_train[i]
}

for (i in 1:N) {
  calc_metrics(Z_test[, i], pred_actual_test[, i], colnames(x)[i])
}

# =====================================================================
# --- 6. FORECASTING & VISUALIZATION ---
# =====================================================================
cat("\n============================================\n")
cat("   FORECAST 7 HARI KEDEPAN                  \n")
cat("============================================\n")

x_centered <- coredata(x)
for(i in 1:N) x_centered[, i] <- x_centered[, i] - mean_train[i]

last_known_Z <- tail(x_centered, 1)[1, ]
forecast_centered <- matrix(NA, nrow = 7, ncol = N)

for (k in 1:7) {
  spat_lag <- (last_known_Z %*% t(weight))[1, ]
  
  next_Z <- numeric(N)
  for (i in 1:N) {
    coefs <- coef(gstar_models[[i]])
    next_Z[i] <- (coefs["X_Time"] * last_known_Z[i]) + (coefs["X_Space"] * spat_lag[i])
  }
  
  forecast_centered[k, ] <- next_Z
  last_known_Z <- next_Z 
}

forecast_actual <- forecast_centered
for(i in 1:N) {
  forecast_actual[, i] <- forecast_actual[, i] + mean_train[i]
}

last_date <- index(tail(x, 1))
future_dates <- seq(last_date + 1, by = "days", length.out = 7)

forecast_xts <- xts(forecast_actual, order.by = future_dates)
colnames(forecast_xts) <- colnames(x)

print(forecast_xts)
plot(forecast_xts, main = "Forecast Suhu 7 Hari ke Depan", legend.loc = "topright")
