# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 02: GSTAR Modeling, Performance, and Forecasting
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================

# --- 1. Load Libraries ---
if (!require("gstar")) install.packages("gstar")
if (!require("xts")) install.packages("xts")

library(gstar)
library(xts)

# --- 2. Data Preparation (Perbaikan Stasioneritas) ---
# ... (kode sebelumnya sama sampai pembentukan objek 'x')

# Lakukan differencing lag 1 agar data stasioner
x_diff <- diff(x, differences = 1)
x_diff <- na.omit(x_diff) # Hapus baris NA di awal akibat differencing

# Train / test split (80:20) MENGGUNAKAN data yang sudah di-diff
s <- round(nrow(x_diff) * 0.8)   
x_train <- x_diff[1:s, ]
x_test  <- x_diff[-c(1:s), ]


# --- 3. Spatial Weight Matrix (Perbaikan Jarak Haversine) ---
if (!require("geosphere")) install.packages("geosphere")
library(geosphere)

# Koordinat (Longitude, Latitude)
coords <- matrix(c(
  112.7297, -7.1927,   # Tanjung Perak
  112.7274, -7.2000,   # Perak I
  112.7870, -7.3799    # Juanda
), ncol = 2, byrow = TRUE)

# Hitung jarak sesungguhnya di bumi (output dalam meter)
dist_mat <- distm(coords, fun = distHaversine)
w        <- 1 / dist_mat
diag(w)  <- 0
weight   <- w / rowSums(w) # Row-normalization

# --- 4. GSTAR Model Fitting ---
# p=1 (Lag 1), d=0 (No differencing since data is stationary)
model <- gstar(x      = x_train,
               weight = weight,
               p      = 1,      
               d      = 0,      
               est    = "OLS")

print(summary(model))

# --- 5. Accuracy Evaluation ---
cat("\n--- Performance Metrics ---\n")
# In-sample performance
performance(model)

# R-squared (in-sample)
fitted_vals <- fitted(model)
cat("\n--- R-squared (In-sample) ---\n")
for (col in colnames(x_train)) {
  actual    <- as.numeric(x_train[, col])
  predicted <- as.numeric(fitted_vals[, col])
  
  # Trim to the shorter length (fitted may drop initial rows due to lags)
  n_min  <- min(length(actual), length(predicted))
  actual    <- tail(actual,    n_min)
  predicted <- tail(predicted, n_min)
  
  ss_res <- sum((actual - predicted)^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  cat(col, ": R² =", round(1 - ss_res / ss_tot, 4), "\n")
}

# Out-of-sample performance (Testing set)
performance(model, x_test)

# R-squared (out-of-sample)
pred_test <- predict(model, n = nrow(x_test))
cat("\n--- R-squared (Out-of-sample) ---\n")
for (col in colnames(x_test)) {
  actual    <- as.numeric(x_test[, col])
  predicted <- as.numeric(pred_test[1:nrow(x_test), col])
  valid     <- complete.cases(actual, predicted)
  ss_res <- sum((actual[valid] - predicted[valid])^2)
  ss_tot <- sum((actual[valid] - mean(actual[valid]))^2)
  cat(col, ": R² =", round(1 - ss_res / ss_tot, 4), "\n")
}

# --- 6. Forecasting & Visualization ---
# Forecast for the next 7 days
forecast_results <- predict(model, n = 7)
print(forecast_results)

# Plots
plot(model, testing = x_test) # Train vs Test
plot(model, n_predict = 7)    # Historical vs Forecast
