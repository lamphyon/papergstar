# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 02: GSTAR Modeling, Performance, and Forecasting
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# Change : d=1 (First-order differencing)
# =====================================================================

# --- 1. Load Libraries ---
if (!require("gstar")) install.packages("gstar")
if (!require("xts"))   install.packages("xts")
library(gstar)
library(xts)

# --- 2. Data Preparation ---
# setwd("path/to/your/project") # Adjust this to your local directory
df <- read.csv(file.choose(), check.names = FALSE)

# Reorder columns to match spatial coordinates
df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak",
                  "Stasiun Meteorologi Perak I",
                  "Stasiun Meteorologi Juanda")]

# Convert to XTS object
x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format = "%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

# --- 3. First-order Differencing (d = 1) ---
# Remove first row (NA from diff) and store last level for back-transformation
x_diff <- diff(x, differences = 1)
x_diff <- x_diff[!is.na(x_diff[, 1]), ]   # drop NA row

# Store the last observed LEVEL values (needed to undo differencing later)
last_levels_all <- x[nrow(x) - nrow(x_diff), ]  # level just before diff series

# --- 4. Train / Test Split (80:20) on Differenced Series ---
s        <- round(nrow(x_diff) * 0.8)
x_train  <- x_diff[1:s, ]
x_test   <- x_diff[-c(1:s), ]

# Keep the corresponding level values to undo differencing at test boundary
# Last level value before the test set begins (for back-transforming test predictions)
last_level_train <- x[s, ]   # level at the last training time point

# --- 5. Spatial Weight Matrix (Inverse-Distance) ---
# Coordinates (lon, lat)
coords <- matrix(c(
  112.7297, -7.1927,   # Tanjung Perak
  112.7274, -7.2000,   # Perak I
  112.7870, -7.3799    # Juanda
), ncol = 2, byrow = TRUE)

dist_mat <- as.matrix(dist(coords, method = "euclidean"))
w        <- 1 / dist_mat
diag(w)  <- 0
weight   <- w / rowSums(w)   # Row-normalisation

# --- 6. GSTAR Model Fitting on Differenced Data ---
# p=1 (Lag-1), d=1 (differencing already applied manually above)
model <- gstar(x      = x_train,
               weight = weight,
               p      = 1,
               d      = 0,       # d=0 here because series is already differenced
               est    = "OLS")

print(summary(model))

# --- 7. Accuracy Evaluation ---
cat("\n--- Performance Metrics (on differenced scale) ---\n")

# In-sample performance (differenced scale)
performance(model)

# R-squared in-sample (differenced scale)
fitted_vals <- fitted(model)
cat("\n--- R-squared In-sample (differenced scale) ---\n")
for (col in colnames(x_train)) {
  actual    <- as.numeric(x_train[, col])
  predicted <- as.numeric(fitted_vals[, col])
  n_min     <- min(length(actual), length(predicted))
  actual    <- tail(actual,    n_min)
  predicted <- tail(predicted, n_min)
  ss_res <- sum((actual - predicted)^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  cat(col, ": R² =", round(1 - ss_res / ss_tot, 4), "\n")
}

# Out-of-sample performance (differenced scale)
performance(model, x_test)

# Predict on test horizon (differenced scale)
pred_test_diff <- predict(model, n = nrow(x_test))

# --- Back-transform test predictions to original (level) scale ---
# Cumulative sum of differenced predictions + last training level
undiff <- function(diff_mat, start_levels) {
  # diff_mat   : matrix/xts of predicted differences (rows = time, cols = locations)
  # start_levels: named vector of last known levels (1 row)
  result <- matrix(NA, nrow = nrow(diff_mat), ncol = ncol(diff_mat))
  colnames(result) <- colnames(diff_mat)
  prev <- as.numeric(start_levels)
  for (i in seq_len(nrow(diff_mat))) {
    curr <- prev + as.numeric(diff_mat[i, ])
    result[i, ] <- curr
    prev <- curr
  }
  result
}

pred_test_level  <- undiff(pred_test_diff[1:nrow(x_test), ], last_level_train)
# Actual levels for test period
x_test_level     <- x[-c(1:s), ]            # original (un-differenced) test levels
# Align: x_test_level has one extra row at position s+1 corresponding to s level;
# actual test observations are rows (s+1):nrow(x) in the original xts
x_test_level_mat <- coredata(x_test_level)

cat("\n--- R-squared Out-of-sample (original / level scale) ---\n")
for (col in colnames(x_test_level)) {
  actual    <- as.numeric(x_test_level_mat[, col])
  predicted <- pred_test_level[, col]
  valid     <- complete.cases(actual, predicted)
  ss_res    <- sum((actual[valid] - predicted[valid])^2)
  ss_tot    <- sum((actual[valid] - mean(actual[valid]))^2)
  cat(col, ": R² =", round(1 - ss_res / ss_tot, 4), "\n")
}

# --- 8. Forecasting Next 7 Days ---
# Predict 7-step-ahead differences
forecast_diff <- predict(model, n = 7)
cat("\n--- 7-Day Forecast (differenced scale) ---\n")
print(forecast_diff)

# Back-transform forecast to original scale
last_level_all  <- x[nrow(x), ]    # very last observed level
forecast_level  <- undiff(forecast_diff, last_level_all)
forecast_dates  <- seq(index(x)[nrow(x)] + 1, by = "day", length.out = 7)
forecast_xts    <- xts(forecast_level, order.by = forecast_dates)
colnames(forecast_xts) <- colnames(x)

cat("\n--- 7-Day Forecast (original / level scale) ---\n")
print(forecast_xts)

# --- 9. Plots ---
plot(model, testing = x_test)   # Train vs Test (differenced scale)
plot(model, n_predict = 7)      # Historical vs Forecast (differenced scale)
