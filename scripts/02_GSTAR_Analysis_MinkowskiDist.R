# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 03: GSTAR Modeling, Performance, and Forecasting with Minkowski Distance
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================

# --- 1. Load Libraries ---
if (!require("gstar")) install.packages("gstar")
if (!require("xts")) install.packages("xts")
library(gstar)
library(xts)

# --- 2. Data Preparation ---
df <- read.csv(file.choose(), check.names = FALSE)

df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak",
                  "Stasiun Meteorologi Perak I",
                  "Stasiun Meteorologi Juanda")]

x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format="%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

s       <- round(nrow(x) * 0.8)
x_train <- x[1:s, ]
x_test  <- x[-c(1:s), ]

# --- 3. Minkowski Distance Function ---
minkowski_dist <- function(coords, p) {
  n     <- nrow(coords)
  d_mat <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        diff <- abs(coords[i, ] - coords[j, ])
        if (is.infinite(p)) {
          d_mat[i, j] <- max(diff)
        } else {
          d_mat[i, j] <- sum(diff^p)^(1 / p)
        }
      }
    }
  }
  return(d_mat)
}

# Coordinates (lon, lat)
coords <- matrix(c(
  112.7297, -7.1927,
  112.7274, -7.2000,
  112.7870, -7.3799
), ncol = 2, byrow = TRUE)

# --- 4. Helper: Compute R-squared ---
compute_r2 <- function(actual, predicted) {
  valid  <- complete.cases(actual, predicted)
  ss_res <- sum((actual[valid] - predicted[valid])^2)
  ss_tot <- sum((actual[valid] - mean(actual[valid]))^2)
  return(1 - ss_res / ss_tot)
}

# --- 5. Helper: Compute RMSE ---
compute_rmse <- function(actual, predicted) {
  valid <- complete.cases(actual, predicted)
  return(sqrt(mean((actual[valid] - predicted[valid])^2)))
}

# --- 6. Helper: Compute MAE ---
compute_mae <- function(actual, predicted) {
  valid <- complete.cases(actual, predicted)
  return(mean(abs(actual[valid] - predicted[valid])))
}

# --- 7. Helper: Compute MAPE ---
compute_mape <- function(actual, predicted) {
  valid <- complete.cases(actual, predicted)
  return(mean(abs((actual[valid] - predicted[valid]) / actual[valid])) * 100)
}

# --- 8. Loop p = 1 to 10: Build Table ---
cat("=============================================================\n")
cat(" Minkowski Distance Comparison: p = 1 to 10\n")
cat("=============================================================\n")

p_values   <- 1:10
stations   <- colnames(x_train)
all_results <- list()

for (p_val in p_values) {

  # Build weight matrix
  dist_mat <- minkowski_dist(coords, p_val)
  w        <- 1 / dist_mat
  diag(w)  <- 0
  weight   <- w / rowSums(w)

  # Fit GSTAR model
  model <- tryCatch(
    gstar(x = x_train, weight = weight, p = 1, d = 0, est = "OLS"),
    error = function(e) NULL
  )

  if (is.null(model)) {
    cat(sprintf("p = %2d : Model failed to converge.\n", p_val))
    next
  }

  # In-sample metrics
  fitted_vals <- fitted(model)
  in_r2   <- numeric(length(stations))
  in_rmse <- numeric(length(stations))
  in_mae  <- numeric(length(stations))
  in_mape <- numeric(length(stations))

  for (k in seq_along(stations)) {
    col       <- stations[k]
    actual    <- as.numeric(x_train[, col])
    predicted <- as.numeric(fitted_vals[, col])
    n_min     <- min(length(actual), length(predicted))
    actual    <- tail(actual,    n_min)
    predicted <- tail(predicted, n_min)

    in_r2[k]   <- compute_r2(actual, predicted)
    in_rmse[k] <- compute_rmse(actual, predicted)
    in_mae[k]  <- compute_mae(actual, predicted)
    in_mape[k] <- compute_mape(actual, predicted)
  }

  # Out-of-sample metrics
  pred_test <- tryCatch(
    predict(model, n = nrow(x_test)),
    error = function(e) NULL
  )

  out_r2   <- numeric(length(stations))
  out_rmse <- numeric(length(stations))
  out_mae  <- numeric(length(stations))
  out_mape <- numeric(length(stations))

  if (!is.null(pred_test)) {
    for (k in seq_along(stations)) {
      col       <- stations[k]
      actual    <- as.numeric(x_test[, col])
      predicted <- as.numeric(pred_test[1:nrow(x_test), col])

      out_r2[k]   <- compute_r2(actual, predicted)
      out_rmse[k] <- compute_rmse(actual, predicted)
      out_mae[k]  <- compute_mae(actual, predicted)
      out_mape[k] <- compute_mape(actual, predicted)
    }
  }

  # Store results
  all_results[[as.character(p_val)]] <- list(
    p        = p_val,
    model    = model,
    in_r2    = in_r2,
    in_rmse  = in_rmse,
    in_mae   = in_mae,
    in_mape  = in_mape,
    out_r2   = out_r2,
    out_rmse = out_rmse,
    out_mae  = out_mae,
    out_mape = out_mape
  )
}

# --- 9. Print Per-Station Tables ---
for (station_idx in seq_along(stations)) {
  st <- stations[station_idx]
  cat(sprintf("\n=============================================================\n"))
  cat(sprintf(" Station: %s\n", st))
  cat(sprintf("=============================================================\n"))
  cat(sprintf("%-4s | %-10s %-8s %-8s %-8s | %-10s %-8s %-8s %-8s\n",
              "p",
              "In-R²", "In-RMSE", "In-MAE", "In-MAPE",
              "Out-R²", "Out-RMSE", "Out-MAE", "Out-MAPE"))
  cat(strrep("-", 80), "\n")

  for (p_val in p_values) {
    r <- all_results[[as.character(p_val)]]
    if (is.null(r)) next
    cat(sprintf("%-4d | %-10.4f %-8.4f %-8.4f %-8.2f%% | %-10.4f %-8.4f %-8.4f %-8.2f%%\n",
                r$p,
                r$in_r2[station_idx],   r$in_rmse[station_idx],
                r$in_mae[station_idx],  r$in_mape[station_idx],
                r$out_r2[station_idx],  r$out_rmse[station_idx],
                r$out_mae[station_idx], r$out_mape[station_idx]))
  }
}

# --- 10. Summary Table: Average Across All Stations ---
cat("\n=============================================================\n")
cat(" Summary: Average Metrics Across All Stations\n")
cat("=============================================================\n")
cat(sprintf("%-4s | %-10s %-8s %-8s %-8s | %-10s %-8s %-8s %-8s\n",
            "p",
            "In-R²", "In-RMSE", "In-MAE", "In-MAPE",
            "Out-R²", "Out-RMSE", "Out-MAE", "Out-MAPE"))
cat(strrep("-", 80), "\n")

best_out_r2 <- -Inf
best_p      <- NA

for (p_val in p_values) {
  r <- all_results[[as.character(p_val)]]
  if (is.null(r)) next

  avg_in_r2    <- mean(r$in_r2)
  avg_in_rmse  <- mean(r$in_rmse)
  avg_in_mae   <- mean(r$in_mae)
  avg_in_mape  <- mean(r$in_mape)
  avg_out_r2   <- mean(r$out_r2)
  avg_out_rmse <- mean(r$out_rmse)
  avg_out_mae  <- mean(r$out_mae)
  avg_out_mape <- mean(r$out_mape)

  cat(sprintf("%-4d | %-10.4f %-8.4f %-8.4f %-8.2f%% | %-10.4f %-8.4f %-8.4f %-8.2f%%\n",
              p_val,
              avg_in_r2,  avg_in_rmse,  avg_in_mae,  avg_in_mape,
              avg_out_r2, avg_out_rmse, avg_out_mae, avg_out_mape))

  if (avg_out_r2 > best_out_r2) {
    best_out_r2 <- avg_out_r2
    best_p      <- p_val
  }
}

cat(strrep("-", 80), "\n")
cat(sprintf(" Best p (highest avg Out-R²): p = %d (Out-R² = %.4f)\n", best_p, best_out_r2))

# --- 11. Final Model with Best p ---
cat("\n=============================================================\n")
cat(sprintf(" Final Model: GSTAR with Minkowski p = %d\n", best_p))
cat("=============================================================\n")

best_dist   <- minkowski_dist(coords, best_p)
best_w      <- 1 / best_dist
diag(best_w) <- 0
best_weight <- best_w / rowSums(best_w)

best_model  <- gstar(x = x_train, weight = best_weight, p = 1, d = 0, est = "OLS")
print(summary(best_model))

# Forecast next 7 days
forecast_results <- predict(best_model, n = 7)
cat("\n--- 7-Day Forecast ---\n")
print(forecast_results)

# Plots
plot(best_model, testing = x_test)
plot(best_model, n_predict = 7)
