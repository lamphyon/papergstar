# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 02: GSTAR Modeling, Performance, and Forecasting
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================


# --- 1. Load Libraries ---
if (!require("gstar")) install.packages("gstar")
if (!require("xts"))   install.packages("xts")
library(gstar)
library(xts)


# --- 2. Data Preparation ---
# setwd("path/to/your/project")   # Adjust to your local directory
df <- read.csv(file.choose(), check.names = FALSE)

# Reorder columns to match spatial coordinates
df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak",
                  "Stasiun Meteorologi Perak I",
                  "Stasiun Meteorologi Juanda")]

# Convert to XTS object
x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format = "%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

# Train / test split (80:20)
s       <- round(nrow(x) * 0.8)
x_train <- x[1:s, ]
x_test  <- x[-c(1:s), ]


# --- 3. Spatial Weight Matrix (Inverse-Distance) ---
# Coordinates: (longitude, latitude)
coords <- matrix(c(
  112.7297, -7.1927,   # Z1 – Tanjung Perak
  112.7274, -7.2000,   # Z2 – Perak I
  112.7870, -7.3799    # Z3 – Juanda
), ncol = 2, byrow = TRUE)

dist_mat <- as.matrix(dist(coords, method = "euclidean"))
w        <- 1 / dist_mat
diag(w)  <- 0
weight   <- w / rowSums(w)   # Row-normalisation


# --- 4. GSTAR Model Fitting ---
# p = 1 (temporal lag 1), d = 0 (no differencing – data already stationary)
model <- gstar(x      = x_train,
               weight = weight,
               p      = 1,
               d      = 0,
               est    = "OLS")

print(summary(model))


# =====================================================================
# --- 5. Extract Z1 / Z2 / Z3 Parameters per GSTAR Equation ----------
# =====================================================================
# The gstar package stores OLS coefficients in model$model$B.
# Row names follow the pattern:
#   psi{lag}{spatial_order}_{location}
# For p = 1, three locations (Z1, Z2, Z3) we get six rows total:
#   psi10_TanjungPerak, psi10_PerakI, psi10_Juanda   <- phi (own-lag)
#   psi11_TanjungPerak, psi11_PerakI, psi11_Juanda   <- theta (spatial-lag)
# =====================================================================

B     <- model$model$B           # coefficient vector (rows = predictors)
p_val <- model$model$p_value      # matching p-values

loc_names <- colnames(x_train)   # c("TanjungPerak", "PerakI", "Juanda")
n_loc     <- length(loc_names)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  GSTAR(1;1) Parameter Estimates – One Equation per Location\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat(sprintf("%-20s  %-12s  %-12s  %-12s  %-12s\n",
            "Parameter", "Estimate", "Std.Error", "t-value", "p-value"))
cat(paste(rep("-", 70), collapse = ""), "\n")

# Derive standard errors from the stored SSE and design matrix
# (mirrors the gstar_est internal calculation)
z_vals <- model$model$z              # stacked response vector
z_hat  <- model$model$Xv %*% B
sse    <- sum((z_vals - z_hat)^2)
std_err <- sqrt(sse * diag(solve(t(model$model$Xv) %*% model$model$Xv)))
t_stat  <- B / std_err

# ---- Helper: print a single-equation block -------------------------
print_equation <- function(eq_name, phi_name, theta_name) {
  cat("\n  Equation:", eq_name, "\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")

  for (row_nm in c(phi_name, theta_name)) {
    idx <- which(rownames(B) == row_nm)
    if (length(idx) == 0) next

    param_label <- if (grepl("psi10", row_nm)) {
      paste0("phi_", gsub("psi10_", "", row_nm))   # own-lag (temporal)
    } else {
      paste0("theta_", gsub("psi11_", "", row_nm)) # spatial-lag
    }

    cat(sprintf("  %-18s  %12.6f  %12.6f  %12.4f  %12.6f\n",
                param_label,
                B[idx, 1],
                std_err[idx],
                t_stat[idx, 1],
                p_val[idx, 1]))
  }
}

# Print one block per location (Z1, Z2, Z3)
for (j in seq_along(loc_names)) {
  loc  <- loc_names[j]
  phi_nm   <- paste0("psi10_", loc)
  theta_nm <- paste0("psi11_", loc)
  print_equation(paste0("Z", j, " (", loc, ")"), phi_nm, theta_nm)
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")

# ---- Compact wide-format table (one row per location) --------------
cat("\n--- Compact Parameter Table (Z1 / Z2 / Z3) ---\n\n")
cat(sprintf("%-6s  %-18s  %-12s  %-12s\n",
            "Loc", "Station", "phi (own-lag)", "theta (spatial)"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (j in seq_along(loc_names)) {
  loc      <- loc_names[j]
  phi_idx  <- which(rownames(B) == paste0("psi10_", loc))
  thet_idx <- which(rownames(B) == paste0("psi11_", loc))

  phi_val  <- if (length(phi_idx))  round(B[phi_idx,  1], 6) else NA
  thet_val <- if (length(thet_idx)) round(B[thet_idx, 1], 6) else NA

  cat(sprintf("Z%-5s  %-18s  %12.6f  %12.6f\n",
              j, loc, phi_val, thet_val))
}
cat("\n")

# ---- Write extracted parameters to CSV for reporting ---------------
param_df <- data.frame(
  Location  = character(),
  Station   = character(),
  Parameter = character(),
  Estimate  = numeric(),
  Std.Error = numeric(),
  t.value   = numeric(),
  p.value   = numeric(),
  stringsAsFactors = FALSE
)

for (j in seq_along(loc_names)) {
  loc <- loc_names[j]
  for (type in c("psi10", "psi11")) {
    row_nm <- paste0(type, "_", loc)
    idx    <- which(rownames(B) == row_nm)
    if (!length(idx)) next
    param_df <- rbind(param_df, data.frame(
      Location  = paste0("Z", j),
      Station   = loc,
      Parameter = ifelse(type == "psi10",
                         paste0("phi_", loc),
                         paste0("theta_", loc)),
      Estimate  = round(B[idx, 1],       6),
      Std.Error = round(std_err[idx],    6),
      t.value   = round(t_stat[idx, 1],  4),
      p.value   = round(p_val[idx, 1],   6),
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(param_df, "GSTAR_parameters_Z1Z2Z3.csv", row.names = FALSE)
cat("Parameter table saved to: GSTAR_parameters_Z1Z2Z3.csv\n\n")


# --- 6. Accuracy Evaluation ---
cat("\n--- In-sample Performance ---\n")
performance(model)

# R-squared (in-sample)
fitted_vals <- fitted(model)
cat("\n--- R-squared (In-sample) ---\n")
for (col in colnames(x_train)) {
  actual    <- as.numeric(x_train[, col])
  predicted <- as.numeric(fitted_vals[, col])
  n_min     <- min(length(actual), length(predicted))
  actual    <- tail(actual,    n_min)
  predicted <- tail(predicted, n_min)
  ss_res    <- sum((actual - predicted)^2,             na.rm = TRUE)
  ss_tot    <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  cat(col, ": R² =", round(1 - ss_res / ss_tot, 4), "\n")
}

cat("\n--- Out-of-sample Performance ---\n")
performance(model, x_test)

# R-squared (out-of-sample)
pred_test <- predict(model, n = nrow(x_test))
cat("\n--- R-squared (Out-of-sample) ---\n")
for (col in colnames(x_test)) {
  actual    <- as.numeric(x_test[, col])
  predicted <- as.numeric(pred_test[1:nrow(x_test), col])
  valid     <- complete.cases(actual, predicted)
  ss_res    <- sum((actual[valid] - predicted[valid])^2)
  ss_tot    <- sum((actual[valid] - mean(actual[valid]))^2)
  cat(col, ": R² =", round(1 - ss_res / ss_tot, 4), "\n")
}


# --- 7. Forecasting & Visualisation ---
forecast_results <- predict(model, n = 7)
cat("\n--- 7-Day Forecast ---\n")
print(forecast_results)

plot(model, testing    = x_test)   # Training vs Test
plot(model, n_predict  = 7)        # Historical vs Forecast
