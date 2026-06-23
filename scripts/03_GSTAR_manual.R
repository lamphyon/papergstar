# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 02: GSTAR Manual (Fixed std_err & p-value)
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================

# --- 1. Load Libraries ---
if (!require("gstar"))      install.packages("gstar")
if (!require("xts"))        install.packages("xts")
if (!require("geosphere"))  install.packages("geosphere")

library(gstar)
library(xts)
library(geosphere)

# --- 2. Data Preparation ---
df <- read.csv(file.choose(), check.names = FALSE)

df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak",
                  "Stasiun Meteorologi Perak I",
                  "Stasiun Meteorologi Juanda")]

x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format = "%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

x_asli <- x

# Differencing lag 1
x_diff <- diff(x, differences = 1)
x_diff <- na.omit(x_diff)

# Train / test split 90:10
s <- round(nrow(x_diff) * 0.9)
x_train <- x_diff[1:s, ]
x_test  <- x_diff[-c(1:s), ]

x_train_asli <- x_asli[1:(s + 1), ]
x_test_asli  <- x_asli[-c(1:(s + 1)), ]

# --- 3. Spatial Weight Matrix (Cross-Correlation) ---
cor_mat <- cor(x_train)
w_cor   <- abs(cor_mat)
diag(w_cor) <- 0
weight  <- w_cor / rowSums(w_cor)

cat("\n--- Matriks Bobot Cross-Correlation ---\n")
print(weight)

# =====================================================================
# --- 4. FUNGSI MANUAL GSTAR (pad_zero diambil dari namespace gstar) ---
# =====================================================================

pad_zero <- gstar:::pad_zero   # ambil helper internal dari package gstar

gstar_manual <- function(x, W, p = 1) {
  x    <- as.matrix(x)
  sp_loc <- colnames(x)
  n_loc  <- ncol(x)

  # Buat matriks prediktor: lag spasial & non-spasial
  xt   <- pad_zero(x, p)       # lag sendiri
  w_xt <- pad_zero(x %*% W, p) # lag spasial (W * x)
  Xv   <- cbind(xt, w_xt)      # gabung jadi design matrix

  # Vektor respons (semua lokasi digabung kolom-by-kolom)
  z <- matrix(x[-seq(p), ], ncol = 1)

  # --- OLS Estimasi ---
  XtX_inv <- solve(t(Xv) %*% Xv)
  B       <- XtX_inv %*% t(Xv) %*% z     # koefisien
  z_hat   <- Xv %*% B                     # fitted values

  # --- Residual & Derajat Bebas ---
  resid <- z - z_hat
  n     <- nrow(z)          # jumlah observasi efektif
  k     <- ncol(Xv)         # jumlah parameter
  df    <- n - k            # derajat bebas

  # --- Std Error yang Benar ---
  # s^2 = SSE / df  (bukan SSE saja seperti di package asli!)
  sse     <- sum(resid^2)
  s2      <- sse / df
  std_err <- sqrt(s2 * diag(XtX_inv))   # sqrt dari varian tiap koef

  # --- t-value & p-value ---
  t_value <- as.vector(B) / std_err
  p_value <- 2 * pt(-abs(t_value), df = df)

  # --- Nama baris koefisien ---
  # Format: psi{lag}{0=autoregressive, 1=spasial}(lokasi)
  nama_ar  <- paste0("psi", rep(1:p, each = n_loc), "0(",
                     rep(sp_loc, p), ")")
  nama_sp  <- paste0("psi", rep(1:p, each = n_loc), "1(",
                     rep(sp_loc, p), ")")
  nama_all <- c(nama_ar, nama_sp)

  # --- Tabel Koefisien ---
  coef_table <- data.frame(
    Estimate   = as.vector(B),
    Std.Err    = std_err,
    t.value    = t_value,
    p.value    = p_value,
    row.names  = nama_all
  )
  colnames(coef_table) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")

  # --- AIC ---
  aic <- n * log(sse / n) + 2 * k

  # --- Fitted values per lokasi ---
  fitted_mat <- matrix(z_hat, ncol = n_loc)
  colnames(fitted_mat) <- sp_loc

  # --- Output ---
  list(
    coefficients = coef_table,
    B            = B,
    std_err      = std_err,
    t_value      = t_value,
    p_value      = p_value,
    fitted       = fitted_mat,
    residuals    = matrix(resid, ncol = n_loc),
    SSE          = sse,
    df           = df,
    s2           = s2,
    AIC          = aic,
    Xv           = Xv,
    z            = z,
    p            = p,
    sp_loc       = sp_loc,
    W            = W
  )
}

# --- 5. Fit Model Manual ---
fit_manual <- gstar_manual(x_train, W = weight, p = 1)

cat("\n============================================\n")
cat("   SUMMARY GSTAR MANUAL (p-value fixed)     \n")
cat("============================================\n\n")
cat("Coefficients:\n")
printCoefmat(fit_manual$coefficients, digits = 4, signif.stars = TRUE)
cat("\nAIC:", round(fit_manual$AIC, 4), "\n")
cat("Degrees of freedom:", fit_manual$df, "\n")
cat("MSE:", round(fit_manual$SSE / nrow(fit_manual$z), 6), "\n")

# --- 6. Fungsi Predict Manual ---
predict_manual <- function(fit, n_ahead = 1) {
  W      <- fit$W
  p      <- fit$p
  B      <- fit$B
  sp_loc <- fit$sp_loc
  n_loc  <- length(sp_loc)

  # Ambil nilai terakhir dari data training sebagai seed
  # fit$z sudah dalam bentuk panjang (n_obs * n_loc) x 1
  # Rekonstruksi matrix observasi terakhir sebanyak p baris
  z_mat    <- matrix(fit$z, ncol = n_loc)
  last_obs <- tail(z_mat, p)  # p baris terakhir sebagai seed

  preds <- matrix(NA, nrow = n_ahead, ncol = n_loc)
  colnames(preds) <- sp_loc

  history <- last_obs  # matriks seed (p x n_loc)

  for (t in 1:n_ahead) {
    # Bangun vektor prediktor untuk 1 langkah ke depan
    x_lag  <- as.vector(t(history))          # lag sendiri (semua lag, semua lokasi)
    wx_lag <- as.vector(t(history %*% W))    # lag spasial
    xv_new <- matrix(c(x_lag, wx_lag), nrow = 1)

    pred_vec <- xv_new %*% B                 # prediksi semua lokasi sekaligus
    preds[t, ] <- as.vector(pred_vec)

    # Geser history: buang baris terlama, tambah prediksi baru
    history <- rbind(history[-1, ], preds[t, ])
  }

  preds
}

# --- 7. Evaluasi In-Sample ---
cat("\n============================================\n")
cat("   EVALUASI PERFORMA PADA SKALA SUHU ASLI   \n")
cat("============================================\n")

calc_metrics <- function(actual, predicted, label) {
  mse    <- mean((actual - predicted)^2, na.rm = TRUE)
  mape   <- mean(abs((actual - predicted) / actual), na.rm = TRUE) * 100
  ss_res <- sum((actual - predicted)^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  r2     <- 1 - (ss_res / ss_tot)
  cat(sprintf("%-15s | MSE: %8.4f | MAPE: %6.2f%% | R2: %6.4f\n",
              label, mse, mape, r2))
}

cat("\n--- Performa In-Sample (Data Training) ---\n")
fitted_diff       <- fit_manual$fitted
aktual_sebelumnya <- coredata(x_train_asli[2:s, ])
fitted_actual     <- aktual_sebelumnya + fitted_diff
actual_train_eval <- coredata(x_train_asli[3:(s + 1), ])

for (i in 1:ncol(x)) {
  calc_metrics(actual_train_eval[, i], fitted_actual[, i], colnames(x)[i])
}

# --- 8. Evaluasi Out-of-Sample ---
cat("\n--- Performa Out-of-Sample (Data Testing) ---\n")
pred_diff_test   <- predict_manual(fit_manual, n_ahead = nrow(x_test))
pred_actual_test <- matrix(NA, nrow = nrow(pred_diff_test), ncol = ncol(x))

for (j in 1:ncol(x)) {
  last_train_val        <- as.numeric(tail(x_train_asli[, j], 1))
  pred_actual_test[, j] <- last_train_val + cumsum(pred_diff_test[, j])
}

actual_test_eval <- coredata(x_test_asli)
for (j in 1:ncol(x)) {
  calc_metrics(actual_test_eval[, j], pred_actual_test[, j], colnames(x)[j])
}

# --- 9. Forecast 7 Hari ke Depan ---
cat("\n============================================\n")
cat("   FORECAST 7 HARI KEDEPAN (SUHU AKTUAL)    \n")
cat("============================================\n\n")

forecast_diff   <- predict_manual(fit_manual, n_ahead = 7)
forecast_actual <- matrix(NA, nrow = 7, ncol = ncol(x))
colnames(forecast_actual) <- colnames(x)

last_date    <- index(tail(x, 1))
future_dates <- seq(last_date + 1, by = "days", length.out = 7)

for (j in 1:ncol(x)) {
  last_known_val        <- as.numeric(tail(x[, j], 1))
  forecast_actual[, j]  <- last_known_val + cumsum(forecast_diff[, j])
}

forecast_xts <- xts(forecast_actual, order.by = future_dates)
print(forecast_xts)
