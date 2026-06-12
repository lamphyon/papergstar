# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 02: GSTAR Parameter Estimation (Manual OLS, per Z equation)
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
#
# Model GSTAR(1;1) untuk tiga lokasi:
#   Z1(t) = a1*Z1(t-1) + b1*w12*Z2(t-1) + b1*w13*Z3(t-1) + e1(t)
#   Z2(t) = a2*Z2(t-1) + b2*w21*Z1(t-1) + b2*w23*Z3(t-1) + e2(t)
#   Z3(t) = a3*Z3(t-1) + b3*w31*Z1(t-1) + b3*w32*Z2(t-1) + e3(t)
#
# Parameter dihitung dengan MKT (OLS): b = (X'X)^-1 * X'Y
# =====================================================================


# --- 1. Load Libraries ---
if (!require("xts")) install.packages("xts")
library(xts)


# --- 2. Data Preparation ---
# setwd("path/to/your/project")   # Sesuaikan dengan direktori lokal Anda
df <- read.csv(file.choose(), check.names = FALSE)

# Urutkan kolom sesuai lokasi Z1, Z2, Z3
df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak",
                  "Stasiun Meteorologi Perak I",
                  "Stasiun Meteorologi Juanda")]

# Konversi ke XTS
x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format = "%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")   # Z1, Z2, Z3

# Train / test split (80:20)
s       <- round(nrow(x) * 0.8)
x_train <- x[1:s, ]
x_test  <- x[-c(1:s), ]

n <- nrow(x_train)


# --- 3. Data Lokasi sebagai Vektor ---
Z1 <- as.numeric(x_train[, "TanjungPerak"])   # Z1 = Tanjung Perak
Z2 <- as.numeric(x_train[, "PerakI"])          # Z2 = Perak I
Z3 <- as.numeric(x_train[, "Juanda"])          # Z3 = Juanda


# --- 4. Spatial Weight Matrix (Inverse-Distance, Row-Normalized) ---
# Koordinat (longitude, latitude)
coords <- matrix(c(
  112.7297, -7.1927,   # Z1 – Tanjung Perak
  112.7274, -7.2000,   # Z2 – Perak I
  112.7870, -7.3799    # Z3 – Juanda
), ncol = 2, byrow = TRUE)

dist_mat <- as.matrix(dist(coords, method = "euclidean"))
w_raw    <- 1 / dist_mat
diag(w_raw) <- 0
W <- w_raw / rowSums(w_raw)   # Bobot ternormalisasi baris

cat("=== Matriks Bobot Lokasi (W) ===\n")
rownames(W) <- colnames(W) <- c("Z1_TanjungPerak", "Z2_PerakI", "Z3_Juanda")
print(round(W, 6))
cat("\n")


# =====================================================================
# --- 5. Pembentukan Lag dan Bobot Spasial ---
# Sesuai gambar: Z_i(t-1) dan bobot spasial W_ij * Z_j(t-1)
# =====================================================================

# Lagged values (t-1)
Z1_lag <- Z1[-n]   # Z1(t-1),  panjang = n-1
Z2_lag <- Z2[-n]   # Z2(t-1)
Z3_lag <- Z3[-n]   # Z3(t-1)

# Nilai saat t (respons), buang t=1
Z1_t <- Z1[-1]
Z2_t <- Z2[-1]
Z3_t <- Z3[-1]

T <- n - 1   # jumlah observasi efektif


# =====================================================================
# --- 6. Estimasi Parameter per Persamaan (MKT / OLS) ---
#
# Untuk setiap persamaan Zi, matriks X disusun dari:
#   - Kolom 1: Zi(t-1)             <- koefisien ai (temporal lag)
#   - Kolom 2: wij*Zj(t-1)        <- koefisien bi (spatial lag ke lokasi j)
#   - Kolom 3: wik*Zk(t-1)        <- koefisien bi (spatial lag ke lokasi k)
#
# b = (X'X)^-1 * X'Y
# =====================================================================

estimate_param <- function(Y, X) {
  b    <- solve(t(X) %*% X) %*% t(X) %*% Y
  Y_hat <- X %*% b
  resid <- Y - Y_hat
  SSE  <- sum(resid^2)
  SST  <- sum((Y - mean(Y))^2)
  R2   <- 1 - SSE / SST
  df_r <- nrow(X) - ncol(X)
  s2   <- SSE / df_r
  se   <- sqrt(s2 * diag(solve(t(X) %*% X)))
  tval <- b[, 1] / se
  pval <- 2 * pt(abs(tval), df = df_r, lower.tail = FALSE)
  list(b = b, se = se, tval = tval, pval = pval,
       R2 = R2, SSE = SSE, Y_hat = Y_hat, resid = resid)
}

# ---- Persamaan Z1: TanjungPerak ----
# Z1(t) = a1*Z1(t-1) + b1*[w12*Z2(t-1) + w13*Z3(t-1)]
X1 <- cbind(
  Z1_lag,
  W[1, 2] * Z2_lag,   # w12 * Z2(t-1)
  W[1, 3] * Z3_lag    # w13 * Z3(t-1)
)
Y1  <- Z1_t
res1 <- estimate_param(Y1, X1)

# ---- Persamaan Z2: Perak I ----
# Z2(t) = a2*Z2(t-1) + b2*[w21*Z1(t-1) + w23*Z3(t-1)]
X2 <- cbind(
  Z2_lag,
  W[2, 1] * Z1_lag,   # w21 * Z1(t-1)
  W[2, 3] * Z3_lag    # w23 * Z3(t-1)
)
Y2  <- Z2_t
res2 <- estimate_param(Y2, X2)

# ---- Persamaan Z3: Juanda ----
# Z3(t) = a3*Z3(t-1) + b3*[w31*Z1(t-1) + w32*Z2(t-1)]
X3 <- cbind(
  Z3_lag,
  W[3, 1] * Z1_lag,   # w31 * Z1(t-1)
  W[3, 2] * Z2_lag    # w32 * Z2(t-1)
)
Y3  <- Z3_t
res3 <- estimate_param(Y3, X3)


# =====================================================================
# --- 7. Tampilkan Hasil Parameter per Persamaan ---
# =====================================================================

print_result <- function(res, eq_label, station, param_names) {
  cat(paste(rep("=", 65), collapse = ""), "\n")
  cat(sprintf("  Persamaan %s  (%s)\n", eq_label, station))
  cat(paste(rep("=", 65), collapse = ""), "\n")
  cat(sprintf("  %-18s %12s %12s %10s %10s\n",
              "Parameter", "Estimate", "Std.Error", "t-value", "p-value"))
  cat(paste(rep("-", 65), collapse = ""), "\n")
  for (i in seq_along(param_names)) {
    sig <- ifelse(res$pval[i] < 0.001, "***",
           ifelse(res$pval[i] < 0.01,  "**",
           ifelse(res$pval[i] < 0.05,  "*",
           ifelse(res$pval[i] < 0.1,   ".", ""))))
    cat(sprintf("  %-18s %12.6f %12.6f %10.4f %10.6f  %s\n",
                param_names[i],
                res$b[i, 1], res$se[i], res$tval[i], res$pval[i], sig))
  }
  cat(sprintf("\n  R-squared (in-sample) : %.6f\n\n", res$R2))
}

cat("\n")
print_result(res1, "Z1", "Tanjung Perak",
             c("a1  [Z1(t-1)]",
               "b1w12 [Z2(t-1)]",
               "b1w13 [Z3(t-1)]"))

print_result(res2, "Z2", "Perak I",
             c("a2  [Z2(t-1)]",
               "b2w21 [Z1(t-1)]",
               "b2w23 [Z3(t-1)]"))

print_result(res3, "Z3", "Juanda",
             c("a3  [Z3(t-1)]",
               "b3w31 [Z1(t-1)]",
               "b3w32 [Z2(t-1)]"))


# =====================================================================
# --- 8. Model GSTAR – Bentuk Persamaan Akhir ---
# =====================================================================
a1 <- res1$b[1, 1]; b1w12 <- res1$b[2, 1]; b1w13 <- res1$b[3, 1]
a2 <- res2$b[1, 1]; b2w21 <- res2$b[2, 1]; b2w23 <- res2$b[3, 1]
a3 <- res3$b[1, 1]; b3w31 <- res3$b[2, 1]; b3w32 <- res3$b[3, 1]

cat(paste(rep("=", 65), collapse = ""), "\n")
cat("  Model GSTAR(1;1) – Persamaan Prediksi\n")
cat(paste(rep("=", 65), collapse = ""), "\n\n")

cat(sprintf(
  "  Z1_hat(t) = %.6f*Z1(t-1) + %.6f*Z2(t-1) + %.6f*Z3(t-1)\n",
  a1, b1w12, b1w13))
cat(sprintf(
  "  Z2_hat(t) = %.6f*Z2(t-1) + %.6f*Z1(t-1) + %.6f*Z3(t-1)\n",
  a2, b2w21, b2w23))
cat(sprintf(
  "  Z3_hat(t) = %.6f*Z3(t-1) + %.6f*Z1(t-1) + %.6f*Z2(t-1)\n\n",
  a3, b3w31, b3w32))


# =====================================================================
# --- 9. Simpan Tabel Parameter ke CSV ---
# =====================================================================
param_df <- data.frame(
  Equation  = c(rep("Z1_TanjungPerak", 3),
                rep("Z2_PerakI", 3),
                rep("Z3_Juanda", 3)),
  Parameter = c("a1 [Z1(t-1)]",  "b1*w12 [Z2(t-1)]", "b1*w13 [Z3(t-1)]",
                "a2 [Z2(t-1)]",  "b2*w21 [Z1(t-1)]", "b2*w23 [Z3(t-1)]",
                "a3 [Z3(t-1)]",  "b3*w31 [Z1(t-1)]", "b3*w32 [Z2(t-1)]"),
  Estimate  = round(c(res1$b, res2$b, res3$b), 8),
  Std.Error = round(c(res1$se, res2$se, res3$se), 8),
  t.value   = round(c(res1$tval, res2$tval, res3$tval), 4),
  p.value   = round(c(res1$pval, res2$pval, res3$pval), 6),
  stringsAsFactors = FALSE
)

write.csv(param_df, "GSTAR_parameters_Z1Z2Z3.csv", row.names = FALSE)
cat("  Parameter table saved  →  GSTAR_parameters_Z1Z2Z3.csv\n\n")


# =====================================================================
# --- 10. Fitted Values & Akurasi (In-sample) ---
# =====================================================================
cat(paste(rep("=", 65), collapse = ""), "\n")
cat("  Akurasi In-sample\n")
cat(paste(rep("=", 65), collapse = ""), "\n")
cat(sprintf("  %-18s %10s %10s %10s\n", "Lokasi", "MSE", "MAPE (%)", "R²"))
cat(paste(rep("-", 55), collapse = ""), "\n")

for (i in 1:3) {
  res  <- list(res1, res2, res3)[[i]]
  loc  <- c("Z1 TanjungPerak", "Z2 PerakI", "Z3 Juanda")[i]
  Y    <- list(Y1, Y2, Y3)[[i]]
  mse  <- mean(res$resid^2)
  mape <- mean(abs(res$resid / Y)) * 100
  cat(sprintf("  %-18s %10.6f %10.4f %10.6f\n", loc, mse, mape, res$R2))
}
cat("\n")


# =====================================================================
# --- 11. Forecasting (Rekursif, h langkah ke depan) ---
# =====================================================================
forecast_gstar <- function(Z1_hist, Z2_hist, Z3_hist,
                           params1, params2, params3, W, h = 7) {
  Z1f <- Z1_hist; Z2f <- Z2_hist; Z3f <- Z3_hist

  for (i in seq_len(h)) {
    z1_prev <- tail(Z1f, 1); z2_prev <- tail(Z2f, 1); z3_prev <- tail(Z3f, 1)

    z1_new <- params1[1]*z1_prev + params1[2]*z2_prev + params1[3]*z3_prev
    z2_new <- params2[1]*z2_prev + params2[2]*z1_prev + params2[3]*z3_prev
    z3_new <- params3[1]*z3_prev + params3[2]*z1_prev + params3[3]*z2_prev

    Z1f <- c(Z1f, z1_new)
    Z2f <- c(Z2f, z2_new)
    Z3f <- c(Z3f, z3_new)
  }

  tail_idx <- (length(Z1f) - h + 1):length(Z1f)
  data.frame(
    h          = 1:h,
    Z1_TanjungPerak = round(Z1f[tail_idx], 4),
    Z2_PerakI       = round(Z2f[tail_idx], 4),
    Z3_Juanda       = round(Z3f[tail_idx], 4)
  )
}

forecast_df <- forecast_gstar(
  Z1, Z2, Z3,
  params1 = c(a1, b1w12, b1w13),
  params2 = c(a2, b2w21, b2w23),
  params3 = c(a3, b3w31, b3w32),
  W = W, h = 7
)

cat(paste(rep("=", 65), collapse = ""), "\n")
cat("  Forecast 7 Hari ke Depan\n")
cat(paste(rep("=", 65), collapse = ""), "\n")
print(forecast_df)
cat("\n")

write.csv(forecast_df, "GSTAR_forecast_7hari.csv", row.names = FALSE)
cat("  Forecast saved  →  GSTAR_forecast_7hari.csv\n")
