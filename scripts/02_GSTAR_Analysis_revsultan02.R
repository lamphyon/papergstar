# =====================================================================
# Project: GSTAR Modeling for Temperature in Surabaya & Sidoarjo
# Script 02: GSTAR Modeling, Performance, and Forecasting
# Authors: Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# =====================================================================

# --- 1. Load Libraries ---
if (!require("gstar")) install.packages("gstar")
if (!require("xts")) install.packages("xts")
if (!require("geosphere")) install.packages("geosphere") # Tambahan untuk Haversine

library(gstar)
library(xts)
library(geosphere)

# --- 2. Data Preparation ---
# (Pastikan objek 'x' yang berisi xts raw data suhu sudah berjalan normal sebelum ini)
# df <- read.csv(file.choose(), check.names = FALSE) 
# x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format="%d-%m-%Y"))
# colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

# Simpan data asli (suhu aktual) sebagai referensi inverse-differencing
x_asli <- x

# Lakukan differencing lag 1 agar data stasioner
x_diff <- diff(x, differences = 1)
x_diff <- na.omit(x_diff) # Hapus baris NA pertama akibat differencing

# Train / test split (80:20) pada data differencing untuk model
s <- round(nrow(x_diff) * 0.8)   
x_train <- x_diff[1:s, ]
x_test  <- x_diff[-c(1:s), ]

# Train / test split pada data asli (bergeser 1 indeks karena baris pertama NA dihapus)
x_train_asli <- x_asli[1:(s+1), ]
x_test_asli  <- x_asli[-c(1:(s+1)), ]

# --- 3. Spatial Weight Matrix (Cross-Correlation) ---
# Menghitung matriks korelasi Pearson dari data training yang sudah stasioner
cor_mat <- cor(x_train)

# Memastikan nilai korelasi absolut, dan diagonalnya 0 (tidak ada korelasi dengan diri sendiri di matriks pembobot)
w_cor <- abs(cor_mat)
diag(w_cor) <- 0

# Normalisasi baris (Row-normalization) agar jumlah tiap baris = 1
weight <- w_cor / rowSums(w_cor)

cat("\n--- Matriks Bobot Cross-Correlation ---\n")
print(weight)

# --- 4. GSTAR Model Fitting ---
model <- gstar(x      = x_train,
               weight = weight,
               p      = 1,   # Lag 1 optimal
               d      = 0,   # Data dimasukkan sudah di-diff
               est    = "OLS")

print(summary(model))

# --- 5. Accuracy Evaluation (Inverse-Differencing) ---
cat("\n============================================\n")
cat("   EVALUASI PERFORMA PADA SKALA SUHU ASLI   \n")
cat("============================================\n")

# Fungsi kustom untuk metrik skala asli
calc_metrics <- function(actual, predicted, label) {
  mse  <- mean((actual - predicted)^2, na.rm = TRUE)
  mape <- mean(abs((actual - predicted) / actual), na.rm = TRUE) * 100
  ss_res <- sum((actual - predicted)^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- 1 - (ss_res / ss_tot)
  cat(sprintf("%-15s | MSE: %6.4f | MAPE: %5.2f%% | R²: %6.4f\n", label, mse, mape, r2))
}

cat("\n--- Performa In-Sample (Data Training) ---\n")
fitted_diff <- fitted(model)

# Mengembalikan fitted value selisih menjadi suhu: Suhu_Hat(t) = Suhu_Asli(t-1) + Diff_Hat(t)
aktual_sebelumnya <- coredata(x_train_asli[2:s, ])
fitted_actual     <- aktual_sebelumnya + coredata(fitted_diff)
actual_train_eval <- coredata(x_train_asli[3:(s+1), ]) 

for (i in 1:ncol(x)) {
  calc_metrics(actual_train_eval[, i], fitted_actual[, i], colnames(x)[i])
}

cat("\n--- Performa Out-of-Sample (Data Testing) ---\n")
pred_diff_test <- predict(model, n = nrow(x_test))
pred_actual_test <- matrix(NA, nrow = nrow(pred_diff_test), ncol = ncol(x))

for(j in 1:ncol(x)) {
  # Suhu_Hat(t+k) = Suhu_Asli_Terakhir_Train + Kumulatif_Prediksi_Selisih
  last_train_val <- as.numeric(tail(x_train_asli[, j], 1))
  pred_actual_test[, j] <- last_train_val + cumsum(as.numeric(pred_diff_test[, j]))
}

actual_test_eval <- coredata(x_test_asli)
for (j in 1:ncol(x)) {
  calc_metrics(actual_test_eval[, j], pred_actual_test[, j], colnames(x)[j])
}

# --- 6. Forecasting & Visualization (Suhu Aktual) ---
cat("\n============================================\n")
cat("   FORECAST 7 HARI KEDEPAN (SUHU AKTUAL)    \n")
cat("============================================\n")

forecast_diff <- predict(model, n = 7)
forecast_actual <- matrix(NA, nrow = 7, ncol = ncol(x))
colnames(forecast_actual) <- colnames(x)

# Buat deret waktu untuk 7 hari ke depan
last_date <- index(tail(x, 1))
future_dates <- seq(last_date + 1, by = "days", length.out = 7)

for(j in 1:ncol(x)) {
  # Patokan awal forecast adalah nilai asli hari terakhir pada seluruh dataset
  last_known_val <- as.numeric(tail(x[, j], 1))
  forecast_actual[, j] <- last_known_val + cumsum(as.numeric(forecast_diff[, j]))
}

# Tampilkan prediksi suhu sesungguhnya dalam bentuk XTS
forecast_xts <- xts(forecast_actual, order.by = future_dates)
print(forecast_xts)
