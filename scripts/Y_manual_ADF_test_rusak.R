# =====================================================================
# --- 1. Load Libraries ---
# =====================================================================
if (!require("xts")) install.packages("xts")
if (!require("geosphere")) install.packages("geosphere") 

library(xts)
library(geosphere)

# =====================================================================
# --- 2. Data Preparation ---
# =====================================================================
# (Bagian ini sama seperti sebelumnya)
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
# --- 2.5. UJI STASIONERITAS (ADF TEST - MANUAL) ---
# =====================================================================
cat("\n============================================\n")
cat("   UJI STASIONERITAS (MANUAL ADF TEST)      \n")
cat("============================================\n")

# Fungsi untuk menghitung ADF statistic menggunakan regresi linear (lm)
manual_adf_stat <- function(y) {
  n <- length(y)
  # Menentukan jumlah lag (k) default (sama seperti standar rumus tseries)
  k <- trunc((n - 1)^(1/3)) 
  
  dy <- diff(y)          # Delta y_t
  yt_1 <- y[1:(n-1)]     # Lag level (y_{t-1})
  tt <- 1:(n-1)          # Time trend (t)
  
  # Indeks data yang digunakan regresi setelah dikurangi oleh jumlah lag
  idx <- (k + 1):(n - 1) 
  
  dy_target <- dy[idx]
  yt_1_target <- yt_1[idx]
  tt_target <- tt[idx]
  
  # Pembuatan matriks lag diff
  if (k > 0) {
    lagged_dy <- matrix(NA, nrow = length(idx), ncol = k)
    for (i in 1:k) {
      lagged_dy[, i] <- dy[(k + 1 - i):(n - 1 - i)]
    }
    # Fit regresi linear dengan lag
    fit <- lm(dy_target ~ yt_1_target + tt_target + lagged_dy)
  } else {
    # Fit regresi linear tanpa lag
    fit <- lm(dy_target ~ yt_1_target + tt_target)
  }
  
  # Ambil nilai t-statistic untuk koefisien y_{t-1} 
  # Baris ke-2 adalah yt_1_target (Baris ke-1 selalu Intercept)
  t_stat <- summary(fit)$coefficients[2, "t value"]
  
  return(t_stat)
}

# Nilai kritis ADF tabel MacKinnon (model dengan konstanta dan trend, sampel besar)
# Alpha = 0.05
cv_5 <- -3.41

# Menguji stasioneritas masing-masing stasiun pada data training
for(i in 1:ncol(x_train_centered)) {
  stasiun <- colnames(x_train_centered)[i]
  
  # Ekstrak data menjadi vektor numerik biasa
  data_uji <- as.numeric(coredata(x_train_centered[, i]))
  
  # Jalankan manual ADF test
  adf_statistic <- manual_adf_stat(data_uji)
  
# =====================================================================
# --- 2.6. PENANGANAN DATA NON-STASIONER (DIFFERENCING) ---
# =====================================================================
cat("\n========================================================\n")
cat("   PENANGANAN DATA NON-STASIONER (DIFFERENCING d=1)      \n")
cat("========================================================\n")

# Copy data untuk hasil yang sudah di-diff
x_train_final <- x_train_centered

for(i in 1:ncol(x_train_centered)) {
  stasiun <- colnames(x_train_centered)[i]
  data_uji <- as.numeric(coredata(x_train_centered[, i]))
  
  # Cek ADF manual kembali
  adf_statistic <- manual_adf_stat(data_uji)
  
  if(adf_statistic >= cv_5) {
    cat(sprintf("\nStasiun %s TIDAK STASIONER. Melakukan differencing...\n", stasiun))
    
    # Melakukan differencing (d=1)
    # Kita menggunakan diff() pada data xts
    diff_data <- diff(x_train_centered[, i], differences = 1)
    
    # Hapus baris pertama yang menjadi NA akibat diff
    x_train_final[, i] <- rbind(0, diff_data[-1, ]) 
    
    cat("Status: Data telah didiferensiasi.\n")
  } else {
    cat(sprintf("\nStasiun %s SUDAH STASIONER. Tidak perlu differencing.\n", stasiun))
  }
}

# Verifikasi hasil akhir
cat("\n--- Verifikasi Stasioneritas Akhir ---\n")
for(i in 1:ncol(x_train_final)) {
  # Ambil data non-zero (setelah diff, baris pertama jadi 0)
  data_check <- as.numeric(x_train_final[, i])
  data_check <- data_check[data_check != 0] 
  
  final_adf <- manual_adf_stat(data_check)
  cat(sprintf("Stasiun %s: ADF = %.4f | Kesimpulan: %s\n", 
              colnames(x_train_final)[i], 
              final_adf, 
              ifelse(final_adf < cv_5, "STASIONER", "TETAP TIDAK STASIONER")))
}
