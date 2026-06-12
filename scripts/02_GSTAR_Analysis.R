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

# --- 2. Data Preparation ---
# setwd("path/to/your/project") # Adjust this to your local directory
# df <- read.csv("data/processed/data_bmkg_combined_clean_v2.csv", check.names = FALSE) 
df <- read.csv(file.choose(), check.names = FALSE)

# Reorder columns to match spatial coordinates
df_suhu <- df[, c("Stasiun Meteorologi Maritim Tanjung Perak", 
                  "Stasiun Meteorologi Perak I", 
                  "Stasiun Meteorologi Juanda")]

# Convert to XTS object
x <- xts(df_suhu, order.by = as.Date(df$TANGGAL, format="%d-%m-%Y"))
colnames(x) <- c("TanjungPerak", "PerakI", "Juanda")

# Train / test split (80:20)
s <- round(nrow(x) * 0.8)   
x_train <- x[1:s, ]
x_test  <- x[-c(1:s), ]

# --- 3. Spatial Weight Matrix (Inverse-Distance) ---
# Coordinates (lon, lat)
coords <- matrix(c(
  112.7297, -7.1927,   # Tanjung Perak
  112.7274, -7.2000,   # Perak I
  112.7870, -7.3799    # Juanda
), ncol = 2, byrow = TRUE)

dist_mat <- as.matrix(dist(coords, method = "euclidean"))
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

# Out-of-sample performance (Testing set)
performance(model, x_test)

# --- 6. Forecasting & Visualization ---
# Forecast for the next 7 days
forecast_results <- predict(model, n = 7)
print(forecast_results)

# Plots
plot(model, testing = x_test) # Train vs Test
plot(model, n_predict = 7)    # Historical vs Forecast
