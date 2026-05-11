# =========================================================================================
# PROJECT   : GSTAR Modeling for Temperature Prediction in Surabaya & Sidoarjo
# SCRIPT 00 : Data Extraction, Transformation, and Loading (ETL)
# AUTHORS   : Willy Dava Nugraha, Farikh Muhammad Fauzan, Abdullah Sultan Barizy
# PURPOSE   : Cleans and merges raw BMKG monthly CSV files into a unified time-series dataset.
# INPUT     : /data/raw/*.csv (Monthly raw data separated by station)
# OUTPUT    : /data/processed/data_tavg_bmkg_202505_202604_clean.csv
# =========================================================================================

# --- 1. Environment Setup ---
# Dynamically check and install required packages
required_packages <- c("reshape2", "zoo")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(zoo))

# --- 2. Data Extraction ---
raw_folder <- "data/raw/"
file_list <- list.files(path = raw_folder, pattern = "*.csv", full.names = TRUE)

if (length(file_list) == 0) {
  stop("CRITICAL ERROR: No CSV files found in the 'data/raw/' directory.")
}

cat(sprintf("[INFO] Discovered %d raw data files. Initiating ETL process...\n", length(file_list)))

# --- 3. Data Transformation (Parsing & Merging) ---
data_gabungan <- data.frame()

for (file in file_list) {
  # Attempt standard comma separation
  temp_df <- read.csv(file, stringsAsFactors = FALSE)
  
  # Fallback to semicolon separation if standard parsing fails (Locale differences)
  if (ncol(temp_df) < 2) {
    temp_df <- read.csv(file, sep = ";", stringsAsFactors = FALSE)
  }
  
  if (ncol(temp_df) < 2) {
    stop(sprintf("CRITICAL ERROR: File %s is corrupted or malformed.", basename(file)))
  }
  
  # Standardize core columns
  temp_df <- temp_df[, 1:2]
  colnames(temp_df) <- c("TANGGAL", "SUHU")
  
  # Extract station name from filename (Expected format: YYYYMM_StationName_*.csv)
  nama_stasiun <- strsplit(basename(file), "_")[[1]][2]
  temp_df$STASIUN <- nama_stasiun
  
  data_gabungan <- rbind(data_gabungan, temp_df)
}

# --- 4. Data Transformation (Cleaning & Type Conversion) ---
# Clean whitespace and stray headers
data_gabungan$TANGGAL <- trimws(data_gabungan$TANGGAL)
data_gabungan <- data_gabungan[tolower(data_gabungan$TANGGAL) != "tanggal", ]
data_gabungan <- data_gabungan[data_gabungan$TANGGAL != "", ]

# Parse dates flexibly accommodating multiple formats
data_gabungan$TANGGAL_DATE <- as.Date(data_gabungan$TANGGAL, 
                                      tryFormats = c("%d-%m-%Y", "%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d", "%d.%m.%Y"))
data_gabungan <- data_gabungan[!is.na(data_gabungan$TANGGAL_DATE), ]

# Standardize decimal separator and convert to numeric
data_gabungan$SUHU <- gsub(",", ".", data_gabungan$SUHU)
data_gabungan$SUHU <- as.numeric(data_gabungan$SUHU)

# --- 5. Data Reshaping (Long to Wide Pivot) ---
df_wide <- dcast(data_gabungan, TANGGAL_DATE ~ STASIUN, value.var = "SUHU", fun.aggregate = mean, na.rm = TRUE)
df_wide <- df_wide[order(df_wide$TANGGAL_DATE), ]

# Map column names to official station definitions
mapping_stasiun <- c(
  "TanjungPerak" = "Stasiun Meteorologi Maritim Tanjung Perak",
  "PerakI"       = "Stasiun Meteorologi Perak I",
  "Juanda"       = "Stasiun Meteorologi Juanda"
)

for (old_name in names(mapping_stasiun)) {
  if (old_name %in% colnames(df_wide)) {
    colnames(df_wide)[colnames(df_wide) == old_name] <- mapping_stasiun[[old_name]]
  }
}

# --- 6. Missing Value Imputation ---
# Apply linear interpolation to fill gaps in time-series data ensuring continuous lag calculation
for (col_idx in 2:ncol(df_wide)) {
  df_wide[, col_idx] <- na.approx(df_wide[, col_idx], rule = 2)
}

# --- 7. Final Formatting & Loading ---
df_wide$TANGGAL <- format(df_wide$TANGGAL_DATE, "%d-%m-%Y")

# Enforce strict column order for downstream GSTAR modeling
df_final <- df_wide[, c("TANGGAL", 
                        "Stasiun Meteorologi Maritim Tanjung Perak", 
                        "Stasiun Meteorologi Perak I", 
                        "Stasiun Meteorologi Juanda")]

output_path <- "data/processed/data_tavg_bmkg_202505_202604_clean.csv"
write.csv(df_final, output_path, row.names = FALSE)

cat("[INFO] === ETL PROCESS COMPLETED SUCCESSFULLY ===\n")
cat(sprintf("[INFO] Total temporal observations: %d days\n", nrow(df_final)))
cat(sprintf("[INFO] Clean dataset exported to: %s\n", output_path))