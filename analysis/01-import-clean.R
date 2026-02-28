# analysis/01-import-clean.R
# Purpose:
# - Import the latest Excel file
# - Clean "納入資料" sheet
# - Create:
#   1) df_long  : long-form measurements (original rows)
#   2) df_patient: one row per patient (Table 1-ready)
# - Save outputs to data/processed/

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

# ========= 0) Paths =========
# If you already moved the file into your Quarto project:
#   thesis/data/raw/project_info_00004(with HR).xlsx
# then use the default path below.
path_default <- "data/raw/project_info_00004(with HR).xlsx"

# Fallback: if you haven't moved it yet and you want to read directly from the uploaded path:
path_uploaded <- "/mnt/data/project_info_00004(with HR).xlsx"

path <- if (file.exists(path_default)) path_default else path_uploaded

if (!file.exists(path)) {
  stop("Excel file not found. Put it into data/raw/ or check the path: ", path)
}

# ========= 1) Read sheet =========
sheet_name <- "納入資料"

df_raw <- read_excel(path, sheet = sheet_name)

# ========= 2) Standardize column names (keep Chinese OK) =========
# Trim whitespace in column names to avoid invisible mismatches
names(df_raw) <- str_trim(names(df_raw))

# ========= 3) Basic cleaning helpers =========
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Convert checkbox-like values to Yes/No
# - Your file uses "☑" to indicate TRUE; empty/NA means FALSE
to_yesno <- function(x) {
  x_chr <- str_trim(as.character(x))
  out <- ifelse(!is.na(x_chr) & x_chr == "☑", "Yes", "No")
  factor(out, levels = c("No", "Yes"))
}

# Parse time robustly (if present)
parse_time <- function(x) {
  # If already POSIXct, keep it
  if (inherits(x, "POSIXct")) return(x)

  x_chr <- str_trim(as.character(x))
  # Try several common formats
  parsed <- suppressWarnings(ymd_hms(x_chr, tz = "Asia/Taipei"))
  if (all(is.na(parsed))) parsed <- suppressWarnings(ymd_hm(x_chr, tz = "Asia/Taipei"))
  if (all(is.na(parsed))) parsed <- suppressWarnings(ymd(x_chr, tz = "Asia/Taipei"))
  parsed
}

# ========= 4) Identify columns =========
# Mandatory: ID
if (!("ID" %in% names(df_raw))) {
  stop("Column 'ID' not found in sheet '納入資料'. Please check the column name in Excel.")
}

# Candidate baseline columns (adjust if your sheet uses different names)
cont_vars <- intersect(c("Age", "Height", "Weight", "BMI", "arm"), names(df_raw))
cat_vars  <- intersect(c("Sex"), names(df_raw))

# Comorbidity columns you mentioned / commonly used
dx_vars <- intersect(c("HTN","DM","Stroke","AMI","CAD","PAD","Af","HF"), names(df_raw))

# Time column (optional)
time_col <- intersect(c("Time", "Datetime", "DateTime", "測量時間", "時間"), names(df_raw))
time_col <- if (length(time_col) >= 1) time_col[[1]] else NA_character_

# ========= 5) Clean long-form data =========
df_long <- df_raw %>%
  mutate(
    # ID: keep as character (prevents numeric -> scientific notation issues)
    ID = as.character(ID),

    # Sex: factor, keep original categories (M/F, 男/女, 0/1 etc.)
    across(all_of(cat_vars), ~ factor(as.character(.x))),

    # Continuous variables to numeric
    across(all_of(cont_vars), to_num),

    # Dx checkbox columns -> Yes/No
    across(all_of(dx_vars), to_yesno)
  )

# Parse time if exists
if (!is.na(time_col)) {
  df_long <- df_long %>%
    mutate(
      !!time_col := parse_time(.data[[time_col]])
    )
}

# ========= 6) Create patient-level dataset (one row per ID) =========
# Rules:
# - Baseline continuous/categorical: first non-missing value per patient
# - Dx: if any row is Yes => patient Yes
df_patient <- df_long %>%
  group_by(ID) %>%
  summarise(
    # categorical baseline
    across(all_of(cat_vars), ~ {
      v <- .x[!is.na(.x) & as.character(.x) != ""]
      if (length(v) == 0) NA else v[[1]]
    }),

    # continuous baseline
    across(all_of(cont_vars), ~ {
      v <- .x[!is.na(.x)]
      if (length(v) == 0) NA_real_ else v[[1]]
    }),

    # Dx: any Yes => Yes
    across(all_of(dx_vars), ~ factor(ifelse(any(.x == "Yes", na.rm = TRUE), "Yes", "No"),
                                     levels = c("No","Yes"))),

    .groups = "drop"
  )

# ========= 7) Quick QA =========
cat("\n===== QA Summary =====\n")
cat("File:", path, "\n")
cat("Sheet:", sheet_name, "\n")
cat("Rows (long):", nrow(df_long), "\n")
cat("Unique patients:", n_distinct(df_long$ID), "\n")
cat("Rows (patient-level):", nrow(df_patient), "\n")

if (!is.na(time_col)) {
  cat("Time column detected:", time_col, "\n")
  cat("Time NA rate:", mean(is.na(df_long[[time_col]])), "\n")
}

# check missingness in baseline
if (length(cont_vars) > 0) {
  miss_cont <- sapply(df_patient[cont_vars], function(x) mean(is.na(x)))
  cat("\nMissing rate (patient-level continuous):\n")
  print(miss_cont)
}

if (length(dx_vars) > 0) {
  yes_rate <- sapply(df_patient[dx_vars], function(x) mean(x == "Yes", na.rm = TRUE))
  cat("\nYes rate (patient-level Dx):\n")
  print(yes_rate)
}

# ========= 8) Save outputs =========
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

saveRDS(df_long,    "data/processed/df_long_included.rds")
saveRDS(df_patient, "data/processed/df_patient_included.rds")

# Optionally also export CSV for quick inspection
write.csv(df_patient, "data/processed/df_patient_included.csv", row.names = FALSE)
write.csv(df_long,    "data/processed/df_long_included.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("- data/processed/df_long_included.rds\n")
cat("- data/processed/df_patient_included.rds\n")
cat("- data/processed/df_patient_included.csv\n")
cat("- data/processed/df_long_included.csv\n")
