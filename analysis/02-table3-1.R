# analysis/02-table3-1.R
# Purpose:
# - Create Table 3-1 (baseline characteristics) as CSV
# - Continuous: mean +/- SD
# - Binary: n (%), percent with 1 decimal
# - Output: outputs/tables/table 3-1.csv

library(dplyr)

# ========= 0) Load cleaned patient-level data =========
df <- readRDS("data/processed/df_patient_included.rds")

# ========= 1) Helpers =========
fmt_mean_sd <- function(x, digits = 2) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s)) return(sprintf(paste0("%.", digits, "f +/- NA"), m))
  sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), m, s)
}

fmt_n_pct <- function(n_event, n_total, digits = 1) {
  if (is.na(n_total) || n_total == 0) return(NA_character_)
  sprintf(paste0("%d (%.", digits, "f%%)"), n_event, 100 * n_event / n_total)
}

# For binary/factor variables, define event as:
# - Yes for Yes/No variables
# - Male for Sex
count_event <- function(x, var_name) {
  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)
  non_miss <- !is.na(x_chr) & x_chr != ""
  n_total <- sum(non_miss)

  if (var_name == "Sex") {
    n_event <- sum(x_low %in% c("m", "male"), na.rm = TRUE)
    return(c(n_event = n_event, n_total = n_total))
  }

  n_event <- sum(x_low %in% c("yes", "y", "1", "true", "t"), na.rm = TRUE)
  if (n_event == 0) {
    x_num <- suppressWarnings(as.numeric(x_chr))
    n_event <- sum(x_num == 1, na.rm = TRUE)
  }
  c(n_event = n_event, n_total = n_total)
}

# ========= 2) Variables =========
required_cont <- c("Age", "Height", "Weight", "BMI")
required_cat <- c("Sex", "HTN", "DM", "CAD", "AMI", "Stroke", "HF", "Af", "PAD")

miss_cont <- setdiff(required_cont, names(df))
miss_cat <- setdiff(required_cat, names(df))
if (length(miss_cont) > 0 || length(miss_cat) > 0) {
  stop(
    "Missing columns in df. Continuous: ", paste(miss_cont, collapse = ", "),
    "; Categorical: ", paste(miss_cat, collapse = ", ")
  )
}

# ========= 3) Build rows =========
rows <- tibble(
  Variable = character(),
  `mean +/- SD or n(%)` = character()
)

add_header <- function(label) {
  tibble(Variable = label, `mean +/- SD or n(%)` = "")
}

add_row <- function(label, value) {
  tibble(Variable = paste0("    ", label), `mean +/- SD or n(%)` = value)
}

# Demographics
sex_cnt <- count_event(df$Sex, "Sex")
rows <- bind_rows(
  rows,
  add_header("Demographics"),
  add_row("Age, years", fmt_mean_sd(df$Age, digits = 2)),
  add_row("Sex", fmt_n_pct(sex_cnt["n_event"], sex_cnt["n_total"], digits = 1))
)

# Anthropometrics
rows <- bind_rows(
  rows,
  add_header("Anthropometrics"),
  add_row("Height, cm", fmt_mean_sd(df$Height, digits = 2)),
  add_row("Weight, kg", fmt_mean_sd(df$Weight, digits = 2)),
  add_row("BMI, kg/m2", fmt_mean_sd(df$BMI, digits = 2))
)

# Comorbidities
htn_cnt <- count_event(df$HTN, "HTN")
dm_cnt <- count_event(df$DM, "DM")
cad_cnt <- count_event(df$CAD, "CAD")
ami_cnt <- count_event(df$AMI, "AMI")
stroke_cnt <- count_event(df$Stroke, "Stroke")
hf_cnt <- count_event(df$HF, "HF")
af_cnt <- count_event(df$Af, "Af")
pad_cnt <- count_event(df$PAD, "PAD")

rows <- bind_rows(
  rows,
  add_header("Comorbidities"),
  add_row("HTN", fmt_n_pct(htn_cnt["n_event"], htn_cnt["n_total"], digits = 1)),
  add_row("DM", fmt_n_pct(dm_cnt["n_event"], dm_cnt["n_total"], digits = 1)),
  add_row("CAD", fmt_n_pct(cad_cnt["n_event"], cad_cnt["n_total"], digits = 1)),
  add_row("AMI", fmt_n_pct(ami_cnt["n_event"], ami_cnt["n_total"], digits = 1)),
  add_row("Stroke", fmt_n_pct(stroke_cnt["n_event"], stroke_cnt["n_total"], digits = 1)),
  add_row("HF", fmt_n_pct(hf_cnt["n_event"], hf_cnt["n_total"], digits = 1)),
  add_row("AF", fmt_n_pct(af_cnt["n_event"], af_cnt["n_total"], digits = 1)),
  add_row("PAD", fmt_n_pct(pad_cnt["n_event"], pad_cnt["n_total"], digits = 1))
)

# ========= 4) Save =========
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(rows, "outputs/tables/table 3-1.csv", row.names = FALSE, fileEncoding = "UTF-8")

message("Saved: outputs/tables/table 3-1.csv")
