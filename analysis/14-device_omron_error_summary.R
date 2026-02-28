# analysis/14-device_omron_error_summary.R
# Purpose:
# - Analyze raw "device analysis" sheet
# - Treat `omron error` == 1 as measurement failure
# - Summarize per-patient Omron error counts
# - Plot:
#   1) Histogram of error counts per patient
#   2) Bar chart of frequency by error-count category
# - Report proportion of patients with >= 1 Omron error

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)

# ---- paths ----
path_xlsx <- "data/raw/project_info_00004(with HR).xlsx"
if (!file.exists(path_xlsx)) {
  stop("Raw file not found: ", path_xlsx)
}

# ---- helpers ----
norm <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[[:space:]_\\-]+", "")
}

# ---- sheet detection ----
sheets <- excel_sheets(path_xlsx)
sheet_norm <- norm(sheets)

target_sheet <- sheets[which(sheet_norm == norm("device analysis"))]
if (length(target_sheet) == 0) {
  target_sheet <- sheets[str_detect(sheet_norm, "device") & str_detect(sheet_norm, "analysis")]
}
if (length(target_sheet) == 0) {
  stop("Cannot find sheet 'device analysis'. Available sheets: ", paste(sheets, collapse = ", "))
}
target_sheet <- target_sheet[1]
message("Using sheet: ", target_sheet)

df <- read_excel(path_xlsx, sheet = target_sheet)
names(df) <- str_trim(names(df))
nm <- names(df)
nm_norm <- norm(nm)

# ---- required columns ----
id_candidates_norm <- c(norm("ID"), norm("patient id"), norm("patient_id"))
id_idx <- which(nm_norm %in% id_candidates_norm)
if (length(id_idx) == 0) {
  id_idx <- which(str_detect(nm_norm, "^id$|patient"))
}
if (length(id_idx) == 0) {
  stop("Cannot detect patient ID column in sheet: ", target_sheet)
}
id_col <- nm[id_idx[1]]

omron_candidates_norm <- c(norm("omron error"), norm("omron_error"), norm("omronerror"))
omron_idx <- which(nm_norm %in% omron_candidates_norm)
if (length(omron_idx) == 0) {
  omron_idx <- which(str_detect(nm_norm, "omron") & str_detect(nm_norm, "error"))
}
if (length(omron_idx) == 0) {
  stop("Cannot detect 'omron error' column in sheet: ", target_sheet)
}
omron_col <- nm[omron_idx[1]]

message("Using ID column: ", id_col)
message("Using Omron error column: ", omron_col)

# ---- clean and summarize ----
to01 <- function(x) {
  x_chr <- str_trim(as.character(x))
  x_num <- suppressWarnings(as.numeric(x_chr))
  # Treat explicit 1 as error; others as 0
  out <- ifelse(!is.na(x_num) & x_num == 1, 1L, 0L)
  out[is.na(out)] <- 0L
  out
}

df2 <- df %>%
  transmute(
    ID = as.character(.data[[id_col]]),
    omron_error = to01(.data[[omron_col]])
  ) %>%
  filter(!is.na(ID), ID != "")

patient_err <- df2 %>%
  group_by(ID) %>%
  summarise(
    n_records = n(),
    omron_error_count = sum(omron_error, na.rm = TRUE),
    has_omron_error = as.integer(omron_error_count >= 1),
    .groups = "drop"
  )

n_patient <- nrow(patient_err)
n_with_error <- sum(patient_err$has_omron_error, na.rm = TRUE)
prop_with_error <- ifelse(n_patient == 0, NA_real_, n_with_error / n_patient)

# ---- outputs ----
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

write.csv(patient_err, "outputs/tables/device_omron_error_by_patient.csv",
  row.names = FALSE, fileEncoding = "UTF-8")

summary_tbl <- tibble(
  n_patient = n_patient,
  n_with_omron_error = n_with_error,
  proportion_with_omron_error = prop_with_error
)
write.csv(summary_tbl, "outputs/tables/device_omron_error_summary.csv",
  row.names = FALSE, fileEncoding = "UTF-8")

# Histogram
p_hist <- ggplot(patient_err, aes(x = omron_error_count)) +
  geom_histogram(binwidth = 1, boundary = -0.5, closed = "right", fill = "#2C7FB8", color = "white") +
  scale_x_continuous(breaks = seq(0, max(patient_err$omron_error_count, na.rm = TRUE), by = 1)) +
  labs(
    title = "Histogram: Omron Error Count per Patient",
    x = "Omron error count",
    y = "Number of patients"
  ) +
  theme_minimal(base_size = 13)

# Bar by count category
bar_df <- patient_err %>%
  count(omron_error_count, name = "n_patient") %>%
  arrange(omron_error_count)

p_bar <- ggplot(bar_df, aes(x = factor(omron_error_count), y = n_patient)) +
  geom_col(fill = "#41AE76") +
  labs(
    title = "Bar: Frequency of Omron Error Count",
    x = "Omron error count",
    y = "Number of patients"
  ) +
  theme_minimal(base_size = 13)

ggsave("outputs/figures/device_omron_error_histogram.png", p_hist, width = 8, height = 5, dpi = 300)
ggsave("outputs/figures/device_omron_error_bar.png", p_bar, width = 8, height = 5, dpi = 300)

cat("\n===== Device Omron Error Summary =====\n")
cat("Sheet:", target_sheet, "\n")
cat("ID column:", id_col, "\n")
cat("Omron error column:", omron_col, "\n")
cat("Patients:", n_patient, "\n")
cat("Patients with >=1 Omron error:", n_with_error, "\n")
cat("Proportion with >=1 Omron error:", round(prop_with_error, 4), "\n")
cat("\nSaved:\n")
cat("- outputs/tables/device_omron_error_by_patient.csv\n")
cat("- outputs/tables/device_omron_error_summary.csv\n")
cat("- outputs/figures/device_omron_error_histogram.png\n")
cat("- outputs/figures/device_omron_error_bar.png\n")
