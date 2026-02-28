# analysis/15-omron_error_descriptive.R
# Descriptive comparison only (no modeling, no p-values):
# 1) Attempt level: Omron error (1) vs non-error (0)
# 2) Patient level: any_Omron_error (1) vs no error (0)
# Output metrics: group summaries + difference + 95% CI + SMD

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

# -----------------------------
# 0) Paths and settings
# -----------------------------
sheet_name_target <- "device analysis"

# If not NULL, force use this patient ID column name.
# Example: patient_id_col_force <- "SN"
patient_id_col_force <- NULL

# Robust path resolution for different working directories
path_candidates <- c(
  "data/raw/project_info_00004(with HR).xlsx",
  file.path("..", "data", "raw", "project_info_00004(with HR).xlsx"),
  file.path("..", "..", "data", "raw", "project_info_00004(with HR).xlsx"),
  file.path("g:/我的雲端硬碟/研究/PhD Thesis", "data/raw/project_info_00004(with HR).xlsx")
)
path_exists <- path_candidates[file.exists(path_candidates)]
if (length(path_exists) == 0) {
  stop(
    "Raw file not found. Checked paths:\n- ",
    paste(path_candidates, collapse = "\n- ")
  )
}
path_xlsx <- path_exists[1]
message("Using raw file: ", normalizePath(path_xlsx, winslash = "/", mustWork = FALSE))

# -----------------------------
# 1) Utilities
# -----------------------------
norm <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[[:space:]_\\-]+", "")
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

to_binary01 <- function(x) {
  x_chr <- str_trim(as.character(x))
  x_num <- suppressWarnings(as.numeric(x_chr))
  out <- ifelse(!is.na(x_num) & x_num == 1, 1L, 0L)
  out[is.na(out)] <- 0L
  out
}

to_yes01 <- function(x) {
  x_chr <- str_to_lower(str_trim(as.character(x)))
  out <- ifelse(
    x_chr %in% c("yes", "y", "1", "true", "t"),
    1L,
    ifelse(x_chr %in% c("no", "n", "0", "false", "f"), 0L, NA_integer_)
  )
  out
}

# Disease-style parser:
# - explicit yes/no/0/1 first
# - if still unresolved: non-empty => 1, empty => 0
to_disease01 <- function(x) {
  x_chr <- str_trim(as.character(x))
  y <- to_yes01(x_chr)
  x_num <- suppressWarnings(as.numeric(x_chr))
  y[is.na(y) & !is.na(x_num) & x_num %in% c(0, 1)] <- as.integer(x_num[is.na(y) & !is.na(x_num) & x_num %in% c(0, 1)])
  y[is.na(y) & !is.na(x_chr) & x_chr != ""] <- 1L
  y[is.na(y) & (is.na(x_chr) | x_chr == "")] <- 0L
  y
}

# Sex parser for SMD (binary coding):
# - male-like => 1, female-like => 0
# - if unknown labels but exactly 2 non-missing levels, map by sorted order.
to_sex01 <- function(x) {
  x_chr <- str_to_lower(str_trim(as.character(x)))
  out <- rep(NA_integer_, length(x_chr))
  out[x_chr %in% c("m", "male", "man", "boy", "1", "yes")] <- 1L
  out[x_chr %in% c("f", "female", "woman", "girl", "0", "no")] <- 0L

  miss <- is.na(out) & !is.na(x_chr) & x_chr != ""
  lv <- sort(unique(x_chr[miss]))
  if (length(lv) == 2) {
    out[miss & x_chr == lv[1]] <- 0L
    out[miss & x_chr == lv[2]] <- 1L
  }
  out
}

parse_binary_var <- function(x, var_name, disease_vars) {
  if (var_name %in% disease_vars) return(to_disease01(x))
  if (var_name == "Sex") return(to_sex01(x))
  y <- to_yes01(x)
  if (all(is.na(y))) {
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    y <- ifelse(x_num %in% c(0, 1), as.integer(x_num), NA_integer_)
  }
  y
}

mean_ci_diff <- function(x1, x0) {
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]
  n1 <- length(x1)
  n0 <- length(x0)
  if (n1 < 2 || n0 < 2) {
    return(c(diff = NA_real_, lower = NA_real_, upper = NA_real_, smd = NA_real_))
  }
  m1 <- mean(x1)
  m0 <- mean(x0)
  s1 <- stats::sd(x1)
  s0 <- stats::sd(x0)
  diff <- m1 - m0
  se <- sqrt((s1^2 / n1) + (s0^2 / n0))
  # Welch-Satterthwaite df
  df <- (s1^2 / n1 + s0^2 / n0)^2 /
    (((s1^2 / n1)^2 / (n1 - 1)) + ((s0^2 / n0)^2 / (n0 - 1)))
  tcrit <- stats::qt(0.975, df = df)
  lower <- diff - tcrit * se
  upper <- diff + tcrit * se
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
  smd <- ifelse(is.finite(pooled_sd) && pooled_sd > 0, diff / pooled_sd, NA_real_)
  c(diff = diff, lower = lower, upper = upper, smd = smd)
}

prop_ci_diff <- function(y1, y0) {
  y1 <- y1[!is.na(y1)]
  y0 <- y0[!is.na(y0)]
  n1 <- length(y1)
  n0 <- length(y0)
  if (n1 == 0 || n0 == 0) {
    return(c(diff = NA_real_, lower = NA_real_, upper = NA_real_, smd = NA_real_))
  }
  p1 <- mean(y1 == 1)
  p0 <- mean(y0 == 1)
  diff <- p1 - p0
  se <- sqrt((p1 * (1 - p1) / n1) + (p0 * (1 - p0) / n0))
  lower <- diff - 1.96 * se
  upper <- diff + 1.96 * se
  # Binary SMD
  denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  smd <- ifelse(is.finite(denom) && denom > 0, diff / denom, NA_real_)
  c(diff = diff, lower = lower, upper = upper, smd = smd)
}

fmt_mean_sd <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s)) return(sprintf(paste0("%.", digits, "f ± NA"), m))
  sprintf(paste0("%.", digits, "f ± %.", digits, "f"), m, s)
}

fmt_n_pct <- function(x01, pct_digits = 2) {
  x01 <- x01[!is.na(x01)]
  n <- length(x01)
  if (n == 0) return(NA_character_)
  e <- sum(x01 == 1, na.rm = TRUE)
  sprintf(paste0("%d (%.", pct_digits, "f%%)"), e, 100 * e / n)
}

fmt_diff_ci <- function(diff, lower, upper, digits = 2) {
  if (!is.finite(diff) || !is.finite(lower) || !is.finite(upper)) return(NA_character_)
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), diff, lower, upper)
}

fmt_smd <- function(x, digits = 2) {
  if (!is.finite(x)) return(NA_character_)
  sprintf(paste0("%.", digits, "f"), x)
}

make_compare_table <- function(dat, grp_col, cont_vars, bin_vars, disease_vars) {
  g <- dat[[grp_col]]
  if (!all(g %in% c(0L, 1L, NA))) stop("Group column must be 0/1.")

  out_cont <- lapply(cont_vars, function(v) {
    x <- to_num(dat[[v]])
    x1 <- x[g == 1]
    x0 <- x[g == 0]
    st <- mean_ci_diff(x1, x0)
    tibble(
      variables = v,
      Error = fmt_mean_sd(x1),
      `No Error` = fmt_mean_sd(x0),
      `Difference (error-No)[95% CI ]` = fmt_diff_ci(st["diff"], st["lower"], st["upper"]),
      SMD = fmt_smd(st["smd"]),
      type = "continuous"
    )
  }) %>% bind_rows()

  out_bin <- lapply(bin_vars, function(v) {
    x_raw <- dat[[v]]
    x <- parse_binary_var(x_raw, v, disease_vars)
    x1 <- x[g == 1]
    x0 <- x[g == 0]
    st <- prop_ci_diff(x1, x0)
    tibble(
      variables = v,
      Error = fmt_n_pct(x1),
      `No Error` = fmt_n_pct(x0),
      `Difference (error-No)[95% CI ]` = fmt_diff_ci(st["diff"], st["lower"], st["upper"]),
      SMD = fmt_smd(st["smd"]),
      type = "binary"
    )
  }) %>% bind_rows()

  bind_rows(out_cont, out_bin) %>%
    select(variables, Error, `No Error`, `Difference (error-No)[95% CI ]`, SMD, type)
}

label_variable <- function(x) {
  dplyr::case_when(
    x == "Age" ~ "Age (years)",
    x == "arm" ~ "Arm circumference (cm)",
    x == "BMI" ~ "BMI (kg/m^2)",
    x %in% c("K-HR", "k-hr", "k_hr", "K_HR", "KHR") ~ "Heart rate (bpm)",
    x == "ASCVD_any" ~ "ASCVD (any)",
    x == "Af" ~ "Atrial fibrillation",
    x == "DM" ~ "Diabetes mellitus",
    x == "HF" ~ "Heart failure",
    x == "HTN" ~ "Hypertension",
    x == "Sex" ~ "Male sex",
    x == "PAD" ~ "Peripheral artery disease",
    TRUE ~ x
  )
}

# -----------------------------
# 2) Load raw sheet
# -----------------------------
sheets <- excel_sheets(path_xlsx)
sheet_norm <- norm(sheets)
sheet_pick <- sheets[which(sheet_norm == norm(sheet_name_target))]
if (length(sheet_pick) == 0) {
  sheet_pick <- sheets[str_detect(sheet_norm, "device") & str_detect(sheet_norm, "analysis")]
}
if (length(sheet_pick) == 0) {
  stop("Cannot find sheet 'device analysis'. Available sheets: ", paste(sheets, collapse = ", "))
}
sheet_pick <- sheet_pick[1]

df <- read_excel(path_xlsx, sheet = sheet_pick)
names(df) <- str_trim(names(df))
nm <- names(df)
nm_norm <- norm(nm)

# -----------------------------
# 3) Detect key columns
# -----------------------------
if (!is.null(patient_id_col_force)) {
  if (!(patient_id_col_force %in% nm)) stop("patient_id_col_force not found: ", patient_id_col_force)
  id_col <- patient_id_col_force
} else {
  id_priority <- c("ID", "SN", "patient id", "patient_id")
  id_col <- NA_character_
  for (cand in id_priority) {
    idx <- which(nm_norm == norm(cand))
    if (length(idx) > 0) {
      id_col <- nm[idx[1]]
      break
    }
  }
  if (is.na(id_col)) {
    idx <- which(str_detect(nm_norm, "^id$|^sn$|patient"))
    if (length(idx) == 0) stop("Cannot detect patient ID column.")
    id_col <- nm[idx[1]]
  }
}

omron_idx <- which(nm_norm %in% c(norm("omron error"), norm("omron_error"), norm("omronerror")))
if (length(omron_idx) == 0) {
  omron_idx <- which(str_detect(nm_norm, "omron") & str_detect(nm_norm, "error"))
}
if (length(omron_idx) == 0) stop("Cannot detect omron error column.")
omron_col <- nm[omron_idx[1]]

# Covariates of interest (fixed as requested)
cont_wanted <- c("Age", "arm", "BMI")
hr_candidates <- c("K-HR", "k-hr", "k_hr")
hr_col <- hr_candidates[hr_candidates %in% nm]
if (length(hr_col) > 0) cont_wanted <- c(cont_wanted, hr_col[1])
cont_vars <- intersect(cont_wanted, nm)

bin_wanted <- c("Sex", "HTN", "DM", "Af", "HF", "PAD", "ASCVD_any")
bin_vars <- intersect(bin_wanted, names(df))
disease_vars <- intersect(c("HTN", "DM", "Af", "HF", "PAD", "ASCVD_any"), names(df))

# Build ASCVD_any from components (robust to checkbox-like values)
need_ascvd <- c("CAD", "AMI", "Stroke", "PAD")
if (all(need_ascvd %in% nm)) {
  df <- df %>%
    mutate(
      ASCVD_any = if_else(
        to_disease01(.data[["CAD"]]) == 1L |
          to_disease01(.data[["AMI"]]) == 1L |
          to_disease01(.data[["Stroke"]]) == 1L |
          to_disease01(.data[["PAD"]]) == 1L,
        "Yes", "No",
        missing = "No"
      )
    )
}
bin_vars <- unique(c(bin_vars, intersect("ASCVD_any", names(df))))
disease_vars <- unique(c(disease_vars, intersect("ASCVD_any", bin_vars)))

# Keep only requested variables in requested order
requested_order <- c("Age", "arm", "BMI", if (length(hr_col) > 0) hr_col[1] else character(0),
                     "Sex", "HTN", "DM", "Af", "HF", "PAD", "ASCVD_any")

# -----------------------------
# 4) Attempt-level dataset
# -----------------------------
df_attempt <- df %>%
  transmute(
    patient_id = as.character(.data[[id_col]]),
    omron_error = to_binary01(.data[[omron_col]]),
    across(all_of(unique(c(cont_vars, bin_vars))), ~ .x)
  ) %>%
  mutate(patient_id = str_trim(patient_id)) %>%
  filter(!is.na(patient_id), patient_id != "")

attempt_tbl <- make_compare_table(
  dat = df_attempt,
  grp_col = "omron_error",
  cont_vars = cont_vars,
  bin_vars = bin_vars,
  disease_vars = disease_vars
) %>%
  mutate(variables = factor(variables, levels = requested_order)) %>%
  arrange(variables) %>%
  mutate(
    variables = as.character(variables),
    Variable = label_variable(variables)
  ) %>%
  select(-type, -variables) %>%
  rename(
    `No error` = `No Error`,
    `Difference (Error - No) [95% CI]` = `Difference (error-No)[95% CI ]`
  ) %>%
  relocate(Variable)

# -----------------------------
# 5) Patient-level dataset (any error)
# -----------------------------
df_patient <- df_attempt %>%
  group_by(patient_id) %>%
  summarise(
    any_omron_error = as.integer(any(omron_error == 1L, na.rm = TRUE)),
    across(all_of(cont_vars), ~ {
      x <- to_num(.x)
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else x[1]
    }),
    across(all_of(bin_vars), ~ {
      v <- dplyr::cur_column()
      x <- parse_binary_var(.x, v, disease_vars)
      if (v == "Sex") {
        x2 <- x[!is.na(x)]
        if (length(x2) == 0) NA_integer_ else x2[1]
      } else {
        if (all(is.na(x))) NA_integer_ else ifelse(any(x == 1, na.rm = TRUE), 1L, 0L)
      }
    }),
    .groups = "drop"
  )

patient_tbl <- make_compare_table(
  dat = df_patient %>% rename(omron_error = any_omron_error),
  grp_col = "omron_error",
  cont_vars = cont_vars,
  bin_vars = bin_vars,
  disease_vars = disease_vars
) %>%
  mutate(variables = factor(variables, levels = requested_order)) %>%
  arrange(variables) %>%
  mutate(
    variables = as.character(variables),
    Variable = label_variable(variables)
  ) %>%
  select(-type, -variables) %>%
  rename(
    `No error` = `No Error`,
    `Difference (Error - No) [95% CI]` = `Difference (error-No)[95% CI ]`
  ) %>%
  relocate(Variable)

# -----------------------------
# 6) Save outputs
# -----------------------------
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(attempt_tbl, "outputs/tables/table3-6omron_error_descripitive_attempt.csv",
  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(patient_tbl, "outputs/tables/table3-7omron_error_descripitive_patent.csv",
  row.names = FALSE, fileEncoding = "UTF-8")

summary_tbl <- tibble(
  sheet_used = sheet_pick,
  patient_id_col = id_col,
  omron_error_col = omron_col,
  n_attempt = nrow(df_attempt),
  n_patient = n_distinct(df_patient$patient_id),
  patient_any_error_n = sum(df_patient$any_omron_error == 1, na.rm = TRUE),
  patient_any_error_prop = mean(df_patient$any_omron_error == 1, na.rm = TRUE)
)
write.csv(summary_tbl, "outputs/tables/omron_error_descriptive_summary.csv",
  row.names = FALSE, fileEncoding = "UTF-8")

cat("\nSaved:\n")
cat("- outputs/tables/table3-6omron_error_descripitive_attempt.csv\n")
cat("- outputs/tables/table3-7omron_error_descripitive_patent.csv\n")
cat("- outputs/tables/omron_error_descriptive_summary.csv\n")
