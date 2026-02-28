# analysis/17-omron_error_patient_logistic.R
# Patient-level logistic regression:
#   outcome: patient_any_error (0/1), any Omron error within patient
#   predictors: Af + HR + arm + age + sex
#   model: logistic (glm) + Firth logistic (logistf, if available)

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

# ---- paths ----
path_xlsx <- "data/raw/project_info_00004(with HR).xlsx"
if (!file.exists(path_xlsx)) stop("Raw file not found: ", path_xlsx)
omron_nonempty_is_one <- TRUE
af_nonempty_is_one <- TRUE

# ---- helpers ----
norm <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[[:space:]_\\-]+", "")
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

to_01 <- function(x, nonempty_is_one = FALSE) {
  x_chr <- str_trim(as.character(x))
  x_low <- str_to_lower(x_chr)
  x_num <- suppressWarnings(as.numeric(x_chr))
  out <- ifelse(
    !is.na(x_num), ifelse(x_num == 1, 1L, 0L),
    ifelse(x_low %in% c("yes", "y", "true", "t"), 1L,
      ifelse(x_low %in% c("no", "n", "false", "f", ""), 0L, NA_integer_)
    )
  )
  if (nonempty_is_one) {
    out[is.na(out) & !is.na(x_chr) & x_chr != ""] <- 1L
    out[is.na(out) & (is.na(x_chr) | x_chr == "")] <- 0L
  }
  out
}

to_yesno <- function(x, nonempty_is_one = FALSE) {
  z <- to_01(x, nonempty_is_one = nonempty_is_one)
  ifelse(z == 1L, "Yes", ifelse(z == 0L, "No", NA_character_))
}

first_non_na <- function(x) {
  idx <- which(!is.na(x) & as.character(x) != "")
  if (length(idx) == 0) NA else x[idx[1]]
}

# ---- 1) sheet detection ----
sheets <- excel_sheets(path_xlsx)
sheet_norm <- norm(sheets)
sheet_pick <- sheets[which(sheet_norm == norm("device analysis"))]
if (length(sheet_pick) == 0) {
  sheet_pick <- sheets[str_detect(sheet_norm, "device") & str_detect(sheet_norm, "analysis")]
}
if (length(sheet_pick) == 0) {
  stop("Cannot find sheet 'device analysis'. Available sheets: ", paste(sheets, collapse = ", "))
}
sheet_pick <- sheet_pick[1]
message("Using sheet: ", sheet_pick)

df <- read_excel(path_xlsx, sheet = sheet_pick)
names(df) <- str_trim(names(df))
nm <- names(df)
nm_norm <- norm(nm)

# ---- 2) detect columns ----
id_col <- NA_character_
id_priority <- c("ID", "SN", "patient id", "patient_id")
for (cand in id_priority) {
  idx <- which(nm_norm == norm(cand))
  if (length(idx) > 0) {
    id_col <- nm[idx[1]]
    break
  }
}
if (is.na(id_col)) stop("Cannot detect ID column in sheet: ", sheet_pick)

omron_idx <- which(nm_norm %in% c(norm("omron error"), norm("omron_error"), norm("omronerror")))
if (length(omron_idx) == 0) {
  omron_idx <- which(str_detect(nm_norm, "omron") & str_detect(nm_norm, "error"))
}
if (length(omron_idx) == 0) stop("Cannot detect omron_error column in sheet: ", sheet_pick)
omron_col <- nm[omron_idx[1]]

hr_idx <- which(nm_norm %in% c(norm("HR"), norm("K-HR"), norm("k_hr"), norm("KHR")))
if (length(hr_idx) == 0) {
  hr_idx <- which(str_detect(nm_norm, "hr"))
}
if (length(hr_idx) == 0) stop("Cannot detect HR column.")
hr_col <- nm[hr_idx[1]]

need_cov <- c("Age", "Sex", "arm", "Af")
miss_cov <- setdiff(need_cov, nm)
if (length(miss_cov) > 0) stop("Missing covariates: ", paste(miss_cov, collapse = ", "))

message("Using ID: ", id_col)
message("Using Omron error: ", omron_col)
message("Using HR: ", hr_col)

# ---- 3) build patient-level dataset ----
df_long <- df %>%
  transmute(
    ID = as.character(.data[[id_col]]),
    omron_error = to_01(.data[[omron_col]], nonempty_is_one = omron_nonempty_is_one),
    Age = to_num(Age),
    Sex = as.character(Sex),
    arm = to_num(arm),
    HR = to_num(.data[[hr_col]]),
    Af = to_yesno(Af, nonempty_is_one = af_nonempty_is_one)
  ) %>%
  filter(!is.na(ID), ID != "")

df_pt_raw <- df_long %>%
  group_by(ID) %>%
  summarise(
    patient_any_error = as.integer(any(omron_error == 1L, na.rm = TRUE)),
    Age = first_non_na(Age),
    Sex = first_non_na(Sex),
    arm = first_non_na(arm),
    HR = first_non_na(HR),
    Af = first_non_na(Af),
    .groups = "drop"
  ) %>%
  mutate(
    Sex = factor(Sex),
    Af = factor(Af, levels = c("No", "Yes"))
  )

missing_df <- data.frame(
  variable = c("patient_any_error", "Age", "Sex", "arm", "HR", "Af"),
  n_missing = c(
    sum(is.na(df_pt_raw$patient_any_error)),
    sum(is.na(df_pt_raw$Age)),
    sum(is.na(df_pt_raw$Sex)),
    sum(is.na(df_pt_raw$arm)),
    sum(is.na(df_pt_raw$HR)),
    sum(is.na(df_pt_raw$Af))
  ),
  stringsAsFactors = FALSE
)

df_pt <- df_pt_raw %>%
  drop_na(patient_any_error, Age, Sex, arm, HR, Af)

if (!all(df_pt$patient_any_error %in% c(0L, 1L))) {
  stop("patient_any_error is not strictly 0/1 after conversion.")
}

# ---- 4) save outputs + model/diagnostics ----
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)

outcome_tab <- table(df_pt$patient_any_error, useNA = "ifany")
if (length(outcome_tab) == 0) {
  outcome_df <- data.frame(
    patient_any_error = character(0),
    n = integer(0),
    stringsAsFactors = FALSE
  )
} else {
  outcome_df <- data.frame(
    patient_any_error = as.character(names(outcome_tab)),
    n = as.integer(outcome_tab),
    stringsAsFactors = FALSE
  )
}
write.csv(
  outcome_df,
  "outputs/tables/omron_error_patient_outcome_distribution.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  missing_df,
  "outputs/tables/omron_error_patient_missing_by_variable.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (n_distinct(df_pt$patient_any_error) < 2) {
  write.csv(
    df_pt,
    "outputs/tables/omron_error_patient_dataset_filtered.csv",
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  sink("outputs/models/omron_error_patient_logistic_summary.txt")
  cat("Patient-level logistic regression for any Omron error\n")
  cat("Sheet used:", sheet_pick, "\n")
  cat("ID column:", id_col, "\n")
  cat("Omron error column:", omron_col, "\n")
  cat("HR column:", hr_col, "\n")
  cat("omron_nonempty_is_one:", omron_nonempty_is_one, "\n")
  cat("af_nonempty_is_one:", af_nonempty_is_one, "\n")
  cat("N patient before drop_na:", nrow(df_pt_raw), "\n")
  cat("N patient:", nrow(df_pt), "\n")
  cat("Missing by variable before drop_na:\n")
  print(missing_df)
  cat("Outcome distribution (patient_any_error):\n")
  print(outcome_tab)
  cat("\nModel skipped: patient_any_error has only one class after filtering.\n")
  cat("Try checking Omron error coding or set omron_nonempty_is_one = TRUE if the source column stores non-empty strings for error.\n")
  sink()

  message("Outcome has one class; model skipped.")
  message("Saved distribution: outputs/tables/omron_error_patient_outcome_distribution.csv")
  message("Saved missing summary: outputs/tables/omron_error_patient_missing_by_variable.csv")
  message("Saved filtered dataset: outputs/tables/omron_error_patient_dataset_filtered.csv")
  message("Saved summary: outputs/models/omron_error_patient_logistic_summary.txt")
} else {
  df_pt <- df_pt %>%
    mutate(
      Age_z = as.numeric(scale(Age)),
      arm_z = as.numeric(scale(arm)),
      HR_z = as.numeric(scale(HR))
    )

  form <- patient_any_error ~ Af + HR_z + arm_z + Age_z + Sex

  # ---- 5) standard logistic ----
  m_glm <- glm(form, data = df_pt, family = binomial(link = "logit"))
  sm_glm <- summary(m_glm)
  cf_glm <- as.data.frame(sm_glm$coefficients)
  cf_glm$term <- rownames(cf_glm)
  rownames(cf_glm) <- NULL

  glm_tbl <- cf_glm %>%
    rename(
      estimate = Estimate,
      se = `Std. Error`,
      p = `Pr(>|z|)`
    ) %>%
    mutate(
      OR = exp(estimate),
      OR_low = exp(estimate - 1.96 * se),
      OR_high = exp(estimate + 1.96 * se),
      Predictor = case_when(
        term == "AfYes" ~ "Atrial fibrillation (Yes vs No)",
        term == "Age_z" ~ "Age (per 1 SD)",
        term == "HR_z" ~ "Heart rate (per 1 SD)",
        term == "arm_z" ~ "Arm circumference (per 1 SD)",
        str_starts(term, "Sex") ~ "Male (ref: Female)",
        term == "(Intercept)" ~ "Intercept",
        TRUE ~ term
      ),
      `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", OR, OR_low, OR_high),
      p = ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
    ) %>%
    filter(Predictor != "Intercept") %>%
    mutate(
      predictor_order = match(
        Predictor,
        c(
          "Atrial fibrillation (Yes vs No)",
          "Age (per 1 SD)",
          "Heart rate (per 1 SD)",
          "Arm circumference (per 1 SD)",
          "Male (ref: Female)"
        )
      )
    ) %>%
    arrange(predictor_order) %>%
    select(Predictor, `OR (95% CI)`, p)

  # ---- 6) Firth logistic (if available) ----
  has_logistf <- requireNamespace("logistf", quietly = TRUE)
  firth_tbl <- NULL
  m_firth <- NULL

  if (has_logistf) {
    m_firth <- logistf::logistf(form, data = df_pt)

    firth_tbl <- data.frame(
      term = names(m_firth$coefficients),
      estimate = as.numeric(m_firth$coefficients),
      se = as.numeric(m_firth$se),
      p = as.numeric(m_firth$prob),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        OR = exp(estimate),
        OR_low = exp(estimate - 1.96 * se),
        OR_high = exp(estimate + 1.96 * se),
        estimate = round(estimate, 4),
        se = round(se, 4),
        OR = round(OR, 4),
        OR_low = round(OR_low, 4),
        OR_high = round(OR_high, 4),
        p = signif(p, 3),
        model = "firth_logistf",
        n_patient = nrow(df_pt),
        n_event = sum(df_pt$patient_any_error == 1L)
      ) %>%
      select(model, term, estimate, se, p, OR, OR_low, OR_high, n_patient, n_event)
  }

  write.csv(
    glm_tbl,
    "outputs/tables/table3-9omron_error_logistic_patient.csv",
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  if (!is.null(firth_tbl)) {
    write.csv(
      firth_tbl,
      "outputs/tables/omron_error_patient_logistic_firth_coef.csv",
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  sink("outputs/models/omron_error_patient_logistic_summary.txt")
  cat("Patient-level logistic regression for any Omron error\n")
  cat("Sheet used:", sheet_pick, "\n")
  cat("ID column:", id_col, "\n")
  cat("Omron error column:", omron_col, "\n")
  cat("HR column:", hr_col, "\n")
  cat("omron_nonempty_is_one:", omron_nonempty_is_one, "\n")
  cat("af_nonempty_is_one:", af_nonempty_is_one, "\n")
  cat("N patient before drop_na:", nrow(df_pt_raw), "\n")
  cat("N patient:", nrow(df_pt), "\n")
  cat("Missing by variable before drop_na:\n")
  print(missing_df)
  cat("N event (patient_any_error=1):", sum(df_pt$patient_any_error == 1L), "\n")
  cat("Outcome distribution (patient_any_error):\n")
  print(outcome_tab)
  cat("\nFormula:", deparse(form), "\n\n")

  cat("== Standard logistic (glm) ==\n")
  print(summary(m_glm))
  cat("\n")

  if (!is.null(m_firth)) {
    cat("== Firth logistic (logistf) ==\n")
    print(summary(m_firth))
    cat("\n")
  } else {
    cat("logistf package not installed; Firth logistic skipped.\n\n")
  }
  sink()

  message("Saved GLM coefficients: outputs/tables/table3-9omron_error_logistic_patient.csv")
  if (!is.null(firth_tbl)) {
    message("Saved Firth coefficients: outputs/tables/omron_error_patient_logistic_firth_coef.csv")
  }
  message("Saved distribution: outputs/tables/omron_error_patient_outcome_distribution.csv")
  message("Saved missing summary: outputs/tables/omron_error_patient_missing_by_variable.csv")
  message("Saved summary: outputs/models/omron_error_patient_logistic_summary.txt")
}
