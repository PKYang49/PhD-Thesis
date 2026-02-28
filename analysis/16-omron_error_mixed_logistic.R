# analysis/16-omron_error_mixed_logistic.R
# Mixed-effects logistic regression:
#   outcome: omron_error (0/1)
#   predictors: Age + Sex + arm + K-HR + Af + PAD
#   random effect: (1 | ID)

library(readxl)
library(dplyr)
library(stringr)
library(lme4)

# ---- paths ----
path_xlsx <- "data/raw/project_info_00004(with HR).xlsx"
if (!file.exists(path_xlsx)) stop("Raw file not found: ", path_xlsx)

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

to_yesno_factor <- function(x) {
  z <- to_01(x, nonempty_is_one = TRUE)
  factor(ifelse(z == 1, "Yes", "No"), levels = c("No", "Yes"))
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

kh_candidates <- c("K-HR", "k-hr", "k_hr", "K_HR", "KHR")
kh_col <- kh_candidates[kh_candidates %in% nm]
if (length(kh_col) == 0) {
  kh_col <- nm[str_detect(nm_norm, "k.*hr|hr.*k")]
}
if (length(kh_col) == 0) stop("Cannot detect K-HR column.")
kh_col <- kh_col[1]

need_cov <- c("Age", "Sex", "arm", "Af", "PAD")
miss_cov <- setdiff(need_cov, nm)
if (length(miss_cov) > 0) stop("Missing covariates: ", paste(miss_cov, collapse = ", "))

message("Using ID: ", id_col)
message("Using outcome: ", omron_col)
message("Using K-HR: ", kh_col)

# ---- 3) build modeling dataset ----
df_m <- df %>%
  transmute(
    ID = as.factor(as.character(.data[[id_col]])),
    omron_error = to_01(.data[[omron_col]], nonempty_is_one = TRUE),
    Age = to_num(Age),
    Sex = as.factor(as.character(Sex)),
    arm = to_num(arm),
    K_HR = to_num(.data[[kh_col]]),
    Af = to_yesno_factor(Af),
    PAD = to_yesno_factor(PAD)
  ) %>%
  filter(!is.na(ID), as.character(ID) != "") %>%
  drop_na(omron_error, Age, Sex, arm, K_HR, Af, PAD)

# Rescale continuous predictors to improve numerical stability
df_m <- df_m %>%
  mutate(
    Age_z = as.numeric(scale(Age)),
    arm_z = as.numeric(scale(arm)),
    K_HR_z = as.numeric(scale(K_HR))
  )

if (!all(df_m$omron_error %in% c(0L, 1L))) {
  stop("omron_error is not strictly 0/1 after conversion.")
}

if (n_distinct(df_m$omron_error) < 2) {
  cat("\nOutcome distribution after filtering:\n")
  print(table(df_m$omron_error, useNA = "ifany"))
  stop("omron_error has only one class after NA filtering; model cannot be fit.")
}

# ---- 4) fit model ----
form <- omron_error ~ Age_z + Sex + arm_z + K_HR_z + Af + PAD + (1 | ID)
m <- glmer(
  form,
  data = df_m,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# ---- 5) outputs ----
sm <- summary(m)
cf <- as.data.frame(sm$coefficients)
cf$term <- rownames(cf)
rownames(cf) <- NULL

coef_tbl <- cf %>%
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
      term == "K_HR_z" ~ "Heart rate (per 1 SD)",
      term == "arm_z" ~ "Arm circumference (per 1 SD)",
      term == "SexM" ~ "Male (ref: Female)",
      term == "PADYes" ~ "Peripheral artery disease (Yes vs No)",
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
        "Male (ref: Female)",
        "Peripheral artery disease (Yes vs No)"
      )
    )
  ) %>%
  arrange(predictor_order) %>%
  select(Predictor, `OR (95% CI)`, p)

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)

write.csv(
  coef_tbl,
  "outputs/tables/table3-8omron_error_logistic_attempt.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

sink("outputs/models/omron_error_mixed_logistic_summary.txt")
cat("Mixed-effects logistic regression for omron_error\n")
cat("Sheet used:", sheet_pick, "\n")
cat("ID column:", id_col, "\n")
cat("Outcome column:", omron_col, "\n")
cat("K-HR column:", kh_col, "\n")
cat("Rows used:", nrow(df_m), "\n")
cat("Unique ID:", n_distinct(df_m$ID), "\n\n")
print(form)
cat("\n")
print(summary(m))
sink()

message("Saved coefficients: outputs/tables/table3-8omron_error_logistic_attempt.csv")
message("Saved summary: outputs/models/omron_error_mixed_logistic_summary.txt")
