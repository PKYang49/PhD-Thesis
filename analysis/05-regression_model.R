# analysis/05-regression_model.R
# Mixed-effects linear models for:
#   O-SYS-C-SYS, K-SYS-C-SYS, O-DIA-C-DIA, K-DIA-C-DIA
# Fixed effects:
#   Age + Sex + Af + HR_c + HTN + DM + ASCVD_any
# Random effect:
#   (1 | ID)
#
# Outputs:
#   1) outputs/tables/regression_mixed_4models_coefficients.csv
#   2) outputs/tables/table3-4_sbp_mixedmodel.csv
#   3) outputs/tables/table 3-5_dbp_mixedmodel.csv
#   4) outputs/models/*_summary.txt

library(dplyr)
library(tidyr)

if (!requireNamespace("lme4", quietly = TRUE)) stop("Need lme4. Run: renv::install('lme4')")
if (!requireNamespace("lmerTest", quietly = TRUE)) stop("Need lmerTest. Run: renv::install('lmerTest')")

library(lme4)
library(lmerTest)

path_long <- "data/processed/df_long_included.rds"
if (!file.exists(path_long)) stop("File not found: ", path_long)

df_long <- readRDS(path_long)

need_bp <- c("O-SYS", "K-SYS", "C-SYS", "O-DIA", "K-DIA", "C-DIA")
miss_bp <- setdiff(need_bp, names(df_long))
if (length(miss_bp) > 0) stop("Missing BP columns: ", paste(miss_bp, collapse = ", "))

need_cov <- c("ID", "Age", "Sex", "Af", "HTN", "DM", "CAD", "AMI", "Stroke", "PAD")
miss_cov <- setdiff(need_cov, names(df_long))
if (length(miss_cov) > 0) stop("Missing covariate columns: ", paste(miss_cov, collapse = ", "))

hr_candidates <- c("K-HR", "k-hr", "K_HR", "k_hr", "KHR", "kHR", "k.hr", "K.hr")
hr_col <- intersect(hr_candidates, names(df_long))
if (length(hr_col) == 0) {
  hr_col <- names(df_long)[grepl("k.*hr|hr.*k", names(df_long), ignore.case = TRUE)]
}
if (length(hr_col) == 0) stop("Cannot detect HR column for HR_c.")
hr_col <- hr_col[1]
message("Using HR column: ", hr_col)

is_yes <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr %in% c("Yes", "Y", "1", "TRUE", "True", "yes", "y", "true", "☑") | x %in% c(1, TRUE)
}

to_yesno_factor <- function(x) {
  factor(ifelse(is_yes(x), "Yes", "No"), levels = c("No", "Yes"))
}

fmt_num <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), x)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    sprintf("%.2f", p)
  )
}

label_key <- function(term) {
  dplyr::case_when(
    term == "Age" ~ "age",
    grepl("^Sex", term) ~ "sex",
    grepl("^Af", term) ~ "af",
    term == "HR_c" ~ "hr",
    grepl("^HTN", term) ~ "htn",
    grepl("^DM", term) ~ "dm",
    grepl("^ASCVD_any", term) ~ "ascvd",
    TRUE ~ NA_character_
  )
}

label_text <- c(
  age = "Age (per year)",
  sex = "Male sex",
  af = "Af",
  hr = "HR (per bpm)",
  htn = "HTN",
  dm = "DM",
  ascvd = "ASCVD"
)

label_order <- c("age", "sex", "af", "hr", "htn", "dm", "ascvd")

df_m <- df_long %>%
  mutate(
    ID = as.factor(as.character(ID)),
    Age = suppressWarnings(as.numeric(as.character(Age))),
    Sex = as.factor(as.character(Sex)),
    Af = to_yesno_factor(Af),
    HTN = to_yesno_factor(HTN),
    DM = to_yesno_factor(DM),
    ASCVD_any = factor(
      ifelse(is_yes(CAD) | is_yes(AMI) | is_yes(Stroke) | is_yes(PAD), "Yes", "No"),
      levels = c("No", "Yes")
    ),
    HR = suppressWarnings(as.numeric(as.character(.data[[hr_col]]))),
    HR_c = HR - mean(HR, na.rm = TRUE),
    Y_OSYS_CSYS = suppressWarnings(as.numeric(as.character(`O-SYS`))) -
      suppressWarnings(as.numeric(as.character(`C-SYS`))),
    Y_KSYS_CSYS = suppressWarnings(as.numeric(as.character(`K-SYS`))) -
      suppressWarnings(as.numeric(as.character(`C-SYS`))),
    Y_ODIA_CDIA = suppressWarnings(as.numeric(as.character(`O-DIA`))) -
      suppressWarnings(as.numeric(as.character(`C-DIA`))),
    Y_KDIA_CDIA = suppressWarnings(as.numeric(as.character(`K-DIA`))) -
      suppressWarnings(as.numeric(as.character(`C-DIA`)))
  )

model_specs <- tibble::tribble(
  ~model_name, ~outcome,
  "O-SYS-C-SYS", "Y_OSYS_CSYS",
  "K-SYS-C-SYS", "Y_KSYS_CSYS",
  "O-DIA-C-DIA", "Y_ODIA_CDIA",
  "K-DIA-C-DIA", "Y_KDIA_CDIA"
)

fit_one_model <- function(data, outcome, model_name) {
  need <- c("ID", outcome, "Age", "Sex", "Af", "HR_c", "HTN", "DM", "ASCVD_any")
  d <- data %>% drop_na(all_of(need))

  form <- as.formula(
    paste0(outcome, " ~ Age + Sex + Af + HR_c + HTN + DM + ASCVD_any + (1|ID)")
  )

  m <- lmer(form, data = d, REML = TRUE)
  sm <- summary(m)

  cf <- as.data.frame(sm$coefficients)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  has_p <- "Pr(>|t|)" %in% names(cf)

  cf <- cf %>%
    transmute(
      model = model_name,
      term = term,
      estimate = Estimate,
      se = `Std. Error`,
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      p = if (has_p) `Pr(>|t|)` else NA_real_,
      n = nrow(d),
      n_id = dplyr::n_distinct(d$ID)
    )

  list(model = m, coef = cf, n = nrow(d), n_id = dplyr::n_distinct(d$ID))
}

build_formatted_table <- function(coef_df, model_left, model_right, left_label, right_label) {
  part <- coef_df %>%
    filter(model %in% c(model_left, model_right)) %>%
    mutate(
      key = label_key(term)
    ) %>%
    filter(!is.na(key)) %>%
    mutate(
      est_ci = paste0(fmt_num(estimate, 2), " (", fmt_num(lower, 2), ", ", fmt_num(upper, 2), ")"),
      p_fmt = fmt_p(p)
    ) %>%
    select(model, key, est_ci, p_fmt)

  left_df <- part %>%
    filter(model == model_left) %>%
    select(key, left_est_ci = est_ci, left_p = p_fmt)

  right_df <- part %>%
    filter(model == model_right) %>%
    select(key, right_est_ci = est_ci, right_p = p_fmt)

  out <- tibble(key = label_order) %>%
    left_join(left_df, by = "key") %>%
    left_join(right_df, by = "key") %>%
    mutate(Variable = unname(label_text[key])) %>%
    select(Variable, left_est_ci, left_p, right_est_ci, right_p)

  out_df <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(out_df) <- c("Variable", left_label, "p", right_label, "p")
  out_df
}

fits <- lapply(seq_len(nrow(model_specs)), function(i) {
  fit_one_model(df_m, model_specs$outcome[i], model_specs$model_name[i])
})

coef_all <- bind_rows(lapply(fits, `[[`, "coef"))

sbp_table <- build_formatted_table(
  coef_df = coef_all,
  model_left = "O-SYS-C-SYS",
  model_right = "K-SYS-C-SYS",
  left_label = "Oscillo-SBP (95% CI)",
  right_label = "Ksens-SBP (95% CI)"
)

dbp_table <- build_formatted_table(
  coef_df = coef_all,
  model_left = "O-DIA-C-DIA",
  model_right = "K-DIA-C-DIA",
  left_label = "Oscillo-DBP (95% CI)",
  right_label = "Ksens-DBP (95% CI)"
)

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)

write.csv(
  coef_all,
  "outputs/tables/regression_mixed_4models_coefficients.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sbp_table,
  "outputs/tables/table3-4_sbp_mixedmodel.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  dbp_table,
  "outputs/tables/table 3-5_dbp_mixedmodel.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

for (i in seq_len(nrow(model_specs))) {
  model_name <- model_specs$model_name[i]
  model_file_tag <- gsub("[^A-Za-z0-9]+", "_", tolower(model_name))
  sink(file.path("outputs/models", paste0(model_file_tag, "_summary.txt")))
  cat("Model:", model_name, "\n")
  cat("Outcome:", model_specs$outcome[i], "\n")
  cat("Rows used:", fits[[i]]$n, "\n")
  cat("Unique ID:", fits[[i]]$n_id, "\n\n")
  print(summary(fits[[i]]$model))
  sink()
}

message("Saved coefficients: outputs/tables/regression_mixed_4models_coefficients.csv")
message("Saved SBP table: outputs/tables/table3-4_sbp_mixedmodel.csv")
message("Saved DBP table: outputs/tables/table 3-5_dbp_mixedmodel.csv")
message("Saved model summaries: outputs/models/*_summary.txt")
