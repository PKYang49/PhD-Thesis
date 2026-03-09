# analysis/18-subgroupdiaabs_hr100.R
# Outcome (DIA):
#   Y_absdiff_dia = abs(O-DIA - C-DIA) - abs(K-DIA - C-DIA)
# Interpretation:
#   Y > 0 => K has smaller absolute error than O (K better)
#   Y < 0 => O has smaller absolute error than K
#
# Subgroups: HTN, DM, other_ascvd, PAD, Af, HR clinical (<100, >=100)

library(dplyr)

if (!requireNamespace("lme4", quietly = TRUE)) stop("Need lme4. Run: renv::install('lme4')")
if (!requireNamespace("lmerTest", quietly = TRUE)) stop("Need lmerTest. Run: renv::install('lmerTest')")
if (!requireNamespace("emmeans", quietly = TRUE)) stop("Need emmeans. Run: renv::install('emmeans')")

library(lme4)
library(lmerTest)
library(emmeans)

path_long <- "data/processed/df_long_included.rds"
path_pat  <- "data/processed/df_patient_included.rds"

df_long <- readRDS(path_long)

need_bp <- c("O-DIA", "K-DIA", "C-DIA")
miss_bp <- setdiff(need_bp, names(df_long))
if (length(miss_bp) > 0) stop("Missing BP columns: ", paste(miss_bp, collapse = ", "))

hr_candidates <- c("k-hr", "K-HR", "k_hr", "K_HR", "KHR", "kHR", "k.hr", "K.hr", "HR_K", "hr_k")
hr_col <- intersect(hr_candidates, names(df_long))
if (length(hr_col) == 0) hr_col <- names(df_long)[grepl("k.*hr|hr.*k", names(df_long), ignore.case = TRUE)]
if (length(hr_col) == 0) stop("Cannot find HR column. Update hr_candidates based on names(df_long).")
hr_col <- hr_col[1]
message("Using HR column: ", hr_col)

arm_candidates <- c("arm", "ARM", "Arm")
arm_col <- intersect(arm_candidates, names(df_long))
if (length(arm_col) == 0) arm_col <- names(df_long)[grepl("^arm$", names(df_long), ignore.case = TRUE)]
arm_col <- if (length(arm_col) > 0) arm_col[1] else NA_character_
message("Using arm column: ", ifelse(is.na(arm_col), "<missing>", arm_col))

need_from_pat <- c("CAD", "AMI", "Stroke", "PAD")
miss_long <- setdiff(need_from_pat, names(df_long))
if (length(miss_long) > 0) {
  message("CAD/AMI/Stroke/PAD not complete in long data. Joining from patient-level...")
  df_pat <- readRDS(path_pat)
  miss_pat <- setdiff(need_from_pat, names(df_pat))
  if (length(miss_pat) > 0) stop("Missing in patient data: ", paste(miss_pat, collapse = ", "))
  df_long <- df_long %>% left_join(df_pat %>% select(ID, all_of(need_from_pat)), by = "ID")
}

is_yes <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr %in% c("Yes", "Y", "1", "TRUE", "True", "yes", "y", "true") | x %in% c(1, TRUE)
}

drop_na_cols <- function(df, cols) {
  keep <- complete.cases(df[, cols, drop = FALSE])
  df[keep, , drop = FALSE]
}

build_formula <- function(outcome, subgroup_var, adjusters) {
  rhs <- c(subgroup_var, adjusters)
  as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1|ID)"))
}

df_m <- df_long %>%
  mutate(
    ID = as.factor(ID),
    `O-DIA` = suppressWarnings(as.numeric(as.character(`O-DIA`))),
    `K-DIA` = suppressWarnings(as.numeric(as.character(`K-DIA`))),
    `C-DIA` = suppressWarnings(as.numeric(as.character(`C-DIA`))),
    Y_absdiff_dia = abs(`O-DIA` - `C-DIA`) - abs(`K-DIA` - `C-DIA`),
    Age = suppressWarnings(as.numeric(as.character(Age))),
    Sex = as.factor(Sex),
    arm = if (!is.na(arm_col)) suppressWarnings(as.numeric(as.character(.data[[arm_col]]))) else NA_real_,
    Af = as.factor(Af),
    HTN = as.factor(HTN),
    DM = as.factor(DM),
    PAD = as.factor(PAD),
    other_ascvd = factor(ifelse(is_yes(CAD) | is_yes(AMI) | is_yes(Stroke), "Yes", "No"), levels = c("No", "Yes")),
    HR = suppressWarnings(as.numeric(as.character(.data[[hr_col]])))
  ) %>%
  mutate(
    HR_c = HR - mean(HR, na.rm = TRUE),
    HR_group = cut(HR, breaks = c(-Inf, 100, Inf), right = FALSE, labels = c("<100", ">=100")),
    HR_group = factor(HR_group, levels = c("<100", ">=100"))
  )

run_subgroup_abs_dia <- function(data, subgroup_var, tag) {
  adjusters <- c("Age", "Sex", "arm", "Af", "HR_c", "HTN", "DM", "other_ascvd", "PAD")
  adjusters <- setdiff(adjusters, subgroup_var)
  if (subgroup_var == "HR_group") adjusters <- setdiff(adjusters, "HR_c")
  adjusters <- adjusters[adjusters %in% names(data) &
    !vapply(data[adjusters], function(x) all(is.na(x)), logical(1))]

  need <- unique(c("ID", "Y_absdiff_dia", subgroup_var, adjusters))
  d <- drop_na_cols(data, need)
  d[[subgroup_var]] <- as.factor(d[[subgroup_var]])

  m <- lmer(build_formula("Y_absdiff_dia", subgroup_var, adjusters), data = d, REML = TRUE)

  emm <- emmeans(m, as.formula(paste0("~ ", subgroup_var)))
  emm_df <- as.data.frame(emm) %>%
    mutate(subgroup = tag, level = as.character(.data[[subgroup_var]])) %>%
    transmute(
      subgroup, level,
      estimate = emmean,
      se = SE,
      df = df,
      lower = emmean + qt(0.025, df = df) * SE,
      upper = emmean + qt(0.975, df = df) * SE
    )

  con_df <- NULL
  if (length(levels(d[[subgroup_var]])) >= 2) {
    con <- contrast(emm, method = "trt.vs.ctrl", ref = 1)
    con_df <- as.data.frame(con) %>%
      mutate(subgroup = tag) %>%
      transmute(
        subgroup, contrast,
        estimate,
        se = SE,
        df,
        lower = estimate + qt(0.025, df = df) * SE,
        upper = estimate + qt(0.975, df = df) * SE,
        p = p.value
      )
  }

  list(means = emm_df, contrasts = con_df)
}

subgroups <- list(
  HTN = "HTN",
  DM = "DM",
  other_ascvd = "other_ascvd",
  PAD = "PAD",
  Af = "Af",
  HR_binary = "HR_group"
)

means_all <- bind_rows(lapply(names(subgroups), function(nm) {
  run_subgroup_abs_dia(df_m, subgroups[[nm]], nm)$means
}))
cons_all <- bind_rows(lapply(names(subgroups), function(nm) {
  run_subgroup_abs_dia(df_m, subgroups[[nm]], nm)$contrasts
}))

hr_adjusters <- c("Age", "Sex", "arm", "Af", "HTN", "DM", "other_ascvd", "PAD")
hr_adjusters <- hr_adjusters[hr_adjusters %in% names(df_m) &
  !vapply(df_m[hr_adjusters], function(x) all(is.na(x)), logical(1))]
d_hr <- drop_na_cols(df_m, c("ID", "Y_absdiff_dia", "HR_c", hr_adjusters))

m_hr <- lmer(build_formula("Y_absdiff_dia", "HR_c", hr_adjusters), data = d_hr, REML = TRUE)
cf_hr <- as.data.frame(summary(m_hr)$coefficients)
cf_hr$term <- rownames(cf_hr)
rownames(cf_hr) <- NULL

hr_row <- cf_hr %>%
  filter(term == "HR_c") %>%
  transmute(
    subgroup = "HR_continuous",
    level = "per 1 bpm",
    estimate = Estimate,
    se = `Std. Error`,
    df = NA_real_,
    lower = Estimate - 1.96 * `Std. Error`,
    upper = Estimate + 1.96 * `Std. Error`,
    p = `Pr(>|t|)`
  )

out_tables_dir <- file.path(getwd(), "outputs", "tables")
out_models_dir <- file.path(getwd(), "outputs", "models")
dir.create(out_tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_models_dir, showWarnings = FALSE, recursive = TRUE)

means_out <- file.path(out_tables_dir, "subgroup_absdiff_dia_means_HR100.csv")
cons_out <- file.path(out_tables_dir, "subgroup_absdiff_dia_contrasts_HR100.csv")
hr_out <- file.path(out_tables_dir, "subgroup_absdiff_dia_HRcontinuous_HR100.csv")
hr_txt_out <- file.path(out_models_dir, "subgroup_absdiff_dia_HRcontinuous_summary_HR100.txt")

write.csv(means_all, means_out, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cons_all, cons_out, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hr_row, hr_out, row.names = FALSE, fileEncoding = "UTF-8")

sink(hr_txt_out)
print(summary(m_hr))
sink()

message("Saved subgroup means: ", means_out)
message("Saved subgroup contrasts: ", cons_out)
message("Saved HR continuous row: ", hr_out)
message("Saved HR continuous summary: ", hr_txt_out)
