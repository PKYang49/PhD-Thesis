# analysis/06-subgroupdiadiff.R
# Outcome (DIA):
#   Y_diff_dia = (O-DIA - C-DIA) - (K-DIA - C-DIA)
# Interpretation:
#   Y > 0 => K is closer to C than O (K better)
#   Y < 0 => O is closer to C than K
#
# Subgroups: HTN, DM, other_ascvd, Af, HR clinical (<60, 60-100, >100)
# Mixed model (random intercept per ID):
#   Y ~ subgroup + (adjusters) + (1|ID)
# Then use emmeans to estimate subgroup-specific mean(Y) and contrasts.

library(dplyr)
library(tidyr)

if (!requireNamespace("lme4", quietly = TRUE)) stop("Need lme4. Run: renv::install('lme4')")
if (!requireNamespace("lmerTest", quietly = TRUE)) stop("Need lmerTest. Run: renv::install('lmerTest')")
if (!requireNamespace("emmeans", quietly = TRUE)) stop("Need emmeans. Run: renv::install('emmeans')")

library(lme4)
library(lmerTest)
library(emmeans)

# ---- paths ----
path_long <- "data/processed/df_long_included.rds"
path_pat  <- "data/processed/df_patient_included.rds"

# ---- 1) load ----
df_long <- readRDS(path_long)

# ---- 2) check BP columns ----
need_bp <- c("O-DIA", "K-DIA", "C-DIA")
miss_bp <- setdiff(need_bp, names(df_long))
if (length(miss_bp) > 0) stop("Missing BP columns: ", paste(miss_bp, collapse = ", "))

# ---- 3) detect HR column ----
hr_candidates <- c("k-hr", "K-HR", "k_hr", "K_HR", "KHR", "kHR", "k.hr", "K.hr", "HR_K", "hr_k")
hr_col <- intersect(hr_candidates, names(df_long))
if (length(hr_col) == 0) hr_col <- names(df_long)[grepl("k.*hr|hr.*k", names(df_long), ignore.case = TRUE)]
if (length(hr_col) == 0) stop("Cannot find HR column (k-hr). Update hr_candidates based on names(df_long).")
hr_col <- hr_col[1]
message("Using HR column: ", hr_col)

# ---- 4) ensure CAD/AMI/Stroke/PAD exist ----
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

# ---- 5) build modeling dataset ----
df_m <- df_long %>%
  mutate(
    ID = as.factor(ID),

    # numeric BP
    `O-DIA` = suppressWarnings(as.numeric(as.character(`O-DIA`))),
    `K-DIA` = suppressWarnings(as.numeric(as.character(`K-DIA`))),
    `C-DIA` = suppressWarnings(as.numeric(as.character(`C-DIA`))),

    # outcome
    Y_diff_dia = (`O-DIA` - `C-DIA`) - (`K-DIA` - `C-DIA`),

    # covariates
    Age = suppressWarnings(as.numeric(as.character(Age))),
    Sex = as.factor(Sex),
    arm = suppressWarnings(as.numeric(as.character(arm))),
    Af = as.factor(Af),
    HTN = as.factor(HTN),
    DM = as.factor(DM),
    PAD = as.factor(PAD),
    other_ascvd = factor(ifelse(is_yes(CAD) | is_yes(AMI) | is_yes(Stroke), "Yes", "No"), levels = c("No", "Yes")),
    HR = suppressWarnings(as.numeric(as.character(.data[[hr_col]])))
  ) %>%
  mutate(
    HR_c = HR - mean(HR, na.rm = TRUE)
  )

# HR clinical groups: <60, 60-100, >100 bpm
df_m <- df_m %>%
  mutate(
    HR_group = cut(
      HR,
      breaks = c(-Inf, 60, 100, Inf),
      right = FALSE,
      labels = c("<60", "60-100", ">100")
    ),
    HR_group = factor(HR_group, levels = c("<60", "60-100", ">100"))
  )

# ---- 6) helper: fit subgroup model + emmeans ----
run_subgroup_dia <- function(data, subgroup_var, tag) {
  adjusters <- c("Age", "Sex", "arm", "Af", "HR_c", "HTN", "DM", "other_ascvd", "PAD")
  adjusters <- setdiff(adjusters, subgroup_var)
  if (subgroup_var == "HR_group") adjusters <- setdiff(adjusters, "HR_c")

  need <- unique(c("ID", "Y_diff_dia", subgroup_var, adjusters))
  d <- data %>% drop_na(all_of(need))
  d[[subgroup_var]] <- as.factor(d[[subgroup_var]])

  form <- as.formula(
    paste0("Y_diff_dia ~ ", subgroup_var, " + ", paste(adjusters, collapse = " + "), " + (1|ID)")
  )
  m <- lmer(form, data = d, REML = TRUE)

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
  lv <- levels(d[[subgroup_var]])
  if (length(lv) >= 2) {
    con <- contrast(emm, method = "trt.vs.ctrl", ref = 1)
    con_df <- as.data.frame(con) %>%
      mutate(subgroup = tag) %>%
      transmute(
        subgroup,
        contrast = contrast,
        estimate = estimate,
        se = SE,
        df = df,
        lower = estimate + qt(0.025, df = df) * SE,
        upper = estimate + qt(0.975, df = df) * SE,
        p = p.value
      )
  }

  list(means = emm_df, contrasts = con_df)
}

# ---- 7) run subgroup analyses ----
subgroups <- list(
  HTN = "HTN",
  DM = "DM",
  other_ascvd = "other_ascvd",
  PAD = "PAD",
  Af = "Af",
  HR_tertiles = "HR_group"
)

means_list <- list()
cons_list <- list()

for (nm in names(subgroups)) {
  v <- subgroups[[nm]]
  out <- run_subgroup_dia(df_m, v, nm)
  means_list[[nm]] <- out$means
  cons_list[[nm]] <- out$contrasts
}

means_all <- bind_rows(means_list)
cons_all <- bind_rows(cons_list)

# ---- 8) HR continuous effect (separate model) ----
d_hr <- df_m %>%
  drop_na(ID, Y_diff_dia, HR_c, Age, Sex, arm, Af, HTN, DM, other_ascvd, PAD)

m_hr <- lmer(
  Y_diff_dia ~ HR_c + Age + Sex + arm + Af + HTN + DM + other_ascvd + PAD + (1|ID),
  data = d_hr, REML = TRUE
)

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

# ---- 9) save outputs ----
out_tables_dir <- file.path(getwd(), "outputs", "tables")
out_models_dir <- file.path(getwd(), "outputs", "models")

if (file.exists(out_tables_dir) && !dir.exists(out_tables_dir)) {
  stop("Path exists but is not a directory: ", out_tables_dir)
}
if (file.exists(out_models_dir) && !dir.exists(out_models_dir)) {
  stop("Path exists but is not a directory: ", out_models_dir)
}

dir.create(out_tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_models_dir, showWarnings = FALSE, recursive = TRUE)

if (file.access(out_tables_dir, 2) != 0) {
  stop("No write permission for directory: ", out_tables_dir)
}
if (file.access(out_models_dir, 2) != 0) {
  stop("No write permission for directory: ", out_models_dir)
}

write.csv(means_all, file.path(out_tables_dir, "subgroup_diff_dia_means.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cons_all, file.path(out_tables_dir, "subgroup_diff_dia_contrasts.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(hr_row, file.path(out_tables_dir, "subgroup_diff_dia_HRcontinuous.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")

sink(file.path(out_models_dir, "subgroup_diff_dia_HRcontinuous_summary.txt"))
print(summary(m_hr))
sink()

message("Saved subgroup means: outputs/tables/subgroup_diff_dia_means.csv")
message("Saved subgroup contrasts: outputs/tables/subgroup_diff_dia_contrasts.csv")
message("Saved HR continuous row: outputs/tables/subgroup_diff_dia_HRcontinuous.csv")
message("Saved HR continuous summary: outputs/models/subgroup_diff_dia_HRcontinuous_summary.txt")
