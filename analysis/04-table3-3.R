# analysis/04-table3-3.R
# Goal:
# 1) For each pair vs C, compute:
#    - mean difference (x - y)
#    - mean absolute difference
#    - SD of differences
# 2) P-values:
#    - SYS: compare |O-SYS - C-SYS| vs |K-SYS - C-SYS| (paired Wilcoxon)
#    - DIA: compare |O-DIA - C-DIA| vs |K-DIA - C-DIA| (paired Wilcoxon)
#
# Output: outputs/tables/diff_summary_and_OvsK_abs_p.csv

library(dplyr)
library(tidyr)

# ---- 1) Load long-form data ----
df <- readRDS("data/processed/df_long_included.rds")

# ---- 2) Confirmed column names ----
cols <- c("O-SYS", "K-SYS", "C-SYS", "O-DIA", "K-DIA", "C-DIA")
missing_cols <- setdiff(cols, names(df))
if (length(missing_cols) > 0) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

# Ensure numeric
df <- df %>%
  mutate(across(all_of(cols), ~ suppressWarnings(as.numeric(as.character(.x)))))

# ---- 3) Helper: summary for one difference vector ----
summ_diff <- function(diff_vec) {
  tibble(
    n = sum(!is.na(diff_vec)),
    mean_diff = mean(diff_vec, na.rm = TRUE),
    mean_abs_diff = mean(abs(diff_vec), na.rm = TRUE),
    sd_diff = sd(diff_vec, na.rm = TRUE)
  )
}

# ---- 4) Compute pairwise differences (x - C) ----
df_diff <- df %>%
  transmute(
    dO_sys = `O-SYS` - `C-SYS`,
    dK_sys = `K-SYS` - `C-SYS`,
    dO_dia = `O-DIA` - `C-DIA`,
    dK_dia = `K-DIA` - `C-DIA`
  )

# ---- 5) Per-pair summaries ----
pair_summary <- bind_rows(
  summ_diff(df_diff$dO_sys) %>% mutate(comparison = "O-SYS vs C-SYS", bp_type = "SYS"),
  summ_diff(df_diff$dK_sys) %>% mutate(comparison = "K-SYS vs C-SYS", bp_type = "SYS"),
  summ_diff(df_diff$dO_dia) %>% mutate(comparison = "O-DIA vs C-DIA", bp_type = "DIA"),
  summ_diff(df_diff$dK_dia) %>% mutate(comparison = "K-DIA vs C-DIA", bp_type = "DIA")
) %>%
  relocate(comparison, bp_type)

# ---- 6) P-values: compare absolute differences O-C vs K-C ----
a_sys <- df_diff %>%
  transmute(aO_sys = abs(dO_sys), aK_sys = abs(dK_sys)) %>%
  drop_na()
p_sys_abs <- if (nrow(a_sys) >= 3) {
  wilcox.test(a_sys$aO_sys, a_sys$aK_sys, paired = TRUE, exact = FALSE, correct = FALSE)$p.value
} else {
  NA_real_
}

a_dia <- df_diff %>%
  transmute(aO_dia = abs(dO_dia), aK_dia = abs(dK_dia)) %>%
  drop_na()
p_dia_abs <- if (nrow(a_dia) >= 3) {
  wilcox.test(a_dia$aO_dia, a_dia$aK_dia, paired = TRUE, exact = FALSE, correct = FALSE)$p.value
} else {
  NA_real_
}

p_by_type <- tibble(
  bp_type = c("SYS", "DIA"),
  p_wilcox_abs = c(p_sys_abs, p_dia_abs)
) %>%
  mutate(p_wilcox_abs = signif(p_wilcox_abs, 3))

# ---- 7) Final output table ----
fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ NA_character_,
    p < 0.001 ~ "< 0.001",
    TRUE      ~ formatC(p, digits = 3, format = "f")
  )
}

out <- pair_summary %>%
  left_join(p_by_type, by = "bp_type") %>%
  mutate(
    Comparison = dplyr::case_when(
      comparison == "O-SYS vs C-SYS" ~ "Oscillo-SYS \u2013 Auscl-SYS",
      comparison == "K-SYS vs C-SYS" ~ "Ksens-SYS \u2013 Auscl-SYS",
      comparison == "O-DIA vs C-DIA" ~ "Oscillo-DIA \u2013 Auscl-DIA",
      comparison == "K-DIA vs C-DIA" ~ "Ksens-DIA \u2013 Auscl-DIA"
    ),
    `Mean of absolute differences` = round(mean_abs_diff, 1),
    `Mean differences`             = round(mean_diff, 1),
    SD                             = round(sd_diff, 1),
    `P value` = dplyr::if_else(
      startsWith(comparison, "K-"),
      fmt_p(p_wilcox_abs),
      NA_character_
    )
  ) %>%
  select(Comparison, `Mean of absolute differences`, `Mean differences`, SD, `P value`)

print(out, n = 50)

dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "outputs/tables/table3-3.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

message("Saved: outputs/tables/table3-3.csv")
