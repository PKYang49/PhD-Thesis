# analysis/03-table3-2_ccc.R
# Table 3-2: Lin's CCC (95% CI)
# Pairs:
# 1) O-SYS vs C-SYS
# 2) K-SYS vs C-SYS
# 3) O-DIA vs C-DIA
# 4) K-DIA vs C-DIA

library(dplyr)
library(tidyr)

if (!requireNamespace("epiR", quietly = TRUE)) {
  stop("Package 'epiR' is not installed. Run: renv::install('epiR')")
}

# ---- 1) Load cleaned long-form data (each measurement = one row) ----
df <- readRDS("data/processed/df_long_included.rds")

# ---- 2) Confirmed BP column names ----
col_map <- list(
  o_sys = "O-SYS",
  k_sys = "K-SYS",
  c_sys = "C-SYS",
  o_dia = "O-DIA",
  k_dia = "K-DIA",
  c_dia = "C-DIA"
)

missing_cols <- setdiff(unlist(col_map), names(df))
if (length(missing_cols) > 0) {
  stop("Missing columns in df: ", paste(missing_cols, collapse = ", "))
}

# ---- 3) Ensure numeric ----
df_num <- df %>%
  mutate(across(all_of(unlist(col_map)), ~ suppressWarnings(as.numeric(as.character(.x)))))

# ---- 4) CCC function for epiR 2.0.90 ----
calc_ccc <- function(data, x, y) {
  d <- data %>% select(all_of(c(x, y))) %>% drop_na()
  n <- nrow(d)

  if (n < 3) {
    return(tibble(n = n, ccc = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  res <- epiR::epi.ccc(d[[x]], d[[y]], ci = "z-transform", conf.level = 0.95)

  tibble(
    n = n,
    ccc = as.numeric(res$rho.c$est[1]),
    ci_low = as.numeric(res$rho.c$lower[1]),
    ci_high = as.numeric(res$rho.c$upper[1])
  )
}

# ---- 5) Pairs (Table 3-2) ----
pairs <- tibble(
  `device comparison` = c("O-SYS vs C-SYS", "K-SYS vs C-SYS", "O-DIA vs C-DIA", "K-DIA vs C-DIA"),
  x = c(col_map$o_sys, col_map$k_sys, col_map$o_dia, col_map$k_dia),
  y = c(col_map$c_sys, col_map$c_sys, col_map$c_dia, col_map$c_dia)
)

# ---- 6) Compute ----
res_tbl <- pairs %>%
  mutate(result = purrr::map2(x, y, ~ calc_ccc(df_num, .x, .y))) %>%
  tidyr::unnest(result) %>%
  mutate(
    `CCC (95% CI)` = ifelse(
      is.na(ccc),
      NA_character_,
      sprintf("%.3f [%.3f, %.3f]", ccc, ci_low, ci_high)
    )
  ) %>%
  select(`device comparison`, `CCC (95% CI)`)

print(res_tbl)

# ---- 7) Save ----
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(res_tbl, "outputs/tables/table3-2_ccc.csv", row.names = FALSE, fileEncoding = "UTF-8")

message("Saved: outputs/tables/table3-2_ccc.csv")
