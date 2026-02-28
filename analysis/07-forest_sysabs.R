# analysis/12-forest_SYS_absdiff_subgroups.R
# Layout like your screenshot:
# left table (labels + mean + CI) + forest plot + p for interaction column
# Subgroups fixed: HTN/DM/ASCVD_any/Af (Yes/No) + HR clinical (<60, 60-100, >100)
#
# Requested tweaks:
# - HR labels: "HR < 60", "HR 60-100", "HR > 100"
# - Forest x-range fixed to -4 ~ +4
# - Under the Forest plot title:
#     + side labeled "Favor Ksens-SBP"
#     - side labeled "Favor Oscillo-SBP"
#
# Inputs expected:
#   outputs/tables/subgroup_absdiff_sys_means.csv
#   outputs/tables/subgroup_absdiff_sys_contrasts.csv

library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(readr)

means_path <- "outputs/tables/subgroup_absdiff_sys_means.csv"
cons_path  <- "outputs/tables/subgroup_absdiff_sys_contrasts.csv"

means <- read_csv(means_path, show_col_types = FALSE)
cons  <- read_csv(cons_path,  show_col_types = FALSE)

# ------------------------------------------------------------
# 1) Normalize subgroup names (your files may still say HR_tertiles)
# ------------------------------------------------------------
means <- means %>%
  mutate(
    subgroup = case_when(
      subgroup %in% c("HTN","DM","ASCVD_any","other_ascvd","PAD","Af") ~ subgroup,
      subgroup %in% c("HR_tertiles","HR_clinical","HR") ~ "HR",
      TRUE ~ subgroup
    )
  )

cons <- cons %>%
  mutate(
    subgroup = case_when(
      subgroup %in% c("HTN","DM","ASCVD_any","other_ascvd","PAD","Af") ~ subgroup,
      subgroup %in% c("HR_tertiles","HR_clinical","HR") ~ "HR",
      TRUE ~ subgroup
    )
  )

keep_sg <- c("HTN","DM","other_ascvd","PAD","Af","HR")
means <- means %>% filter(subgroup %in% keep_sg)
cons  <- cons  %>% filter(subgroup %in% keep_sg)

# ------------------------------------------------------------
# 2) Keep original sign (negative => K better) to match your axis labels request
#   If you wanted positive to favor Ksens, you'd set flip_sign=TRUE and swap labels.
# ------------------------------------------------------------
flip_sign <- TRUE
if (flip_sign) {
  means <- means %>%
    mutate(estimate = -estimate,
           lower0 = lower, upper0 = upper,
           lower = -upper0, upper = -lower0) %>%
    select(-lower0, -upper0)

  cons <- cons %>%
    mutate(estimate = -estimate,
           lower0 = lower, upper0 = upper,
           lower = -upper0, upper = -lower0) %>%
    select(-lower0, -upper0)
}

# ------------------------------------------------------------
# 3) p for interaction (one per subgroup; show only on last row in block)
# ------------------------------------------------------------
p_int <- cons %>%
  group_by(subgroup) %>%
  summarise(p_int = min(p, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 4) Labels
# ------------------------------------------------------------
subgroup_print <- function(sg) {
  sg <- as.character(sg)
  dplyr::if_else(sg %in% c("ASCVD_any", "other_ascvd"), "ASCVD", sg)
}

level_print <- function(sg, level) {
  level <- as.character(level)
  if (sg %in% c("HTN","DM","ASCVD_any","other_ascvd","PAD","Af")) {
    if (level %in% c("Yes","Y","1",1,TRUE)) return(paste0("with ", subgroup_print(sg)))
    if (level %in% c("No","N","0",0,FALSE))  return(paste0("without ", subgroup_print(sg)))
  }
  if (sg == "HR") {
    # enforce requested labels
    if (level %in% c("<60","< 60")) return("HR < 60")
    if (level %in% c("60-100","60 - 100","60–100")) return("HR 60-100")
    if (level %in% c(">100","> 100")) return("HR > 100")
    return(paste0("HR ", level))
  }
  paste0(subgroup_print(sg), ": ", level)
}

df <- means %>%
  mutate(
    subgroup2  = subgroup_print(subgroup),
    level_disp = mapply(level_print, subgroup, level, SIMPLIFY = TRUE),
    mean_txt   = sprintf("%.2f", estimate),
    low_txt    = sprintf("%.2f", lower),
    high_txt   = sprintf("%.2f", upper)
  ) %>%
  left_join(p_int, by = c("subgroup" = "subgroup"))

# Add a parent row for HR so children can be shown indented beneath it
hr_parent <- df %>%
  filter(subgroup2 == "HR") %>%
  slice(1) %>%
  mutate(
    level = NA,
    level_disp = "HR",
    mean_txt = "",
    low_txt = "",
    high_txt = "",
    estimate = NA_real_,
    lower = NA_real_,
    upper = NA_real_
  )

df <- bind_rows(df, hr_parent)

# ------------------------------------------------------------
# 5) Fixed ordering: HTN, DM, ASCVD, Af, HR
#    Within subgroup: with -> without; HR: <60 -> 60-100 -> >100
# ------------------------------------------------------------
df <- df %>%
  mutate(
    subgroup2 = factor(subgroup2, levels = c("HTN","DM","ASCVD","PAD","Af","HR")),
    level_rank = case_when(
      subgroup2 == "HR" & level_disp == "HR" ~ 0L,
      str_detect(level_disp, "^with") ~ 1L,
      str_detect(level_disp, "^without") ~ 2L,
      level_disp == "HR < 60" ~ 1L,
      level_disp == "HR 60-100" ~ 2L,
      level_disp == "HR > 100" ~ 3L,
      TRUE ~ 99L
    )
  ) %>%
  arrange(subgroup2, level_rank) %>%
  mutate(row_id = row_number(),
         y = rev(row_id))

# Show p only on last row of each subgroup block
df <- df %>%
  group_by(subgroup2) %>%
  mutate(show_p = row_number() == n()) %>%
  ungroup() %>%
  mutate(
    level_disp_plot = case_when(
      as.character(subgroup2) == "HR" & level_disp == "HR" ~ "HR",
      as.character(subgroup2) == "HR" ~ paste0("   ", str_remove(level_disp, "^HR\\s*")),
      TRUE ~ level_disp
    ),
    p_txt = ifelse(show_p,
                   ifelse(is.na(p_int), "", format.pval(p_int, digits = 3, eps = 1e-3)),
                   "")
  )

# ------------------------------------------------------------
# 6) Header row
# ------------------------------------------------------------
header <- tibble(
  y = max(df$y) + 1.25,
  mean_txt = "mean\ndifference",
  low_txt  = "95% CI Low",
  high_txt = "95% CI high",
  p_txt    = "p for interaction"
)

# ------------------------------------------------------------
# 7) Forest fixed range -5 ~ +5
# ------------------------------------------------------------
x_lim <- c(-5, 5)

# ------------------------------------------------------------
# 8) Left table panel
# ------------------------------------------------------------
left_tbl <- ggplot() +
  geom_text(data = df, aes(x = 0,    y = y, label = level_disp_plot), hjust = 0, size = 6) +
  geom_text(data = df, aes(x = 1.05, y = y, label = mean_txt),   hjust = 1, size = 6) +
  geom_text(data = df, aes(x = 1.60, y = y, label = low_txt),    hjust = 1, size = 6) +
  geom_text(data = df, aes(x = 2.15, y = y, label = high_txt),   hjust = 1, size = 6) +
  geom_text(data = header, aes(x = 1.05, y = y, label = mean_txt), hjust = 1, fontface = "bold", size = 5) +
  geom_text(data = header, aes(x = 1.60, y = y, label = low_txt),  hjust = 1, fontface = "bold", size = 5) +
  geom_text(data = header, aes(x = 2.15, y = y, label = high_txt), hjust = 1, fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(-0.05, 2.25), ylim = c(0.5, max(df$y) + 2.0), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 0, 10, 10))

# ------------------------------------------------------------
# 9) Forest panel + favor labels
# ------------------------------------------------------------
df_forest <- df %>% filter(!is.na(estimate), !is.na(lower), !is.na(upper))

forest <- ggplot(df_forest, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.0, orientation = "y") +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  labs(title = "Forest plot", x = NULL, y = NULL) +
  # Favor labels under the title: right side = Favor Ksens-SBP, left side = Favor Oscillo-SBP
  annotate("text",
           x = -2.0, y = max(df$y) + 1.35,
           label = "Favor Oscillo-SBP", size = 6, hjust = 0.5) +
  annotate("text",
           x =  2.0, y = max(df$y) + 1.35,
           label = "Favor Ksens-SBP", size = 6, hjust = 0.5) +
  coord_cartesian(xlim = x_lim, ylim = c(0.5, max(df$y) + 2.0), clip = "off") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 24),
    plot.margin = margin(10, 0, 10, 0)
  )

# ------------------------------------------------------------
# 10) Right p-value panel
# ------------------------------------------------------------
right_p <- ggplot() +
  geom_text(data = df, aes(x = 1, y = y, label = p_txt), hjust = 1, size = 6) +
  geom_text(data = header, aes(x = 1, y = y, label = p_txt), hjust = 1, fontface = "bold", size = 6) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.5, max(df$y) + 2.0), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 0))

# ------------------------------------------------------------
# 11) Combine and save
# ------------------------------------------------------------
p <- left_tbl + forest + right_p + plot_layout(widths = c(1.35, 1.70, 0.60))

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("outputs/figures/forest_SYS_absdiff_subgroups.png", p, width = 16, height = 8, dpi = 300)

p
