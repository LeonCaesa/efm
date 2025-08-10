# ---- Packages ----
need <- c("tidyverse", "ggh4x", "readr")
to_install <- setdiff(need, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggh4x)
  library(readr)
})

# If running in RStudio, set working directory to this script's folder (safe no-op otherwise)
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)), silent = TRUE)
}

# ---- Inputs ----
file_dir  <- "/projectnb/dmfgrp/efm/CovResult1209/"
file_name <- "total_466_20.RData"
family_list <- c("binomial", "negbinom(20)", "poisson", "quasipoisson")

# ---- Load & Combine (robust typing) ----
agg_df <- tibble()

for (fam in family_list) {
  rdata_path <- file.path(file_dir, fam, file_name)
  loaded <- load(rdata_path)               # expects object `error_matrix`
  if (!("error_matrix" %in% loaded || exists("error_matrix"))) {
    stop("`error_matrix` not found in: ", rdata_path)
  }
  df <- as_tibble(error_matrix, .name_repair = "minimal")
  names(df)[1:5] <- c("error", "error_type", "d", "esti_method", "repeat_idx")
  df$family <- fam

  df <- type_convert(df, col_types = cols(
    error       = col_double(),
    error_type  = col_character(),
    d           = col_integer(),
    esti_method = col_character(),
    repeat_idx  = col_integer(),
    family      = col_character()
  )) |>
    mutate(
      error       = suppressWarnings(as.numeric(error)),
      d           = suppressWarnings(as.integer(d)),
      repeat_idx  = suppressWarnings(as.integer(repeat_idx)),
      error_type  = as.character(error_type),
      esti_method = as.character(esti_method),
      family      = as.character(family)
    )

  agg_df <- bind_rows(agg_df, df)
  rm(error_matrix)
}

# ---- Clean labels ----
agg_df <- agg_df |>
  mutate(
    esti_method = if_else(esti_method == "fagqem", "efm(em)", esti_method),
    error_type  = recode(error_type, "l2entrophy" = "l2entropy"),
    d           = factor(d, levels = c(66, 166, 216, 266, 316, 366, 416, 466)),
    family      = factor(family, levels = c("binomial", "negbinom(20)", "poisson", "quasipoisson")),
    esti_method = factor(esti_method, levels = c("efm(em)", "naive"))
  )

plot_df <- agg_df |>
  filter(error_type %in% c("l2entropy", "l2frobenius", "l2normalized"),
         as.numeric(as.character(d)) > 16)

# Facet label for binomial
fam_labels <- c(
  "binomial"      = "Binomial (logit)",
  "negbinom(20)"  = "Negbin (\u03C6 = 20)",  # φ = 20
  "poisson"       = "Poisson",
  "quasipoisson"  = "Quasipoisson"
)

# ---- Colors (match ggplot2 default hues used in your original) ----
cols <- c("efm(em)" = "#F8766D",  # reddish
          "naive"   = "#00BFC4")  # teal

# Short k/M labeler (robust for faceting)
label_k <- function(accuracy = 1) {
  force(accuracy)
  function(x) {
    if (!length(x)) return(character())
    m <- max(abs(x), na.rm = TRUE)
    if (!is.finite(m)) return(rep("", length(x)))
    if (m >= 1e6) paste0(scales::number(x/1e6, accuracy = accuracy), "M")
    else if (m >= 1e3) paste0(scales::number(x/1e3, accuracy = accuracy), "k")
    else scales::number(x, accuracy = accuracy)
  }
}

# ---- Plot (capitalized titles/labels) ----
OUT_WIDTH_IN  <- 8
OUT_HEIGHT_IN <- 4.6
BASE_PT <- 11

p <- ggplot(plot_df, aes(x = d, y = as.numeric(error))) +
  geom_boxplot(
    aes(color = esti_method, fill = esti_method),
    width = 0.7,
    position = position_dodge2(width = 0.75, preserve = "single"),
    outlier.size = 0.7, alpha = 0.35, size = 0.3
  ) +
  stat_summary(
    aes(group = esti_method, shape = esti_method, color = esti_method),
    fun = median, geom = "point",
    position = position_dodge(width = 0.75),
    size = 1.6, show.legend = TRUE
  ) +
  ggh4x::facet_nested(
    'Error Type' + error_type ~ 'EFM Family' + family,   # <- capitalized strip headers
    scales = "free_y", independent = "y",
    labeller = labeller(family = fam_labels)
  ) +
  labs(
    x = "D",                               # <- capitalized axis title
    y = "Covariance Error",                # <- capitalized axis title
    fill  = "Estimator",                   # <- capitalized legend title
    color = "Estimator",
    shape = "Estimator"
    # subtitle/title can be added later and will already be capitalized by you
    # e.g., title = "EFM Family", subtitle = "Simulation Results"
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values  = cols) +
  scale_shape_manual(values = c("efm(em)" = 16, "naive" = 17)) +
  scale_x_discrete(breaks = c("66","166","266","366","466")) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(5),
    labels = label_k(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  theme_bw(base_size = BASE_PT) +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(face = "plain"),
    strip.text.y = element_text(face = "plain"),
    strip.text.y.right = element_text(margin = margin(l = 8, r = 4)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.25),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    panel.spacing.x = unit(0.12, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    plot.margin = margin(t = 5, r = 16, b = 8, l = 3)
  )

print(p)

ggsave("/projectnb/dmfgrp/efm/figures/cov_error.png", p,
       width = OUT_WIDTH_IN, height = OUT_HEIGHT_IN, units = "in", dpi = 300)
