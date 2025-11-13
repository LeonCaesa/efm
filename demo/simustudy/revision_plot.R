# =========================================================
# revision_plot_all_in_one.R — one-call plot generator
# Generates:
#   1) EFMOptiComparep5/10/512.png
#   2) EFMOptiComparep512LargeV.png
#   3) EFMOptiComparep512Largeq.png
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ---------- shared style ----------
mk_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      legend.position = "bottom",
      text         = element_text(face = "plain", family = "Arial"),
      strip.text.x = element_text(size = base_size, face = "plain", family = "Arial"),
      axis.title   = element_text(size = base_size + 1, face = "plain", family = "Arial"),
      axis.text    = element_text(size = base_size, face = "plain", family = "Arial"),
      legend.title = element_text(size = base_size, face = "plain", family = "Arial"),
      legend.text  = element_text(size = base_size, face = "plain", family = "Arial")
    )
}
algo_scale_all <- function() {
  scale_colour_manual(
    values = c(em="#F8766D", lapl="#7CAE00", ps="#00BFC4", sml="#C77CFF"),
    breaks = c("em","lapl","ps","sml"), name="algo"
  )
}
size_shapes <- function() {
  scale_shape_manual(
    values = c(`50`=17, `300`=16, `500`=15),
    breaks = c("50","300","500"), name="size"
  )
}

# =========================================================
# TYPE 1: EFMOptiComparep[5,10,512].png
# =========================================================
build_type1 <- function(dir_type1, d_values = c(5,10,512),
                        n = 512, q = 2,
                        families = c("poisson","binomial","negbinom"),
                        algos = c("ps","sml","lapl","em"),
                        sizes = c(50,300,500), T_epochs = 25) {

  out <- tibble(loss=numeric(), Model=character(), size=numeric(), algo=character(),
                time=numeric(), d=numeric(), q=numeric(), comp_time=numeric())

  for (d in d_values) for (fam in families) for (alg in algos) {
    s_idx <- if (alg %in% c("ps","sml")) seq_along(sizes) else 1
    for (i in s_idx) {
      f <- paste0(
        dir_type1,
        paste(alg, fam, paste0("s", sizes[i]),
              paste0("d", d), paste0("q", q), paste0("T", T_epochs),
              sep = "_"),
        ".RData"
      )
      try({
        load(f)  # -> efm_result
        it <- length(efm_result$like_list)
        efm_result$efm_time <- as.double(efm_result$efm_time) +
          as.double(efm_result$eval_time) -
          as.numeric(efm_result$eval_time, units = "secs")
        out <- bind_rows(out, tibble(
          loss      = as.numeric(efm_result$like_list),
          Model     = fam,
          size      = sizes[i],
          algo      = alg,
          time      = seq_len(it),
          d         = d,
          q         = q,
          comp_time = as.numeric(efm_result$efm_time)/it
        ))
      }, silent = TRUE)
    }
  }
  list(data = out, n = n)
}

# NOTE: binomial facet label uses (logit)
# REPLACE your existing plot_type1_for_d() with this
plot_type1_for_d <- function(dat, n, d_val, x_max = NULL, x_by = NULL, base_size = 11) {
  fam_labels <- c(
    binomial = "Binomial(logit)",
    negbinom = "Negbinom(log)",
    poisson  = "Poisson(log)"
  )

  df <- dat %>%
    filter(d == d_val,
           algo %in% c("ps","sml","lapl","em"),
           Model %in% c("poisson","binomial","negbinom"),
           size %in% c(50,300,500)) %>%
    mutate(
      Model = factor(Model, levels = c("binomial","negbinom","poisson")),
      algo  = factor(algo, levels = c("em","lapl","ps","sml")),
      size  = factor(size, levels = c(50,300,500), labels = c("50","300","500"))
    )

  p <- ggplot(df, aes(x = as.numeric(time), y = log(as.numeric(loss)/n))) +
    geom_point(aes(shape = size, colour = algo), alpha = 0.85, stroke = 0.4, size = 1.8) +
    algo_scale_all() + size_shapes() + mk_theme(base_size) +
    xlab("Adam Steps") + ylab("Avged Negative likelihood") +
    facet_wrap(~ Model, scales = "free", labeller = labeller(Model = fam_labels))

  # Axis handling: only enforce limits/ticks if provided; otherwise let ggplot decide
  if (!is.null(x_max) || !is.null(x_by)) {
    breaks <- if (!is.null(x_by) && !is.null(x_max)) {
      seq(0, x_max, by = x_by)
    } else if (!is.null(x_by)) {
      seq(0, max(df$time, na.rm = TRUE), by = x_by)
    } else {
      waiver()
    }
    limits <- if (!is.null(x_max)) c(0, x_max) else NULL

    p <- p + scale_x_continuous(limits = limits, breaks = breaks,
                                expand = expansion(mult = c(0, 0.02)))
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
  }

  p
}


# =========================================================
# TYPE 2: EFMOptiComparep512LargeV.png (large-p time traces)
# =========================================================
build_type2 <- function(dir_type2, n = 512, d = 512, T_epochs = 10) {
  fam_labels <- c("Poisson(log)", "Gamma(log)", "Negbiom(log)")  # index 1..3 in files
  out <- tibble(loss=numeric(), Model=character(), size_lbl=character(),
                algo=character(), time=numeric())

  for (fam_idx in 1:3) for (s in 1:4) {
    f <- paste0(dir_type2, paste(fam_idx, s, n, d, T_epochs, sep = "_"), ".RData")
    try({
      load(f)
      if (s == 4) {
        it  <- length(lapl_result$like_list)
        tms <- cumsum(seq_len(it) * lapl_result$lapl_time / it)
        tmp <- tibble(loss = as.numeric(lapl_result$like_list),
                      Model = fam_labels[fam_idx], size_lbl = "50",
                      algo  = "lapl", time = as.numeric(tms))
      } else {
        it  <- length(ps_result$like_list)
        tms <- cumsum(seq_len(it) * ps_result$ps_time / it)
        tmp <- tibble(loss = as.numeric(ps_result$like_list),
                      Model = fam_labels[fam_idx],
                      size_lbl = c("50","300","500")[s],
                      algo  = "ps", time = as.numeric(tms))
      }
      out <- bind_rows(out, tmp)
    }, silent = TRUE)
  }
  list(data = out, n = n)
}

plot_type2 <- function(dat, n, base_size = 11) {
  dat %>%
    filter(time <= 10000) %>%
    mutate(Model = factor(Model, levels = c("Gamma(log)","Negbiom(log)","Poisson(log)")),
           algo  = factor(algo, levels = c("lapl","ps")),
           size_lbl = factor(size_lbl, levels = c("50","300","500"))) %>%
    ggplot(aes(x = as.numeric(time), y = log(as.numeric(loss)/n))) +
    geom_point(aes(colour = algo, shape = size_lbl), alpha = 0.85, stroke = 0.4) +
    scale_colour_manual(values = c(lapl="#F8766D", ps="#00BFC4"), name = "algo") +
    scale_shape_manual(values = c(`50`=15, `300`=17, `500`=1), name = "size") +
    mk_theme(base_size) + xlab("Time (s)") + ylab("Avged Negative likelihood") +
    facet_wrap(~ Model, scales = "free")
}

# =========================================================
# TYPE 3: EFMOptiComparep512Largeq.png (binomial, d=512, varying q)
# =========================================================
build_type3 <- function(dir_type3,
                        n = 512, d = 512, q_list = c(6,8,12),
                        algos = c("ps","sml","lapl","em"),
                        T_epochs = 25) {

  out <- tibble(loss = numeric(), size = character(),
                algo = character(), time = numeric(), q = numeric())

  for (alg in algos) for (qq in q_list) {
    pat <- sprintf("^%s_binomial.*_d%s_q%s_T%s\\.RData$", alg, d, qq, T_epochs)
    files <- list.files(dir_type3, pattern = pat, full.names = TRUE)

    if (length(files) == 0) {
      message(sprintf("Skipping: no files for alg=%s, q=%s in %s", alg, qq, dir_type3))
      next
    }

    for (f in files) {
      size_match <- str_match(basename(f), "_s(\\d+)")
      size_lbl   <- ifelse(is.na(size_match[,2]), "50", size_match[,2])

      load(f)  # -> efm_result
      it <- length(efm_result$like_list)

      out <- bind_rows(out, tibble(
        loss = as.numeric(efm_result$like_list),
        size = size_lbl,
        algo = tolower(alg),
        time = seq_len(it),
        q    = qq
      ))
    }
  }
  list(data = out, n = n)
}

plot_type3 <- function(dat, n, base_size = 11) {
  dat %>%
    mutate(
      algo = factor(algo, levels = c("em","lapl","ps","sml")),
      size = factor(size, levels = c("50","300","500"), labels = c("50","300","500")),
      q_lab = factor(q, levels = sort(unique(q)), labels = paste0("Rank = ", sort(unique(q))))
    ) %>%
    ggplot(aes(x = as.numeric(time), y = log(as.numeric(loss)/n))) +
    geom_point(aes(shape = size, colour = algo), alpha = 0.85, stroke = 0.4, size = 1.8) +
    scale_colour_manual(
      values = c(em="#F8766D", lapl="#7CAE00", ps="#00BFC4", sml="#C77CFF"),
      breaks = c("em","lapl","ps","sml"), name="algo"
    ) +
    scale_shape_manual(values = c(`50`=17, `300`=16, `500`=15),
                       breaks = c("50","300","500"), name="size") +
    mk_theme(base_size) +
    xlab("Adam Steps") + ylab("Avged Negative likelihood") +
    facet_wrap(~ q_lab, scales = "free")
}

# =========================================================
# ONE-CALL WRAPPER
# =========================================================
render_revision_plots <- function(
    # ---- Default paths (edit if needed) ----
    dir_type1 = "results/OptiResult0108/",
    dir_type2 = "results/LargeP2/",
    dir_type3 = "results/OptiResult0118_2025/",
    out_dir   = "figures/",
    which     = c(1, 2, 3),            # any subset of {1,2,3}
    # ---- Type 1 options ----
    d_values_type1 = c(5, 10, 512),
    n_type1 = 512, q_type1 = 2, sizes_type1 = c(50,300,500), T1_epochs = 25,
    x_max_type1 = NULL, x_by_type1 = 10,    # set x_max_type1=80/100 to show full axis
    title_template_type1 = NULL,            # e.g., "EFMOptiCompare p=%s"
    # ---- Type 2 options ----
    n_type2 = 512, d_type2 = 512, T2_epochs = 10,
    # ---- Type 3 options ----
    n_type3 = 512, d_type3 = 512, q_list_type3 = c(6,8,12), T3_epochs = 25,
    # ---- IO & styling ----
    width = 8, height = 4, dpi = 300, overwrite = TRUE, return_plots = TRUE,
    base_size = 11
) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  safe_ggsave <- function(path, plot) {
    if (!overwrite && file.exists(path)) {
      message("Skip (exists): ", path)
      return(invisible(FALSE))
    }
    tryCatch({
      ggplot2::ggsave(path, plot, width = width, height = height, units = "in", dpi = dpi)
      message("Wrote: ", path)
      TRUE
    }, error = function(e) { warning("Failed to save ", path, ": ", conditionMessage(e)); FALSE })
  }

  plots <- list()

  # ----- TYPE 1 -----
  if (1 %in% which) {
    t1 <- build_type1(
      dir_type1, d_values = d_values_type1,
      n = n_type1, q = q_type1, sizes = sizes_type1, T_epochs = T1_epochs
    )
    type1_plots <- lapply(d_values_type1, function(dv) {
      p <- plot_type1_for_d(t1$data, t1$n, d_val = dv,
                            x_max = x_max_type1, x_by = x_by_type1,
                            base_size = base_size)
      if (!is.null(title_template_type1)) {
        p <- p + ggplot2::labs(title = sprintf(title_template_type1, dv))
      }
      fname <- file.path(out_dir, sprintf("EFMOptiComparep%s.png", dv))
      safe_ggsave(fname, p)
      p
    })
    names(type1_plots) <- paste0("p", d_values_type1)
    plots$type1 <- type1_plots
  }

  # ----- TYPE 2 -----
  if (2 %in% which) {
    t2 <- build_type2(dir_type2, n = n_type2, d = d_type2, T_epochs = T2_epochs)
    p2 <- plot_type2(t2$data, t2$n, base_size = base_size)
    safe_ggsave(file.path(out_dir, "EFMOptiComparep512LargeV.png"), p2)
    plots$type2 <- p2
  }

  # ----- TYPE 3 -----
  if (3 %in% which) {
    t3 <- build_type3(dir_type3, n = n_type3, d = d_type3,
                      q_list = q_list_type3, T_epochs = T3_epochs)
    p3 <- plot_type3(t3$data, t3$n, base_size = base_size)
    safe_ggsave(file.path(out_dir, "EFMOptiComparep512Largeq.png"), p3)
    plots$type3 <- p3
  }

  invisible(if (return_plots) plots else NULL)
}

# ---------------- One-line usage (edit as needed) ----------------
# Show x-axis to 100 with ticks every 20, no top title, base font ≥10pt.
render_revision_plots(
  x_max_type1 = NULL,
  x_by_type1 = 20,
  title_template_type1 = NULL,  # e.g., "EFMOptiCompare p=%s" to add a title
  base_size = 11                # >=10pt for print
)
