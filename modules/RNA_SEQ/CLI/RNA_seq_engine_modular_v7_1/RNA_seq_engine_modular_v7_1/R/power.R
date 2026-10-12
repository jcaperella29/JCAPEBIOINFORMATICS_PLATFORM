# Extracted from RNA_SEQ_CLI_V7_1_STAT_CLEARANCE_CANDIDATE_2026_08_10
# Modularization only: statistical behavior intentionally unchanged.

run_power_analysis <- function(pheno, phenotype_column, effect_size) {
  groups <- as.factor(pheno[[phenotype_column]])
  group_sizes <- table(groups)
  k <- length(group_sizes)
  sig <- 0.05
  note <- "Generic effect-size power reference only; not RNA-seq-specific, gene-wise, or FDR-aware power."
  if (k == 2) {
    n1 <- as.numeric(group_sizes[1])
    n2 <- as.numeric(group_sizes[2])
    power_res <- pwr::pwr.t2n.test(n1 = n1, n2 = n2, d = effect_size, sig.level = sig)
    data.frame(Analysis = "Generic effect-size power reference", Test = "t-test", `Effect Size (d)` = effect_size, Groups = paste(names(group_sizes), collapse = ", "), n1 = n1, n2 = n2, `Significance Level` = sig, `Estimated Power` = round(power_res$power, 3), Note = note, check.names = FALSE)
  } else if (k > 2) {
    total_n <- sum(group_sizes)
    power_res <- pwr::pwr.anova.test(k = k, n = total_n / k, f = effect_size, sig.level = sig)
    data.frame(Analysis = "Generic effect-size power reference", Test = "ANOVA", `Effect Size (f)` = effect_size, Groups = paste(names(group_sizes), collapse = ", "), `Samples per Group` = round(total_n / k), `Significance Level` = sig, `Estimated Power` = round(power_res$power, 3), Note = note, check.names = FALSE)
  } else {
    data.frame(Analysis = "Generic effect-size power reference", Message = "Not enough groups for power analysis", Note = note)
  }
}

make_power_curve_plot <- function(effect_size, test_type, n_min, n_max) {
  n_seq <- seq(n_min, n_max)
  sig <- 0.05
  power_vals <- sapply(n_seq, function(n) {
    if (test_type == "ttest") pwr::pwr.t.test(n = n, d = effect_size, sig.level = sig, type = "two.sample")$power
    else pwr::pwr.anova.test(k = 3, n = n, f = effect_size, sig.level = sig)$power
  })
  df <- data.frame(SampleSize = n_seq, Power = power_vals)
  plot_ly(df, x = ~SampleSize, y = ~Power, type = "scatter", mode = "lines+markers", line = list(color = "#00cc99", width = 3)) %>%
    layout(title = "Generic Effect-size Power Reference (Not RNA-seq-specific)", xaxis = list(title = "Sample Size per Group"), yaxis = list(title = "Power", range = c(0, 1)),
           shapes = list(list(type = "line", x0 = min(n_seq), x1 = max(n_seq), y0 = 0.8, y1 = 0.8, line = list(dash = "dash", color = "red"))))
}

