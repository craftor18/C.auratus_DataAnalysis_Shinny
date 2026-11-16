# ---- p-value helpers --------------------------------------------------------

fmt_p <- function(p, digits = 3) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 1e-300) return("< 1e-300")
    if (x < 1e-3) return(formatC(x, format = "e", digits = 2))
    formatC(x, format = "f", digits = digits)
  }, character(1))
}

p_stars <- function(p) {
  out <- cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"))
  out <- as.character(out)
  out[is.na(out)] <- "NA"
  out
}

P_ADJUST_CHOICES <- c(
  "Benjamini–Hochberg (FDR)" = "BH",
  "Holm" = "holm",
  "Bonferroni" = "bonferroni",
  "不校正" = "none"
)

# ---- genotype coding --------------------------------------------------------

MISSING_GENOTYPES <- c("", "NA", ".", "./.", ".|.", "-", "--", "N/N", "NN")

# Phased (0|1 vs 1|0) and unphased calls of the same alleles are the same genotype
# for association testing, so alleles are sorted after splitting.
normalize_genotype <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% MISSING_GENOTYPES] <- NA
  parts <- strsplit(x, "[/|]")
  vapply(parts, function(a) {
    if (length(a) == 1 && !is.na(a) && nchar(a) == 2) a <- strsplit(a, "")[[1]]
    if (length(a) != 2 || anyNA(a) || any(a %in% c(".", "N", ""))) return(NA_character_)
    paste(sort(a), collapse = "/")
  }, character(1))
}

genotype_info <- function(g) {
  g <- normalize_genotype(g)
  obs <- g[!is.na(g)]
  if (!length(obs)) stop("没有可识别的基因型。请使用 0/1、0|1、A/G 或 AG 等格式。")
  alleles_all <- unlist(strsplit(obs, "/"))
  tab <- sort(table(alleles_all), decreasing = TRUE)
  alleles <- names(tab)
  if (length(alleles) > 2) {
    stop(sprintf("检测到 %d 个等位基因（%s）。本工具只支持双等位 SNP。", length(alleles), paste(alleles, collapse = ", ")))
  }
  ref <- if (all(grepl("^[0-9]+$", alleles)) && "0" %in% alleles) "0" else alleles[1]
  alt <- setdiff(alleles, ref)
  alt <- if (length(alt)) alt else NA_character_

  hom_ref <- paste(ref, ref, sep = "/")
  het <- if (is.na(alt)) NA_character_ else paste(sort(c(ref, alt)), collapse = "/")
  hom_alt <- if (is.na(alt)) NA_character_ else paste(alt, alt, sep = "/")
  levels <- c(hom_ref, het, hom_alt)
  levels <- levels[!is.na(levels)]

  dosage <- vapply(strsplit(g, "/"), function(p) {
    if (anyNA(p) || is.na(alt)) return(if (anyNA(p)) NA_real_ else 0)
    sum(p == alt)
  }, numeric(1))

  n_hom_ref <- sum(g == hom_ref, na.rm = TRUE)
  n_het <- if (is.na(het)) 0L else sum(g == het, na.rm = TRUE)
  n_hom_alt <- if (is.na(hom_alt)) 0L else sum(g == hom_alt, na.rm = TRUE)
  n <- n_hom_ref + n_het + n_hom_alt
  alt_freq <- (n_het + 2 * n_hom_alt) / (2 * n)

  list(
    genotype = factor(g, levels = levels),
    dosage = dosage,
    ref = ref,
    alt = alt,
    counts = c(hom_ref = n_hom_ref, het = n_het, hom_alt = n_hom_alt),
    n_called = n,
    n_missing = sum(is.na(g)),
    alt_freq = alt_freq,
    maf = min(alt_freq, 1 - alt_freq),
    minor_allele = if (alt_freq <= 0.5) alt else ref,
    hwe_p = hwe_exact(n_het, n_hom_ref, n_hom_alt)
  )
}

# Exact HWE test (Wigginton, Cutler & Abecasis 2005, AJHG 76:887).
hwe_exact <- function(obs_hets, obs_hom1, obs_hom2) {
  n <- obs_hets + obs_hom1 + obs_hom2
  if (n == 0) return(NA_real_)
  obs_homc <- max(obs_hom1, obs_hom2)
  obs_homr <- min(obs_hom1, obs_hom2)
  rare <- 2 * obs_homr + obs_hets
  if (rare == 0) return(1)
  probs <- numeric(rare + 1)
  mid <- floor(rare * (2 * n - rare) / (2 * n))
  if ((rare %% 2) != (mid %% 2)) mid <- mid + 1
  probs[mid + 1] <- 1

  hets <- mid
  homr <- (rare - mid) / 2
  homc <- n - hets - homr
  while (hets >= 2) {
    probs[hets - 1] <- probs[hets + 1] * hets * (hets - 1) / (4 * (homr + 1) * (homc + 1))
    hets <- hets - 2
    homr <- homr + 1
    homc <- homc + 1
  }
  hets <- mid
  homr <- (rare - mid) / 2
  homc <- n - hets - homr
  while (hets <= rare - 2) {
    probs[hets + 3] <- probs[hets + 1] * 4 * homr * homc / ((hets + 2) * (hets + 1))
    hets <- hets + 2
    homr <- homr - 1
    homc <- homc - 1
  }
  probs <- probs / sum(probs)
  min(1, sum(probs[probs <= probs[obs_hets + 1] * (1 + 1e-7)]))
}

# ---- association models -----------------------------------------------------

# Additive model: y ~ allele dosage (+ covariates). Partial R2 is the variance
# explained by the SNP beyond the covariates.
assoc_additive <- function(y, dosage, covariates = NULL) {
  d <- data.frame(.y = y, .dosage = dosage)
  if (!is.null(covariates) && ncol(covariates)) d <- cbind(d, covariates)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 3 || length(unique(d$.dosage)) < 2) return(NULL)
  cov_names <- setdiff(names(d), c(".y", ".dosage"))
  for (cv in cov_names) {
    if (is.character(d[[cv]])) d[[cv]] <- factor(d[[cv]])
    if (is.factor(d[[cv]]) && nlevels(droplevels(d[[cv]])) < 2) d[[cv]] <- NULL
  }
  cov_names <- setdiff(names(d), c(".y", ".dosage"))
  rhs_reduced <- if (length(cov_names)) paste(sprintf("`%s`", cov_names), collapse = " + ") else "1"
  fit_full <- stats::lm(stats::as.formula(paste(".y ~ .dosage +", rhs_reduced)), data = d)
  fit_red <- stats::lm(stats::as.formula(paste(".y ~", rhs_reduced)), data = d)
  co <- summary(fit_full)$coefficients
  if (!".dosage" %in% rownames(co)) return(NULL)
  rss_full <- sum(stats::residuals(fit_full)^2)
  rss_red <- sum(stats::residuals(fit_red)^2)
  data.frame(
    model = if (length(cov_names)) paste0("加性模型（协变量：", paste(cov_names, collapse = ", "), "）") else "加性模型",
    n = nrow(d),
    beta = co[".dosage", 1],
    se = co[".dosage", 2],
    statistic = co[".dosage", 3],
    p = co[".dosage", 4],
    r2 = if (rss_red > 0) (rss_red - rss_full) / rss_red else NA_real_,
    check.names = FALSE
  )
}

omnibus_tests <- function(y, g) {
  d <- data.frame(y = y, g = droplevels(as.factor(g)))
  d <- d[stats::complete.cases(d), ]
  if (nlevels(d$g) < 2) return(NULL)
  n_per <- table(d$g)
  rows <- list()
  if (all(n_per >= 2)) {
    a <- summary(stats::aov(y ~ g, data = d))[[1]]
    rows$anova <- data.frame(test = "单因素 ANOVA", statistic = a[1, "F value"],
                             df = sprintf("%d, %d", a[1, "Df"], a[2, "Df"]), p = a[1, "Pr(>F)"])
    w <- stats::oneway.test(y ~ g, data = d, var.equal = FALSE)
    rows$welch <- data.frame(test = "Welch ANOVA（方差不齐）", statistic = unname(w$statistic),
                             df = sprintf("%.0f, %.1f", w$parameter[1], w$parameter[2]), p = w$p.value)
  }
  k <- stats::kruskal.test(y ~ g, data = d)
  rows$kw <- data.frame(test = "Kruskal–Wallis", statistic = unname(k$statistic),
                        df = as.character(k$parameter), p = k$p.value)
  do.call(rbind, unname(rows))
}

PAIRWISE_CHOICES <- c(
  "Tukey HSD（参数，自带校正）" = "tukey",
  "Welch t 检验" = "welch",
  "Dunn 检验（非参数）" = "dunn",
  "Wilcoxon 秩和检验" = "wilcox"
)

pairwise_tests <- function(y, g, method = "tukey", p_adjust = "BH") {
  d <- data.frame(y = y, g = droplevels(as.factor(g)))
  d <- d[stats::complete.cases(d), ]
  lv <- levels(d$g)
  if (length(lv) < 2) return(NULL)
  pairs <- utils::combn(lv, 2, simplify = FALSE)
  n_of <- table(d$g)
  mean_of <- tapply(d$y, d$g, mean)
  med_of <- tapply(d$y, d$g, stats::median)

  res <- data.frame(
    group1 = vapply(pairs, `[`, "", 1),
    group2 = vapply(pairs, `[`, "", 2),
    stringsAsFactors = FALSE
  )
  res$n1 <- as.integer(n_of[res$group1])
  res$n2 <- as.integer(n_of[res$group2])
  res$diff <- if (method %in% c("dunn", "wilcox")) {
    unname(med_of[res$group2] - med_of[res$group1])
  } else {
    unname(mean_of[res$group2] - mean_of[res$group1])
  }

  if (method == "tukey") {
    ok <- all(n_of >= 2)
    if (ok) {
      tk <- stats::TukeyHSD(stats::aov(y ~ g, data = d))$g
      key <- paste(res$group2, res$group1, sep = "-")
      key_rev <- paste(res$group1, res$group2, sep = "-")
      idx <- ifelse(key %in% rownames(tk), match(key, rownames(tk)), match(key_rev, rownames(tk)))
      res$statistic <- NA_real_
      res$p <- tk[idx, "p adj"]
      res$p.adj <- res$p
    } else {
      res$statistic <- NA_real_
      res$p <- NA_real_
      res$p.adj <- NA_real_
    }
    res$method <- "Tukey HSD"
  } else if (method == "dunn") {
    dn <- dunn_test(d$y, d$g)
    res$statistic <- dn$z[match(paste(res$group1, res$group2), paste(dn$group1, dn$group2))]
    res$p <- dn$p[match(paste(res$group1, res$group2), paste(dn$group1, dn$group2))]
    res$p.adj <- stats::p.adjust(res$p, method = p_adjust)
    res$method <- "Dunn"
  } else {
    stat_p <- t(vapply(seq_len(nrow(res)), function(i) {
      a <- d$y[d$g == res$group1[i]]
      b <- d$y[d$g == res$group2[i]]
      if (length(a) < 2 || length(b) < 2) return(c(NA_real_, NA_real_))
      tt <- if (method == "welch") {
        stats::t.test(b, a, var.equal = FALSE)
      } else {
        suppressWarnings(stats::wilcox.test(b, a, exact = FALSE))
      }
      c(unname(tt$statistic), tt$p.value)
    }, numeric(2)))
    res$statistic <- stat_p[, 1]
    res$p <- stat_p[, 2]
    res$p.adj <- stats::p.adjust(res$p, method = p_adjust)
    res$method <- if (method == "welch") "Welch t" else "Wilcoxon"
  }
  res$signif <- p_stars(res$p.adj)
  res
}

# Dunn (1964) z test on Kruskal–Wallis ranks with tie correction.
dunn_test <- function(y, g) {
  g <- droplevels(as.factor(g))
  N <- length(y)
  r <- rank(y)
  ties <- table(y)
  tie_corr <- sum(ties^3 - ties) / (12 * (N - 1))
  lv <- levels(g)
  mr <- tapply(r, g, mean)
  ni <- table(g)
  pairs <- utils::combn(lv, 2, simplify = FALSE)
  out <- lapply(pairs, function(p) {
    se <- sqrt((N * (N + 1) / 12 - tie_corr) * (1 / ni[[p[1]]] + 1 / ni[[p[2]]]))
    z <- (mr[[p[2]]] - mr[[p[1]]]) / se
    data.frame(group1 = p[1], group2 = p[2], z = z, p = 2 * stats::pnorm(-abs(z)))
  })
  do.call(rbind, out)
}

bracket_positions <- function(pw, y, step_frac = 0.08, only_signif = FALSE) {
  if (is.null(pw) || !nrow(pw)) return(pw)
  if (only_signif) pw <- pw[!is.na(pw$p.adj) & pw$p.adj < 0.05, , drop = FALSE]
  if (!nrow(pw)) return(pw)
  rng <- range(y, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- abs(rng[2]) + 1
  pw$y.position <- rng[2] + span * step_frac * seq_len(nrow(pw))
  pw
}

group_summary <- function(y, g) {
  d <- data.frame(y = y, g = as.factor(g))
  d <- d[!is.na(d$y) & !is.na(d$g), ]
  d |>
    group_by(组别 = g) |>
    summarise(
      n = n(),
      均值 = mean(y),
      标准差 = stats::sd(y),
      标准误 = stats::sd(y) / sqrt(n()),
      中位数 = stats::median(y),
      Q1 = stats::quantile(y, 0.25),
      Q3 = stats::quantile(y, 0.75),
      `CV%` = 100 * stats::sd(y) / mean(y),
      .groups = "drop"
    )
}

handle_outliers <- function(df, value_col, group_col, method = "none") {
  df$.outlier <- FALSE
  if (method == "none") return(df)
  split_idx <- split(seq_len(nrow(df)), df[[group_col]])
  for (idx in split_idx) {
    v <- df[[value_col]][idx]
    if (method == "iqr") {
      q <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE)
      iqr <- q[2] - q[1]
      df$.outlier[idx] <- !is.na(v) & (v < q[1] - 1.5 * iqr | v > q[2] + 1.5 * iqr)
    } else if (method == "winsor") {
      q <- stats::quantile(v, c(0.05, 0.95), na.rm = TRUE)
      df[[value_col]][idx] <- pmin(pmax(v, q[1]), q[2])
    }
  }
  if (method == "iqr") df <- df[!df$.outlier, , drop = FALSE]
  df
}

OUTLIER_CHOICES <- c(
  "保留全部数据（推荐）" = "none",
  "组内 1.5×IQR 剔除" = "iqr",
  "组内 5–95% 缩尾 (winsorize)" = "winsor"
)

# ---- expression -------------------------------------------------------------

sample_group <- function(samples, rule = "prefix") {
  if (rule == "prefix") sub("[_.-][^_.-]*$", "", samples) else gsub("[0-9]+$", "", samples)
}

expression_long <- function(mat_df, genes, groups, scale = "raw", group_rule = "prefix") {
  sample_cols <- setdiff(names(mat_df), "gene_symbol")
  grp <- sample_group(sample_cols, group_rule)
  keep <- grp %in% groups
  sub <- mat_df[mat_df$gene_symbol %in% genes, c("gene_symbol", sample_cols[keep]), drop = FALSE]
  long <- tidyr::pivot_longer(sub, -gene_symbol, names_to = "sample", values_to = "value")
  long$value <- suppressWarnings(as.numeric(long$value))
  long$group <- factor(sample_group(long$sample, group_rule), levels = groups)
  long$gene_symbol <- factor(long$gene_symbol, levels = genes[genes %in% sub$gene_symbol])
  long$expr <- if (scale == "raw") log2(long$value + 1) else long$value
  long
}

expression_tests <- function(long, method = "welch", p_adjust = "BH") {
  groups <- levels(long$group)
  out <- lapply(split(long, long$gene_symbol), function(d) {
    if (!nrow(d)) return(NULL)
    means <- tapply(d$expr, d$group, mean, na.rm = TRUE)
    row <- data.frame(gene = as.character(d$gene_symbol[1]))
    for (gname in groups) row[[paste0("mean_", gname)]] <- unname(means[gname])
    if (length(groups) == 2) {
      a <- d$expr[d$group == groups[1]]
      b <- d$expr[d$group == groups[2]]
      row$log2FC <- mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE)
      row$p <- tryCatch({
        if (method == "welch") stats::t.test(b, a)$p.value else suppressWarnings(stats::wilcox.test(b, a, exact = FALSE))$p.value
      }, error = function(e) NA_real_)
      row$test <- if (method == "welch") "Welch t" else "Wilcoxon"
    } else {
      row$log2FC <- NA_real_
      row$p <- tryCatch({
        if (method == "welch") stats::oneway.test(expr ~ group, data = d)$p.value else stats::kruskal.test(expr ~ group, data = d)$p.value
      }, error = function(e) NA_real_)
      row$test <- if (method == "welch") "Welch ANOVA" else "Kruskal–Wallis"
    }
    row
  })
  res <- do.call(rbind, out)
  if (is.null(res)) return(NULL)
  res$padj <- stats::p.adjust(res$p, method = p_adjust)
  res$signif <- p_stars(res$padj)
  res
}
