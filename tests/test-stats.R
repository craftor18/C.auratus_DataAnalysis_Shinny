# Run from the app directory: Rscript tests/test-stats.R  (or shiny::runTests())
app_dir <- if (file.exists("app.R")) "." else ".."
old <- setwd(app_dir)
on.exit(setwd(old), add = TRUE)
for (f in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) source(f)

check <- function(desc, ok) {
  if (!isTRUE(ok)) stop("FAILED: ", desc, call. = FALSE)
  cat("ok -", desc, "\n")
}

# Genotype normalisation: phase must not create separate genotype classes.
g <- normalize_genotype(c("0|1", "1|0", "0/0", "1/1", "./.", "AG", "G/A", NA, ""))
check("phased hets collapse to 0/1", identical(g[1:4], c("0/1", "0/1", "0/0", "1/1")))
check("missing codes become NA", all(is.na(g[c(5, 8, 9)])))
check("letter genotypes are sorted", identical(g[6:7], c("A/G", "A/G")))

gi <- genotype_info(c("0|0", "0|1", "1|0", "1|1", "0|0", NA))
check("REF is 0 for numeric coding", gi$ref == "0" && gi$alt == "1")
check("dosage counts ALT copies", identical(gi$dosage, c(0, 1, 1, 2, 0, NA)))
check("genotype counts", identical(unname(gi$counts), c(2L, 2L, 1L)))
check("ALT frequency", abs(gi$alt_freq - 0.4) < 1e-12)

gl <- genotype_info(c("A/A", "A/A", "A/G", "G/G", "A/A"))
check("major allele is REF for letters", gl$ref == "A" && gl$alt == "G")
check("multi-allelic sites are rejected", inherits(try(genotype_info(c("A/C", "G/T")), silent = TRUE), "try-error"))

# HWE exact test: perfect equilibrium vs. total heterozygote deficit.
check("HWE p ~ 1 at equilibrium", hwe_exact(50, 25, 25) > 0.8)
check("HWE p tiny with no hets", hwe_exact(0, 50, 50) < 1e-20)
check("HWE p in [0, 1]", all(vapply(list(c(10, 80, 10), c(57, 14, 29)), function(x) {
  p <- hwe_exact(x[1], x[2], x[3]); p >= 0 && p <= 1
}, logical(1))))

# Additive model recovers a known effect.
set.seed(1)
dos <- sample(0:2, 300, replace = TRUE, prob = c(.49, .42, .09))
y <- 10 + 1.5 * dos + rnorm(300)
a <- assoc_additive(y, dos)
check("additive beta recovered", abs(a$beta - 1.5) < 0.25 && a$p < 1e-20)
sex <- sample(c("F", "M"), 300, replace = TRUE)
a2 <- assoc_additive(y + 3 * (sex == "M"), dos, data.frame(sex = sex))
check("covariate-adjusted beta recovered", abs(a2$beta - 1.5) < 0.25)

# Pairwise tests and adjustment.
grp <- factor(rep(c("a", "b", "c"), each = 20))
val <- c(rnorm(20, 0), rnorm(20, 0), rnorm(20, 3))
for (m in c("tukey", "welch", "dunn", "wilcox")) {
  pw <- pairwise_tests(val, grp, m, "BH")
  check(paste(m, "returns 3 pairs"), nrow(pw) == 3)
  check(paste(m, "detects shifted group"), all(pw$p.adj[pw$group2 == "c"] < 0.01))
  check(paste(m, "adjusted p >= raw p"), all(pw$p.adj >= pw$p - 1e-12))
}
om <- omnibus_tests(val, grp)
check("omnibus has ANOVA, Welch, KW", nrow(om) == 3 && all(om$p < 1e-6))

# Dunn's z equals manual computation without ties.
dn <- dunn_test(c(1, 2, 3, 4, 5, 6), factor(c("x", "x", "x", "y", "y", "y")))
z_manual <- (5 - 2) / sqrt((6 * 7 / 12) * (1 / 3 + 1 / 3))
check("dunn z matches formula", abs(dn$z - z_manual) < 1e-10)

# Expression helpers.
m <- data.frame(gene_symbol = c("g1", "g2"), A_1 = c(1, 10), A_2 = c(3, 12), B_1 = c(15, 10), B_2 = c(17, 11))
check("group prefix parsing", identical(sample_group(c("HF1_3", "WJ_10", "A_1")), c("HF1", "WJ", "A")))
l <- expression_long(m, c("g1", "g2"), c("A", "B"), scale = "raw")
check("log2(x+1) transform", abs(l$expr[1] - log2(2)) < 1e-12)
et <- expression_tests(l, "welch", "BH")
check("log2FC direction is B - A", et$log2FC[et$gene == "g1"] > 2)

# Outliers and examples.
od <- handle_outliers(data.frame(y = c(1, 2, 2, 3, 100), g = "a"), "y", "g", "iqr")
check("IQR removes extreme value", nrow(od) == 4)
check("p_stars thresholds", identical(p_stars(c(1e-5, 0.0005, 0.005, 0.03, 0.2)), c("****", "***", "**", "*", "ns")))
for (k in names(EXAMPLE_FILES)) check(paste("example loads:", k), nrow(load_example(k)) > 0)

cat("\nAll tests passed.\n")
