pkgs <- c(
  "shiny", "bslib", "ggplot2", "dplyr", "tidyr", "readr", "DT", "patchwork",
  "ggrepel", "ggpubr", "scales", "showtext", "markdown"
)
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need, repos = getOption("repos", "https://cloud.r-project.org"))
invisible(lapply(pkgs, requireNamespace))
cat("All packages available:", paste(pkgs, collapse = ", "), "\n")
