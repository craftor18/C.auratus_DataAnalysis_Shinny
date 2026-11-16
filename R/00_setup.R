suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(DT)
  library(patchwork)
  library(ggrepel)
  library(ggpubr)
  library(scales)
  library(showtext)
})

options(shiny.maxRequestSize = 100 * 1024^2)

`%||%` <- function(x, y) if (is.null(x)) y else x

if (file.exists("fonts/times.ttf")) {
  sysfonts::font_add("Times New Roman", "fonts/times.ttf")
}
showtext_auto()
SCREEN_DPI <- 96
showtext_opts(dpi = SCREEN_DPI)

EXAMPLE_FILES <- list(
  genotype_phenotype = "data/examples/genotype_phenotype.tsv",
  population_phenotype = "data/examples/population_phenotype.csv",
  locuszoom = "data/examples/locuszoom_lg34.csv",
  gene_annotation = "data/examples/gene_annotation.csv",
  expression = "data/examples/expression_tpm.csv"
)

FONT_CHOICES <- c("无衬线 (Sans)" = "sans", "衬线 (Serif)" = "serif", "Times New Roman" = "Times New Roman")
