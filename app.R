# app.R

# 1. Load libraries
library(shiny)
library(shinythemes)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(stats)
library(FSA)
library(rlang)
library(grid)
library(colourpicker)
library(stringr)
library(reticulate)
library(DT)
library(showtext)
library(patchwork)
library(ggrepel)
library(scales)

# 2. Source modules
source("R/mod_instructions.R")
source("R/mod_gwas.R")
source("R/mod_expression.R")
source("R/mod_phenotype.R")
source("R/mod_locuszoom.R")

# 3. Setup fonts
font_add("Times New Roman", "./fonts/times.ttf")
showtext_auto()

# 4. Define UI
ui <- navbarPage(
  theme = shinytheme("flatly"),
  title = "星海之约ZYの❥Shiny",
  
  # Custom CSS
  tags$head(
    tags$style(HTML("
      body { font-family: 'Times New Roman', sans-serif; font-size: 14px; }
      .navbar-default .navbar-brand, .navbar-default .navbar-nav > li > a { font-weight: 600; }
      .nav-pills > li > a { font-weight: 500; }
      .nav-pills > li.active > a { background-color: #337ab7; color: #fff; }
      .well { background-color: #f9f9f9; border: 1px solid #eee; box-shadow: 0 2px 3px rgba(0,0,0,0.05); }
      table.table-hover > tbody > tr:hover { background-color: #f6f6f6; }
      .btn.btn-primary { margin-top: 5px; font-weight: 500; }
      .btn.btn-primary.btn-block { width: 100%; }
      .shiny-input-container { margin-bottom: 15px; }
    "))
  ),
  
  # Module UIs
  instructions_ui("instructions"),
  gwas_ui("gwas"),
  expression_ui("expression"),
  phenotype_ui("phenotype"),
  locuszoom_ui("locuszoom")
)

# 5. Define Server
server <- function(input, output, session) {
  # Module Servers
  instructions_server("instructions")
  gwas_server("gwas")
  expression_server("expression")
  phenotype_server("phenotype")
  locuszoom_server("locuszoom")
}

# 6. Run Application
shinyApp(ui, server)
