# C. auratus GWAS downstream analysis workbench.
# Modules and helpers live in R/ and are sourced automatically by shiny::runApp().

ui <- page_navbar(
  id = "main_nav",
  title = tags$span(class = "brand", tags$span(class = "brand-mark", icon("dna")), "Auratus", tags$b("GWAS")),
  window_title = "Auratus GWAS · 鲫鱼关联分析工作台",
  theme = app_theme(),
  fillable = FALSE,
  header = tagList(
    tags$head(
      tags$link(rel = "stylesheet", href = "styles.css"),
      tags$script(src = "app.js"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1")
    ),
    useBusyIndicators(spinners = TRUE, pulse = TRUE)
  ),
  nav_panel("首页", value = "home", icon = icon("house"), home_ui("home")),
  nav_panel("群体表型", value = "pheno", icon = icon("fish"), phenotype_ui("pheno")),
  nav_panel("区域定位", value = "locus", icon = icon("magnifying-glass-chart"), locuszoom_ui("locus")),
  nav_panel("等位效应", value = "geno", icon = icon("dna"), genopheno_ui("geno")),
  nav_panel("候选表达", value = "expr", icon = icon("chart-simple"), expression_ui("expr")),
  nav_panel("指南", value = "guide", icon = icon("book-open"), guide_ui("guide")),
  nav_spacer(),
  nav_item(input_dark_mode(id = "dark_mode", mode = "light")),
  nav_item(tags$a(href = "https://github.com/craftor18/C.auratus_DataAnalysis_Shinny", target = "_blank",
                  class = "nav-link", title = "GitHub", icon("github")))
)

server <- function(input, output, session) {
  observeEvent(input$goto_tab, nav_select("main_nav", input$goto_tab))
  phenotype_server("pheno")
  locuszoom_server("locus")
  genopheno_server("geno")
  expression_server("expr")
}

shinyApp(ui, server)
