# R/mod_instructions.R

instructions_ui <- function(id) {
  ns <- NS(id)

  tabPanel("Instructions",
           # 在这里放置说明文档
           uiOutput(ns("instructions_ui"))
  )
}

instructions_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$instructions_ui <- renderUI({
        includeMarkdown("instructions.md")
      })
    }
  )
}
