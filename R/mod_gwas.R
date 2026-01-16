# R/mod_gwas.R

gwas_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "GWAS表型统计分析",
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("load_sample1"), "Load Sample Data", class = "btn btn-info btn-block"),
        fileInput(ns("gwas_file1"), "上传 TSV 文件", accept = c(".tsv"), multiple = TRUE),
        uiOutput(ns("file_select_ui1")),
        selectInput(ns("trait1"), "选择表型变量", choices = NULL),
        selectInput(ns("test_method1"), "选择统计检验方法",
                    choices = c(
                      "T检验（两两比较）" = "t.test",
                      "ANOVA + Tukey HSD" = "anova",
                      "Kruskal-Wallis + Dunn's Test" = "kruskal"
                    )
        ),
        uiOutput(ns("pairwise_choices1")),
        radioButtons(ns("sig_format1"), "显著性标记格式",
                     choices = c("仅 * 号" = "stars", "完整 p 值" = "pval")),
        radioButtons(ns("color_mode1"), "选择颜色方式",
                     choices = c("Built-in Palette" = "palette", "Custom Colors" = "custom"),
                     selected = "palette"),
        uiOutput(ns("color_ui1")),
        selectInput(ns("plot_format1"), "选择图形格式", choices = c("png", "pdf"), selected = "png"),
        numericInput(ns("plot_dpi1"), "设置 DPI (仅适用于 PNG)", value = 300, min = 72, max = 600),
        actionButton(ns("analyze1"), "分析数据", class = "btn btn-primary btn-block"),
        downloadButton(ns("downloadTable1"), "下载统计表", class = "btn btn-success"),
        downloadButton(ns("downloadPlot1"), "下载小提琴图", class = "btn btn-success")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("数据概览", tableOutput(ns("summaryTable1"))),
          tabPanel("统计分析", tableOutput(ns("anovaTable1"))),
          tabPanel("两两比较", tableOutput(ns("pairwiseTable1"))),
          tabPanel("小提琴图", plotOutput(ns("violinPlot1")))
        )
      )
    )
  )
}

gwas_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data_list_rv1 <- reactiveVal(NULL)

      observeEvent(input$gwas_file1, {
        req(input$gwas_file1)
        files <- input$gwas_file1
        lst <- lapply(seq_len(nrow(files)), function(i) {
          df <- read_tsv(files$datapath[i], show_col_types = FALSE)
          df <- df %>% mutate(Genotype = gsub("\\|", "/", Genotype))
          df
        })
        names(lst) <- files$name
        data_list_rv1(lst)
      })

      observeEvent(input$load_sample1, {
        df <- read_tsv("data/gwas_sample.tsv", show_col_types = FALSE)
        df <- df %>% mutate(Genotype = gsub("\\|", "/", Genotype))
        lst <- list("gwas_sample.tsv" = df)
        data_list_rv1(lst)
      })

      data_list1 <- reactive({
        req(data_list_rv1())
        data_list_rv1()
      })

      output$file_select_ui1 <- renderUI({
        req(data_list1())
        choices <- names(data_list1())
        selectInput(session$ns("selected_file1"), "选择文件", choices = choices, selected = choices[1])
      })

      current_data1 <- reactive({
        req(data_list1(), input$selected_file1)
        data_list1()[[input$selected_file1]]
      })

      observe({
        req(current_data1())
        updateSelectInput(session, "trait1", choices = names(current_data1())[4:9])
      })

      output$pairwise_choices1 <- renderUI({
        req(current_data1(), input$test_method1 == "t.test")
        df <- current_data1()
        genotypes <- unique(df$Genotype)
        pairwise_comparisons <- combn(genotypes, 2, simplify = FALSE)
        pairwise_labels <- sapply(pairwise_comparisons, function(p) paste(p[1], "vs", p[2]))
        checkboxGroupInput(session$ns("selected_pairs1"), "选择两两比较组", choices = pairwise_labels, selected = pairwise_labels)
      })

      output$color_ui1 <- renderUI({
        if (input$color_mode1 == "palette") {
          selectInput(session$ns("color_palette1"), "选择调色板",
                      choices = c("Default", "Set1", "Set2", "Dark2", "Pastel1", "Pastel2"),
                      selected = "Default")
        } else {
          req(current_data1())
          genotypes <- sort(unique(current_data1()$Genotype))
          n <- length(genotypes)
          defaultColors <- rainbow(n)
          color_inputs <- lapply(seq_along(genotypes), function(i) {
            colourInput(
              inputId = session$ns(paste0("custom_color1_", make.names(genotypes[i]))),
              label = genotypes[i],
              value = defaultColors[i]
            )
          })
          do.call(tagList, color_inputs)
        }
      })

      output$summaryTable1 <- renderTable({
        req(current_data1())
        df <- current_data1()
        summary_df <- df %>%
          summarise(across(where(is.numeric), list(
            Mean = mean, SD = sd, Min = min, Max = max
          )))
        summary_df
      }, rownames = TRUE)

      analysis_result_cn1 <- eventReactive(input$analyze1, {
        req(current_data1(), input$trait1, input$test_method1)
        df <- current_data1()
        if (input$test_method1 == "t.test") {
          return(data.frame(提示 = "T 检验不适用于整体方差分析，请查看‘两两比较’"))
        } else if (input$test_method1 == "anova") {
          res <- aov(as.formula(paste(input$trait1, "~ Genotype")), data = df)
          anova_res <- summary(res)[[1]]
          data.frame(
            检验方法 = "ANOVA",
            F值 = anova_res["F value", "F value"],
            自由度 = anova_res["Df", "Df"],
            p值 = anova_res["Pr(>F)", "Pr(>F)"]
          )
        } else if (input$test_method1 == "kruskal") {
          res <- kruskal.test(as.formula(paste(input$trait1, "~ Genotype")), data = df)
          data.frame(
            检验方法 = "Kruskal-Wallis",
            统计量 = res$statistic,
            自由度 = res$parameter,
            p值 = res$p.value
          )
        }
      })

      output$anovaTable1 <- renderTable({
        analysis_result_cn1()
      }, rownames = TRUE)

      pairwise_result_cn1 <- eventReactive(input$analyze1, {
        req(current_data1(), input$trait1, input$test_method1)
        if (input$test_method1 != "t.test") return(NULL)
        df <- current_data1()
        selected_pairs <- input$selected_pairs1
        if (is.null(selected_pairs) || length(selected_pairs) == 0) {
          return(data.frame(提示 = "请选择至少一个两两比较组"))
        }
        results <- lapply(selected_pairs, function(pair) {
          groups <- unlist(strsplit(pair, " vs "))
          subset_data <- df %>% filter(Genotype %in% groups)
          res <- t.test(as.formula(paste(input$trait1, "~ Genotype")), data = subset_data)
          data.frame(比较组 = pair, 统计量 = res$statistic, 自由度 = res$parameter, p值 = res$p.value)
        })
        do.call(rbind, results)
      })

      output$pairwiseTable1 <- renderTable({
        if (input$test_method1 == "t.test") {
          pairwise_result_cn1()
        } else {
          data.frame(提示 = "仅 T 检验支持两两比较")
        }
      }, rownames = TRUE)

      violin_plot1 <- eventReactive(input$analyze1, {
        req(current_data1(), input$trait1)
        df <- current_data1()

        p <- ggplot(df, aes(x = Genotype, y = !!sym(input$trait1), fill = Genotype)) +
          geom_violin(trim = FALSE, alpha = 0.7, color = "black",
                      size = 1.2, scale = "width",
                      draw_quantiles = c(0.25, 0.5, 0.75)) +
          labs(title = paste("Genotype to", input$trait1, "'s impact"), x = "Genotype", y = input$trait1) +
          theme_minimal(base_size = 16) +
          theme(
            axis.line.x.bottom = element_line(color = "black", linewidth = 0.8),
            axis.line.y.left   = element_line(color = "black", linewidth = 0.8),
            axis.ticks.length.x = unit(-0.15, "cm"),
            axis.ticks.length.y = unit(-0.15, "cm"),
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),
            axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")),
            panel.grid = element_blank()
          ) +
          coord_cartesian(clip = "off")

        if (input$color_mode1 == "custom") {
          req(current_data1())
          genotypes <- sort(unique(current_data1()$Genotype))
          colors <- sapply(genotypes, function(g) {
            input[[paste0("custom_color1_", make.names(g))]]
          })
          p <- p + scale_fill_manual(values = colors)

        } else if (input$color_mode1 == "palette") {
          if (input$color_palette1 != "Default") {
            p <- p + scale_fill_brewer(palette = input$color_palette1)
          }
        }

        if (input$test_method1 == "t.test") {
          req(input$selected_pairs1)
          comparisons <- lapply(input$selected_pairs1, function(pair) unlist(strsplit(pair, " vs ")))
          p <- p + stat_compare_means(
            comparisons = comparisons,
            label = ifelse(input$sig_format1 == "stars", "p.signif", "p.format")
          )
        }
        p
      })

      output$violinPlot1 <- renderPlot({
        print(violin_plot1())
      }, height = 600, width = 800)

      output$downloadPlot1 <- downloadHandler(
        filename = function() {
          ext <- ifelse(input$plot_format1 == "png", "png", "pdf")
          paste("ViolinPlot_", input$trait1, ".", ext, sep = "")
        },
        content = function(file) {
          if (input$plot_format1 == "pdf") {
            ggsave(file, plot = violin_plot1(), device = "pdf", width = 8, height = 6)
          } else {
            ggsave(file, plot = violin_plot1(), width = 8, height = 6, dpi = input$plot_dpi1)
          }
        }
      )

      output$downloadTable1 <- downloadHandler(
        filename = function() { paste("Statistics_", input$test_method1, ".csv", sep = "") },
        content = function(file) {
          req(data_list1(), input$trait1, input$test_method1)
          combined_table <- do.call(rbind, lapply(names(data_list1()), function(fname) {
            df <- data_list1()[[fname]]
            if (input$test_method1 == "t.test") {
              selected_pairs <- input$selected_pairs1
              if (is.null(selected_pairs) || length(selected_pairs) == 0) {
                return(data.frame(File = fname, Message = "No pairwise comparisons selected"))
              }
              results <- lapply(selected_pairs, function(pair) {
                groups <- unlist(strsplit(pair, " vs "))
                subset_data <- df %>% filter(Genotype %in% groups)
                res <- t.test(as.formula(paste(input$trait1, "~ Genotype")), data = subset_data)
                data.frame(File = fname, Comparison = pair,
                           Statistic = res$statistic,
                           DF = res$parameter,
                           p_value = res$p.value)
              })
              do.call(rbind, results)
            } else if (input$test_method1 == "anova") {
              res <- aov(as.formula(paste(input$trait1, "~ Genotype")), data = df)
              anova_res <- summary(res)[[1]]
              data.frame(File = fname,
                         Test_Method = "ANOVA",
                         F_value = anova_res["F value", "F value"],
                         DF = anova_res["Df", "Df"],
                         p_value = anova_res["Pr(>F)", "Pr(>F)"])
            } else if (input$test_method1 == "kruskal") {
              res <- kruskal.test(as.formula(paste(input$trait1, "~ Genotype")), data = df)
              data.frame(File = fname,
                         Test_Method = "Kruskal-Wallis",
                         Statistic = res$statistic,
                         DF = res$parameter,
                         p_value = res$p.value)
            }
          }))
          write.csv(combined_table, file, row.names = FALSE, fileEncoding = "UTF-8")
        }
      )
    }
  )
}
