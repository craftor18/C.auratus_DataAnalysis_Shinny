# R/mod_expression.R

expression_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "基因表达可视化",
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("load_sample2"), "Load Sample Data", class = "btn btn-info btn-block"),
        fileInput(ns("upload_data2"), "上传 CSV 文件（若不上传则使用内置数据）", accept = c(".csv")),
        textOutput(ns("upload_status2")),
        textInput(ns("gene_input2"), "请输入基因名称（多个基因用逗号或空格分隔）：", value = ""),
        textInput(ns("group1_2"), "匹配数据列前缀 1：", value = "HF2"),
        textInput(ns("group2_2"), "匹配数据列前缀 2：", value = "WJ"),
        selectInput(ns("test_method2"), "选择检验方法：",
                    choices = c("t-test", "one-tailed t-test", "Mann-Whitney U", "Kolmogorov-Smirnov"),
                    selected = "t-test"),
        textInput(ns("y_axis_label2"), "Y轴标签标题：", value = "Log2 Normalized Expression"),
        selectInput(ns("transform_method2"), "数据变换方式：",
                    choices = c("无" = "none", "Z-score 标准化" = "zscore", "平方根变换" = "sqrt"),
                    selected = "none"),
        selectInput(ns("plot_type2"), "选择图形类型：",
                    choices = c("小提琴图" = "violin", "箱线图" = "boxplot", "小提琴图+箱线图" = "both"),
                    selected = "violin"),
        selectInput(ns("sig_display2"), "显著性显示方式：",
                    choices = c("仅显示符号" = "symbol", "显示完整信息" = "full"),
                    selected = "full"),
        checkboxInput(ns("adjust_outliers2"), "调整极值", value = FALSE),
        actionButton(ns("generate_plot2"), "生成图表", class = "btn btn-primary btn-block"),

        h4("图例与配色设置"),
        fluidRow(
          column(6, textInput(ns("legend_label1_2"), "图例标签 1：", value = "HF2")),
          column(6, textInput(ns("legend_label2_2"), "图例标签 2：", value = "WJ"))
        ),
        selectInput(ns("legend_font2"), "图例字体：",
                    choices = c("Times New Roman", "Arial", "Helvetica", "Courier New", "Georgia"),
                    selected = "Times New Roman"),
        fluidRow(
          column(6, colourInput(ns("group1_color2"), "组1颜色：", value = "#FF6347")),
          column(6, colourInput(ns("group2_color2"), "组2颜色：", value = "#87CEFA"))
        ),
        actionButton(ns("update_legend2"), "更新图例", class = "btn btn-info btn-block"),

        h4("保存图片选项"),
        numericInput(ns("save_width2"), "宽度（英寸）：", value = 6, min = 1),
        numericInput(ns("save_height2"), "高度（英寸）：", value = 4, min = 1),
        numericInput(ns("save_dpi2"), "DPI：", value = 300, min = 72),
        selectInput(ns("save_format2"), "保存格式：", choices = c("png", "pdf"), selected = "png"),
        downloadButton(ns("download_plot2"), "下载图片", class = "btn btn-success")
      ),
      mainPanel(
        plotOutput(ns("violin_plot2"), height = "600px")
      )
    )
  )
}

expression_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      data_reactive2 <- reactiveVal(NULL)
      ggplot_obj2 <- reactiveVal(NULL)

      observeEvent(input$load_sample2, {
        df <- read.csv("data/expression_sample.csv", stringsAsFactors = FALSE)
        data_reactive2(df)
        updateTextInput(session, "gene_input2", value = "GeneA, GeneB, GeneC")
        output$upload_status2 <- renderText("Sample expression data loaded!")
      })

      DEFAULT_DATA_PATH2 <- file.path("data", "dataset.csv")
      load_default_data2 <- function() {
        if (file.exists(DEFAULT_DATA_PATH2)) {
          tryCatch({
            read.csv(DEFAULT_DATA_PATH2, stringsAsFactors = FALSE)
          }, error = function(e) {
            message("Error loading default dataset: ", e)
            NULL
          })
        } else {
          message("Default dataset not found at: ", DEFAULT_DATA_PATH2)
          NULL
        }
      }

      observeEvent(input$upload_data2, {
        req(input$upload_data2)
        file <- input$upload_data2
        df <- tryCatch({
          read.csv(file$datapath, stringsAsFactors = FALSE)
        }, error = function(e) NULL)
        if (is.null(df)) {
          output$upload_status2 <- renderText("上传文件处理错误。")
          data_reactive2(NULL)
        } else if (!("gene_symbol" %in% names(df))) {
          output$upload_status2 <- renderText("上传文件中必须包含 'gene_symbol' 列。")
          data_reactive2(NULL)
        } else {
          output$upload_status2 <- renderText(paste("文件", file$name, "上传成功！"))
          data_reactive2(df)
        }
      })

      observe({
        if (is.null(input$upload_data2) && is.null(data_reactive2())) {
          df_default <- load_default_data2()
          if (!is.null(df_default)) {
            output$upload_status2 <- renderText("使用内置数据集。")
            data_reactive2(df_default)
          } else {
            output$upload_status2 <- renderText("未上传文件且未找到内置数据集。")
          }
        }
      })

      generate_ggplot2 <- function() {
        df <- data_reactive2()
        req(df)
        gene_input <- input$gene_input2
        gene_list <- str_split(gene_input, pattern = "\\s*,\\s*|\\s+")[[1]]
        gene_list <- gene_list[gene_list != ""]
        if (length(gene_list) == 0) {
          showNotification("请至少输入一个基因名称。", type = "error")
          return(NULL)
        }
        gene_data <- df %>% filter(gene_symbol %in% gene_list)
        if (nrow(gene_data) == 0) {
          showNotification("未找到匹配的基因数据。", type = "error")
          return(NULL)
        }

        prefix1 <- input$group1_2
        prefix2 <- input$group2_2
        group1_cols <- grep(paste0("^", prefix1), names(gene_data), value = TRUE)
        group2_cols <- grep(paste0("^", prefix2), names(gene_data), value = TRUE)

        data_group1 <- gene_data %>%
          select(gene_symbol, all_of(group1_cols)) %>%
          pivot_longer(cols = group1_cols, names_to = "sample", values_to = "Expression") %>%
          mutate(Group = prefix1)

        data_group2 <- gene_data %>%
          select(gene_symbol, all_of(group2_cols)) %>%
          pivot_longer(cols = group2_cols, names_to = "sample", values_to = "Expression") %>%
          mutate(Group = prefix2)

        data_long <- bind_rows(data_group1, data_group2)
        data_long$Expression <- as.numeric(data_long$Expression)
        data_long <- data_long %>% mutate(Expression = log2(Expression + 1))

        if (input$adjust_outliers2) {
          data_long <- data_long %>%
            group_by(gene_symbol, Group) %>%
            mutate(
              lower = quantile(Expression, 0.25, na.rm = TRUE),
              upper = quantile(Expression, 0.75, na.rm = TRUE),
              iqr = upper - lower,
              lower_bound = lower - 1.5 * iqr,
              upper_bound = upper + 1.5 * iqr,
              Expression = ifelse(Expression < lower_bound, lower_bound,
                                  ifelse(Expression > upper_bound, upper_bound, Expression))
            ) %>%
            ungroup() %>%
            select(-lower, -upper, -iqr, -lower_bound, -upper_bound)
        }

        test_method <- input$test_method2
        p_values_df <- data_long %>% group_by(gene_symbol) %>%
          summarise(
            p_val = if (test_method == "t-test") {
              t.test(Expression[Group == prefix1],
                     Expression[Group == prefix2],
                     var.equal = FALSE)$p.value
            } else if (test_method == "one-tailed t-test") {
              t.test(Expression[Group == prefix1],
                     Expression[Group == prefix2],
                     alternative = "greater",
                     var.equal = FALSE)$p.value
            } else if (test_method == "Mann-Whitney U") {
              wilcox.test(Expression[Group == prefix1],
                          Expression[Group == prefix2],
                          alternative = "greater")$p.value
            } else if (test_method == "Kolmogorov-Smirnov") {
              ks.test(Expression[Group == prefix1],
                      Expression[Group == prefix2])$p.value
            }
          ) %>%
          mutate(significance = case_when(
            p_val < 0.001 ~ "***",
            p_val < 0.01  ~ "**",
            p_val < 0.05  ~ "*",
            TRUE          ~ "ns"
          ))

        if (input$transform_method2 == "zscore") {
          data_long <- data_long %>%
            group_by(gene_symbol) %>%
            mutate(Expression_plot = (Expression - mean(Expression, na.rm = TRUE)) /
                     sd(Expression, na.rm = TRUE)) %>%
            ungroup()
        } else if (input$transform_method2 == "sqrt") {
          data_long <- data_long %>% mutate(Expression_plot = sqrt(Expression))
        } else {
          data_long <- data_long %>% mutate(Expression_plot = Expression)
        }

        data_long <- data_long %>%
          mutate(
            gene_numeric = as.numeric(factor(gene_symbol, levels = gene_list)),
            dodge_offset = ifelse(Group == prefix1, -0.2, 0.2),
            x_dodge = gene_numeric + dodge_offset
          )

        lab1 <- if (!is.null(input$legend_label1_2) && input$legend_label1_2 != "") input$legend_label1_2 else prefix1
        lab2 <- if (!is.null(input$legend_label2_2) && input$legend_label2_2 != "") input$legend_label2_2 else prefix2

        p <- ggplot(data_long, aes(x = x_dodge, y = Expression_plot, fill = Group,
                                   group = interaction(gene_symbol, Group))) +
          labs(x = "Gene", y = input$y_axis_label2) +
          scale_fill_manual(
            values = setNames(c(input$group1_color2, input$group2_color2), c(prefix1, prefix2)),
            breaks = c(prefix1, prefix2),
            labels = c(lab1, lab2)
          ) +
          scale_x_continuous(breaks = 1:length(gene_list), labels = gene_list) +
          theme_classic(base_family = input$legend_font2) +
          theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white"),
            legend.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12)
          )

        if (input$plot_type2 == "violin") {
          p <- p + geom_violin(trim = FALSE, width = 0.3) +
            stat_summary(fun = "mean", geom = "point", shape = 23, size = 3, fill = "white")
        } else if (input$plot_type2 == "boxplot") {
          p <- p + geom_boxplot(width = 0.2)
        } else if (input$plot_type2 == "both") {
          p <- p + geom_violin(trim = FALSE, width = 0.3, alpha = 0.5) +
            geom_boxplot(width = 0.15)
        }

        annotation_df <- data_long %>%
          group_by(gene_symbol) %>%
          summarise(
            x_center = first(gene_numeric),
            y_max = max(Expression_plot, na.rm = TRUE)
          ) %>%
          left_join(p_values_df, by = "gene_symbol")

        if (input$sig_display2 == "symbol") {
          annotation_df$label <- annotation_df$significance
        } else {
          annotation_df$label <- paste0(
            "p = ", formatC(annotation_df$p_val, format = "e", digits = 2),
            " ", annotation_df$significance
          )
        }

        p <- p + geom_text(
          data = annotation_df,
          aes(
            x = x_center,
            y = y_max + 0.05 * (max(data_long$Expression_plot, na.rm = TRUE)),
            label = label
          ),
          vjust = 0, size = 3, inherit.aes = FALSE
        ) +
          expand_limits(y = max(data_long$Expression_plot, na.rm = TRUE) * 1.1)

        p
      }

      observeEvent(input$generate_plot2, {
        req(input$gene_input2)
        req(data_reactive2())
        updateTextInput(session, "legend_label1_2", value = input$group1_2)
        updateTextInput(session, "legend_label2_2", value = input$group2_2)

        p <- generate_ggplot2()
        if (is.null(p)) return()

        ggplot_obj2(p)
        output$violin_plot2 <- renderPlot({ p })
      })

      observeEvent(input$update_legend2, {
        req(ggplot_obj2())
        p <- generate_ggplot2()
        if (is.null(p)) return()

        ggplot_obj2(p)
        output$violin_plot2 <- renderPlot({ p })
      })

      output$download_plot2 <- downloadHandler(
        filename = function() {
          paste0("gene_expression_plot.", input$save_format2)
        },
        content = function(file) {
          req(ggplot_obj2())
          if (input$save_format2 == "pdf") {
            ggsave(
              filename = file,
              plot     = ggplot_obj2(),
              device   = "pdf",
              width    = input$save_width2,
              height   = input$save_height2,
              dpi      = input$save_dpi2
            )
          } else {
            ggsave(
              filename = file,
              plot     = ggplot_obj2(),
              device   = "png",
              width    = input$save_width2,
              height   = input$save_height2,
              dpi      = input$save_dpi2
            )
          }
        }
      )
    }
  )
}
