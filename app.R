###############################################################################
# 将四个独立的 Shiny 功能合并为单一应用示例（仅增强UI美观，逻辑代码原封不动）
###############################################################################

# --------------------- 新增：美学库与自定义CSS等 ----------------------------- #
library(shiny)
library(shinythemes)  # 新增：可选漂亮主题
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(stats)
library(FSA)          # Dunn’s Test
library(rlang)        # 用于 !!sym()
library(grid)         # 用于 unit()
library(colourpicker) # 用于交互式颜色选择
library(plotly)
library(stringr)
library(reticulate)
library(DT)
library(showtext)
font_add("Times New Roman", "./fonts/times.ttf")
showtext_auto()

library(patchwork)
library(ggrepel)
library(scales)


###############################################################################
# UI：在 fluidPage 中套用 shinytheme，并添加一些简单的自定义 CSS
###############################################################################
ui <- fluidPage(
  # 应用一个漂亮的主题（可自由更换为cerulean、superhero、united等）
  theme = shinytheme("flatly"),
  
  # 在页面头部添加一些自定义CSS，用于微调表格、按钮、标题等外观
  tags$head(
    tags$style(HTML("
      /* 让整体排版稍微紧凑一些 */
      body {
        font-family: 'Times New Roman', sans-serif;
        font-size: 14px;
      }

      /* 提升标题栏的视觉层次感 */
      .navbar-default .navbar-brand, 
      .navbar-default .navbar-nav > li > a {
        font-weight: 600;
      }

      /* 优化tabPanel显示 */
      .nav-pills > li > a {
        font-weight: 500;
      }
      .nav-pills > li.active > a {
        background-color: #337ab7;
        color: #fff;
      }

      /* 让 sidebarPanel 看起来更像卡片 */
      .well {
        background-color: #f9f9f9;
        border: 1px solid #eee;
        box-shadow: 0 2px 3px rgba(0,0,0,0.05);
      }

      /* 表格加上细网格线和 hover 效果 */
      table.table-hover > tbody > tr:hover {
        background-color: #f6f6f6;
      }

      /* 让 actionButton 显示为实心的 primary 按钮 */
      .btn.btn-primary {
        margin-top: 5px;
        font-weight: 500;
      }
      .btn.btn-primary.btn-block {
        width: 100%;
      }

      /* 提高一些 panel 等容器的可读性 */
      .shiny-input-container {
        margin-bottom: 15px;
      }
    "))
  ),
  
  # 页标题
  titlePanel("星海之约ZYの❥Shiny"),
  
  # 使用 type = 'pills' 让标签页更扁平化
  tabsetPanel(type = "pills",
              #=========================================================================#
              # 1. GWAS 表型统计分析工具
              #=========================================================================#
              tabPanel(
                "GWAS表型统计分析",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("gwas_file1", "上传 TSV 文件", accept = c(".tsv"), multiple = TRUE),
                    uiOutput("file_select_ui1"),
                    selectInput("trait1", "选择表型变量", choices = NULL),
                    selectInput("test_method1", "选择统计检验方法",
                                choices = c(
                                  "T检验（两两比较）" = "t.test",
                                  "ANOVA + Tukey HSD" = "anova",
                                  "Kruskal-Wallis + Dunn's Test" = "kruskal"
                                )
                    ),
                    uiOutput("pairwise_choices1"),
                    radioButtons("sig_format1", "显著性标记格式",
                                 choices = c("仅 * 号" = "stars", "完整 p 值" = "pval")),
                    radioButtons("color_mode1", "选择颜色方式",
                                 choices = c("Built-in Palette" = "palette", "Custom Colors" = "custom"),
                                 selected = "palette"),
                    uiOutput("color_ui1"),
                    selectInput("plot_format1", "选择图形格式", choices = c("png", "pdf"), selected = "png"),
                    numericInput("plot_dpi1", "设置 DPI (仅适用于 PNG)", value = 300, min = 72, max = 600),
                    actionButton("analyze1", "分析数据", class = "btn btn-primary btn-block"),
                    downloadButton("downloadTable1", "下载统计表", class = "btn btn-success"),
                    downloadButton("downloadPlot1", "下载小提琴图", class = "btn btn-success")
                  ),
                  mainPanel(
                    tabsetPanel(
                      tabPanel("数据概览", tableOutput("summaryTable1")),
                      tabPanel("统计分析", tableOutput("anovaTable1")),
                      tabPanel("两两比较", tableOutput("pairwiseTable1")),
                      tabPanel("小提琴图", plotOutput("violinPlot1"))
                    )
                  )
                )
              ),
              
              #=========================================================================#
              # 2. Gene Expression Plot (plotly)
              #=========================================================================#
              tabPanel(
                "基因表达可视化",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("upload_data2", "上传 CSV 文件（若不上传则使用内置数据）", accept = c(".csv")),
                    textOutput("upload_status2"),
                    textInput("gene_input2", "请输入基因名称（多个基因用逗号或空格分隔）：", value = ""),
                    textInput("group1_2", "匹配数据列前缀 1：", value = "HF2"),
                    textInput("group2_2", "匹配数据列前缀 2：", value = "WJ"),
                    selectInput("test_method2", "选择检验方法：",
                                choices = c("t-test", "one-tailed t-test", "Mann-Whitney U", "Kolmogorov-Smirnov"),
                                selected = "t-test"),
                    textInput("y_axis_label2", "Y轴标签标题：", value = "Log2 Normalized Expression"),
                    selectInput("transform_method2", "数据变换方式：",
                                choices = c("无" = "none", "Z-score 标准化" = "zscore", "平方根变换" = "sqrt"),
                                selected = "none"),
                    selectInput("plot_type2", "选择图形类型：",
                                choices = c("小提琴图" = "violin", "箱线图" = "boxplot", "小提琴图+箱线图" = "both"),
                                selected = "violin"),
                    selectInput("sig_display2", "显著性显示方式：",
                                choices = c("仅显示符号" = "symbol", "显示完整信息" = "full"),
                                selected = "full"),
                    checkboxInput("adjust_outliers2", "调整极值", value = FALSE),
                    actionButton("generate_plot2", "生成图表", class = "btn btn-primary btn-block"),
                    
                    h4("图例与配色设置"),
                    fluidRow(
                      column(6, textInput("legend_label1_2", "图例标签 1：", value = "HF2")),
                      column(6, textInput("legend_label2_2", "图例标签 2：", value = "WJ"))
                    ),
                    selectInput("legend_font2", "图例字体：",
                                choices = c("Times New Roman", "Arial", "Helvetica", "Courier New", "Georgia"),
                                selected = "Times New Roman"),
                    fluidRow(
                      column(6, colourInput("group1_color2", "组1颜色：", value = "#FF6347")),
                      column(6, colourInput("group2_color2", "组2颜色：", value = "#87CEFA"))
                    ),
                    actionButton("update_legend2", "更新图例", class = "btn btn-info btn-block"),
                    
                    h4("保存图片选项"),
                    numericInput("save_width2", "宽度（英寸）：", value = 6, min = 1),
                    numericInput("save_height2", "高度（英寸）：", value = 4, min = 1),
                    numericInput("save_dpi2", "DPI：", value = 300, min = 72),
                    selectInput("save_format2", "保存格式：", choices = c("png", "pdf"), selected = "png"),
                    downloadButton("download_plot2", "下载图片", class = "btn btn-success")
                  ),
                  mainPanel(
                    plotlyOutput("violin_plot2", height = "600px")
                  )
                )
              ),
              
              #=========================================================================#
              # 3. 表型数据可视化
              #=========================================================================#
              tabPanel(
                "表型数据可视化",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("file3", "Upload CSV File", accept = c(".csv")),
                    selectInput("x_var3", "Select Categorical Variable", choices = NULL),
                    selectInput("y_var3", "Select Phenotype Variable(s)", choices = NULL, multiple = TRUE),
                    selectInput("stat_test3", "Select Statistical Test",
                                choices = c("t-test" = "t.test", "ANOVA" = "anova", "Kruskal-Wallis" = "kruskal.test")),
                    checkboxInput("show_points3", "Show Points", TRUE),
                    checkboxInput("show_lines3", "显示显著性标记及分组连线", TRUE),
                    checkboxInput("show_ns3", "显示不显著标记及分组连线", TRUE),
                    selectInput("selected_comparisons3", "选择比较组对", choices = NULL, multiple = TRUE),
                    
                    h4("Customize Colors"),
                    uiOutput("color_selectors3"),
                    
                    h4("Font Settings"),
                    numericInput("font_size3", "Font Size", value = 16, min = 8, max = 30),
                    
                    h4("Save Image Settings"),
                    numericInput("subplot_width3", "Subplot Width (inch)", value = 6, min = 4, max = 20),
                    numericInput("subplot_height3", "Subplot Height (inch)", value = 4, min = 3, max = 15),
                    numericInput("img_dpi3", "Resolution (dpi)", value = 300, min = 72, max = 600),
                    radioButtons("img_format3", "Format", choices = c("PNG", "PDF")),
                    
                    h4("Layout Settings"),
                    numericInput("layout_ncol3", "Number of columns (0 for auto)", value = 0, min = 0, max = 10),
                    numericInput("layout_nrow3", "Number of rows (0 for auto)", value = 0, min = 0, max = 10),
                    
                    downloadButton("save_plot3", "Download Image", class = "btn btn-success")
                  ),
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Data Table", DTOutput("table3")),
                      tabPanel("Plot", plotOutput("boxplot3"))
                    )
                  )
                )
              ),
              
              #=========================================================================#
              # 4. LocusZoom绘图
              #=========================================================================#
              tabPanel(
                "LocusZoom绘图",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("gwas_file4", "Upload GWAS File (e.g., 47_16447785.csv):", accept = c(".csv")),
                    fileInput("gene_file4", "Upload Gene Track File (e.g., genes.up.csv):", accept = c(".csv")),
                    textInput("lead_snp4", "Enter Lead SNP ID (e.g., '47_16447785'):",
                              placeholder = "Enter SNP ID here"),
                    hr(),
                    h4("Advanced Options"),
                    numericInput("legend_square_size4", "Legend Square Size", value = 15, min = 1),
                    numericInput("lead_snp_size4", "Lead SNP Diamond Size", value = 5, min = 1),
                    numericInput("display_width4", "Display Plot Width (inches)", value = 10, min = 1),
                    numericInput("display_height4", "Display Plot Height (inches)", value = 6, min = 1),
                    numericInput("display_dpi4", "Display Plot DPI", value = 300, min = 72),
                    numericInput("save_width4", "Save Image Width (inches)", value = 3.2, min = 0.1),
                    numericInput("save_height4", "Save Image Height (inches)", value = 1.8, min = 0.1),
                    numericInput("save_dpi4", "Save Image DPI", value = 300, min = 72),
                    textInput("font_family4", "Font Family", value = "sans"),
                    actionButton("update_plot4", "Generate Plot", class = "btn btn-primary btn-block")
                  ),
                  mainPanel(
                    plotOutput("combined_plot4", width = "100%", height = "800px"),
                    downloadButton("download_pdf4", "Download PDF", class = "btn btn-success"),
                    downloadButton("download_png4", "Download PNG", class = "btn btn-success")
                  )
                )
              )
  )
)

###############################################################################
# SERVER：保持与原逻辑完全一致，未对内部任何处理做改动
###############################################################################
server <- function(input, output, session) {
  
  #---------------------------------------------------------------------------
  # 1. GWAS表型统计分析工具（原逻辑不变）
  #---------------------------------------------------------------------------
  data_list1 <- reactive({
    req(input$gwas_file1)
    files <- input$gwas_file1
    lst <- lapply(seq_len(nrow(files)), function(i) {
      df <- read_tsv(files$datapath[i], show_col_types = FALSE)
      # 替换 Genotype 中的竖线 | 为 /，与原始示例保持一致
      df <- df %>% mutate(Genotype = gsub("\\|", "/", Genotype))
      df
    })
    names(lst) <- files$name
    lst
  })
  
  output$file_select_ui1 <- renderUI({
    req(input$gwas_file1)
    choices <- input$gwas_file1$name
    selectInput("selected_file1", "选择文件", choices = choices, selected = choices[1])
  })
  
  current_data1 <- reactive({
    req(data_list1(), input$selected_file1)
    data_list1()[[input$selected_file1]]
  })
  
  observe({
    req(current_data1())
    # 假设只取第4到第9列作为“表型变量”
    updateSelectInput(session, "trait1", choices = names(current_data1())[4:9])
  })
  
  output$pairwise_choices1 <- renderUI({
    req(current_data1(), input$test_method1 == "t.test")
    df <- current_data1()
    genotypes <- unique(df$Genotype)
    pairwise_comparisons <- combn(genotypes, 2, simplify = FALSE)
    pairwise_labels <- sapply(pairwise_comparisons, function(p) paste(p[1], "vs", p[2]))
    checkboxGroupInput("selected_pairs1", "选择两两比较组", choices = pairwise_labels, selected = pairwise_labels)
  })
  
  output$color_ui1 <- renderUI({
    if (input$color_mode1 == "palette") {
      selectInput("color_palette1", "选择调色板",
                  choices = c("Default", "Set1", "Set2", "Dark2", "Pastel1", "Pastel2"),
                  selected = "Default")
    } else {
      req(current_data1())
      genotypes <- sort(unique(current_data1()$Genotype))
      n <- length(genotypes)
      defaultColors <- rainbow(n)
      color_inputs <- lapply(seq_along(genotypes), function(i) {
        colourInput(
          inputId = paste0("custom_color1_", make.names(genotypes[i])),
          label = genotypes[i],
          value = defaultColors[i]
        )
      })
      do.call(tagList, color_inputs)
    }
  })
  
  # 数据概览
  output$summaryTable1 <- renderTable({
    req(current_data1())
    df <- current_data1()
    summary_df <- df %>%
      summarise(across(where(is.numeric), list(
        Mean = mean, SD = sd, Min = min, Max = max
      )))
    summary_df
  }, rownames = TRUE)
  
  # 统计检验结果（中文，仅展示）
  analysis_result_cn1 <- reactive({
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
  
  # 英文版本（便于下载导表等）
  analysis_result_eng1 <- reactive({
    req(current_data1(), input$trait1, input$test_method1)
    df <- current_data1()
    if (input$test_method1 == "t.test") {
      return(data.frame(Message = "T-test is not applicable for overall variance analysis. Please refer to the 'Pairwise Comparison' tab."))
    } else if (input$test_method1 == "anova") {
      res <- aov(as.formula(paste(input$trait1, "~ Genotype")), data = df)
      anova_res <- summary(res)[[1]]
      data.frame(
        Test_Method = "ANOVA",
        F_value = anova_res["F value", "F value"],
        DF = anova_res["Df", "Df"],
        p_value = anova_res["Pr(>F)", "Pr(>F)"]
      )
    } else if (input$test_method1 == "kruskal") {
      res <- kruskal.test(as.formula(paste(input$trait1, "~ Genotype")), data = df)
      data.frame(
        Test_Method = "Kruskal-Wallis",
        Statistic = res$statistic,
        DF = res$parameter,
        p_value = res$p.value
      )
    }
  })
  
  output$anovaTable1 <- renderTable({
    analysis_result_cn1()
  }, rownames = TRUE)
  
  # 两两比较结果
  pairwise_result_cn1 <- reactive({
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
  
  # 小提琴图
  violin_plot1 <- reactive({
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
    
    # 自定义颜色
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
    
    # 若为 T 检验，添加显著性标记
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
  
  # 小提琴图下载
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
  
  # 统计表下载
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
  
  #---------------------------------------------------------------------------
  # 2. Gene Expression Plot (plotly) （原逻辑不变）
  #---------------------------------------------------------------------------
  data_reactive2 <- reactiveVal(NULL)
  ggplot_obj2 <- reactiveVal(NULL)
  plotly_obj2 <- reactiveVal(NULL)
  
  # 如需默认数据集，调整路径；此处仅示例
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
    # 若未上传文件且默认路径有内置数据则读入
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
    # log2 转换
    data_long <- data_long %>% mutate(Expression = log2(Expression + 1))
    
    # 调整极值
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
    
    # 计算检验 P 值
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
    
    # 其它变换方式
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
    # 更新图例文字
    updateTextInput(session, "legend_label1_2", value = input$group1_2)
    updateTextInput(session, "legend_label2_2", value = input$group2_2)
    
    p <- generate_ggplot2()
    if (is.null(p)) return()
    
    ggplot_obj2(p)
    pltly <- ggplotly(p, source = "legend") %>% layout(legend = list(itemclick = "toggleothers"))
    plotly_obj2(pltly)
    output$violin_plot2 <- renderPlotly({ pltly })
  })
  
  observeEvent(input$update_legend2, {
    req(ggplot_obj2())
    p <- generate_ggplot2()
    if (is.null(p)) return()
    
    ggplot_obj2(p)
    pltly <- ggplotly(p, source = "legend") %>% layout(legend = list(itemclick = "toggleothers"))
    
    # 更新 legend 名称
    for (i in seq_along(pltly$x$data)) {
      if (!is.null(pltly$x$data[[i]]$name)) {
        if (pltly$x$data[[i]]$name == input$group1_2) {
          pltly$x$data[[i]]$name <- input$legend_label1_2
        }
        if (pltly$x$data[[i]]$name == input$group2_2) {
          pltly$x$data[[i]]$name <- input$legend_label2_2
        }
      }
    }
    
    plotly_obj2(pltly)
    output$violin_plot2 <- renderPlotly({ pltly })
  })
  
  # 下载图片
  output$download_plot2 <- downloadHandler(
    filename = function() {
      paste0("gene_expression_plot.", input$save_format2)
    },
    content = function(file) {
      req(ggplot_obj2())  # 注意要对 ggplot_obj2() 进行保存
      # 宽高与 DPI 依然用你在UI里输入的数值
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
  
  #---------------------------------------------------------------------------
  # 3. 表型数据可视化（原逻辑不变）
  #---------------------------------------------------------------------------
  rv3 <- reactiveValues(comparisons = NULL)
  
  data3 <- reactive({
    req(input$file3)
    read.csv(input$file3$datapath, stringsAsFactors = FALSE)
  })
  
  observeEvent(data3(), {
    df <- data3()
    # 找出非数值型列
    cat_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    # 找出数值型列
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    updateSelectInput(session, "x_var3", choices = cat_cols)
    updateSelectInput(session, "y_var3", choices = numeric_cols)
  })
  
  observeEvent({
    data3()
    input$x_var3
  }, {
    df <- data3()
    if (is.null(input$x_var3) || !(input$x_var3 %in% colnames(df))) return()
    
    groups <- levels(as.factor(df[[input$x_var3]]))
    if (length(groups) >= 2) {
      comps <- combn(groups, 2, simplify = FALSE)
      comp_names <- sapply(comps, function(pair) paste(pair, collapse = " vs "))
      rv3$comparisons <- setNames(comps, comp_names)
      updateSelectInput(session, "selected_comparisons3", choices = comp_names)
    } else {
      updateSelectInput(session, "selected_comparisons3", choices = NULL)
    }
  })
  
  output$color_selectors3 <- renderUI({
    req(data3(), input$x_var3)
    df <- data3()
    if(! input$x_var3 %in% colnames(df)) {
      return(NULL)
    }
    groups <- levels(as.factor(df[[input$x_var3]]))
    if(length(groups) == 0) return(NULL)
    
    color_inputs <- lapply(groups, function(group) {
      colourInput(
        inputId = paste0("color3_", group),
        label = paste("Color for", group),
        value = sample(grDevices::colors(), 1)
      )
    })
    do.call(tagList, color_inputs)
  })
  
  get_color_mapping3 <- reactive({
    req(data3(), input$x_var3)
    df <- data3()
    groups <- levels(as.factor(df[[input$x_var3]]))
    setNames(
      sapply(groups, function(g) {
        color <- input[[paste0("color3_", g)]]
        if (is.null(color) || color == "") {
          return(sample(grDevices::colors(), 1))
        } else {
          return(color)
        }
      }, simplify = TRUE),
      groups
    )
  })
  
  create_plot3 <- function(current_y) {
    df <- data3()
    req(input$x_var3 %in% colnames(df))
    
    df[[input$x_var3]] <- as.factor(df[[input$x_var3]])
    
    # winsorize
    df_adjusted <- do.call(rbind, lapply(split(df, df[[input$x_var3]]), function(subdf) {
      if(nrow(subdf) > 0){
        low <- quantile(subdf[[current_y]], 0.05, na.rm = TRUE)
        high <- quantile(subdf[[current_y]], 0.95, na.rm = TRUE)
        subdf[[current_y]] <- pmin(pmax(subdf[[current_y]], low), high)
      }
      subdf
    }))
    
    colors <- get_color_mapping3()
    
    p <- ggplot(df_adjusted, aes_string(x = input$x_var3, y = current_y, fill = input$x_var3)) +
      scale_fill_manual(values = colors) +
      theme_classic(base_size = input$font_size3) +
      theme(
        text = element_text(family = "Times New Roman", color = "black", size = input$font_size3),
        axis.title = element_text(color = "black", size = input$font_size3),
        axis.text = element_text(color = "black", size = input$font_size3),
        plot.title = element_text(color = "black", size = input$font_size3 + 2),
        legend.title = element_text(color = "black", size = input$font_size3),
        legend.text = element_text(color = "black", size = input$font_size3)
      ) +
      geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA, lwd = 1, coef = 0)
    
    # 绘制误差棒
    error_df_top <- df_adjusted %>%
      group_by(gr = .data[[input$x_var3]]) %>%
      summarise(
        Q3 = quantile(.data[[current_y]], 0.75, na.rm = TRUE),
        se = sd(.data[[current_y]], na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(
        ymin = Q3,
        ymax = Q3 + se
      )
    
    error_df_bottom <- df_adjusted %>%
      group_by(gr = .data[[input$x_var3]]) %>%
      summarise(
        Q1 = quantile(.data[[current_y]], 0.25, na.rm = TRUE),
        se = sd(.data[[current_y]], na.rm = TRUE) / sqrt(n())
      ) %>%
      mutate(
        ymin = Q1 - se,
        ymax = Q1
      )
    
    error_df <- bind_rows(error_df_top, error_df_bottom)
    error_df$gr <- factor(error_df$gr, levels = levels(df_adjusted[[input$x_var3]]))
    
    p <- p + geom_errorbar(
      data = error_df,
      aes(x = gr, ymin = ymin, ymax = ymax),
      width = 0.15, color = "black", lwd = 1, inherit.aes = FALSE
    )
    
    if (input$show_points3) {
      p <- p + geom_jitter(aes_string(color = input$x_var3), width = 0.2, alpha = 0.5) +
        scale_color_manual(values = colors)
    }
    
    if (input$show_lines3) {
      if (!is.null(input$selected_comparisons3) && length(input$selected_comparisons3) > 0) {
        comps_to_use <- lapply(input$selected_comparisons3, function(name) rv3$comparisons[[name]])
      } else {
        x_levels <- levels(df_adjusted[[input$x_var3]])
        comps_to_use <- combn(x_levels, 2, simplify = FALSE)
      }
      
      if (input$stat_test3 == "t.test") {
        if (input$show_ns3) {
          p <- p + stat_compare_means(
            comparisons = comps_to_use,
            method = "t.test",
            label = "p.signif",
            hide.ns = FALSE
          )
        } else {
          sig_comparisons <- list()
          for (comp in comps_to_use) {
            grp1 <- comp[[1]]
            grp2 <- comp[[2]]
            data1 <- df_adjusted[df_adjusted[[input$x_var3]] == grp1, current_y]
            data2 <- df_adjusted[df_adjusted[[input$x_var3]] == grp2, current_y]
            ttest_res <- t.test(data1, data2)
            if (ttest_res$p.value < 0.05) {
              sig_comparisons <- c(sig_comparisons, list(comp))
            }
          }
          if (length(sig_comparisons) > 0) {
            p <- p + stat_compare_means(
              comparisons = sig_comparisons,
              method = "t.test",
              label = "p.signif",
              hide.ns = FALSE
            )
          }
        }
      } else if (input$stat_test3 == "anova") {
        formula <- as.formula(paste(current_y, "~", input$x_var3))
        aov_test <- aov(formula, data = df_adjusted)
        aov_summary <- summary(aov_test)
        pval <- aov_summary[[1]][["Pr(>F)"]][1]
        sig_label <- ifelse(
          pval < 0.0001, "****",
          ifelse(pval < 0.001, "***",
                 ifelse(pval < 0.01, "**",
                        ifelse(pval < 0.05, "*", "ns")))
        )
        x_levels <- levels(df_adjusted[[input$x_var3]])
        x_center <- (length(x_levels) + 1) / 2
        y_max <- max(df_adjusted[[current_y]], na.rm = TRUE)
        if (input$show_ns3 || sig_label != "ns") {
          p <- p + annotate(
            "text", x = x_center, y = y_max * 1.1,
            label = paste("ANOVA:", sig_label, "\np =", format(pval, digits = 3)),
            size = input$font_size3, color = "black"
          )
        }
      } else if (input$stat_test3 == "kruskal.test") {
        formula <- as.formula(paste(current_y, "~", input$x_var3))
        kt_test <- kruskal.test(formula, data = df_adjusted)
        pval <- kt_test$p.value
        sig_label <- ifelse(
          pval < 0.0001, "****",
          ifelse(pval < 0.001, "***",
                 ifelse(pval < 0.01, "**",
                        ifelse(pval < 0.05, "*", "ns")))
        )
        x_levels <- levels(df_adjusted[[input$x_var3]])
        x_center <- (length(x_levels) + 1) / 2
        y_max <- max(df_adjusted[[current_y]], na.rm = TRUE)
        if (input$show_ns3 || sig_label != "ns") {
          p <- p + annotate(
            "text", x = x_center, y = y_max * 1.1,
            label = paste("Kruskal:", sig_label, "\np =", format(pval, digits = 3)),
            size = input$font_size3, color = "black"
          )
        }
      }
    }
    
    p <- p + ggtitle(current_y) +
      theme(legend.position = "none", aspect.ratio = 1)
    
    p
  }
  
  combined_plot3 <- reactive({
    req(input$y_var3)
    plots <- lapply(input$y_var3, function(var) create_plot3(var))
    n <- length(plots)
    
    if (!is.null(input$layout_ncol3) && input$layout_ncol3 > 0) {
      ncol <- input$layout_ncol3
    } else {
      ncol <- if(n == 1) 1 else 2
    }
    if (!is.null(input$layout_nrow3) && input$layout_nrow3 > 0) {
      nrow <- input$layout_nrow3
    } else {
      nrow <- ceiling(n / ncol)
    }
    
    ggarrange(
      plotlist = plots,
      labels = LETTERS[1:n],
      ncol = ncol,
      nrow = nrow
    )
  })
  
  output$boxplot3 <- renderPlot({
    print(combined_plot3())
  },
  height = function(){
    req(input$y_var3)
    n <- length(input$y_var3)
    if (!is.null(input$layout_ncol3) && input$layout_ncol3 > 0) {
      ncol <- input$layout_ncol3
    } else {
      ncol <- if(n == 1) 1 else 2
    }
    if (!is.null(input$layout_nrow3) && input$layout_nrow3 > 0) {
      nrow <- input$layout_nrow3
    } else {
      nrow <- ceiling(n / ncol)
    }
    input$subplot_height3 * nrow * 96
  },
  width = function(){
    req(input$y_var3)
    n <- length(input$y_var3)
    if (!is.null(input$layout_ncol3) && input$layout_ncol3 > 0) {
      ncol <- input$layout_ncol3
    } else {
      ncol <- if(n == 1) 1 else 2
    }
    input$subplot_width3 * ncol * 96
  })
  
  output$table3 <- renderDT({
    req(data3())
    datatable(data3())
  })
  
  output$save_plot3 <- downloadHandler(
    filename = function() {
      paste0("plot.", tolower(input$img_format3))
    },
    content = function(file) {
      req(input$y_var3)
      n <- length(input$y_var3)
      if (!is.null(input$layout_ncol3) && input$layout_ncol3 > 0) {
        ncol <- input$layout_ncol3
      } else {
        ncol <- if(n == 1) 1 else 2
      }
      if (!is.null(input$layout_nrow3) && input$layout_nrow3 > 0) {
        nrow <- input$layout_nrow3
      } else {
        nrow <- ceiling(n / ncol)
      }
      overall_width <- input$subplot_width3 * ncol
      overall_height <- input$subplot_height3 * nrow
      
      combined <- combined_plot3()
      
      if (input$img_format3 == "PNG") {
        ggsave(file, plot = combined, width = overall_width, height = overall_height, dpi = input$img_dpi3)
      } else {
        ggsave(file, plot = combined, width = overall_width, height = overall_height, device = cairo_pdf)
      }
    }
  )
  
  #---------------------------------------------------------------------------
  # 4. LocusZoom Plot Generator（原逻辑不变）
  #---------------------------------------------------------------------------
  observeEvent(input$update_plot4, {
    req(input$gwas_file4)
    req(input$gene_file4)
    gwas_data <- read.csv(input$gwas_file4$datapath)
    gene_data <- read.csv(input$gene_file4$datapath)
    
    req(input$lead_snp4)
    target_snp <- input$lead_snp4
    highlight_snp <- gwas_data %>% filter(SNP == target_snp)
    
    if (nrow(highlight_snp) == 0) {
      showModal(modalDialog(
        title = "Error",
        paste("Lead SNP ID", target_snp, "not found in the GWAS file."),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    
    target_chr <- as.numeric(highlight_snp$chr)
    target_pos <- as.numeric(highlight_snp$pos)
    target_start <- min(gwas_data$pos)
    target_end <- max(gwas_data$pos)
    nudge_val <- (target_end - target_start) * 0.01
    
    gene_track_data <- gene_data %>%
      filter(chr == target_chr, start <= target_end, stop >= target_start) %>%
      filter(name != "unknown") %>%
      mutate(label_y = ifelse(row_number() %% 2 == 0, 1, 0))
    
    threshold_value <- -log10(5e-8)
    
    custom_colors <- c(
      ">0.8" = "#d73027",
      ">0.6" = "#fc8d59",
      ">0.4" = "#fee090",
      ">0.2" = "#91bfdb",
      "<0.2" = "#4575b4"
    )
    
    gwas_plot <- ggplot() +
      geom_point(
        data = gwas_data %>% filter(R2 <= 0.2),
        aes(x = pos, y = pmin(-log10(P), threshold_value), color = "<0.2"),
        size = 2.5
      ) +
      geom_point(
        data = gwas_data %>% filter(R2 > 0.2 & R2 <= 0.4),
        aes(x = pos, y = pmin(-log10(P), threshold_value), color = ">0.2"),
        size = 2.5
      ) +
      geom_point(
        data = gwas_data %>% filter(R2 > 0.4 & R2 <= 0.6),
        aes(x = pos, y = pmin(-log10(P), threshold_value), color = ">0.4"),
        size = 2.5
      ) +
      geom_point(
        data = gwas_data %>% filter(R2 > 0.6 & R2 <= 0.8),
        aes(x = pos, y = pmin(-log10(P), threshold_value), color = ">0.6"),
        size = 2.5
      ) +
      geom_point(
        data = gwas_data %>% filter(R2 > 0.8),
        aes(x = pos, y = -log10(P), color = ">0.8"),
        size = 2.5
      ) +
      # Lead SNP 菱形标记
      geom_point(
        data = highlight_snp,
        aes(x = pos, y = -log10(P)),
        shape = 23, size = input$lead_snp_size4, fill = "red", color = "black", stroke = 1
      ) +
      # 将 P 值信息贴在 Lead SNP 右侧
      geom_text(
        data = highlight_snp,
        aes(
          x = pos, y = -log10(P),
          label = paste0("LG", target_chr, ":", target_pos, " (P = ",
                         formatC(P, format = "e", digits = 2), ")")
        ),
        nudge_x = nudge_val,
        hjust = 0,
        size = 4, fontface = "bold", color = "black"
      ) +
      scale_color_manual(
        values = custom_colors,
        name = "R2",
        labels = c(">0.8", ">0.6", ">0.4", ">0.2", "<0.2"),
        breaks = c(">0.8", ">0.6", ">0.4", ">0.2", "<0.2"),
        guide = guide_legend(override.aes = list(shape = 15, size = input$legend_square_size4))
      ) +
      theme_classic(base_size = 14, base_family = ifelse(input$font_family4 == "", "sans", input$font_family4)) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        axis.line = element_line(color = "black", size = 0.8),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "right",
        legend.title = element_text(
          size = 12, face = "bold",
          family = ifelse(input$font_family4 == "", "sans", input$font_family4),
          color = "black"
        ),
        legend.text = element_text(
          size = 10,
          family = ifelse(input$font_family4 == "", "sans", input$font_family4),
          color = "black"
        ),
        plot.title = element_text(
          hjust = 0.5, size = 16, face = "bold",
          family = ifelse(input$font_family4 == "", "sans", input$font_family4),
          color = "black"
        )
      ) +
      labs(
        x = "Chromosome Position (bp)",
        y = expression(-log[10](P)),
        title = paste("LocusZoom Plot for LG", target_chr, ":", target_pos)
      )
    
    if (any(gwas_data$R2 > 0.8)) {
      gwas_plot <- gwas_plot +
        geom_hline(yintercept = threshold_value, linetype = "dashed", color = "black", size = 0.8)
    }
    
    # Gene Track
    gene_track_plot <- ggplot(gene_track_data, aes(xmin = start, xmax = stop, y = name)) +
      geom_segment(aes(x = start, xend = stop, y = name, yend = name),
                   color = "blue", size = 1.5) +
      geom_text_repel(
        aes(x = (start + stop) / 2, y = name, label = name),
        size = 4, hjust = 0.5, color = "black",
        nudge_y = 0.2, direction = "y",
        segment.color = "grey50", segment.size = 0.5,
        max.overlaps = Inf
      ) +
      scale_y_discrete(expand = c(0.2, 0)) +
      scale_x_continuous(
        breaks = function(x) {
          default_breaks <- pretty_breaks()(x)
          gene_start <- min(gene_track_data$start)
          if(any(abs(default_breaks - gene_start) < 1e-6)) {
            default_breaks
          } else {
            sort(unique(c(default_breaks, gene_start)))
          }
        },
        labels = comma,
        guide = guide_axis(n.dodge = 2, check.overlap = TRUE)
      ) +
      theme_classic(base_size = 14, base_family = ifelse(input$font_family4 == "", "sans", input$font_family4)) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        axis.line.x = element_line(color = "black", size = 0.8),
        axis.text.x = element_text(color = "black")
      ) +
      labs(x = "Chromosome Position (bp)", y = NULL)
    
    combined_plot4 <- gwas_plot / gene_track_plot + plot_layout(heights = c(3, 1))
    
    output$combined_plot4 <- renderPlot({
      combined_plot4
    },
    width = function() { input$display_width4 * input$display_dpi4 },
    height = function() { input$display_height4 * input$display_dpi4 }
    )
    
    # 下载 PDF
    output$download_pdf4 <- downloadHandler(
      filename = function() { "LocusZoom_Plot.pdf" },
      content = function(file) {
        ggsave(
          filename = file,
          plot = combined_plot4,
          width = input$save_width4  * 2,
          height = input$save_height4 * 2,
          scale = 0.5,
          dpi = input$save_dpi4,
          device = cairo_pdf,
          family = ifelse(input$font_family4 == "", "sans", input$font_family4)
        )
      }
    )
    
    # 下载 PNG
    output$download_png4 <- downloadHandler(
      filename = function() { "LocusZoom_Plot.png" },
      content = function(file) {
        ggsave(
          filename = file,
          plot = combined_plot4,
          width = input$save_width4  * 2,
          height = input$save_height4 * 2,
          scale = 0.5,
          dpi = input$save_dpi4
        )
      }
    )
  })
}

# 运行应用
shinyApp(ui, server)
