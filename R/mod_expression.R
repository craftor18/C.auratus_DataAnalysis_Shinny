DEFAULT_GENES <- c("POMC", "prl", "GHR", "IGF1", "MSTN", "EXOC1")

expression_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330,
      sidebar_section(
        "表达矩阵", "table-cells",
        data_source_ui(ns("src"), example_label = "示例：24 样本 × 14,058 基因"),
        helpText("第一列 gene_symbol，其余列为样本；样本名用前缀表示分组，如 HF1_3、WJ_10。"),
        radioButtons(ns("scale"), "数值尺度",
                     choices = c("TPM / FPKM / 标准化计数 → log₂(x+1)" = "raw", "已经是 log 尺度" = "log"))
      ),
      sidebar_section(
        "分组与基因", "object-group",
        radioButtons(ns("group_rule"), "分组规则", inline = TRUE,
                     choices = c("最后一个 _ 之前" = "prefix", "去掉末尾数字" = "digits")),
        selectizeInput(ns("groups"), "比较组（按顺序，第 1 组为对照）", choices = NULL, multiple = TRUE),
        selectizeInput(ns("genes"), "候选基因", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "输入基因名搜索，可多选", maxItems = 24))
      ),
      sidebar_section(
        "统计", "calculator",
        radioButtons(ns("method"), "检验", inline = TRUE,
                     choices = c("Welch t / ANOVA" = "welch", "Wilcoxon / Kruskal" = "wilcox")),
        selectInput(ns("p_adjust"), "多重检验校正（跨基因）", choices = P_ADJUST_CHOICES, selected = "BH")
      ),
      sidebar_section(
        "图形", "palette",
        radioButtons(ns("plot_type"), "样式", inline = TRUE, selected = "box",
                     choices = c("小提琴" = "violin", "箱线" = "box", "均值±SE" = "meanse")),
        helpText("每组样本少于 10 个时，小提琴图的核密度估计并不可靠，建议用箱线图加散点。"),
        checkboxInput(ns("free_y"), "各基因独立纵轴", TRUE),
        sliderInput(ns("ncol"), "每行基因数", min = 1, max = 6, value = 3, step = 1),
        selectInput(ns("palette"), "配色", choices = PALETTE_CHOICES, selected = "npg"),
        textInput(ns("custom_colors"), "自定义颜色（可选）", placeholder = "#E64B35, #4DBBD5"),
        layout_columns(
          col_widths = c(7, 5),
          selectInput(ns("font"), "字体", choices = FONT_CHOICES),
          numericInput(ns("base_size"), "字号", 12, min = 8, max = 24)
        )
      ),
      sidebar_section("导出", "file-export", download_ui(ns("dl"), width = 9, height = 6))
    ),
    page_header("STEP 4", "候选基因表达", "比较候选基因在不同组之间的表达：log₂ 尺度检验、FDR 校正和 Z-score 热图。", "chart-simple"),
    uiOutput(ns("boxes")),
    navset_card_underline(
      id = ns("tabs"),
      full_screen = TRUE,
      nav_panel("分组分布", value = "dist", icon = icon("chart-column"), uiOutput(ns("plot_ui"))),
      nav_panel("热图", value = "heat", icon = icon("border-all"), uiOutput(ns("heat_ui"))),
      nav_panel("统计结果", icon = icon("square-root-variable"), DTOutput(ns("tbl_stats"))),
      nav_panel("长表数据", icon = icon("table"), DTOutput(ns("tbl_long")))
    ),
    method_card(
      p(strong("对数转换："), "转录组表达量呈右偏、方差随均值增大，统一转为 log₂(x+1) 后再检验。",
        "log₂FC = 实验组均值 − 对照组均值（log₂ 尺度），+1 表示表达量高约 2 倍。"),
      p(strong("多重检验："), "同时检验多个基因时，按所选方法（默认 BH-FDR）跨基因校正。"),
      p(strong("适用范围："), "本页用于验证少量候选基因。全转录组差异表达请对原始 counts 使用 DESeq2、edgeR 或 limma-voom，",
        "它们会处理文库大小和离散度。"),
      p(strong("热图："), "每个基因在样本间做 Z-score 标准化，展示相对高低，不同基因之间的颜色不代表绝对表达量。")
    )
  )
}

expression_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    raw <- data_source_server("src", "expression")

    mat <- reactive({
      df <- raw()
      gcol <- first_or(names(df)[tolower(names(df)) %in% c("gene_symbol", "gene", "symbol", "gene_id", "id")], names(df)[1])
      names(df)[names(df) == gcol] <- "gene_symbol"
      df$gene_symbol <- as.character(df$gene_symbol)
      num <- vapply(df, is.numeric, logical(1))
      df <- df[, c("gene_symbol", names(df)[num & names(df) != "gene_symbol"]), drop = FALSE]
      validate(need(ncol(df) >= 3, "至少需要两个数值样本列。"))
      df[!duplicated(df$gene_symbol) & !is.na(df$gene_symbol), ]
    })

    observeEvent(list(mat(), input$group_rule), {
      m <- mat()
      grp <- unique(sample_group(setdiff(names(m), "gene_symbol"), input$group_rule))
      updateSelectizeInput(session, "groups", choices = grp, selected = grp)
      genes <- m$gene_symbol
      sel <- intersect(DEFAULT_GENES, genes)
      if (!length(sel)) sel <- head(genes, 4)
      updateSelectizeInput(session, "genes", choices = genes, selected = sel, server = TRUE)
    })

    long <- reactive({
      req(length(input$groups) >= 2, length(input$genes) >= 1)
      l <- expression_long(mat(), input$genes, input$groups, input$scale, input$group_rule)
      validate(need(nrow(l) > 0, "所选基因在矩阵中不存在。"))
      l
    })

    stats_tbl <- reactive(expression_tests(long(), input$method, input$p_adjust))

    output$boxes <- renderUI({
      l <- long()
      s <- stats_tbl()
      n_sig <- sum(s$padj < 0.05, na.rm = TRUE)
      layout_columns(
        fill = FALSE,
        class = "stat-row",
        stat_box("样本数", length(unique(l$sample)), "vials", "primary",
                 paste(sprintf("%s: %d", levels(l$group), as.integer(table(l$group[!duplicated(l$sample)]))), collapse = " · ")),
        stat_box("候选基因", nlevels(l$gene_symbol), "dna", "info", sprintf("矩阵共 %s 个基因", format(nrow(mat()), big.mark = ","))),
        stat_box("显著基因 (adj < 0.05)", n_sig, "star", if (n_sig > 0) "warning" else "secondary",
                 sprintf("%s，%s 校正", unique(s$test)[1], input$p_adjust))
      )
    })

    plot_obj <- reactive({
      l <- long()
      s <- stats_tbl()
      fam <- input$font %||% "sans"
      lv <- levels(l$group)
      cols <- setNames(get_palette(input$palette, length(lv), input$custom_colors), lv)

      p <- ggplot(l, aes(group, expr, fill = group, colour = group))
      if (input$plot_type == "violin") {
        p <- p + geom_violin(alpha = 0.3, colour = NA, scale = "width", width = 0.85, trim = TRUE) +
          geom_boxplot(width = 0.16, fill = "white", colour = "#1E293B", outlier.shape = NA, linewidth = 0.4)
      } else if (input$plot_type == "box") {
        p <- p + geom_boxplot(width = 0.6, alpha = 0.35, colour = "#1E293B", outlier.shape = NA, linewidth = 0.4)
      }
      p <- p + geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7, shape = 16)
      if (input$plot_type == "meanse") {
        p <- p + stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, colour = "#1E293B") +
          stat_summary(fun = mean, geom = "point", shape = 21, size = 3.4, colour = "#1E293B")
      }

      ann <- data.frame(gene_symbol = factor(s$gene, levels = levels(l$gene_symbol)),
                        label = ifelse(is.na(s$log2FC),
                                       sprintf("adj.P = %s %s", fmt_p(s$padj), s$signif),
                                       sprintf("log2FC = %+.2f\nadj.P = %s %s", s$log2FC, fmt_p(s$padj), s$signif)))
      p +
        geom_text(data = ann, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
                  hjust = -0.05, vjust = 1.25, size = input$base_size / 3.9, family = fam, lineheight = 0.95,
                  colour = "#334155") +
        facet_wrap(~gene_symbol, ncol = input$ncol, scales = if (isTRUE(input$free_y)) "free_y" else "fixed") +
        scale_fill_manual(values = cols, name = "Group") +
        scale_colour_manual(values = cols, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0.06, 0.32))) +
        labs(x = NULL, y = if (input$scale == "raw") expression(log[2]("expression + 1")) else "Expression (log scale)",
             title = "Candidate gene expression",
             caption = sprintf("%s test; %s-adjusted across genes", unique(s$test)[1], input$p_adjust)) +
        theme_pub(input$base_size, fam) +
        theme(legend.position = "top", axis.text.x = element_text(angle = 0))
    })

    plot_height <- reactive({
      n <- nlevels(long()$gene_symbol)
      max(380, ceiling(n / input$ncol) * 290 + 90)
    })
    output$plot_ui <- renderUI(plotOutput(session$ns("plot"), height = paste0(plot_height(), "px")))
    output$plot <- renderPlot(plot_obj(), res = SCREEN_DPI)

    heat_obj <- reactive({
      l <- long()
      fam <- input$font %||% "sans"
      h <- l |>
        group_by(gene_symbol) |>
        mutate(z = (expr - mean(expr, na.rm = TRUE)) / stats::sd(expr, na.rm = TRUE)) |>
        ungroup()
      h$z[!is.finite(h$z)] <- 0
      ord <- unique(h$sample[order(h$group, h$sample)])
      h$sample <- factor(h$sample, levels = ord)
      h$gene_symbol <- factor(h$gene_symbol, levels = rev(levels(l$gene_symbol)))
      ggplot(h, aes(sample, gene_symbol, fill = z)) +
        geom_tile(colour = "white", linewidth = 0.6) +
        facet_grid(. ~ group, scales = "free_x", space = "free_x") +
        scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                             limits = c(-2.5, 2.5), oob = squish, name = "Z-score") +
        labs(x = NULL, y = NULL, title = "Row-scaled expression (Z-score)") +
        theme_pub(input$base_size, fam) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = rel(0.75)),
              axis.text.y = element_text(face = "italic"),
              axis.line = element_blank(), axis.ticks = element_blank(),
              panel.grid.major.y = element_blank(), panel.spacing.x = unit(4, "pt"))
    })
    output$heat_ui <- renderUI({
      n <- nlevels(long()$gene_symbol)
      plotOutput(session$ns("heat"), height = paste0(max(300, n * 34 + 170), "px"))
    })
    output$heat <- renderPlot(heat_obj(), res = SCREEN_DPI)

    current_plot <- reactive(if (identical(input$tabs, "heat")) heat_obj() else plot_obj())
    download_server("dl", current_plot, function() if (identical(input$tabs, "heat")) "expression_heatmap" else "expression_candidates")

    output$tbl_stats <- renderDT({
      s <- stats_tbl()
      dt_table(s, p_cols = c("p", "padj"))
    })
    output$tbl_long <- renderDT({
      l <- long()
      dt_table(data.frame(基因 = l$gene_symbol, 样本 = l$sample, 组 = l$group, 原始值 = l$value, 分析值 = l$expr),
               page_length = 12)
    })
  })
}
