group_plot <- function(d, title, ylab, pw, cols, opts) {
  fam <- opts$font %||% "sans"
  n_lab <- table(d$g)
  p <- ggplot(d, aes(g, y, fill = g, colour = g))
  if (opts$plot_type == "violin") {
    p <- p + geom_violin(alpha = 0.28, colour = NA, trim = TRUE, scale = "width", width = 0.85) +
      geom_boxplot(width = 0.16, fill = "white", colour = "#1E293B", outlier.shape = NA, linewidth = 0.4)
  } else if (opts$plot_type == "box") {
    p <- p + geom_boxplot(width = 0.6, alpha = 0.35, colour = "#1E293B", outlier.shape = NA, linewidth = 0.4)
  }
  if (isTRUE(opts$show_points)) p <- p + geom_jitter(width = 0.14, height = 0, size = 1.3, alpha = 0.55, shape = 16)
  if (opts$plot_type == "meanse") {
    p <- p + stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, colour = "#1E293B") +
      stat_summary(fun = mean, geom = "point", shape = 21, size = 3.4, colour = "#1E293B")
  } else {
    p <- p + stat_summary(fun = mean, geom = "point", shape = 23, size = 2.2, fill = "white", colour = "#1E293B")
  }
  if (isTRUE(opts$show_pairs) && !is.null(pw)) {
    bp <- bracket_positions(pw, d$y, step_frac = 0.075, only_signif = isTRUE(opts$only_signif))
    if (!is.null(bp) && nrow(bp)) {
      bp$label <- if (opts$label_style == "stars") bp$signif else paste0("p = ", fmt_p(bp$p.adj))
      p <- p + ggpubr::stat_pvalue_manual(bp, label = "label", inherit.aes = FALSE, tip.length = 0.008,
                                          size = opts$base_size / 3.8, bracket.size = 0.35, family = fam)
    }
  }
  p +
    scale_fill_manual(values = cols, guide = "none") +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_discrete(labels = function(l) sprintf("%s\n(n=%d)", l, as.integer(n_lab[l]))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.06))) +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub(opts$base_size, fam)
}

phenotype_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330,
      sidebar_section(
        "数据", "database",
        data_source_ui(ns("src"), example_label = "示例：4 群体 400 尾形态学"),
        helpText("每行一个个体：至少一个分组列 + 数值性状列。")
      ),
      sidebar_section(
        "变量", "sliders",
        selectInput(ns("group"), "分组变量", choices = NULL),
        selectizeInput(ns("traits"), "性状（可多选）", choices = NULL, multiple = TRUE),
        selectInput(ns("ratio_by"), "异速校正：性状 ÷", choices = c("不校正" = "")),
        layout_columns(
          col_widths = c(6, 6),
          selectInput(ns("filter_var"), "子集变量", choices = c("全部个体" = "")),
          selectizeInput(ns("filter_vals"), "保留", choices = NULL, multiple = TRUE)
        ),
        selectInput(ns("outliers"), "异常值处理", choices = OUTLIER_CHOICES)
      ),
      sidebar_section(
        "统计", "calculator",
        selectInput(ns("pw_method"), "两两比较", choices = PAIRWISE_CHOICES, selected = "tukey"),
        selectInput(ns("p_adjust"), "多重检验校正", choices = P_ADJUST_CHOICES, selected = "BH")
      ),
      sidebar_section(
        "图形", "palette",
        radioButtons(ns("plot_type"), "样式", inline = TRUE,
                     choices = c("小提琴" = "violin", "箱线" = "box", "均值±SE" = "meanse")),
        checkboxInput(ns("show_points"), "显示个体散点", TRUE),
        checkboxInput(ns("show_pairs"), "显示两两比较括号", TRUE),
        checkboxInput(ns("only_signif"), "只标注显著（p.adj < 0.05）", TRUE),
        radioButtons(ns("label_style"), "显著性标签", inline = TRUE, choices = c("星号" = "stars", "p 值" = "p")),
        sliderInput(ns("ncol"), "每行子图数", min = 1, max = 4, value = 2, step = 1),
        selectInput(ns("palette"), "配色", choices = PALETTE_CHOICES, selected = "okabe_ito"),
        textInput(ns("custom_colors"), "自定义颜色（可选）", placeholder = "#0072B2, #E69F00, ..."),
        layout_columns(
          col_widths = c(7, 5),
          selectInput(ns("font"), "字体", choices = FONT_CHOICES),
          numericInput(ns("base_size"), "字号", 12, min = 8, max = 24)
        )
      ),
      sidebar_section("导出", "file-export", download_ui(ns("dl"), width = 9, height = 7))
    ),
    page_header("STEP 1", "群体表型比较", "先描述各群体（或处理组）的表型分布与差异，再做关联分析。", "fish"),
    uiOutput(ns("boxes")),
    navset_card_underline(
      id = ns("tabs"),
      full_screen = TRUE,
      nav_panel("分布图", value = "dist", icon = icon("chart-column"), uiOutput(ns("plot_ui"))),
      nav_panel("描述统计", value = "desc", icon = icon("list"), DTOutput(ns("tbl_desc"))),
      nav_panel("检验结果", value = "test", icon = icon("square-root-variable"),
                h6(class = "table-title", "总体检验"), DTOutput(ns("tbl_omni")),
                h6(class = "table-title", "两两比较"), DTOutput(ns("tbl_pw"))),
      nav_panel("性状相关", value = "cor", icon = icon("diagram-project"), plotOutput(ns("cor_plot"), height = "520px")),
      nav_panel("数据", value = "data", icon = icon("table"), DTOutput(ns("tbl_data")))
    ),
    method_card(
      p(strong("总体检验："), "同时给出经典 ANOVA、Welch ANOVA（不假设方差齐）和 Kruskal–Wallis（不假设正态）。",
        "三者结论一致时更可靠。"),
      p(strong("两两比较："), "比较 k 个组需要做 k(k−1)/2 次检验，必须做多重校正。Tukey HSD 自带族错误率控制，",
        "其他方法按所选方式（默认 BH-FDR）校正。"),
      p(strong("异速生长："), "体重、体高等形态性状随体长增大而增大。比较体型差异时，",
        "可用「性状 ÷ 体长」消除个体大小影响，更严格的做法是以体长为协变量做 ANCOVA。"),
      p(strong("混杂因素："), "示例中各群体的性别比例差异很大（如 WR 全为雌性），建议用「子集」按性别分开比较。"),
      p(strong("异常值："), "默认保留全部数据。剔除或缩尾会改变结果，务必在论文方法中注明。")
    )
  )
}

phenotype_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    raw <- data_source_server("src", "population_phenotype")

    observeEvent(raw(), {
      df <- raw()
      cats <- categorical_columns(df)
      nums <- numeric_columns(df)
      nums <- nums[!grepl("^(iid|id|sample|no)$", nums, ignore.case = TRUE)]
      updateSelectInput(session, "group", choices = cats,
                        selected = first_or(intersect(c("FID", "Population", "Group", "Treatment"), cats), first_or(cats)))
      sel <- intersect(c("Body_weight", "Body_length", "Body_height", "Head_length"), nums)
      if (!length(sel)) sel <- head(nums, 2)
      updateSelectizeInput(session, "traits", choices = nums, selected = sel)
      updateSelectInput(session, "ratio_by", choices = c("不校正" = "", nums))
      updateSelectInput(session, "filter_var", choices = c("全部个体" = "", cats))
    })

    observeEvent(input$filter_var, {
      df <- raw()
      if (nzchar(input$filter_var %||% "") && input$filter_var %in% names(df)) {
        v <- sort(unique(stats::na.omit(as.character(df[[input$filter_var]]))))
        updateSelectizeInput(session, "filter_vals", choices = v, selected = v)
      } else {
        updateSelectizeInput(session, "filter_vals", choices = character(0), selected = character(0))
      }
    })

    filtered <- reactive({
      df <- raw()
      req(input$group %in% names(df))
      if (nzchar(input$filter_var %||% "") && length(input$filter_vals)) {
        df <- df[as.character(df[[input$filter_var]]) %in% input$filter_vals, , drop = FALSE]
      }
      validate(need(nrow(df) > 0, "筛选后没有个体。"))
      df
    })

    trait_list <- reactive({
      df <- filtered()
      traits <- intersect(input$traits, names(df))
      validate(need(length(traits) >= 1, "请至少选择一个性状。"))
      ratio <- input$ratio_by %||% ""
      lapply(setNames(traits, traits), function(tr) {
        y <- df[[tr]]
        label <- tr
        if (nzchar(ratio) && ratio != tr && ratio %in% names(df)) {
          y <- y / df[[ratio]]
          label <- sprintf("%s / %s", tr, ratio)
        }
        d <- data.frame(y = y, g = factor(df[[input$group]]))
        d <- d[!is.na(d$y) & !is.na(d$g), ]
        n_before <- nrow(d)
        d <- handle_outliers(d, "y", "g", input$outliers)
        d$g <- droplevels(d$g)
        list(
          trait = tr, label = label, d = d, n_removed = n_before - nrow(d),
          omni = omnibus_tests(d$y, d$g),
          pw = pairwise_tests(d$y, d$g, input$pw_method, input$p_adjust)
        )
      })
    })

    output$boxes <- renderUI({
      df <- filtered()
      tl <- trait_list()
      g <- factor(df[[input$group]])
      kw_sig <- sum(vapply(tl, function(x) {
        o <- x$omni
        any(o$test == "Kruskal–Wallis" & o$p < 0.05)
      }, logical(1)))
      removed <- sum(vapply(tl, `[[`, numeric(1), "n_removed"))
      layout_columns(
        fill = FALSE,
        class = "stat-row",
        stat_box("个体数", nrow(df), "fish", "primary",
                 if (nzchar(input$filter_var %||% "")) sprintf("子集：%s = %s", input$filter_var, paste(input$filter_vals, collapse = "/")) else "全部个体"),
        stat_box("分组", nlevels(g), "layer-group", "info", paste(levels(g), collapse = " · ")),
        stat_box("组间差异显著的性状", sprintf("%d / %d", kw_sig, length(tl)), "star", if (kw_sig) "warning" else "secondary",
                 "Kruskal–Wallis p < 0.05"),
        stat_box("剔除的异常值", removed, "filter", if (removed) "danger" else "success",
                 names(OUTLIER_CHOICES)[OUTLIER_CHOICES == input$outliers])
      )
    })

    plot_obj <- reactive({
      tl <- trait_list()
      lv <- levels(factor(filtered()[[input$group]]))
      cols <- setNames(get_palette(input$palette, length(lv), input$custom_colors), lv)
      opts <- reactiveValuesToList(input)[c("plot_type", "show_points", "show_pairs", "only_signif",
                                            "label_style", "base_size", "font")]
      plots <- lapply(tl, function(x) group_plot(x$d, x$label, x$label, x$pw, cols, opts))
      wrap_plots(plots, ncol = min(input$ncol, length(plots))) +
        plot_annotation(
          tag_levels = "A",
          caption = sprintf("Pairwise: %s; adjustment: %s", unique(tl[[1]]$pw$method)[1], input$p_adjust),
          theme = theme(plot.caption = element_text(colour = "#64748B", hjust = 0))
        ) &
        theme(plot.tag = element_text(face = "bold"))
    })

    output$plot_ui <- renderUI({
      n <- length(trait_list())
      rows <- ceiling(n / min(input$ncol, n))
      plotOutput(session$ns("plot"), height = paste0(max(420, rows * 400), "px"))
    })
    output$plot <- renderPlot(plot_obj(), res = SCREEN_DPI)

    cor_obj <- reactive({
      df <- filtered()
      nums <- intersect(numeric_columns(df), names(df))
      nums <- nums[!grepl("^(iid|id|sample|no)$", nums, ignore.case = TRUE)]
      validate(need(length(nums) >= 2, "至少需要两个数值性状。"))
      cm <- stats::cor(df[, nums], use = "pairwise.complete.obs", method = "pearson")
      cl <- as.data.frame(as.table(cm))
      names(cl) <- c("a", "b", "r")
      cl$a <- factor(cl$a, levels = nums)
      cl$b <- factor(cl$b, levels = rev(nums))
      ggplot(cl, aes(a, b, fill = r)) +
        geom_tile(colour = "white", linewidth = 0.8) +
        geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.6), size = 3.6, family = input$font %||% "sans") +
        scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#0F172A"), guide = "none") +
        scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", limits = c(-1, 1), name = "Pearson r") +
        coord_equal() +
        labs(x = NULL, y = NULL, title = "Trait correlation") +
        theme_pub(input$base_size, input$font %||% "sans") +
        theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(),
              axis.ticks = element_blank(), panel.grid.major.y = element_blank())
    })
    output$cor_plot <- renderPlot(cor_obj(), res = SCREEN_DPI)

    current_plot <- reactive(if (identical(input$tabs, "cor")) cor_obj() else plot_obj())
    download_server("dl", current_plot, function() if (identical(input$tabs, "cor")) "trait_correlation" else "phenotype_groups")

    output$tbl_desc <- renderDT({
      tl <- trait_list()
      t <- do.call(rbind, lapply(tl, function(x) cbind(性状 = x$label, group_summary(x$d$y, x$d$g))))
      dt_table(t, page_length = 16)
    })
    output$tbl_omni <- renderDT({
      t <- do.call(rbind, lapply(tl <- trait_list(), function(x) cbind(性状 = x$label, x$omni)))
      names(t) <- c("性状", "检验", "统计量", "自由度", "p")
      dt_table(t, p_cols = "p", page_length = 12)
    })
    output$tbl_pw <- renderDT({
      t <- do.call(rbind, lapply(trait_list(), function(x) {
        cbind(性状 = x$label, x$pw[, c("group1", "group2", "n1", "n2", "diff", "p", "p.adj", "signif", "method")])
      }))
      names(t) <- c("性状", "组 1", "组 2", "n1", "n2", "差值 (2−1)", "p", "p.adj", "显著性", "方法")
      dt_table(t, p_cols = c("p", "p.adj"), page_length = 12)
    })
    output$tbl_data <- renderDT(dt_table(filtered(), page_length = 12))
  })
}
