first_or <- function(x, default = NULL) {
  x <- x[!is.na(x)]
  if (length(x)) x[[1]] else default
}

genopheno_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330,
      sidebar_section(
        "数据", "database",
        data_source_ui(ns("src"), example_label = "示例：鲫鱼 4 群体 180 尾"),
        helpText("每行一个个体：基因型列 + 一个或多个数值表型列。基因型支持 0/1、0|1、A/G、AG。")
      ),
      sidebar_section(
        "变量", "sliders",
        selectInput(ns("geno_col"), "基因型列", choices = NULL),
        selectInput(ns("trait"), "表型（数量性状）", choices = NULL),
        selectizeInput(ns("covars"), "协变量（加性模型）", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "如群体、性别")),
        helpText("群体与性别会造成分层混杂，建议作为协变量纳入。")
      ),
      sidebar_section(
        "统计", "calculator",
        selectInput(ns("pw_method"), "两两比较", choices = PAIRWISE_CHOICES, selected = "welch"),
        selectInput(ns("p_adjust"), "多重检验校正", choices = P_ADJUST_CHOICES, selected = "BH")
      ),
      sidebar_section(
        "图形", "palette",
        radioButtons(ns("plot_type"), "样式", inline = TRUE,
                     choices = c("小提琴" = "violin", "箱线" = "box", "均值±SE" = "meanse")),
        checkboxInput(ns("show_points"), "显示个体散点", TRUE),
        checkboxInput(ns("show_pairs"), "显示两两比较括号", TRUE),
        checkboxInput(ns("only_signif"), "只标注显著（p.adj < 0.05）", FALSE),
        radioButtons(ns("label_style"), "显著性标签", inline = TRUE,
                     choices = c("p 值" = "p", "星号" = "stars")),
        selectInput(ns("palette"), "配色", choices = PALETTE_CHOICES),
        textInput(ns("custom_colors"), "自定义颜色（可选）", placeholder = "#0072B2, #E69F00, #009E73"),
        layout_columns(
          col_widths = c(7, 5),
          selectInput(ns("font"), "字体", choices = FONT_CHOICES),
          numericInput(ns("base_size"), "字号", 13, min = 8, max = 24)
        )
      ),
      sidebar_section("导出", "file-export", download_ui(ns("dl"), width = 5, height = 5))
    ),
    page_header("STEP 3", "基因型–表型关联", "验证候选 SNP 对数量性状的等位基因效应：加性模型、HWE、两两比较。", "dna"),
    uiOutput(ns("boxes")),
    navset_card_underline(
      full_screen = TRUE,
      nav_panel("分布图", icon = icon("chart-column"), plotOutput(ns("plot"), height = "540px")),
      nav_panel(
        "关联检验", icon = icon("square-root-variable"),
        h6(class = "table-title", "加性模型（allele dosage 回归）"),
        DTOutput(ns("tbl_additive")),
        h6(class = "table-title", "总体检验（基因型作为分类变量）"),
        DTOutput(ns("tbl_omni")),
        h6(class = "table-title", "两两比较"),
        DTOutput(ns("tbl_pw"))
      ),
      nav_panel(
        "基因型频率", icon = icon("chart-pie"),
        layout_columns(
          col_widths = c(5, 7),
          DTOutput(ns("tbl_freq")),
          plotOutput(ns("freq_plot"), height = "320px")
        )
      ),
      nav_panel("数据", icon = icon("table"), DTOutput(ns("tbl_data")))
    ),
    method_card(
      p(strong("等位基因计数："), "相位基因型 0|1 与 1|0 表示同一个杂合子，统一转成无相位的 0/1；",
        "按等位基因频率确定参考 (REF) 与替代 (ALT) 等位基因，剂量 = ALT 拷贝数 (0/1/2)。"),
      p(strong("加性模型："), "表型 ~ 剂量 + 协变量。β 是每多一个 ALT 拷贝的表型变化量，",
        "偏 R² 是 SNP 在协变量之外额外解释的方差比例。这与 GEMMA、PLINK 等 GWAS 软件使用的模型一致。"),
      p(strong("Hardy–Weinberg 平衡："), "使用 Wigginton 等 (2005) 的精确检验。在对照群体中显著偏离 (p < 1e-6)",
        "通常提示分型错误；但多群体混合样本本身也会偏离 HWE（Wahlund 效应）。"),
      p(strong("两两比较："), "Tukey HSD 自带多重校正；Welch t、Dunn、Wilcoxon 的 p 值按所选方法校正。",
        "某基因型样本很少（如 n < 5）时结果不稳定，应谨慎解读。")
    )
  )
}

genopheno_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    raw <- data_source_server("src", "genotype_phenotype")

    observeEvent(raw(), {
      df <- raw()
      chars <- names(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
      gt_guess <- first_or(chars[grepl("geno|^gt$|snp", chars, ignore.case = TRUE)], first_or(chars))
      updateSelectInput(session, "geno_col", choices = chars, selected = gt_guess)
      nums <- numeric_columns(df)
      nums <- nums[!grepl("^(iid|id|sample|fid|no)$", nums, ignore.case = TRUE)]
      updateSelectInput(session, "trait", choices = nums, selected = first_or(nums))
      cov_choices <- setdiff(c(categorical_columns(df), nums), c(gt_guess, "Sample", "IID"))
      updateSelectizeInput(session, "covars", choices = cov_choices,
                           selected = intersect(c("FID", "Population", "Gender", "Sex"), cov_choices))
    })

    an <- reactive({
      df <- raw()
      req(input$geno_col %in% names(df), input$trait %in% names(df))
      validate(need(is.numeric(df[[input$trait]]), "所选表型列不是数值型。"))
      gi <- tryCatch(genotype_info(df[[input$geno_col]]),
                     error = function(e) validate(need(FALSE, conditionMessage(e))))
      y <- df[[input$trait]]
      covs <- setdiff(input$covars, c(input$trait, input$geno_col))
      cov_df <- if (length(covs)) df[, covs, drop = FALSE] else NULL
      d <- data.frame(genotype = gi$genotype, dosage = gi$dosage, y = y)
      d <- d[!is.na(d$genotype) & !is.na(d$y), ]
      d$genotype <- droplevels(d$genotype)
      validate(need(nlevels(d$genotype) >= 2, "该位点在当前样本中只有一种基因型，无法比较。"))
      list(
        gi = gi,
        d = d,
        additive = rbind(assoc_additive(y, gi$dosage, NULL), assoc_additive(y, gi$dosage, cov_df)),
        omni = omnibus_tests(d$y, d$genotype),
        pw = pairwise_tests(d$y, d$genotype, input$pw_method, input$p_adjust)
      )
    })

    output$boxes <- renderUI({
      a <- an()
      gi <- a$gi
      add <- a$additive[nrow(a$additive), ]
      hwe_theme <- if (!is.na(gi$hwe_p) && gi$hwe_p < 1e-6) "danger" else "success"
      layout_columns(
        fill = FALSE,
        class = "stat-row",
        stat_box("有效样本", sprintf("%d", gi$n_called), "users", "primary",
                 sprintf("缺失 %d · REF %s / ALT %s", gi$n_missing, gi$ref, gi$alt)),
        stat_box("次等位基因频率 MAF", sprintf("%.3f", gi$maf), "percent", "info",
                 sprintf("次等位基因：%s", gi$minor_allele)),
        stat_box("HWE 精确检验", fmt_p(gi$hwe_p), "scale-balanced", hwe_theme,
                 sprintf("%d / %d / %d（纯合 REF / 杂合 / 纯合 ALT）", gi$counts[1], gi$counts[2], gi$counts[3])),
        stat_box("加性效应 β", if (is.null(add) || !nrow(add)) "NA" else sprintf("%+.3g", add$beta), "arrow-trend-up",
                 if (!is.null(add) && nrow(add) && add$p < 0.05) "warning" else "secondary",
                 if (is.null(add) || !nrow(add)) "无法拟合" else sprintf("p = %s · 偏 R² = %.1f%%", fmt_p(add$p), 100 * add$r2))
      )
    })

    plot_obj <- reactive({
      a <- an()
      d <- a$d
      lv <- levels(d$genotype)
      cols <- setNames(get_palette(input$palette, length(lv), input$custom_colors), lv)
      n_lab <- table(d$genotype)
      fam <- input$font %||% "sans"

      p <- ggplot(d, aes(genotype, y, fill = genotype, colour = genotype))
      if (input$plot_type == "violin") {
        p <- p +
          geom_violin(alpha = 0.28, colour = NA, trim = TRUE, scale = "width", width = 0.8) +
          geom_boxplot(width = 0.16, fill = "white", colour = "#1E293B", outlier.shape = NA, linewidth = 0.45)
      } else if (input$plot_type == "box") {
        p <- p + geom_boxplot(width = 0.55, alpha = 0.35, colour = "#1E293B", outlier.shape = NA, linewidth = 0.45)
      }
      if (isTRUE(input$show_points)) {
        p <- p + geom_jitter(width = 0.12, height = 0, size = 1.7, alpha = 0.6, shape = 16)
      }
      if (input$plot_type == "meanse") {
        p <- p +
          stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.14, colour = "#1E293B", linewidth = 0.6) +
          stat_summary(fun = mean, geom = "point", size = 4, shape = 21, colour = "#1E293B")
      } else {
        p <- p + stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", colour = "#1E293B")
      }

      if (isTRUE(input$show_pairs) && !is.null(a$pw)) {
        bp <- bracket_positions(a$pw, d$y, only_signif = isTRUE(input$only_signif))
        if (!is.null(bp) && nrow(bp)) {
          bp$label <- if (input$label_style == "stars") bp$signif else paste0("p = ", fmt_p(bp$p.adj))
          p <- p + ggpubr::stat_pvalue_manual(bp, label = "label", inherit.aes = FALSE, tip.length = 0.01,
                                              size = input$base_size / 3.6, bracket.size = 0.4, family = fam)
        }
      }

      add <- a$additive[nrow(a$additive), ]
      subtitle <- if (!is.null(add) && nrow(add)) {
        sprintf("Additive model: beta = %+.3g ± %.2g (SE), p = %s, partial R² = %.1f%%  |  MAF = %.3f, HWE p = %s",
                add$beta, add$se, fmt_p(add$p), 100 * add$r2, a$gi$maf, fmt_p(a$gi$hwe_p))
      } else NULL

      p +
        scale_fill_manual(values = cols, guide = "none") +
        scale_colour_manual(values = cols, guide = "none") +
        scale_x_discrete(labels = function(l) sprintf("%s\nn = %d", l, as.integer(n_lab[l]))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
        labs(x = sprintf("Genotype (%s)", input$geno_col), y = input$trait,
             title = sprintf("%s by genotype", input$trait), subtitle = subtitle,
             caption = sprintf("Pairwise: %s; adjustment: %s", unique(a$pw$method)[1], input$p_adjust)) +
        theme_pub(input$base_size, fam)
    })

    output$plot <- renderPlot(plot_obj(), res = SCREEN_DPI)
    download_server("dl", plot_obj, function() paste0("genotype_", input$trait))

    output$tbl_additive <- renderDT({
      t <- an()$additive
      validate(need(!is.null(t), "无法拟合加性模型。"))
      names(t) <- c("模型", "n", "β", "SE", "t", "p", "偏 R²")
      dt_table(t, p_cols = "p", page_length = 5, options = list(dom = "t"))
    })
    output$tbl_omni <- renderDT({
      t <- an()$omni
      names(t) <- c("检验", "统计量", "自由度", "p")
      dt_table(t, p_cols = "p", options = list(dom = "t"))
    })
    output$tbl_pw <- renderDT({
      t <- an()$pw
      t <- t[, c("group1", "group2", "n1", "n2", "diff", "statistic", "p", "p.adj", "signif", "method")]
      names(t) <- c("组 1", "组 2", "n1", "n2", "差值 (2−1)", "统计量", "p", "p.adj", "显著性", "方法")
      dt_table(t, p_cols = c("p", "p.adj"), options = list(dom = "t"))
    })

    freq_df <- reactive({
      gi <- an()$gi
      n <- gi$n_called
      q <- gi$alt_freq
      p <- 1 - q
      labs <- c(paste(gi$ref, gi$ref, sep = "/"),
                paste(sort(c(gi$ref, gi$alt)), collapse = "/"),
                paste(gi$alt, gi$alt, sep = "/"))
      data.frame(
        基因型 = labs,
        类型 = c("纯合 REF", "杂合", "纯合 ALT"),
        观测数 = as.integer(gi$counts),
        观测频率 = as.numeric(gi$counts) / n,
        HWE期望数 = n * c(p^2, 2 * p * q, q^2)
      )
    })
    output$tbl_freq <- renderDT(dt_table(freq_df(), options = list(dom = "t")))
    output$freq_plot <- renderPlot({
      f <- freq_df()
      fl <- tidyr::pivot_longer(f[, c("基因型", "观测数", "HWE期望数")], -基因型, names_to = "type", values_to = "n")
      fl$type <- factor(fl$type, c("观测数", "HWE期望数"), c("Observed", "Expected (HWE)"))
      fl$基因型 <- factor(fl$基因型, f$基因型)
      ggplot(fl, aes(基因型, n, fill = type)) +
        geom_col(position = position_dodge(0.75), width = 0.7) +
        geom_text(aes(label = round(n, 1)), position = position_dodge(0.75), vjust = -0.4, size = 3.4) +
        scale_fill_manual(values = c("#0E7490", "#CBD5E1"), name = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
        labs(x = NULL, y = "Count") +
        theme_pub(12) + theme(legend.position = "top")
    }, res = SCREEN_DPI)

    output$tbl_data <- renderDT(dt_table(raw(), page_length = 12))
  })
}
