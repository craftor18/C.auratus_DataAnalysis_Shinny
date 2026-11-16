ld_bin <- function(r2) {
  b <- cut(r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), right = FALSE,
           labels = c("0.0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"))
  b <- as.character(b)
  b[is.na(b)] <- "NA"
  factor(b, levels = names(LD_COLORS))
}

# Greedy row packing so that gene bodies and their labels never overlap.
pack_genes <- function(g, span_bp) {
  if (!nrow(g)) return(g)
  g <- g[order(g$start), ]
  half_label <- nchar(g$label) * span_bp * 0.0065
  mid <- (g$start + g$stop) / 2
  eff_start <- pmin(g$start, mid - half_label)
  eff_end <- pmax(g$stop, mid + half_label)
  pad <- span_bp * 0.01
  ends <- numeric(0)
  g$row <- integer(nrow(g))
  for (i in seq_len(nrow(g))) {
    r <- which(ends < eff_start[i] - pad)[1]
    if (is.na(r)) {
      ends <- c(ends, eff_end[i])
      r <- length(ends)
    } else {
      ends[r] <- eff_end[i]
    }
    g$row[i] <- r
  }
  g
}

locuszoom_ui <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330,
      sidebar_section(
        "GWAS 结果", "file-lines",
        data_source_ui(ns("gwas"), example_label = "示例：LG34 区域"),
        helpText("列：SNP、chr、pos、P，可选 R2（与先导 SNP 的 LD r²）。")
      ),
      sidebar_section(
        "基因注释", "book-open",
        data_source_ui(ns("genes"), example_label = "示例：鲫鱼基因组注释"),
        helpText("列：chr、start、stop、name，可选 ID、strand、comment。")
      ),
      sidebar_section(
        "区域", "crosshairs",
        selectizeInput(ns("lead"), "先导 SNP", choices = NULL),
        helpText(class = "text-warning-emphasis",
                 "R² 必须是各 SNP 与所选先导 SNP 之间的 LD。更换先导 SNP 后，需要重新计算 LD。"),
        numericInput(ns("flank"), "两侧窗口 (kb，0 = 整条染色体)", value = 100, min = 0, step = 10),
        checkboxInput(ns("hide_unknown"), "隐藏未命名基因 (unknown)", FALSE)
      ),
      sidebar_section(
        "阈值与标注", "ruler-horizontal",
        layout_columns(
          col_widths = c(6, 6),
          numericInput(ns("gw_log"), "显著 −log₁₀P", value = 7.3, min = 0, max = 50, step = 0.1),
          numericInput(ns("sug_log"), "提示 −log₁₀P", value = 5, min = 0, max = 50, step = 0.1)
        ),
        sliderInput(ns("n_label"), "另外标注 Top SNP 数", min = 0, max = 15, value = 3, step = 1),
        layout_columns(
          col_widths = c(7, 5),
          selectInput(ns("font"), "字体", choices = FONT_CHOICES),
          numericInput(ns("base_size"), "字号", 13, min = 8, max = 24)
        ),
        sliderInput(ns("pt_size"), "点大小", min = 0.5, max = 4, value = 2, step = 0.25)
      ),
      sidebar_section("导出", "file-export", download_ui(ns("dl"), width = 8, height = 6))
    ),
    page_header("STEP 2", "LocusZoom 区域关联图", "在显著位点周围同时展示关联强度、LD 结构和基因，用于圈定候选基因。", "magnifying-glass-chart"),
    uiOutput(ns("boxes")),
    navset_card_underline(
      full_screen = TRUE,
      nav_panel("区域图", icon = icon("chart-line"), plotOutput(ns("plot"), height = "640px")),
      nav_panel("候选基因", icon = icon("dna"),
                p(class = "text-muted small mt-2", "按与先导 SNP 的距离排序。距离为 0 表示 SNP 落在基因内部。"),
                DTOutput(ns("tbl_genes"))),
      nav_panel("区域 SNP", icon = icon("list-ol"), DTOutput(ns("tbl_snps")))
    ),
    method_card(
      p(strong("−log₁₀(P)："), "纵轴是关联显著性。图中不截断 P 值，最显著位点按真实高度显示。"),
      p(strong("LD 着色："), "颜色表示每个 SNP 与先导 SNP（紫色菱形）的 r²，分级配色与 LocusZoom 一致。",
        "与先导 SNP 高度连锁 (r² > 0.8) 的信号通常代表同一个关联，不应看作独立位点。"),
      p(strong("阈值："), "5×10⁻⁸ 是人类 GWAS 常用的全基因组阈值。在鱼类等 LD 结构不同的物种中，",
        "也常按有效标记数做 Bonferroni 校正 (0.05 / m)，可在左侧修改。"),
      p(strong("候选基因："), "先看覆盖先导 SNP 的基因，再看高 LD 区间内的基因；还需结合表达（STEP 4）和功能注释综合判断。")
    )
  )
}

locuszoom_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    gwas_raw <- data_source_server("gwas", "locuszoom")
    gw_p <- reactive(10^-(input$gw_log %||% 7.3))
    sug_p <- reactive(10^-(input$sug_log %||% 5))
    genes_raw <- data_source_server("genes", "gene_annotation")

    gwas <- reactive({
      df <- gwas_raw()
      m <- match_columns(df, GWAS_COLUMNS)
      validate(need(!anyNA(m[c("chr", "pos", "P")]), "GWAS 文件需要 chr、pos、P 列（SNP、R2 可选）。"))
      out <- data.frame(
        SNP = if (!is.na(m[["SNP"]])) as.character(df[[m[["SNP"]]]]) else paste(df[[m[["chr"]]]], df[[m[["pos"]]]], sep = ":"),
        chr = clean_chr(df[[m[["chr"]]]]),
        pos = suppressWarnings(as.numeric(df[[m[["pos"]]]])),
        P = suppressWarnings(as.numeric(df[[m[["P"]]]])),
        R2 = if (!is.na(m[["R2"]])) suppressWarnings(as.numeric(df[[m[["R2"]]]])) else NA_real_,
        stringsAsFactors = FALSE
      )
      out <- out[is.finite(out$pos) & is.finite(out$P) & out$P > 0 & out$P <= 1, ]
      validate(need(nrow(out) > 0, "没有有效的 SNP 记录（P 必须在 (0, 1] 之间）。"))
      out$logp <- -log10(out$P)
      out[order(out$P), ]
    })

    observeEvent(gwas(), {
      g <- gwas()
      updateSelectizeInput(session, "lead", choices = g$SNP, selected = g$SNP[1], server = TRUE)
    })

    genes <- reactive({
      df <- genes_raw()
      m <- match_columns(df, GENE_COLUMNS)
      validate(need(!anyNA(m[c("chr", "start", "stop", "name")]), "注释文件需要 chr、start、stop、name 列。"))
      data.frame(
        chr = clean_chr(df[[m[["chr"]]]]),
        start = as.numeric(df[[m[["start"]]]]),
        stop = as.numeric(df[[m[["stop"]]]]),
        name = as.character(df[[m[["name"]]]]),
        ID = if (!is.na(m[["ID"]])) as.character(df[[m[["ID"]]]]) else NA_character_,
        strand = if (!is.na(m[["strand"]])) as.character(df[[m[["strand"]]]]) else NA_character_,
        comment = if (!is.na(m[["comment"]])) as.character(df[[m[["comment"]]]]) else NA_character_,
        stringsAsFactors = FALSE
      )
    })

    region <- reactive({
      g <- gwas()
      req(input$lead)
      lead <- g[g$SNP == input$lead, ][1, ]
      validate(need(!is.na(lead$SNP), "所选先导 SNP 不在 GWAS 文件中。"))
      flank <- (input$flank %||% 0) * 1000
      w <- g[g$chr == lead$chr, ]
      if (flank > 0) w <- w[w$pos >= lead$pos - flank & w$pos <= lead$pos + flank, ]
      xlim <- if (flank > 0) c(lead$pos - flank, lead$pos + flank) else range(w$pos)
      w$ld <- ld_bin(w$R2)
      w$is_lead <- w$SNP == lead$SNP

      gn <- genes()
      gn <- gn[gn$chr == lead$chr & gn$stop >= xlim[1] & gn$start <= xlim[2], ]
      if (isTRUE(input$hide_unknown)) gn <- gn[!tolower(gn$name) %in% c("unknown", "na", ""), ]
      gn$label <- ifelse(tolower(gn$name) %in% c("unknown", "na", "") & !is.na(gn$ID), gn$ID, gn$name)
      gn$distance_kb <- ifelse(gn$start <= lead$pos & gn$stop >= lead$pos, 0,
                               pmin(abs(gn$start - lead$pos), abs(gn$stop - lead$pos)) / 1000)
      gn$relation <- ifelse(gn$distance_kb == 0, "包含先导 SNP",
                            ifelse(gn$stop < lead$pos, "左侧（坐标较小）", "右侧（坐标较大）"))
      gn <- pack_genes(gn, diff(xlim))
      list(lead = lead, w = w, genes = gn, xlim = xlim)
    })

    output$boxes <- renderUI({
      r <- region()
      n_gw <- sum(r$w$P < gw_p(), na.rm = TRUE)
      n_ld <- sum(r$w$R2 >= 0.8 & !r$w$is_lead, na.rm = TRUE)
      host <- r$genes$label[r$genes$distance_kb == 0]
      layout_columns(
        fill = FALSE,
        class = "stat-row",
        stat_box("先导 SNP", r$lead$SNP, "location-crosshairs", "primary",
                 sprintf("LG%s : %s bp", r$lead$chr, format(r$lead$pos, big.mark = ","))),
        stat_box("先导 P 值", fmt_p(r$lead$P), "bolt", if (r$lead$P < gw_p()) "danger" else "warning",
                 sprintf("−log₁₀P = %.2f", r$lead$logp)),
        stat_box("窗口内 SNP", format(nrow(r$w), big.mark = ","), "braille", "info",
                 sprintf("%d 个达全基因组显著 · %d 个 r² ≥ 0.8", n_gw, n_ld)),
        stat_box("窗口内基因", nrow(r$genes), "dna", "success",
                 if (length(host)) paste("SNP 位于", paste(host, collapse = ", "), "内") else "SNP 位于基因间区")
      )
    })

    plot_obj <- reactive({
      r <- region()
      fam <- input$font %||% "sans"
      bs <- input$base_size %||% 13
      w <- r$w[order(r$w$is_lead, ifelse(is.na(r$w$R2), -1, r$w$R2)), ]
      others <- w[!w$is_lead, ]
      lead <- w[w$is_lead, ]
      to_label <- head(others[order(others$P), ], input$n_label)
      xl <- r$xlim / 1e6
      ymax <- max(w$logp, -log10(gw_p()), na.rm = TRUE) * 1.12

      top <- ggplot(w, aes(pos / 1e6, logp)) +
        geom_vline(xintercept = r$lead$pos / 1e6, colour = LEAD_COLOR, linetype = "dotted", linewidth = 0.4) +
        geom_hline(yintercept = -log10(gw_p()), colour = "#B91C1C", linetype = "dashed", linewidth = 0.45) +
        geom_hline(yintercept = -log10(sug_p()), colour = "#64748B", linetype = "dotted", linewidth = 0.45) +
        geom_point(data = others, aes(fill = ld), shape = 21, colour = "#334155", stroke = 0.15,
                   size = input$pt_size, alpha = 0.92) +
        geom_point(data = lead, shape = 23, fill = LEAD_COLOR, colour = "black", size = input$pt_size * 2.3, stroke = 0.6) +
        ggrepel::geom_label_repel(
          data = lead, aes(label = sprintf("%s\nP = %s", SNP, fmt_p(P))),
          family = fam, size = bs / 3.9, label.size = 0.2, fill = "white", nudge_y = ymax * 0.05,
          box.padding = 0.5, segment.colour = "#475569", min.segment.length = 0
        ) +
        ggrepel::geom_text_repel(
          data = to_label, aes(label = SNP), family = fam, size = bs / 4.6, colour = "#334155",
          box.padding = 0.35, segment.colour = "#94A3B8", min.segment.length = 0, max.overlaps = Inf
        ) +
        scale_fill_manual(values = LD_COLORS, drop = FALSE, name = expression(r^2),
                          breaks = if (anyNA(w$R2)) names(LD_COLORS) else setdiff(names(LD_COLORS), "NA"),
                          guide = guide_legend(override.aes = list(size = 4), reverse = FALSE)) +
        scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
        coord_cartesian(xlim = xl) +
        labs(x = NULL, y = expression(-log[10](italic(P))),
             title = sprintf("Regional association on LG%s", r$lead$chr),
             subtitle = sprintf("Lead SNP %s  |  window %s–%s Mb  |  dashed: -log10(P) = %s, dotted: -log10(P) = %s",
                                r$lead$SNP, format(round(xl[1], 3), nsmall = 3), format(round(xl[2], 3), nsmall = 3),
                                input$gw_log, input$sug_log)) +
        theme_pub(bs, fam) +
        theme(
          legend.position = "inside",
          legend.position.inside = c(0.995, 0.99),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = "#E2E8F0"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()
        )

      gn <- r$genes
      n_rows <- max(1, if (nrow(gn)) max(gn$row) else 1)
      bottom <- ggplot() +
        geom_vline(xintercept = r$lead$pos / 1e6, colour = LEAD_COLOR, linetype = "dotted", linewidth = 0.4)
      if (nrow(gn)) {
        gn$host <- gn$distance_kb == 0
        bottom <- bottom +
          geom_segment(data = gn, aes(x = start / 1e6, xend = stop / 1e6, y = -row, yend = -row),
                       colour = "#94A3B8", linewidth = 0.4) +
          geom_rect(data = gn, aes(xmin = start / 1e6, xmax = stop / 1e6, ymin = -row - 0.16, ymax = -row + 0.16, fill = host),
                    colour = NA) +
          geom_text(data = gn, aes(x = (start + stop) / 2e6, y = -row - 0.42, label = label, colour = host),
                    family = fam, fontface = "italic", size = bs / 4.4, vjust = 1)
        if (!all(is.na(gn$strand))) {
          st <- gn[gn$strand %in% c("+", "-"), ]
          bottom <- bottom + geom_segment(
            data = st,
            aes(x = ifelse(strand == "+", start, stop) / 1e6, xend = ifelse(strand == "+", stop, start) / 1e6,
                y = -row, yend = -row),
            arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"), colour = "#334155", linewidth = 0.3
          )
        }
      } else {
        bottom <- bottom + annotate("text", x = mean(xl), y = -1, label = "No annotated genes in window",
                                    colour = "#64748B", family = fam)
      }
      bottom <- bottom +
        scale_fill_manual(values = c(`TRUE` = LEAD_COLOR, `FALSE` = "#0E7490"), guide = "none") +
        scale_colour_manual(values = c(`TRUE` = LEAD_COLOR, `FALSE` = "#1E293B"), guide = "none") +
        scale_y_continuous(limits = c(-n_rows - 0.9, -0.5), expand = expansion(0)) +
        scale_x_continuous(labels = label_number(accuracy = 0.01)) +
        coord_cartesian(xlim = xl) +
        labs(x = sprintf("Position on LG%s (Mb)", r$lead$chr), y = "Genes") +
        theme_pub(bs, fam) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
              panel.grid.major.y = element_blank())

      top / bottom + plot_layout(heights = c(3.2, max(0.9, n_rows * 0.55)))
    })

    output$plot <- renderPlot(plot_obj(), res = SCREEN_DPI)
    download_server("dl", plot_obj, function() paste0("locuszoom_", input$lead))

    output$tbl_genes <- renderDT({
      gn <- region()$genes
      validate(need(nrow(gn) > 0, "窗口内没有注释基因。"))
      gn <- gn[order(gn$distance_kb), ]
      t <- data.frame(基因 = gn$name, ID = gn$ID, 起始 = gn$start, 终止 = gn$stop,
                      `长度 (kb)` = (gn$stop - gn$start) / 1000, `距先导 SNP (kb)` = gn$distance_kb,
                      位置关系 = gn$relation, 功能注释 = gn$comment, check.names = FALSE)
      t$起始 <- as.integer(t$起始)
      t$终止 <- as.integer(t$终止)
      dt_table(t, num_cols = c("长度 (kb)", "距先导 SNP (kb)"), digits = 1)
    })

    output$tbl_snps <- renderDT({
      w <- region()$w
      w <- w[order(w$P), c("SNP", "chr", "pos", "P", "logp", "R2")]
      w$pos <- as.integer(w$pos)
      names(w) <- c("SNP", "染色体", "位置", "P", "−log10P", "r²")
      dt_table(w, p_cols = "P", num_cols = c("−log10P", "r²"))
    })
  })
}
