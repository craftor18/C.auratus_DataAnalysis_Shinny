helix <- function(n = 16) {
  div(
    class = "helix", `aria-hidden` = "true",
    lapply(seq_len(n), function(i) {
      div(class = "helix-row", style = sprintf("--i:%d", i),
          span(class = "helix-dot a"), span(class = "helix-bar"), span(class = "helix-dot b"))
    })
  )
}

counter <- function(value, label, suffix = "") {
  div(class = "hero-stat reveal",
      div(class = "hero-stat-value", span(class = "count-up", `data-target` = value, "0"), suffix),
      div(class = "hero-stat-label", label))
}

step_card <- function(step, target, icon_name, title, text, tags_vec) {
  tags$a(
    class = "step-card reveal", href = "#", `data-nav` = target,
    div(class = "step-card-top", span(class = "step-num", step), span(class = "step-icon", icon(icon_name))),
    h4(title),
    p(text),
    div(class = "step-tags", lapply(tags_vec, function(t) span(class = "chip", t))),
    span(class = "step-go", "打开", icon("arrow-right"))
  )
}

home_ui <- function(id) {
  n_snp <- nrow(load_example("locuszoom"))
  n_gene <- nrow(load_example("gene_annotation"))
  n_expr <- nrow(load_example("expression"))
  n_fish <- nrow(load_example("population_phenotype")) + nrow(load_example("genotype_phenotype"))

  div(
    class = "home",
    tags$section(
      class = "hero",
      div(class = "hero-bg", div(class = "blob b1"), div(class = "blob b2"), div(class = "blob b3")),
      div(
        class = "hero-inner",
        div(
          class = "hero-copy",
          span(class = "eyebrow", icon("fish"), " Carassius auratus · GWAS 下游分析"),
          h1(class = "hero-title", "从关联信号", tags$br(), span(class = "grad-text", "到候选基因")),
          p(class = "hero-lead",
            "面向鲫鱼群体遗传研究的交互式工作台：表型描述、区域关联定位、等位基因效应验证、候选基因表达，",
            "按 GWAS 的分析逻辑串成一条完整流程。"),
          div(
            class = "hero-actions",
            tags$a(class = "btn btn-primary btn-lg", href = "#", `data-nav` = "pheno", icon("play"), " 开始分析"),
            tags$a(class = "btn btn-outline-secondary btn-lg", href = "#", `data-nav` = "guide", icon("book"), " 阅读指南")
          )
        ),
        div(class = "hero-visual", helix())
      ),
      div(
        class = "hero-stats",
        counter(n_fish, "尾个体表型"),
        counter(n_snp, "个区域 SNP"),
        counter(n_gene, "条基因注释"),
        counter(n_expr, "个表达基因")
      )
    ),
    tags$section(
      class = "home-section",
      div(class = "section-head reveal",
          span(class = "eyebrow", "WORKFLOW"),
          h2("四步分析流程"),
          p("每一步的输出都是下一步的输入：先看表型差异，再定位关联区域，然后验证 SNP 效应，最后看候选基因的表达。")),
      div(
        class = "step-grid",
        step_card("01", "pheno", "fish", "群体表型比较",
                  "比较群体或处理组之间的形态性状，支持异速校正、子集筛选和多重比较校正。",
                  c("ANOVA", "Kruskal–Wallis", "Tukey / Dunn")),
        step_card("02", "locus", "magnifying-glass-chart", "LocusZoom 区域定位",
                  "以先导 SNP 为中心，展示 −log₁₀P、LD r² 与基因结构，并按距离列出候选基因。",
                  c("LD r²", "全基因组阈值", "候选基因表")),
        step_card("03", "geno", "dna", "基因型–表型关联",
                  "统一相位基因型、计算 MAF 与 HWE，用加性模型估计等位基因效应 β。",
                  c("加性模型", "HWE 精确检验", "协变量")),
        step_card("04", "expr", "chart-simple", "候选基因表达",
                  "在 log₂ 尺度比较候选基因的组间表达，跨基因 FDR 校正，并绘制 Z-score 热图。",
                  c("log₂FC", "BH-FDR", "热图"))
      )
    ),
    tags$section(
      class = "home-section",
      div(class = "section-head reveal",
          span(class = "eyebrow", "PRINCIPLES"),
          h2("统计上站得住，图可以直接投稿")),
      div(
        class = "feature-grid",
        div(class = "feature reveal", icon("shield-halved"),
            h5("统计严谨"), p("每张图都标明所用检验和校正方法；默认不剔除异常值，不截断 P 值。")),
        div(class = "feature reveal", icon("dna"),
            h5("遗传学正确"), p("0|1 与 1|0 合并为同一杂合子；加性剂量编码；LD 着色与 LocusZoom 一致。")),
        div(class = "feature reveal", icon("eye"),
            h5("色盲友好"), p("默认 Okabe-Ito 配色，另有 NPG、Lancet、JCO 等期刊风格可选。")),
        div(class = "feature reveal", icon("file-export"),
            h5("出版级导出"), p("PDF 矢量、PNG 或 TIFF (LZW)，尺寸和 DPI 可调，表格可导出 CSV。"))
      )
    )
  )
}

guide_ui <- function(id) {
  div(
    class = "guide-wrap",
    page_header("GUIDE", "使用指南", "数据格式、分析原理与常见问题。", "book-open"),
    card(class = "guide-card", card_body(includeMarkdown("docs/guide.md")))
  )
}
