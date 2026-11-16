read_table_auto <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  delim <- switch(ext, tsv = "\t", txt = "\t", csv = ",", NULL)
  if (is.null(delim)) {
    first <- readLines(path, n = 1, warn = FALSE)
    delim <- if (grepl("\t", first)) "\t" else ","
  }
  df <- readr::read_delim(path, delim = delim, show_col_types = FALSE, progress = FALSE,
                          locale = readr::locale(encoding = "UTF-8"), trim_ws = TRUE)
  as.data.frame(df, check.names = FALSE)
}

.example_cache <- new.env(parent = emptyenv())

load_example <- function(key) {
  if (is.null(.example_cache[[key]])) {
    .example_cache[[key]] <- read_table_auto(EXAMPLE_FILES[[key]])
  }
  .example_cache[[key]]
}

# Maps common column spellings to canonical names; returns NULL when a required one is missing.
match_columns <- function(df, spec) {
  low <- tolower(names(df))
  out <- list()
  for (canon in names(spec)) {
    hit <- which(low %in% tolower(spec[[canon]]))
    out[[canon]] <- if (length(hit)) names(df)[hit[1]] else NA_character_
  }
  unlist(out)
}

GWAS_COLUMNS <- list(
  SNP = c("SNP", "rs", "rsid", "id", "marker", "variant", "snp_id"),
  chr = c("chr", "chrom", "chromosome", "lg", "#chrom", "chr_id"),
  pos = c("pos", "bp", "position", "ps", "base_pair_location", "pos_bp"),
  P = c("p", "pval", "p_value", "pvalue", "p_wald", "p_lrt", "p.value"),
  R2 = c("r2", "ld", "rsq", "r^2", "ld_r2")
)

GENE_COLUMNS <- list(
  chr = c("chr", "chrom", "chromosome", "seqname", "lg"),
  start = c("start", "begin", "txstart", "gene_start"),
  stop = c("stop", "end", "txend", "gene_end"),
  name = c("name", "gene", "symbol", "gene_name", "gene_symbol"),
  ID = c("id", "gene_id", "locus_tag"),
  strand = c("strand"),
  comment = c("comment", "description", "annotation", "function")
)

clean_chr <- function(x) sub("^(chr|lg|chromosome)", "", tolower(trimws(as.character(x))))

categorical_columns <- function(df, max_levels = 50) {
  names(df)[vapply(df, function(x) {
    (is.character(x) || is.factor(x) || is.logical(x)) && length(unique(stats::na.omit(x))) <= max_levels
  }, logical(1))]
}

numeric_columns <- function(df) {
  names(df)[vapply(df, is.numeric, logical(1))]
}

# ---- data source module -----------------------------------------------------

data_source_ui <- function(id, accept = c(".csv", ".tsv", ".txt"), example_label = "使用示例数据",
                           upload_label = "上传文件") {
  ns <- NS(id)
  tagList(
    radioButtons(ns("source"), NULL,
                 choices = c(setNames("example", example_label), setNames("upload", upload_label)),
                 inline = TRUE),
    conditionalPanel(
      sprintf("input['%s'] == 'upload'", ns("source")),
      fileInput(ns("file"), NULL, accept = accept, buttonLabel = "选择文件", placeholder = "CSV / TSV")
    )
  )
}

data_source_server <- function(id, example_key) {
  moduleServer(id, function(input, output, session) {
    reactive({
      if (identical(input$source, "upload")) {
        req(input$file)
        tryCatch(
          read_table_auto(input$file$datapath, input$file$name),
          error = function(e) validate(need(FALSE, paste("文件读取失败：", conditionMessage(e))))
        )
      } else {
        load_example(example_key)
      }
    })
  })
}

# ---- figure export module ---------------------------------------------------

download_ui <- function(id, width = 7, height = 5) {
  ns <- NS(id)
  tagList(
    layout_columns(
      col_widths = c(6, 6),
      numericInput(ns("width"), "宽度 (inch)", value = width, min = 1, max = 30, step = 0.5),
      numericInput(ns("height"), "高度 (inch)", value = height, min = 1, max = 30, step = 0.5)
    ),
    layout_columns(
      col_widths = c(6, 6),
      selectInput(ns("format"), "格式", choices = c("PDF（矢量）" = "pdf", "PNG" = "png", "TIFF（LZW）" = "tiff")),
      numericInput(ns("dpi"), "DPI", value = 300, min = 72, max = 1200, step = 50)
    ),
    downloadButton(ns("download"), "导出图片", class = "btn-primary w-100", icon = icon("download"))
  )
}

download_server <- function(id, plot_r, basename) {
  moduleServer(id, function(input, output, session) {
    output$download <- downloadHandler(
      filename = function() {
        nm <- if (is.function(basename)) basename() else basename
        sprintf("%s_%s.%s", make.names(nm), format(Sys.Date(), "%Y%m%d"), input$format)
      },
      content = function(file) {
        p <- plot_r()
        req(p)
        dpi <- input$dpi
        showtext_opts(dpi = if (input$format == "pdf") 72 else dpi)
        on.exit(showtext_opts(dpi = SCREEN_DPI), add = TRUE)
        args <- list(filename = file, plot = p, width = input$width, height = input$height,
                     units = "in", dpi = dpi, bg = "white")
        if (input$format == "pdf") {
          args$device <- grDevices::cairo_pdf
        } else if (input$format == "tiff") {
          args$device <- "tiff"
          args$compression <- "lzw"
        } else {
          args$device <- "png"
        }
        do.call(ggplot2::ggsave, args)
      }
    )
  })
}

table_download_server <- function(id, data_r, basename) {
  moduleServer(id, function(input, output, session) {
    output$csv <- downloadHandler(
      filename = function() sprintf("%s_%s.csv", basename, format(Sys.Date(), "%Y%m%d")),
      content = function(file) readr::write_excel_csv(data_r(), file)
    )
  })
}

table_download_ui <- function(id, label = "导出 CSV") {
  downloadButton(NS(id, "csv"), label, class = "btn-outline-secondary btn-sm", icon = icon("table"))
}
