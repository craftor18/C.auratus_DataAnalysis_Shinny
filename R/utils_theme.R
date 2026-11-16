app_theme <- function() {
  bs_theme(
    version = 5,
    preset = "shiny",
    primary = "#0E7490",
    secondary = "#64748B",
    success = "#059669",
    info = "#0284C7",
    warning = "#D97706",
    danger = "#DC2626",
    bg = "#FFFFFF",
    fg = "#0F172A",
    base_font = font_collection(
      "Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto",
      "PingFang SC", "Hiragino Sans GB", "Microsoft YaHei", "Noto Sans CJK SC", "sans-serif"
    ),
    heading_font = font_collection(
      "Inter", "system-ui", "-apple-system", "Segoe UI",
      "PingFang SC", "Microsoft YaHei", "Noto Sans CJK SC", "sans-serif"
    ),
    code_font = font_collection("JetBrains Mono", "SFMono-Regular", "Menlo", "Consolas", "monospace"),
    "border-radius" = "0.75rem",
    "card-border-radius" = "1rem",
    "navbar-padding-y" = "0.6rem"
  )
}

# Okabe-Ito is colour-blind safe; NPG / Lancet / JCO follow common journal styles.
PALETTES <- list(
  okabe_ito = c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999"),
  npg = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"),
  lancet = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"),
  jco = c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67", "#8F7700", "#3B3B3B"),
  aqua = c("#0E7490", "#F59E0B", "#7C3AED", "#10B981", "#EF4444", "#3B82F6", "#EC4899", "#84CC16")
)

PALETTE_CHOICES <- c(
  "Okabe-Ito（色盲友好）" = "okabe_ito",
  "Nature (NPG)" = "npg",
  "Lancet" = "lancet",
  "JCO" = "jco",
  "Aqua（应用主题）" = "aqua",
  "Viridis" = "viridis"
)

get_palette <- function(name, n, custom = NULL) {
  custom_cols <- parse_custom_colors(custom)
  if (length(custom_cols) >= n) return(custom_cols[seq_len(n)])
  if (identical(name, "viridis")) return(viridis_pal(end = 0.9)(n))
  base <- PALETTES[[name]] %||% PALETTES$okabe_ito
  if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

parse_custom_colors <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return(character(0))
  cols <- trimws(unlist(strsplit(x, "[,，;；\\s]+")))
  cols <- cols[nzchar(cols)]
  ok <- vapply(cols, function(cl) !inherits(try(grDevices::col2rgb(cl), silent = TRUE), "try-error"), logical(1))
  cols[ok]
}

LD_COLORS <- c(
  "0.8–1.0" = "#D43F3A",
  "0.6–0.8" = "#EEA236",
  "0.4–0.6" = "#5CB85C",
  "0.2–0.4" = "#46B8DA",
  "0.0–0.2" = "#357EBD",
  "NA" = "#B8B8B8"
)
LEAD_COLOR <- "#9632B8"

theme_pub <- function(base_size = 13, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(colour = "#475569", size = rel(0.85), margin = margin(b = 8)),
      plot.caption = element_text(colour = "#64748B", size = rel(0.7), hjust = 0),
      plot.title.position = "plot",
      axis.line = element_line(colour = "#1E293B", linewidth = 0.5),
      axis.ticks = element_line(colour = "#1E293B", linewidth = 0.5),
      axis.text = element_text(colour = "#1E293B"),
      axis.title = element_text(colour = "#0F172A"),
      legend.title = element_text(face = "bold", size = rel(0.85)),
      legend.text = element_text(size = rel(0.8)),
      legend.key.size = unit(0.9, "lines"),
      strip.background = element_rect(fill = "#F1F5F9", colour = NA),
      strip.text = element_text(face = "bold", colour = "#0F172A"),
      panel.grid.major.y = element_line(colour = "#EEF2F6", linewidth = 0.4),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(10, 14, 10, 10)
    )
}
