page_header <- function(step, title, subtitle, icon_name) {
  div(
    class = "page-header",
    div(class = "page-header-icon", icon(icon_name)),
    div(
      span(class = "page-step", step),
      h2(class = "page-title", title),
      p(class = "page-subtitle", subtitle)
    )
  )
}

method_card <- function(..., title = "方法与原理") {
  accordion(
    open = FALSE,
    class = "method-accordion",
    accordion_panel(title, icon = icon("flask"), div(class = "method-body", ...))
  )
}

stat_box <- function(title, value, icon_name, theme = "primary", subtitle = NULL) {
  value_box(
    title = title,
    value = value,
    showcase = icon(icon_name),
    theme = theme,
    p(class = "small mb-0", subtitle),
    class = "stat-box"
  )
}

plot_card <- function(..., full_screen = TRUE) {
  card(full_screen = full_screen, class = "plot-card", card_body(class = "p-2", ...))
}

dt_table <- function(df, p_cols = character(0), num_cols = NULL, digits = 3, page_length = 10, ...) {
  if (is.null(num_cols)) num_cols <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], p_cols)
  num_cols <- setdiff(num_cols, names(df)[vapply(df, is.integer, logical(1))])
  dt <- DT::datatable(
    df,
    rownames = FALSE,
    class = "compact hover row-border",
    extensions = "Buttons",
    options = list(
      pageLength = page_length,
      scrollX = TRUE,
      dom = "Bfrtip",
      buttons = list("copy", list(extend = "csv", charset = "utf-8", bom = TRUE)),
      language = list(
        search = "搜索：",
        info = "第 _START_–_END_ 条，共 _TOTAL_ 条",
        infoEmpty = "无数据",
        zeroRecords = "没有匹配的记录",
        paginate = list(previous = "上一页", `next` = "下一页")
      )
    ),
    ...
  )
  if (length(num_cols)) dt <- DT::formatRound(dt, num_cols, digits)
  if (length(p_cols)) dt <- DT::formatSignif(dt, p_cols, 3)
  dt
}

sidebar_section <- function(title, icon_name, ...) {
  div(class = "sidebar-section", div(class = "sidebar-section-title", icon(icon_name), title), ...)
}
