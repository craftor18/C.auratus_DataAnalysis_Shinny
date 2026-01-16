# R/mod_phenotype.R

phenotype_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "表型数据可视化",
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("load_sample3"), "Load Sample Data", class = "btn btn-info btn-block"),
        fileInput(ns("file3"), "Upload CSV File", accept = c(".csv")),
        selectInput(ns("x_var3"), "Select Categorical Variable", choices = NULL),
        selectInput(ns("y_var3"), "Select Phenotype Variable(s)", choices = NULL, multiple = TRUE),
        selectInput(ns("stat_test3"), "Select Statistical Test",
                    choices = c("t-test" = "t.test", "ANOVA" = "anova", "Kruskal-Wallis" = "kruskal.test")),
        checkboxInput(ns("show_points3"), "Show Points", TRUE),
        checkboxInput(ns("show_lines3"), "显示显著性标记及分组连线", TRUE),
        checkboxInput(ns("show_ns3"), "显示不显著标记及分组连线", TRUE),
        selectInput(ns("selected_comparisons3"), "选择比较组对", choices = NULL, multiple = TRUE),

        h4("Customize Colors"),
        uiOutput(ns("color_selectors3")),

        h4("Font Settings"),
        numericInput(ns("font_size3"), "Font Size", value = 16, min = 8, max = 30),

        h4("Save Image Settings"),
        numericInput(ns("subplot_width3"), "Subplot Width (inch)", value = 6, min = 4, max = 20),
        numericInput(ns("subplot_height3"), "Subplot Height (inch)", value = 4, min = 3, max = 15),
        numericInput(ns("img_dpi3"), "Resolution (dpi)", value = 300, min = 72, max = 600),
        radioButtons(ns("img_format3"), "Format", choices = c("PNG", "PDF")),

        h4("Layout Settings"),
        numericInput(ns("layout_ncol3"), "Number of columns (0 for auto)", value = 0, min = 0, max = 10),
        numericInput(ns("layout_nrow3"), "Number of rows (0 for auto)", value = 0, min = 0, max = 10),

        downloadButton(ns("save_plot3"), "Download Image", class = "btn btn-success")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Data Table", DTOutput(ns("table3"))),
          tabPanel("Plot", plotOutput(ns("boxplot3")))
        )
      )
    )
  )
}

phenotype_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      rv3 <- reactiveValues(comparisons = NULL, data = NULL)

      observeEvent(input$file3, {
        req(input$file3)
        rv3$data <- read.csv(input$file3$datapath, stringsAsFactors = FALSE)
      })

      observeEvent(input$load_sample3, {
        rv3$data <- read.csv("data/phenotype_sample.csv", stringsAsFactors = FALSE)
      })

      data3 <- reactive({
        req(rv3$data)
        rv3$data
      })

      observeEvent(data3(), {
        df <- data3()
        cat_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
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
            inputId = ns(paste0("color3_", group)),
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
    }
  )
}
