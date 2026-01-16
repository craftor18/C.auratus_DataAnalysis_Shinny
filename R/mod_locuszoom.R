# R/mod_locuszoom.R

locuszoom_ui <- function(id) {
  ns <- NS(id)

  tabPanel(
    "LocusZoom绘图",
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("load_sample4"), "Load Sample Data", class = "btn btn-info btn-block"),
        fileInput(ns("gwas_file4"), "Upload GWAS File (e.g., 47_16447785.csv):", accept = c(".csv")),
        fileInput(ns("gene_file4"), "Upload Gene Track File (e.g., genes.up.csv):", accept = c(".csv")),
        textInput(ns("lead_snp4"), "Enter Lead SNP ID (e.g., '47_16447785'):",
                  placeholder = "Enter SNP ID here"),
        hr(),
        h4("Advanced Options"),
        numericInput(ns("legend_square_size4"), "Legend Square Size", value = 15, min = 1),
        numericInput(ns("lead_snp_size4"), "Lead SNP Diamond Size", value = 5, min = 1),
        numericInput(ns("display_width4"), "Display Plot Width (inches)", value = 10, min = 1),
        numericInput(ns("display_height4"), "Display Plot Height (inches)", value = 6, min = 1),
        numericInput(ns("display_dpi4"), "Display Plot DPI", value = 300, min = 72),
        numericInput(ns("save_width4"), "Save Image Width (inches)", value = 3.2, min = 0.1),
        numericInput(ns("save_height4"), "Save Image Height (inches)", value = 1.8, min = 0.1),
        numericInput(ns("save_dpi4"), "Save Image DPI", value = 300, min = 72),
        textInput(ns("font_family4"), "Font Family", value = "sans"),
        actionButton(ns("update_plot4"), "Generate Plot", class = "btn btn-primary btn-block")
      ),
      mainPanel(
        plotOutput(ns("combined_plot4"), width = "100%", height = "800px"),
        downloadButton(ns("download_pdf4"), "Download PDF", class = "btn btn-success"),
        downloadButton(ns("download_png4"), "Download PNG", class = "btn btn-success")
      )
    )
  )
}

locuszoom_server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      gwas_data4 <- reactiveVal()
      gene_data4 <- reactiveVal()

      observeEvent(input$gwas_file4, {
        req(input$gwas_file4)
        gwas_data4(read.csv(input$gwas_file4$datapath))
      })

      observeEvent(input$gene_file4, {
        req(input$gene_file4)
        gene_data4(read.csv(input$gene_file4$datapath))
      })

      observeEvent(input$load_sample4, {
        gwas_data4(read.csv("data/locuszoom_gwas_sample.csv"))
        gene_data4(read.csv("data/locuszoom_genes_sample.csv"))
        updateTextInput(session, "lead_snp4", value = "47_16447785")
      })

      plot_reactive <- eventReactive(input$update_plot4, {
        req(gwas_data4(), gene_data4(), input$lead_snp4)
        gwas_data <- gwas_data4()
        gene_data <- gene_data4()

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
          geom_point(
            data = highlight_snp,
            aes(x = pos, y = -log10(P)),
            shape = 23, size = input$lead_snp_size4, fill = "red", color = "black", stroke = 1
          ) +
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

        gwas_plot / gene_track_plot + plot_layout(heights = c(3, 1))
      })

      output$combined_plot4 <- renderPlot({
        plot_reactive()
      },
      width = function() { input$display_width4 * input$display_dpi4 },
      height = function() { input$display_height4 * input$display_dpi4 }
      )

      output$download_pdf4 <- downloadHandler(
        filename = function() { "LocusZoom_Plot.pdf" },
        content = function(file) {
          ggsave(
            filename = file,
            plot = plot_reactive(),
            width = input$save_width4  * 2,
            height = input$save_height4 * 2,
            scale = 0.5,
            dpi = input$save_dpi4,
            device = cairo_pdf,
            family = ifelse(input$font_family4 == "", "sans", input$font_family4)
          )
        }
      )

      output$download_png4 <- downloadHandler(
        filename = function() { "LocusZoom_Plot.png" },
        content = function(file) {
          ggsave(
            filename = file,
            plot = plot_reactive(),
            width = input$save_width4  * 2,
            height = input$save_height4 * 2,
            scale = 0.5,
            dpi = input$save_dpi4
          )
        }
      )
    }
  )
}
