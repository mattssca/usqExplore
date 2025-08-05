library(shiny)
library(plotly)
library(bslib)
library(DT)
library(shinyWidgets)
library(RColorBrewer)
library(pheatmap)
library(later)
library(colourpicker)
library(survival)
library(survminer)

# Load required data
load("uroscanseq_meta.Rdata")
metadata = uroscanseq_meta
load("expr_met_uroscanseq.Rdata")
expr_mat = as.matrix(uroscanseq_expr)
load("sample_order.Rdata")
names(sample_order) <- colnames(expr_mat)[sample_order]
load("genes_to_plot_df.Rdata")
load("exp_genes.Rdata")
load("var_categories.Rdata")

lund_colors <- c(
  Uro = "#3cb44b",
  GU = "#4363d8",
  BaSq = "#CD2626",
  Mes = "#f58231",
  ScNE = "#A020F0",
  UroA = "#3cb44b",
  UroB = "#8B1A1A",
  UroC = "#006400"
)

custom_colorscale <- list(
  list(0, "green"),
  list(0.5, "black"),
  list(1, "red")
)

# --- Forest Plot Signature Names (global scope) ---
forest_signature_names <- c(
  "Proliferation Score", "Mol. grade (WHO1999)", "Mol. grade (WHO2022)", "Progression Score", 
  "Immune 141_UP", "NK Cells", "T Cells", "CD8+ T Cells", "Cytotoxicity Score", "B Cells", 
  "Myeloid DCs", "Monocytic Lineage", "Macrophages", "M2 Macrophages", "Neutrophils", 
  "Stromal 141_UP", "Fibroblasts", "Endothelial Cells", "Smooth Muscle"
)
forest_metadata_names <- c(
  "Proliferation.Score", "Mol.Grade.WHO.1999.Score", "Mol.Grade.WHO.2022.Score", "Progression.Score", 
  "Immune.141.Up", "NK.Cells", "T.Cells", "T.Cells.CD8", "Cytotoxicity.Score", "B.Cells", 
  "Myeloid.DCs", "Monocytic.Lineage", "Macrophages", "M2.Macrophage", "Neutrophils", 
  "Stromal.141.Up", "Fibroblasts", "Endothlial.Cells", "Smooth.Muscle"
)

#' Bin Numeric Variables (internal)
int_bin_numeric_variables <- function(this_data = NULL, num_bins = NULL){
  numeric_columns = sapply(this_data, is.numeric)
  bin_column = function(column, num_bins){
    range_min = min(column, na.rm = TRUE)
    range_max = max(column, na.rm = TRUE)
    binned_column = cut(column, 
                        breaks = seq(range_min, range_max, length.out = num_bins + 1), 
                        labels = FALSE, 
                        include.lowest = TRUE)
    return(binned_column)
  }
  this_data[numeric_columns] <- lapply(this_data[numeric_columns], bin_column, num_bins = num_bins)
  return(this_data)
}

ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty"),
  tags$head(
    tags$style(HTML("
      .annotation-legends-row {
        display: flex;
        flex-wrap: wrap;
        gap: 30px;
        align-items: flex-start;
      }
      .annotation-legend-item {
        min-width: 120px;
        margin-right: 20px;
      }
      .main-title { font-size: 2em; font-weight: bold; margin-bottom: 10px; }
      .plot-container { background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 8px #ccc; }
      .filter-title { font-size: 1.3em; font-weight: bold; margin-bottom: 10px; }
      .after-plot-text { margin-top: 15px; font-size: 1.1em; }
      .section-title { font-size: 1.2em; font-weight: bold; margin-top: 15px; margin-bottom: 8px; }
      /* Custom green style for file input Browse button */
      .btn-file, .fileinput-button {
        background-color: #dd6262 !important;
        border-color: #dd6262 !important;
        color: #fff !important;
      }
    "))
  ),
  fluidRow(
    column(
      width = 12,
      tags$img(src = "logo.png", height = "100px", style = "margin-bottom:10px;"),
      # Removed helpText from here; moved to filter box below
    )
  ),
  tabsetPanel(
    tabPanel("Welcome",
      fluidRow(
        column(12,
          tags$div(style = "margin-top:20px; font-size:1.15em;",
            tags$p("Welcome to the UROSCANSEQ Explorer! This app allows you to interactively explore gene expression and clinical metadata from the UROSCANSEQ dataset."),
            tags$p("Start by filtering samples in the Metadata tab. All downstream analyses will use the selected samples."),
            tags$ul(
              tags$li("Metadata tab: Filter samples and view metadata table."),
              tags$li("Exploration Plots: Create interactive plots based on selected samples."),
              tags$li("Heatmap: Visualize gene expression heatmaps."),
              tags$li("Survival: Kaplan-Meier survival analysis."),
              tags$li("LundTax Signature Forest Plot: Cox model forest plots for gene signatures.")
            ),
            tags$p(
              "For more information and source code, visit our ",
              tags$a(href = "https://github.com/LundBladderCancerGroup", "GitHub", target = "_blank"),
              "."
            )
          )
        )
      )
    ),
    tabPanel("Metadata",
      fluidRow(
        column(5,
          tags$div(class = "main-title", "Metadata"),
          tags$div(style = "margin-bottom:8px; font-size:1.15em; color:#2c3e50;",
            tags$p(
              "You can subset your data in two ways:"
            ),
            tags$ul(
              tags$li("Apply one or more filters below to interactively select samples based on metadata. Each filter will update the sample selection dynamically."),
              tags$li("Alternatively, if you already have a list of sample IDs, upload a txt file (one sample per line) to subset directly to those samples.")
            ),
            tags$p("If no filters are added and no file is uploaded, all samples are included.")
          ),
          fileInput("sample_id_file", "Upload sample ID list", accept = ".txt"),
          tags$div(id = "filter_ui_container"),
          actionBttn("add_filter", "Add Filter", style = "simple", color = "success"),
          tags$div(style = "height:30px"),
          tags$div(style = "height:18px"),
          tags$div(style = "height:18px"),
          downloadBttn("download_ids", "Download Sample IDs", style = "unite", color = "success")
        ),
        column(7,
          tags$div(class = "section-title", "Metadata Table"),
          DT::dataTableOutput("table")
        )
      )
    ),
    tabPanel("Exploration Plots",
      fluidRow(
        column(5,
          tags$div(class = "main-title", "Exploration Plots"),
          tags$div(style = "margin-bottom:8px; font-size:1.08em; color:#2c3e50;",
            "Create interactive plots based on selected samples and metadata. Choose variables and plot type to explore distributions and relationships in the filtered dataset."
          ),
          tags$div(
            style = "margin-bottom:8px;",
            strong("How to construct a plot:"),
            tags$ul(
              tags$li("Select an X variable (required) to define the main grouping or axis for your plot."),
              tags$li("Optionally select a Y variable for boxplots or scatterplots."),
              tags$li("Optionally select a color variable to further group or color samples."),
              tags$li("Choose a plot type below. The plot will update automatically as you change selections."),
              tags$li("You can also sort/rank by the X variable and show proportions using the switches below.")
            )
          ),
          pickerInput("xvar", "X variable", choices = var_categories, selected = "Prediction.5.Class", 
                      options = list(`live-search` = TRUE)),
          pickerInput("yvar", "Y variable (optional)", choices = c("None", var_categories), selected = "None", 
                      options = list(`live-search` = TRUE)),
          pickerInput("color", "Color by (optional)", choices = c("None", var_categories), selected = "Prediction.5.Class", 
                      options = list(`live-search` = TRUE)),
          radioGroupButtons("plot_type", "Plot type",
                            choices = c("Histogram" = "hist", "Barplot" = "bar", "Boxplot" = "box", "Scatterplot" = "scatter"),
                            selected = "hist", justified = TRUE
          ),
          prettySwitch("sort_x", "Sort/Ranks by X variable", value = FALSE, status = "info"),
          prettySwitch("show_proportion", "Show proportions", value = FALSE)
        ),
        column(7,
          div(class = "plot-container", plotlyOutput("plot", height = "565px")),
          div(class = "after-plot-text", textOutput("n_samples"))
        )
      )
    ),
    tabPanel("Heatmap",
      tags$div(class = "main-title", "Heatmap"),
      tags$div(style = "margin-bottom:8px; font-size:1.08em; color:#2c3e50;",
        "Visualize gene expression heatmaps for selected genes and samples. Add annotation tracks and meta-genes to explore sample-level information and gene signatures in the filtered dataset."
      ),
      tags$div(
        style = "margin-bottom:8px;",
        strong("How to construct a heatmap:"),
        tags$ul(
          tags$li("Upload a custom gene list or select a predefined signature."),
          tags$li("Optionally define meta-genes for annotation tracks."),
          tags$li("Add annotation variables to display sample-level information."),
          tags$li("Click 'Plot Heatmap' to generate the heatmap for selected samples and genes.")
        )
      ),
      tags$div(class = "section-title", "Genes"),
      fileInput("gene_file", "Upload gene list (txt, one gene per line)", accept = ".txt"),
      pickerInput(
        "gene_signature",
        "Select predefined gene list (signature)",
        choices = c("Custom Genes", unique(genes_to_plot_df$signature)),
        multiple = FALSE,
        selected = "Early_CC",
        options = list(`live-search` = TRUE)
      ),
      uiOutput("gene_select_ui"),
      tags$hr(),
      tags$div(class = "section-title", "Define Meta-genes"),
      uiOutput("metagene_ui"),
      actionBttn("add_metagene", "Add Meta-gene", color = "success", style = "simple"),
      tags$hr(),
      tags$div(class = "section-title", "Add Annotation Track"),
      pickerInput(
        "heatmap_ann_vars",
        "Add annotation variable(s) for heatmap (in addition to Subtype 5 Class)",
        choices = setdiff(names(metadata), "Prediction.5.Class"),
        multiple = TRUE,
        options = list(`live-search` = TRUE)
      ),
      uiOutput("ann_color_pickers"),
      tags$hr(),
      fluidRow(
        column(
          width = 6,
          actionBttn("plot_heatmap", "Plot Heatmap", width = "100%", style = "simple", color = "success"),
          tags$div(style = "height:20px") 
        ),
        column(
          width = 6,
          div(style = "display: flex; align-items: center; height: 100%; justify-content: flex-end;",
              prettySwitch("sort_heatmap_samples", "Sort samples by Proliferation Score", value = TRUE, status = "info")
          )
        )
      ),
      plotlyOutput("heatmap_plot"),
      tags$div(style = "height:20px"),
      uiOutput("heatmap_legends")
    ),
    tabPanel("Gene Expression Outcomes",
      fluidRow(
        column(5,
          tags$div(class = "main-title", "Gene Expression Outcomes"),
          tags$div(style = "margin-bottom:12px; font-size:1.08em; color:#2c3e50;",
            "This tab allows you to perform survival analysis using gene expression or meta-gene data. Kaplan-Meier curves are generated for selected genes or meta-genes, stratified by expression level and clinical endpoint. All analyses use the currently filtered sample set."
          ),
          tags$div(
            style = "margin-bottom:8px;",
            strong("How to construct a survival plot:"),
            tags$ul(
              tags$li("Select a gene or meta-gene to stratify samples."),
              tags$li("Choose a clinical endpoint and cutoff method for grouping."),
              tags$li("Optionally filter samples by subtype."),
              tags$li("The plot and risk table will update automatically for selected samples.")
            )
          ),
          tags$div(class = "section-title", "Kaplan-Meier Plot"),
          selectizeInput("km_gene", "Select gene or meta-gene",
                      choices = exp_genes, selected = NULL, options = list(server = TRUE)),
          selectInput(
            "km_endpoint", "Select clinical endpoint",
            choices = c(
              "Progression",
              "Clinical Progression",
              "Clinical Progression During Follow-up",
              "BCG Any",
              "BCG Adequate"
            ),
            selected = "Clinical Progression"
          ),
          pickerInput(
            "km_filter_subtypes",
            "Filter samples by subtype",
            choices = levels(metadata$Prediction.5.Class),
            selected = levels(metadata$Prediction.5.Class),
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `live-search` = TRUE)
          ),
          radioButtons("km_cut", "Cutoff method",
                       choices = c("Median" = "median", "Tertiles" = "tertile", "Quartiles" = "quartile"),
                       inline = TRUE
          ),
          checkboxInput("km_show_confint", "Show confidence intervals", value = TRUE)
        ),
        column(7,
          plotOutput("km_plot", height = "700px")
        )
      ),
      fluidRow(
        column(12,
          # The risk table is part of the ggplot output in survminer::ggsurvplot, so this is a placeholder for future extension.
          # If you want a separate risk table, you can add DT::dataTableOutput("km_risk_table") here.
        )
      )
    ),
    tabPanel("Risk by Gene Signature",
      fluidRow(
        column(5,
          tags$div(class = "main-title", "Risk by Gene Signature"),
          tags$div(style = "margin-bottom:12px; font-size:1.08em; color:#2c3e50;",
            "This tab displays risk analysis by predefined gene signature using Cox proportional hazards models. Hazard ratios and confidence intervals are calculated for selected gene signatures and clinical endpoints, using the currently filtered sample set."
          ),
          tags$div(
            style = "margin-bottom:8px;",
            strong("How to construct a forest plot:"),
            tags$ul(
              tags$li("Select a clinical endpoint for survival analysis."),
              tags$li("Filter samples by subtype as needed."),
              tags$li("Choose one or more gene signatures to display."),
              tags$li("Adjust the number of bins for numeric signatures to control grouping granularity."),
              tags$li("The plot and results table will update automatically for selected samples.")
            )
          ),
          tags$div(class = "section-title", "Risk by Gene Signature (Cox Model)"),
          selectInput(
            "forest_endpoint", "Select clinical endpoint",
            choices = c(
              "Progression",
              "Clinical Progression",
              "Clinical Progression During Follow-up",
              "BCG Any",
              "BCG Adequate"
            ),
            selected = "Clinical Progression"
          ),
          pickerInput(
            "forest_filter_subtypes",
            "Filter samples by subtype",
            choices = levels(metadata$Prediction.5.Class),
            selected = levels(metadata$Prediction.5.Class),
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `live-search` = TRUE)
          ),
          pickerInput(
            "forest_signature_select",
            "Select signature(s) to show",
            choices = c("All", forest_signature_names),
            selected = "All",
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `live-search` = TRUE)
          ),
          sliderInput(
            "forest_binsize",
            "Number of bins for numeric signatures",
            min = 2, max = 20, value = 10, step = 1
          )
        ),
        column(7,
          plotOutput("forest_plot", height = "700px")
        )
      ),
      fluidRow(
        column(12,
          DT::dataTableOutput("forest_table"),
          downloadButton("download_forest_table", "Download Table", style = "margin-top:10px;")
        )
      )
    ),
    tabPanel("Correlation Explorer",
      fluidRow(
        column(5,
          tags$div(class = "main-title", "Correlation Explorer"),
          tags$div(style = "margin-bottom:8px; font-size:1.08em; color:#2c3e50;",
            "Explore relationships between clinical metadata variables using scatterplots. Select any two numeric metadata variables to visualize their correlation."
          ),
          tags$div(
            style = "margin-bottom:8px;",
            strong("How to construct a scatter plot:"),
            tags$ul(
              tags$li("Select X and Y variables from the available numeric metadata fields."),
              tags$li("Click 'Plot Scatter' to visualize their relationship."),
              tags$li("A regression line is added for quick assessment of correlation.")
            )
          ),
          selectInput("corr_scatter_x", "X variable", choices = names(metadata)[sapply(metadata, is.numeric)]),
          selectInput("corr_scatter_y", "Y variable", choices = names(metadata)[sapply(metadata, is.numeric)]),
          actionBttn("plot_corr_scatter", "Plot Scatter", style = "simple", color = "success")
        ),
        column(7,
          plotOutput("corr_scatter_plot", height = "500px")
        )
      )
    ),
    tabPanel("External Data Integration",
      fluidRow(
        column(12,
          tags$div(class = "main-title", "External Data Integration"),
          tags$div(style = "margin-top:20px; font-size:1.15em; color:#2c3e50;",
            "Coming soon..."
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # --- Publication/Export Logic ---
  output$download_metadata_csv <- downloadHandler(
    filename = function() {
      paste0("filtered_metadata_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(data_filtered(), file, row.names = FALSE)
    }
  )

  output$preview_metadata_table <- DT::renderDataTable({
    data_filtered()
  }, options = list(pageLength = 10))

  output$download_main_plot <- downloadHandler(
    filename = function() {
      paste0("main_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      # Save the main plot (from Exploration Plots tab) as PNG
      png(file, width = 900, height = 600)
      # Try to re-create the plot using the same logic as output$plot
      df <- data_filtered()
      x <- input$xvar
      y <- if (input$yvar != "None") input$yvar else NULL
      color <- if (input$color != "None") input$color else NULL
      show_prop <- isTRUE(input$show_proportion)
      color_map <- if (!is.null(color) && color %in% c("Prediction.5.Class", "Prediction.7.Class")) lund_colors else NULL
      hide_x_ticks <- identical(x, "Sample.ID")
      # Only support hist/bar/box/scatter for PNG export
      if (input$plot_type == "hist") {
        hist(df[[x]], main = paste("Histogram of", x), xlab = x, col = "#73bac9")
      } else if (input$plot_type == "bar") {
        tab <- as.data.frame(table(df[[x]], if (!is.null(color)) df[[color]] else NULL))
        colnames(tab)[1:2] <- c("x", "color")
        tab$Freq <- as.numeric(tab$Freq)
        tab <- tab[tab$Freq > 0, ]
        barplot(tab$Freq, names.arg = tab$x, col = "#73bac9", main = paste("Barplot of", x))
      } else if (input$plot_type == "box" && !is.null(y)) {
        boxplot(df[[y]] ~ df[[x]], col = "#73bac9", main = paste("Boxplot of", y, "by", x), xlab = x, ylab = y)
      } else if (input$plot_type == "scatter" && !is.null(y)) {
        plot(df[[x]], df[[y]], xlab = x, ylab = y, pch = 19, col = "#73bac9", main = paste("Scatterplot of", x, "vs", y))
        abline(lm(df[[y]] ~ df[[x]]), col = "red")
      }
      dev.off()
    }
  )

  output$download_main_table <- downloadHandler(
    filename = function() {
      paste0("main_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # Export the main results table (from Exploration Plots tab)
      write.csv(data_filtered(), file, row.names = FALSE)
    }
  )
  filter_count <- reactiveVal(0)
  filter_ids <- reactiveVal(character(0))

  # Fix Add Filter button logic
  observeEvent(input$add_filter, {
    # Store current filter input values before adding new filter
    ids_prev <- isolate(filter_ids())
    filter_values <- lapply(ids_prev, function(fid) {
      list(
        cat = input[[paste0(fid, "_cat")]],
        var = input[[paste0(fid, "_var")]],
        val = input[[paste0(fid, "_val")]]
      )
    })
    names(filter_values) <- ids_prev

    count <- filter_count() + 1
    filter_count(count)
    id <- paste0("filter_", count)
    ids <- c(filter_ids(), id)
    filter_ids(ids)
    insertUI(
      selector = "#filter_ui_container",
      where = "beforeEnd",
      ui = tags$div(id = id,
        fluidRow(
          column(5,
            pickerInput(paste0(id, "_cat"), "Category", choices = names(var_categories), width = "100%", options = list(`live-search` = TRUE)),
            uiOutput(paste0(id, "_var_ui")),
            uiOutput(paste0(id, "_val_ui")),
            actionBttn(paste0(id, "_remove"), "Remove", style = "simple", color = "danger", size = "md"),
            tags$div(style = "height:18px")
          )
        )
      )
    )
    # Restore previous filter selections for all filters after UI is rendered
    later::later(function() {
      for (i in seq_along(ids_prev)) {
        fid <- ids_prev[i]
        vals <- filter_values[[fid]]
        if (!is.null(vals$cat)) updatePickerInput(session, paste0(fid, "_cat"), selected = vals$cat)
        # Wait for var UI to be rendered before updating var and val
        later::later(function() {
          if (!is.null(vals$var)) updatePickerInput(session, paste0(fid, "_var"), selected = vals$var)
          later::later(function() {
            if (!is.null(vals$val)) {
              if (!is.null(vals$var) && is.numeric(metadata[[vals$var]])) {
                updateSliderInput(session, paste0(fid, "_val"), value = vals$val)
              } else {
                updatePickerInput(session, paste0(fid, "_val"), selected = vals$val)
              }
            }
          }, delay = 0.1)
        }, delay = 0.1)
      }
    }, delay = 0.2)
  })

  # --- Correlation Explorer Logic ---
  observe({
    req(input$corr_scatter_x, input$corr_scatter_y)
    output$corr_scatter_plot <- renderPlot({
      x <- input$corr_scatter_x
      y <- input$corr_scatter_y
      df <- metadata
      if (is.null(x) || is.null(y) || !(x %in% names(df)) || !(y %in% names(df))) return(NULL)
      xvals <- df[[x]]
      yvals <- df[[y]]
      # Remove NAs
      valid <- complete.cases(xvals, yvals)
      xvals <- xvals[valid]
      yvals <- yvals[valid]
      plot(xvals, yvals, xlab = x, ylab = y, pch = 19, col = "#73bac9", main = paste("Scatterplot of", x, "vs", y))
      fit <- lm(yvals ~ xvals)
      abline(fit, col = "red")
      # Pearson correlation
      cor_test <- cor.test(xvals, yvals)
      r <- cor_test$estimate
      pval <- cor_test$p.value
      # Confidence interval for R
      ci_r <- cor_test$conf.int
      # Regression line equation
      intercept <- coef(fit)[1]
      slope <- coef(fit)[2]
      eqn <- paste0("y = ", formatC(slope, digits = 3, format = "f"), "x + ", formatC(intercept, digits = 3, format = "f"))
      # Compose legend text
      legend_text <- c(
        paste0("R = ", formatC(r, digits = 3, format = "f")),
        paste0("p = ", formatC(pval, digits = 3, format = "e")),
        paste0("95% CI for R: [", formatC(ci_r[1], digits = 3, format = "f"), ", ", formatC(ci_r[2], digits = 3, format = "f"), "]"),
        paste0("Regression: ", eqn)
      )
      legend("topright", legend = legend_text, bty = "n", text.col = c("#d35741", "#4363d8", "#499f8b", "#73bac9"), cex = 1.1)
    })
  })
  # ...existing code...
  observe({
    lapply(filter_ids(), function(id) {
      output[[paste0(id, "_var_ui")]] <- renderUI({
        cat_input <- input[[paste0(id, "_cat")]]
        if (is.null(cat_input)) return(NULL)
        pickerInput(paste0(id, "_var"), "Variable", choices = var_categories[[cat_input]], width = "100%",
                    options = list(`live-search` = TRUE))
      })
    })
  })
  
  observe({
    lapply(filter_ids(), function(id) {
      observeEvent(input[[paste0(id, "_remove")]], {
        removeUI(selector = paste0("#", id))
        filter_ids(setdiff(filter_ids(), id))
      }, ignoreInit = TRUE)
    })
  })
  
  observe({
    lapply(filter_ids(), function(id) {
      output[[paste0(id, "_val_ui")]] <- renderUI({
        var_input <- input[[paste0(id, "_var")]]
        if (is.null(var_input)) return(NULL)
        vals <- metadata[[var_input]]
        if (is.numeric(vals)) {
          rng <- range(vals, na.rm = TRUE)
          sliderInput(
            paste0(id, "_val"), "Range",
            min = rng[1], max = rng[2], value = rng, step = diff(rng)/100,
            width = "100%"
          )
        } else {
          pickerInput(paste0(id, "_val"), "Value(s)", choices = unique(vals), multiple = TRUE, width = "100%",
                      options = list(`live-search` = TRUE))
        }
      })
    })
  })
  
data_filtered <- reactive({
  df <- metadata
  # If a sample ID file is uploaded, subset to those sample IDs first
  if (!is.null(input$sample_id_file)) {
    ids <- tryCatch({
      readLines(input$sample_id_file$datapath)
    }, error = function(e) NULL)
    ids <- trimws(ids)
    ids <- ids[ids != ""]
    df <- df[df$Sample.ID %in% ids, , drop = FALSE]
  }
  for (id in filter_ids()) {
    var <- input[[paste0(id, "_var")]]
    val <- input[[paste0(id, "_val")]]
    if (!is.null(var) && !is.null(val)) {
      if (is.numeric(metadata[[var]])) {
        df <- df[df[[var]] >= val[1] & df[[var]] <= val[2], , drop = FALSE]
      } else if (length(val) > 0) {
        df <- df[df[[var]] %in% val, , drop = FALSE]
      }
    }
  }
  df
})
  
  output$n_samples <- renderText({
    n <- nrow(data_filtered())
    paste("Samples after filtering:", n)
  })
  
  output$plot <- renderPlotly({
    df <- data_filtered()
    x <- input$xvar
    y <- if (input$yvar != "None") input$yvar else NULL
    color <- if (input$color != "None") input$color else NULL
    show_prop <- isTRUE(input$show_proportion)
    color_map <- if (!is.null(color) && color %in% c("Prediction.5.Class", "Prediction.7.Class")) lund_colors else NULL
    
    if (input$sort_x && !is.null(x) && !is.null(y)) {
      ord <- order(df[[y]], na.last = TRUE)
      df <- df[ord, , drop = FALSE]
      df[[x]] <- factor(df[[x]], levels = unique(df[[x]]))
    }
    
    hide_x_ticks <- identical(x, "Sample.ID")
    
    if (input$plot_type == "hist") {
      p <- plot_ly(
        df, x = ~get(x),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "histogram",
        histnorm = if (show_prop) "probability" else ""
      ) %>%
        layout(
          yaxis = list(title = if (show_prop) "Proportion" else "Count"),
          xaxis = list(showticklabels = !hide_x_ticks)
        )
      p
    } else if (input$plot_type == "bar") {
      tab <- as.data.frame(table(df[[x]], if (!is.null(color)) df[[color]] else NULL))
      colnames(tab)[1:2] <- c("x", "color")
      tab$Freq <- as.numeric(tab$Freq)
      tab <- tab[tab$Freq > 0, ]
      if (show_prop) {
        total <- tapply(tab$Freq, tab$x, sum)
        tab$yval <- tab$Freq / total[tab$x]
        ylab <- "Proportion"
      } else {
        tab$yval <- tab$Freq
        ylab <- "Count"
      }
      p <- plot_ly(
        tab, x = ~x, y = ~yval, type = "bar",
        color = if (!is.null(color)) ~color else NULL,
        colors = color_map
      ) %>%
        layout(
          yaxis = list(title = ylab),
          xaxis = list(showticklabels = !hide_x_ticks)
        )
      p
    } else if (input$plot_type == "box") {
      plot_ly(
        df, x = ~get(x), y = ~get(y),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "box"
      ) %>%
        layout(xaxis = list(showticklabels = !hide_x_ticks))
    } else if (input$plot_type == "scatter") {
      plot_ly(
        df, x = ~get(x), y = ~get(y),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "scatter", mode = "markers"
      ) %>%
        layout(xaxis = list(showticklabels = !hide_x_ticks))
    }
  })
  
  output$table <- DT::renderDataTable({
    data_filtered()
  }, options = list(pageLength = 10))
  
  output$download_ids <- downloadHandler(
    filename = function() {
      paste0("sample_ids_", Sys.Date(), ".txt")
    },
    content = function(file) {
      writeLines(as.character(data_filtered()$Sample.ID), file)
    }
  )
  
  # --- Meta-gene logic ---
  metagenes <- reactiveVal(list())
  metagene_exprs <- reactiveVal(list())
  
  # Helper to get current meta-gene input values
  get_current_metagene_inputs <- function(mg_list, input) {
    out <- list()
    for (mg_id in names(mg_list)) {
      out[[mg_id]] <- list(
        name = input[[paste0(mg_id, "_name")]],
        genes = input[[paste0(mg_id, "_genes")]],
        filter = input[[paste0(mg_id, "_filter")]]
      )
    }
    out
  }
  
  observeEvent(input$add_metagene, {
    mg_list <- metagenes()
    current_inputs <- get_current_metagene_inputs(mg_list, input)
    mg_id <- paste0("mg_", length(mg_list) + 1)
    mg_list[[mg_id]] <- list(name = NULL, genes = NULL)
    metagenes(mg_list)
    
    # Delay restoring inputs until UI is updated
    isolate({
      later::later(function() {
        for (id in names(current_inputs)) {
          vals <- current_inputs[[id]]
          updateTextInput(session, paste0(id, "_name"), value = vals$name)
          updateTextInput(session, paste0(id, "_filter"), value = vals$filter)
          updateSelectizeInput(session, paste0(id, "_genes"), selected = vals$genes, server = TRUE)
        }
      }, delay = 0.1)
    })
  })
  
  # Meta-gene UI with filter
  output$metagene_ui <- renderUI({
    mg_list <- metagenes()
    genes <- exp_genes
    mg_ui <- lapply(seq_along(mg_list), function(i) {
      mg_id <- names(mg_list)[i]
      fluidRow(
        column(3, textInput(paste0(mg_id, "_filter"), "Filter genes")),
        column(3, textInput(paste0(mg_id, "_name"), "Meta-gene Name")),
        column(4, selectizeInput(paste0(mg_id, "_genes"), "Genes", choices = NULL, multiple = TRUE, options = list(server = TRUE))),
        column(2, actionButton(paste0(mg_id, "_generate"), "Generate", class = "btn-primary btn-sm", style = "margin-top:25px;"))
      )
    })
    do.call(tagList, mg_ui)
  })
  
  # Server-side filtering for meta-gene gene selection
  observe({
    mg_list <- metagenes()
    for (mg_id in names(mg_list)) {
      observe({
        filter_val <- input[[paste0(mg_id, "_filter")]]
        filtered_genes <- exp_genes
        if (!is.null(filter_val) && nzchar(filter_val)) {
          filtered_genes <- grep(filter_val, exp_genes, value = TRUE, ignore.case = TRUE)
        }
        updateSelectizeInput(session, paste0(mg_id, "_genes"), choices = filtered_genes, server = TRUE)
      })
    }
  })
  
  observe({
    mg_list <- metagenes()
    for (mg_id in names(mg_list)) {
      observeEvent(input[[paste0(mg_id, "_generate")]], {
        name <- input[[paste0(mg_id, "_name")]]
        genes <- input[[paste0(mg_id, "_genes")]]
        if (!is.null(name) && !is.null(genes) && length(genes) > 0) {
          valid_genes <- genes[genes %in% rownames(expr_mat)]
          if (length(valid_genes) > 0) {
            expr <- colMeans(expr_mat[valid_genes, , drop = FALSE], na.rm = TRUE)
            mg_exprs <- metagene_exprs()
            mg_exprs[[name]] <- expr
            metagene_exprs(mg_exprs)
            showNotification(paste("Meta-gene", name, "generated!"), type = "message")
          }
        }
      }, ignoreInit = TRUE)
    }
  })
  
  # --- HEATMAP TAB LOGIC ---
  # Read uploaded gene list if present
  uploaded_genes <- reactive({
    req(input$gene_file)
    genes <- readLines(input$gene_file$datapath)
    genes <- trimws(genes)
    genes <- genes[genes != ""]
    genes
  })
  
  # Determine selected genes for heatmap
  selected_genes_heatmap <- reactive({
    genes <- exp_genes
    # Priority: uploaded file > predefined signature > manual selection
    if (!is.null(input$gene_file)) {
      up_genes <- uploaded_genes()
      selected_genes <- up_genes[up_genes %in% genes]
      if (length(selected_genes) > 0) return(selected_genes)
    }
    if (!is.null(input$gene_signature) && input$gene_signature != "" && input$gene_signature != "Custom Genes") {
      sig_genes <- genes_to_plot_df$gene[genes_to_plot_df$signature == input$gene_signature]
      selected_genes <- sig_genes[sig_genes %in% genes]
      if (length(selected_genes) > 0) return(selected_genes)
    }
    # Fallback to manual selection
    if (!is.null(input$genes)) {
      selected_genes <- input$genes[input$genes %in% genes]
      if (length(selected_genes) > 0) return(selected_genes)
    }
    head(genes, 10)
  })
  
  # Heatmap gene selection UI with filter
  output$gene_select_ui <- renderUI({
    genes <- exp_genes
    if (!is.null(input$gene_signature) && input$gene_signature == "Custom Genes") {
      selectizeInput("genes", "Select Genes for Heatmap (manual entry)", choices = genes, multiple = TRUE, selected = NULL, options = list(server = TRUE))
    }
  })
  
  # Annotation color pickers UI
  output$ann_color_pickers <- renderUI({
    req(input$heatmap_ann_vars)
    out <- list()
    for (var in input$heatmap_ann_vars) {
      # Get unique levels for this variable
      vals <- metadata[[var]]
      if (is.numeric(vals)) next # skip numeric
      if (is.factor(vals)) {
        uniq <- levels(vals)
      } else {
        uniq <- sort(unique(as.character(vals)))
      }
      out[[var]] <- tags$div(
        tags$label(paste("Pick colors for", var)),
        lapply(seq_along(uniq), function(i) {
          colourpicker::colourInput(
            inputId = paste0("color_", var, "_", uniq[i]),
            label = uniq[i],
            value = RColorBrewer::brewer.pal(max(3, length(uniq)), "Set2")[i %% 8 + 1]
          )
        }),
        tags$hr(style="margin:5px 0;")
      )
    }
    if (length(out) > 0) out
  })
  
  observeEvent(input$plot_heatmap, {
    req(expr_mat)
    sample_ids <- data_filtered()$Sample.ID
    valid_samples <- sample_ids[sample_ids %in% colnames(expr_mat)]
    if (length(valid_samples) == 0) {
      output$heatmap_plot <- renderPlotly({ plotly_empty() })
      output$heatmap_legends <- renderUI({ NULL })
      return()
    }
    ann <- data_filtered()
    ann <- ann[ann$Sample.ID %in% valid_samples, ]
    ann <- ann[order(ann$Prediction.5.Class), ]
    
    # --- Sorting logic ---
    if (!is.null(input$sort_heatmap_samples) && input$sort_heatmap_samples && exists("sample_order")) {
      ann$sample_order_val <- sample_order[as.character(ann$Sample.ID)]
      ann <- ann[order(ann$Prediction.5.Class, ann$sample_order_val, na.last = TRUE), ]
    }
    sorted_samples <- ann$Sample.ID
    mat <- expr_mat[selected_genes_heatmap(), sorted_samples, drop = FALSE]
    
    # --- Add meta-genes as annotation tracks ---
    mg_exprs <- metagene_exprs()
    ann_tracks <- list()
    legends <- list()
    subtype_levels <- levels(metadata$Prediction.5.Class)
    subtype_colors <- lund_colors[subtype_levels]
    # Always include Prediction.5.Class annotation track
    ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
      z = matrix(as.numeric(factor(ann$Prediction.5.Class, levels = subtype_levels)), nrow = 1),
      x = sorted_samples,
      y = "Subtype 5 Class",
      type = "heatmap",
      showscale = FALSE,
      colors = subtype_colors,
      hoverinfo = "text",
      text = paste("Subtype:", ann$Prediction.5.Class)
    )
    legends[[length(legends) + 1]] <- tags$div(
      tags$b("Subtype 5 Class:"),
      tags$ul(
        lapply(subtype_levels, function(lvl) {
          tags$li(
            tags$span(style = paste0("display:inline-block;width:15px;height:15px;background:", lund_colors[lvl], ";margin-right:5px;"), ""),
            lvl
          )
        })
      )
    )
    # Add meta-gene annotation tracks
    if (length(mg_exprs) > 0) {
      for (mg_name in names(mg_exprs)) {
        mg_vec <- mg_exprs[[mg_name]][sorted_samples]
        ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
          z = matrix(mg_vec, nrow = 1),
          x = sorted_samples,
          y = mg_name,
          type = "heatmap",
          colorscale = custom_colorscale,
          showscale = TRUE,
          colorbar = list(title = mg_name),
          hoverinfo = "text",
          text = paste(mg_name, ":", signif(mg_vec, 3))
        )
        legends[[length(legends) + 1]] <- tags$div(
          tags$b(mg_name, ": Meta-gene (mean expression of selected genes)")
        )
      }
    }
    # Add user-selected annotation tracks
    ann_vars <- input$heatmap_ann_vars
    if (!is.null(ann_vars) && length(ann_vars) > 0) {
      for (i in seq_along(ann_vars)) {
        var <- ann_vars[i]
        vals <- as.character(ann[[var]])
        # Assign colors (categorical: use user color for each level, numeric: Viridis)
   if (is.numeric(ann[[var]])) {
  # Calculate quantiles from the filtered data
  vals_num <- ann[[var]]
  q5 <- quantile(vals_num, 0.05, na.rm = TRUE)
  q50 <- median(vals_num, na.rm = TRUE)
  q95 <- quantile(vals_num, 0.95, na.rm = TRUE)
  # Normalize for Plotly colorscale (0-1)
  rng <- range(c(q5, q50, q95))
  norm <- function(x) if (diff(rng) == 0) 0.5 else (x - rng[1]) / diff(rng)
  colorscale <- list(
    list(0, "blue"),
    list(norm(q5), "blue"),
    list(norm(q50), "white"),
    list(norm(q95), "red"),
    list(1, "red")
  )
  z <- matrix(vals_num, nrow = 1)
  show_scale <- FALSE
  legend_html <- NULL
} else {
          uniq <- sort(unique(vals))
          n_uniq <- length(uniq)
          # Get user colors for each level
          pal <- sapply(uniq, function(u) {
            input_id <- paste0("color_", var, "_", u)
            col <- input[[input_id]]
            if (is.null(col) || !nzchar(col)) RColorBrewer::brewer.pal(max(3, n_uniq), "Set2")[which(uniq == u) %% 8 + 1] else col
          }, USE.NAMES = TRUE)
          col_map <- setNames(pal, uniq)
          colorscale <- lapply(seq_along(uniq), function(j) list((j-1)/(n_uniq-1), col_map[uniq[j]]))
          z <- matrix(as.numeric(factor(vals, levels = uniq)), nrow = 1)
          show_scale <- FALSE
          legend_html <- tags$div(
            tags$b(var, ":"),
            tags$ul(
              lapply(uniq, function(u) {
                tags$li(
                  tags$span(style = paste0("display:inline-block;width:15px;height:15px;background:", col_map[u], ";margin-right:5px;"), ""),
                  u
                )
              })
            )
          )
        }
        ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
          z = z,
          x = sorted_samples,
          y = var,
          type = "heatmap",
          showscale = show_scale,
          colorscale = colorscale,
          hoverinfo = "text",
          text = paste(var, ":", vals)
        )
        if (!is.null(legend_html)) {
          legends[[length(legends) + 1]] <- legend_html
        }
      }
    }
    
    # Main heatmap
    heatmap_main <- plot_ly(
      z = as.matrix(mat),
      x = colnames(mat),
      y = rownames(mat),
      type = "heatmap",
      colorscale = custom_colorscale,
      zmin = -2,
      zmax = 2,
      showscale = TRUE,
      colorbar = list(title = NULL)
    )
    n_ann <- length(ann_tracks)
    ann_height <- 0.15
    main_height <- 1 - ann_height * n_ann
    if (main_height < 0.1) main_height <- 0.1
    
    output$heatmap_plot <- renderPlotly({
      subplot(
        c(ann_tracks, list(heatmap_main)),
        nrows = n_ann + 1,
        heights = c(rep(ann_height, n_ann), main_height),
        shareX = TRUE,
        titleY = TRUE
      ) %>% layout(
        margin = list(t = 40),
        yaxis = list(title = "", ticks = "", showticklabels = FALSE),
        xaxis = list(showticklabels = FALSE, ticks = ""),
        yaxis2 = list(title = NULL)
      )
    })
    output$heatmap_legends <- renderUI({
      tags$div(
        class = "annotation-legends-row",
        lapply(legends, function(lg) {
          tags$div(class = "annotation-legend-item", lg)
        })
      )
    })
  })

  # --- SURVIVAL TAB LOGIC ---
  prev_gene_choices <- reactiveVal(NULL)

  observe({
    gene_choices <- c(names(metagene_exprs()), exp_genes)
    current_sel <- input$km_gene
    valid_sel <- current_sel[current_sel %in% gene_choices]
    if (!identical(prev_gene_choices(), gene_choices)) {
      updateSelectizeInput(session, "km_gene", choices = gene_choices, selected = if (length(valid_sel) > 0) valid_sel else NULL)
      prev_gene_choices(gene_choices)
    }
  })

  # --- KM Survival Plot logic ---
  output$km_plot <- renderPlot({
    gene_selected <- input$km_gene
    cut_method <- input$km_cut
    endpoint <- input$km_endpoint
    filter_subtypes <- input$km_filter_subtypes
    metagenes <- metagene_exprs()
    genes_all <- exp_genes
    req(gene_selected, cut_method, endpoint, filter_subtypes)
    df <- data_filtered()
    gene <- gene_selected

    # Filter samples by selected subtypes
    df <- df[df$Prediction.5.Class %in% filter_subtypes, , drop = FALSE]

    # Endpoint mapping (adjust column names as needed)
    endpoint_map <- list(
      "Progression" = list(time_col = "Progression.Time", event_col = "Progression.Event"),
      "Clinical Progression" = list(time_col = "Progression.Clinical.Time", event_col = "Progression.Clinical.Event"),
      "Clinical Progression During Follow-up" = list(time_col = "Progression.During.Follow.up.Time", event_col = "Progression.During.Follow.up.Event"),
      "BCG Any" = list(time_col = "BCG.Any.Time", event_col = "BCG.Any.Event"),
      "BCG Adequate" = list(time_col = "BCG.Adequate.Time", event_col = "BCG.Adequate.Event")
    )
    ep <- endpoint_map[[endpoint]]
    if (is.null(ep)) return(NULL)

    # Filter out samples with NA in the event column
    df <- df[!is.na(df[[ep$event_col]]), , drop = FALSE]

    valid_samples <- df$Sample.ID[df$Sample.ID %in% colnames(expr_mat)]
    if (gene %in% names(metagenes)) {
      expr <- metagenes[[gene]][valid_samples]
    } else if (gene %in% rownames(expr_mat)) {
      expr <- expr_mat[gene, valid_samples]
    } else {
      return(NULL)
    }
    df <- df[df$Sample.ID %in% valid_samples, , drop = FALSE]

    surv_time <- suppressWarnings(as.numeric(as.character(df[[ep$time_col]])))
    event_col <- df[[ep$event_col]]
    if (is.factor(event_col) || is.character(event_col)) {
      surv_event <- as.numeric(as.character(event_col))
    } else {
      surv_event <- as.numeric(event_col)
    }

    valid <- !is.na(expr) & !is.na(surv_time) & !is.na(surv_event)
    if (sum(valid) == 0) {
      showNotification("No valid samples for selected gene and endpoint.", type = "error")
      return(NULL)
    }
    expr <- expr[valid]
    surv_time <- surv_time[valid]
    surv_event <- surv_event[valid]
    subtype <- droplevels(df$Prediction.5.Class[valid])

    cut_method <- input$km_cut
    group <- NULL
    if (cut_method == "median") {
      group <- ifelse(expr > median(expr, na.rm = TRUE), "High", "Low")
    } else if (cut_method == "tertile") {
      q <- quantile(expr, probs = c(1/3, 2/3), na.rm = TRUE)
      group <- cut(expr, breaks = c(-Inf, q, Inf), labels = c("Low", "Mid", "High"))
    } else if (cut_method == "quartile") {
      q <- quantile(expr, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      group <- cut(expr, breaks = c(-Inf, q, Inf), labels = c("Q1", "Q2", "Q3", "Q4"))
    } else {
      showNotification("Invalid cutoff method.", type = "error")
      return(NULL)
    }

    plot_df <- data.frame(
      surv_time = surv_time,
      surv_event = surv_event,
      group = group,
      subtype = subtype
    )

    # Always group by gene expression cutoff
    fit <- survival::survfit(survival::Surv(surv_time, surv_event) ~ group, data = plot_df)
    n_total <- length(surv_event)
    n_events <- sum(surv_event == 1, na.rm = TRUE)
    n_no_events <- sum(surv_event == 0, na.rm = TRUE)
    # Dynamic plot title based on selected subtypes
    if (length(filter_subtypes) == length(levels(df$Prediction.5.Class))) {
      plot_title <- "Kaplan-Meier: All Subtypes"
      subset_info <- "All subtypes included"
    } else {
      plot_title <- paste("Kaplan-Meier:", paste(filter_subtypes, collapse = ", "))
      subset_info <- paste("Subset:", paste(filter_subtypes, collapse = ", "))
    }
    sample_info <- paste0("Samples: ", n_total, " (", n_events, " events, ", n_no_events, " censored)")
    subtitle_text <- paste("Gene:", gene_selected, "| Endpoint:", endpoint, "\n", sample_info)

group_palette <- ggsci::pal_npg()(4)

    ggsurv <- survminer::ggsurvplot(
      fit,
      data = plot_df,
      risk.table = TRUE,
      risk.table.col = "strata", # color numbers by group
      risk.table.y.text = FALSE,  # remove y-axis labels
      pval = TRUE,
      pval.coord = c(0.6, 0.6),
      legend.title = "Group",
      ggtheme = theme_light() + theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"),
        plot.caption = element_text(hjust = 0.5, size = 12, face = "italic"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = "top"
      ),
      surv.median.line = "hv",
      conf.int = input$km_show_confint,
      palette = group_palette,
      title = plot_title,
      subtitle = subtitle_text,
      xlab = NULL, # Remove x-axis label from main plot
      ylab = "Survival Probability",
      risk.table.x.text = TRUE, # Show x-axis label below risk table
      risk.table.xlab = "Time (days)", # Add x-axis label only below risk table
    )
    ggsurv$plot <- ggsurv$plot + coord_cartesian(ylim = c(0.5, 1))
    print(ggsurv)
  }, res = 120)
  
  # --- Forest Plot Logic ---
  forest_results <- reactive({
    req(input$forest_endpoint, input$forest_filter_subtypes)
    df <- data_filtered()
    df <- df[df$Prediction.5.Class %in% input$forest_filter_subtypes, , drop = FALSE]
    valid_samples <- df$Sample.ID[df$Sample.ID %in% colnames(expr_mat)]
    df <- df[df$Sample.ID %in% valid_samples, , drop = FALSE]

    endpoint_map <- list(
      "Progression" = list(time_col = "Progression.Time", event_col = "Progression.Event"),
      "Clinical Progression" = list(time_col = "Progression.Clinical.Time", event_col = "Progression.Clinical.Event"),
      "Clinical Progression During Follow-up" = list(time_col = "Progression.During.Follow.up.Time", event_col = "Progression.During.Follow.up.Event"),
      "BCG Any" = list(time_col = "BCG.Any.Time", event_col = "BCG.Any.Event"),
      "BCG Adequate" = list(time_col = "BCG.Adequate.Time", event_col = "BCG.Adequate.Event")
    )
    ep <- endpoint_map[[input$forest_endpoint]]
    if (is.null(ep)) return(NULL)

    surv_time <- suppressWarnings(as.numeric(as.character(df[[ep$time_col]])))
    event_col <- df[[ep$event_col]]
    if (is.factor(event_col) || is.character(event_col)) {
      surv_event <- as.numeric(as.character(event_col))
    } else {
      surv_event <- as.numeric(event_col)
    }
    valid <- !is.na(surv_time) & !is.na(surv_event)
    df <- df[valid, , drop = FALSE]
    surv_time <- surv_time[valid]
    surv_event <- surv_event[valid]

    selected_sigs <- input$forest_signature_select
    if (is.null(selected_sigs) || "All" %in% selected_sigs) {
      sig_indices <- seq_along(forest_metadata_names)
    } else {
      sig_indices <- which(forest_signature_names %in% selected_sigs)
    }

    binsize <- input$forest_binsize
    df_binned <- df
    for (i in sig_indices) {
      sig_col <- forest_metadata_names[i]
      if (sig_col %in% names(df_binned) && is.numeric(df_binned[[sig_col]])) {
        df_binned[[sig_col]] <- int_bin_numeric_variables(data.frame(x = df_binned[[sig_col]]), num_bins = binsize)[[1]]
      }
    }

    results <- lapply(sig_indices, function(i) {
      sig_col <- forest_metadata_names[i]
      sig_label <- forest_signature_names[i]
      if (!sig_col %in% names(df_binned)) {
        cat("[Forest Plot] Signature dropped (missing column):", sig_label, "\n")
        return(NULL)
      }
      sig_vec <- df_binned[[sig_col]]
      if (all(is.na(sig_vec))) {
        cat("[Forest Plot] Signature dropped (all NA):", sig_label, "\n")
        return(NULL)
      }
      # Diagnostic: print sample count for this signature
      cat("[Forest Plot] Signature:", sig_label, "| Sample count:", length(sig_vec), "\n")
      tryCatch({
        cox <- coxph(Surv(surv_time, surv_event) ~ sig_vec)
        s <- summary(cox)
        hr <- s$coefficients[1, "exp(coef)"]
        lower <- s$conf.int[1, "lower .95"]
        upper <- s$conf.int[1, "upper .95"]
        pval <- s$coefficients[1, "Pr(>|z|)"]
        data.frame(Signature = sig_label, HR = hr, Lower = lower, Upper = upper, P = pval)
      }, error = function(e) {
        cat("[Forest Plot] Signature dropped (Cox error):", sig_label, "| Error:", conditionMessage(e), "\n")
        return(NULL)
      })
    })
    forest_df <- do.call(rbind, results)
    # Bonferroni adjustment: multiply p by number of tests (signatures actually tested)
    n_tests <- length(results[!sapply(results, is.null)])
    if (!is.null(forest_df) && nrow(forest_df) > 0) {
      forest_df$Bonferroni_P <- pmin(forest_df$P * n_tests, 1)
    }
    forest_df
  })

  output$forest_table <- DT::renderDataTable({
    forest_df <- forest_results()
    if (is.null(forest_df)) return(NULL)
    num_cols <- c("HR", "Lower", "Upper")
    for (col in num_cols) {
      if (col %in% names(forest_df)) {
        forest_df[[col]] <- formatC(forest_df[[col]], format = "f", digits = 2)
      }
    }
    # Format p-values: scientific notation for small values, 2 digits for others
    format_pval <- function(p) {
      pnum <- suppressWarnings(as.numeric(p))
      if (is.na(pnum)) return(p)
      if (pnum < 0.01) {
        formatC(pnum, format = "e", digits = 2)
      } else {
        formatC(pnum, format = "f", digits = 2)
      }
    }
    if ("P" %in% names(forest_df)) {
      forest_df$P <- vapply(forest_df$P, format_pval, character(1))
    }
    if ("Bonferroni_P" %in% names(forest_df)) {
      forest_df$Bonferroni_P <- vapply(forest_df$Bonferroni_P, format_pval, character(1))
    }
    forest_df
  }, options = list(pageLength = 18))

  output$download_forest_table <- downloadHandler(
    filename = function() {
      paste0("forest_table_", Sys.Date(), ".txt")
    },
    content = function(file) {
      forest_df <- forest_results()
      # Output only raw numeric values with all digits
      num_cols <- c("HR", "Lower", "Upper", "P", "Bonferroni_P")
      for (col in num_cols) {
        if (col %in% names(forest_df)) {
          forest_df[[col]] <- suppressWarnings(as.numeric(forest_df[[col]]))
        }
      }
      write.table(forest_df, file, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  )

  output$forest_plot <- renderPlot({
    forest_df <- forest_results()
    req(forest_df)
    library(ggplot2)
    forest_df$color <- "black"
    forest_df$color[forest_df$P < 0.05] <- "red"
    forest_df$color[forest_df$Bonferroni_P < 0.05] <- "blue"
    all_subtypes <- length(input$forest_filter_subtypes) == length(levels(metadata$Prediction.5.Class))
    # Use sample count from the filtered data
    df <- data_filtered()
    df <- df[df$Prediction.5.Class %in% input$forest_filter_subtypes, , drop = FALSE]
    valid_samples <- df$Sample.ID[df$Sample.ID %in% colnames(expr_mat)]
    df <- df[df$Sample.ID %in% valid_samples, , drop = FALSE]
    endpoint_map <- list(
      "Progression" = list(time_col = "Progression.Time", event_col = "Progression.Event"),
      "Clinical Progression" = list(time_col = "Progression.Clinical.Time", event_col = "Progression.Clinical.Event"),
      "Clinical Progression During Follow-up" = list(time_col = "Progression.During.Follow.up.Time", event_col = "Progression.During.Follow.up.Event"),
      "BCG Any" = list(time_col = "BCG.Any.Time", event_col = "BCG.Any.Event"),
      "BCG Adequate" = list(time_col = "BCG.Adequate.Time", event_col = "BCG.Adequate.Event")
    )
    ep <- endpoint_map[[input$forest_endpoint]]
    # Filter out samples with NA in the event column (and time column)
    if (!is.null(ep)) {
      df <- df[!is.na(df[[ep$event_col]]) & !is.na(df[[ep$time_col]]), , drop = FALSE]
    }
    event_col <- if (!is.null(ep)) df[[ep$event_col]] else NULL
    sample_count <- nrow(df)
    n_events <- if (!is.null(event_col)) sum(event_col == 1, na.rm = TRUE) else NA
    n_censored <- if (!is.null(event_col)) sum(event_col == 0, na.rm = TRUE) else NA
    plot_title <- if (all_subtypes) "Forest Plot: All Subtypes" else paste("Forest Plot:", paste(input$forest_filter_subtypes, collapse = ", "))
    subtitle_text <- paste0(
      "Endpoint: ", input$forest_endpoint,
      " | Samples: ", sample_count,
      " (", n_events, " events, ", n_censored, " censored)"
    )

    legend_labels <- c("p ≥ 0.05", "p < 0.05", "Bonferroni p < 0.05")
    legend_colors <- c("#d35741", "#73bac9", "#499f8b")
    forest_df$legend_group <- factor(
      ifelse(forest_df$Bonferroni_P < 0.05, legend_labels[3],
        ifelse(forest_df$P < 0.05, legend_labels[2], legend_labels[1])
      ),
      levels = legend_labels
    )
    forest_df$legend_color <- as.character(forest_df$color)
    ggplot(forest_df, aes(x = reorder(Signature, HR), y = HR, color = legend_group)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = Lower, ymax = Upper, color = legend_group), width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      coord_flip() +
      scale_color_manual(
        name = "Significance: ",
        values = setNames(legend_colors, legend_labels)
      ) +
      labs(x = "", y = "Hazard Ratio (HR)", title = plot_title, subtitle = subtitle_text) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"),
        legend.position = "top",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 9)
      )
  }, res = 120, height = 700, width = 900)
}

shinyApp(ui, server)