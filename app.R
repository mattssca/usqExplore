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
    "))
  ),
  fluidRow(
    column(
      width = 12,
      div("UROSCANSEQ Explorer", class = "main-title"),
      helpText("Use filters to subset your data. Select variables and plot type below.")
    )
  ),
  sidebarLayout(
    sidebarPanel(width = 5,
                 wellPanel(
                   div("Metadata Filters", class = "filter-title"),
                   actionBttn("add_filter", "Add Filter", style = "jelly", color = "default"),
                   tags$div(id = "filter_ui_container")
                 ),
                 hr(),
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
                 prettySwitch("show_proportion", "Show proportions", value = FALSE),
                 downloadBttn("download_ids", "Download Sample IDs", style = "gradient", color = "primary")
    ),
    mainPanel(width = 7,
              tabsetPanel(
                tabPanel("Plot",
                         div(class = "plot-container", plotlyOutput("plot")),
                         div(class = "after-plot-text", textOutput("n_samples"))
                ),
                tabPanel("Table",
                         DT::dataTableOutput("table")
                ),
                tabPanel("Heatmap",
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
                         actionButton("add_metagene", "Add Meta-gene"),
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
                             actionButton("plot_heatmap", "Plot Heatmap", width = "100%"),
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
                tabPanel("Survival",
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
                           selected = "Progression"
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
                         checkboxInput("km_show_confint", "Show confidence intervals", value = TRUE),
                         plotOutput("km_plot", height = "700px")
                )
              )
    )
  )
)

server <- function(input, output, session) {
  filter_count <- reactiveVal(0)
  filter_ids <- reactiveVal(character(0))
  
  observeEvent(input$add_filter, {
    # Save all current filter input values BEFORE UI changes
    current_ids <- filter_ids()
    current_vals <- lapply(current_ids, function(id) {
      list(
        cat = input[[paste0(id, "_cat")]],
        var = input[[paste0(id, "_var")]],
        val = input[[paste0(id, "_val")]]
      )
    })
    names(current_vals) <- current_ids
    
    # Insert new filter UI
    id <- paste0("filter_", filter_count() + 1)
    insertUI(
      selector = "#filter_ui_container",
      ui = tags$div(
        id = id,
        fluidRow(
          column(3, pickerInput(paste0(id, "_cat"), "Category", choices = names(var_categories), width = "100%")),
          column(5, uiOutput(paste0(id, "_var_ui"))),
          column(4, uiOutput(paste0(id, "_val_ui"))),
          column(2, actionBttn(paste0(id, "_remove"), "Remove", style = "simple", color = "danger", size = "sm"))
        )
      ),
      immediate = TRUE
    )
    filter_count(filter_count() + 1)
    filter_ids(c(filter_ids(), id))
    
    # Restore previous filter input values after UI update
    later::later(function() {
      for (fid in names(current_vals)) {
        vals <- current_vals[[fid]]
        # Restore cat and var immediately
        if (!is.null(vals$cat)) {
          updatePickerInput(session, paste0(fid, "_cat"), selected = vals$cat)
        }
        if (!is.null(vals$var)) {
          updatePickerInput(session, paste0(fid, "_var"), selected = vals$var)
        }
      }
      # Restore val after a further delay to ensure UI is rendered
      later::later(function() {
        for (fid in names(current_vals)) {
          vals <- current_vals[[fid]]
          if (!is.null(vals$val)) {
            # Check if the UI is a slider or picker
            if (is.numeric(vals$val) && length(vals$val) == 2) {
              updateSliderInput(session, paste0(fid, "_val"), value = vals$val)
            } else {
              updatePickerInput(session, paste0(fid, "_val"), selected = vals$val)
            }
          }
        }
      }, delay = 0.1)
    }, delay = 0.1)
  })
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
  }, options = list(pageLength = 7))
  
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

    group_palette <- RColorBrewer::brewer.pal(4, "Set1")

    ggsurv <- survminer::ggsurvplot(
      fit,
      data = plot_df,
      risk.table = TRUE,
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
      xlab = "Time (days)",
      ylab = "Survival Probability"
    )
    ggsurv$plot <- ggsurv$plot + coord_cartesian(ylim = c(0.5, 1))
    print(ggsurv)
  }, res = 120)
}

shinyApp(ui, server)