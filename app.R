library(shiny)
library(plotly)
library(bslib)
library(DT)
library(shinyWidgets)
library(RColorBrewer)
library(pheatmap)
library(later)
library(colourpicker)

# Load metadata
load("uroscanseq_meta.Rdata")
metadata = uroscanseq_meta

# Load expression data (should create expr_mat)
load("expr_met_uroscanseq.Rdata")
expr_mat = uroscanseq_expr

# Load sample order vector
load("sample_order.Rdata")
names(sample_order) <- colnames(expr_mat)[sample_order]

# Load predefined gene lists
load("genes_to_plot_df.Rdata") # Should create genes_to_plot_df with columns: gene, signature

# Load gene names for fast selectize
load("exp_genes.Rdata") # exp_genes should be a character vector of gene names

# Define your color palette
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

# Set levels for Predictions_5classes
if ("subtype_5_class" %in% names(metadata)) {
  metadata$subtype_5_class <- factor(
    metadata$subtype_5_class,
    levels = c("Uro", "GU", "BaSq", "Mes", "ScNE")
  )
}

# Set levels for Predictions_7classes
if ("subtype_7_class" %in% names(metadata)) {
  metadata$subtype_7_class <- factor(
    metadata$subtype_7_class,
    levels = c("UroA", "UroB", "UroC", "GU", "BaSq", "Mes", "ScNE")
  )
}

# Variable categories
var_categories <- list(
  `Clinical Data` = c("clin_prog_tat1_event", "clin_prog_tat1_time", "clin_prog_pat_tat1_event", "clin_prog_pat_tat1_time", "clin_clinprog_fu_tat1_event", "clin_clinprog_fu_tat1_time", "bcg_any_event", "bcg_any_time", "bcg_adequate_event", "bcg_adequate_time", "gem_mmc_event", "gem_mmc_time", "time_to_primary_cystectomy", "progression_event", "time_to_progression", "primary_cystectomy_event_nmi", "all_cystectomy_event_nmi", "time_to_cystectomy_nmi", "recurrence_event", "time_to_recurrence", "primary_recurrence", "palliative", "pdd_turb_returb", "urine_cytology_pre_turb", "tnm", "primary_cystectomy", "neoadj_induction", "preop_chemo_type", "n_chemo_doses", "returb", "pad_returb_t_stage", "prostatic_urethra", "lvi", "t1_invasion_depth", "variant_histology", "bcg_instillations", "adjuvant_mmc", "only_palliative_treatment", "ypt", "ypn", "tma_grade_who16_as"),
  `Cohorts` = c("sample_id", "cohort", "seq_batch", "cohort_group"),
  `Descriptive` = c("age", "eau_risk_category", "eau_score", "hospital_turb", "gender", "smoker", "death_other_cause", "death_bladder_cancer"),
  `EAU Factors` = c("eau_over70", "eau_tumor_status", "eau_n_tumors", "eau_tumor_diam", "eau_stage", "eau_cis", "eau_who73", "eau_var_hist", "eau_prostatic_urethra", "eau_lvi"),
  `LundTax Signatures` = c("proliferation_score", "progression_score", "progression_risk", "mol_grade_2022_score", "mol_grade_2022", "mol_grade_1999_score", "mol_grade_1999", "stromal141_up", "immune141_up", "b_cells", "t_cells", "cd8_t_cells", "nk_cells", "t_cells", "neutrophils", "monocytic_lineage", "macrophages", "m2_macrophage", "myeloid_dc", "cytotoxicity_score", "endothelial_cells", "fibroblasts", "smooth_muscle"),
  `Pathology` = c("multi_hist"),
  `Sequencing Metrics` = c("dna_ngul", "rna_ngul", "rin"),
  `Sample Sets` = c("set_719", "set_676", "set_hq_572", "set_hq_572_index_uc_533", "set_lq_147", "set_lq_147_index_uc_129", "set_rna_tma", "set_eau_risk", "set_clin_prog_tat1", "set_clin_prog_ta", "set_clin_prog_t1", "set_clin_prog_pat_tat1", "set_clin_prog_pat_ta", "set_clin_prog_pat_t1", "set_clin_clinprog_fu_tat1", "set_clin_clinprog_fu_ta", "set_clin_clinprog_fu_t1", "set_bcg_any", "set_bcg_adequate", "set_gem_mmc"),
  `Subtype Predictions` = c("subtype_5_class", "subtype_7_class", "Uro_score", "UroA_score", "UroB_score", "UroC_score", "GU_score", "BaSq_score", "Mes_score", "ScNE_score"),
  `Tumor Information` = c("stage_simp", "stage", "node", "met", "grade", "n_tumors_cat", "tumor_size"),
  `Quality Control` = c("qc_evaluation", "sample_category", "category_group")
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
                 pickerInput("xvar", "X variable", choices = var_categories, selected = "subtype_5_class", 
                             options = list(`live-search` = TRUE)),
                 pickerInput("yvar", "Y variable (optional)", choices = c("None", var_categories), selected = "None", 
                             options = list(`live-search` = TRUE)),
                 pickerInput("color", "Color by (optional)", choices = c("None", var_categories), selected = "subtype_5_class", 
                             options = list(`live-search` = TRUE)),
                 radioGroupButtons("plot_type", "Plot type",
                                   choices = c("Histogram" = "hist", "Barplot" = "bar", "Boxplot" = "box", "Scatterplot" = "scatter"),
                                   selected = "hist", justified = TRUE
                 ),
                 prettySwitch("sort_x", "Sort/Ranks by X variable", value = FALSE, status = "info"),
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
                           choices = setdiff(names(metadata), "subtype_5_class"),
                           multiple = TRUE,
                           options = list(`live-search` = TRUE)
                         ),
                         uiOutput("ann_color_pickers"), # <-- Add this line
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
                )
              )
    )
  )
)

server <- function(input, output, session) {
  filter_count <- reactiveVal(0)
  filter_ids <- reactiveVal(character(0))
  
  observeEvent(input$add_filter, {
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
    
    if (input$sort_x && !is.null(x) && is.numeric(df[[x]])) {
      df <- df[order(df[[x]], decreasing = TRUE), , drop = FALSE]
      if (input$plot_type %in% c("bar", "box", "scatter")) {
        df[[x]] <- factor(df[[x]], levels = unique(df[[x]]))
      }
    }
    
    color_map <- if (!is.null(color) && color %in% c("subtype_5_class", "subtype_7_class")) lund_colors else NULL
    
    if (input$plot_type == "hist") {
      plot_ly(
        df, x = ~get(x),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "histogram"
      )
    } else if (input$plot_type == "bar") {
      plot_ly(
        df, x = ~get(x),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "bar"
      )
    } else if (input$plot_type == "box") {
      plot_ly(
        df, x = ~get(x), y = ~get(y),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "box"
      )
    } else if (input$plot_type == "scatter") {
      plot_ly(
        df, x = ~get(x), y = ~get(y),
        color = if (!is.null(color)) ~get(color) else NULL,
        colors = color_map,
        type = "scatter", mode = "markers"
      )
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
      writeLines(as.character(data_filtered()$sample_id), file)
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
    sample_ids <- data_filtered()$sample_id
    valid_samples <- sample_ids[sample_ids %in% colnames(expr_mat)]
    if (length(valid_samples) == 0) {
      output$heatmap_plot <- renderPlotly({ plotly_empty() })
      output$heatmap_legends <- renderUI({ NULL })
      return()
    }
    ann <- data_filtered()
    ann <- ann[ann$sample_id %in% valid_samples, ]
    ann <- ann[order(ann$subtype_5_class), ]
    
    # --- Sorting logic ---
    if (!is.null(input$sort_heatmap_samples) && input$sort_heatmap_samples && exists("sample_order")) {
      ann$sample_order_val <- sample_order[as.character(ann$sample_id)]
      ann <- ann[order(ann$subtype_5_class, ann$sample_order_val, na.last = TRUE), ]
    }
    sorted_samples <- ann$sample_id
    mat <- expr_mat[selected_genes_heatmap(), sorted_samples, drop = FALSE]
    
    # --- Add meta-genes as annotation tracks ---
    mg_exprs <- metagene_exprs()
    ann_tracks <- list()
    legends <- list()
    subtype_levels <- levels(metadata$subtype_5_class)
    subtype_colors <- lund_colors[subtype_levels]
    # Always include subtype_5_class annotation track
    ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
      z = matrix(as.numeric(factor(ann$subtype_5_class, levels = subtype_levels)), nrow = 1),
      x = sorted_samples,
      y = "Subtype 5 Class",
      type = "heatmap",
      showscale = FALSE,
      colors = subtype_colors,
      hoverinfo = "text",
      text = paste("Subtype:", ann$subtype_5_class)
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
}

shinyApp(ui, server)