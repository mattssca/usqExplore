heatmap_server <- function(input, output, session, data_filtered) {
  
  # Meta-gene reactive values
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
  
  # Add meta-gene
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
    genes <- expr_gene_names
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
        filtered_genes <- expr_gene_names
        if (!is.null(filter_val) && nzchar(filter_val)) {
          filtered_genes <- grep(filter_val, expr_gene_names, value = TRUE, ignore.case = TRUE)
        }
        updateSelectizeInput(session, paste0(mg_id, "_genes"), choices = filtered_genes, server = TRUE)
      })
    }
  })
  
  # Meta-gene generation
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
    genes <- expr_gene_names
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
    genes <- expr_gene_names
    if (!is.null(input$gene_signature) && input$gene_signature == "Custom Genes") {
      selectizeInput("genes", "Select Genes for Heatmap (manual entry)", choices = genes, multiple = TRUE, selected = NULL, options = list(server = TRUE))
    }
  })
  
  # Annotation color pickers UI
  output$ann_color_pickers <- renderUI({
    req(input$heatmap_ann_vars)
    out <- list()
    for (var in input$heatmap_ann_vars) {
      # Hide color picker for predefined subtypes
      if (var %in% c("subtype_5_class", "subtype_7_class")) next
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
  
  # Main heatmap plotting
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
    subtype_levels <- levels(metadata$subtype_5_class)
    ann$subtype_5_class <- factor(ann$subtype_5_class, levels = subtype_levels)
    ann <- ann[order(ann$subtype_5_class), ]
    
    # Sorting logic
    if (!is.null(input$sort_heatmap_samples) && input$sort_heatmap_samples) {
      # Use the proliferation_score column from ann
      sample_order <- ann$proliferation_score
      names(sample_order) <- ann$sample_id
      ann$sample_order_val <- sample_order[as.character(ann$sample_id)]
      ann <- ann[order(ann$subtype_5_class, ann$sample_order_val, na.last = TRUE), ]
    }
    sorted_samples <- ann$sample_id
    mat <- expr_mat[selected_genes_heatmap(), sorted_samples, drop = FALSE]
    
    # Add meta-genes as annotation tracks
    mg_exprs <- metagene_exprs()
    ann_tracks <- list()
    legends <- list()
    
    # Toggle between 5-class and 7-class subtype annotation tracks
    if (isTRUE(input$main_heatmap_subtype)) {
      # 7-class annotation track
      subtype7_levels <- levels(metadata$subtype_7_class)
      present7_levels <- subtype7_levels[subtype7_levels %in% ann$subtype_7_class]
      ann$subtype_7_class <- factor(ann$subtype_7_class, levels = present7_levels)
      ann <- ann[order(ann$subtype_7_class), ]
      subtype7_colorscale <- lapply(seq_along(present7_levels), function(i) {
        list((i-1)/(length(present7_levels)-1), lund_colors[present7_levels[i]])
      })
      ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
        z = matrix(as.numeric(ann$subtype_7_class), nrow = 1),
        x = sorted_samples,
        y = "Subtype 7 Class",
        type = "heatmap",
        showscale = FALSE,
        colorscale = subtype7_colorscale,
        hoverinfo = "text",
        text = paste("Subtype:", ann$subtype_7_class)
      )
      legends[[length(legends) + 1]] <- tags$div(
        tags$b("Subtype 7 Class:"),
        tags$ul(
          lapply(present7_levels, function(lvl) {
            tags$li(
              tags$span(style = paste0("display:inline-block;width:15px;height:15px;background:", lund_colors[lvl], ";margin-right:5px;"), ""),
              lvl
            )
          })
        )
      )
    } else {
      # 5-class annotation track
      subtype_levels <- levels(metadata$subtype_5_class)
      present_levels <- subtype_levels[subtype_levels %in% ann$subtype_5_class]
      ann$subtype_5_class <- factor(ann$subtype_5_class, levels = present_levels)
      ann <- ann[order(ann$subtype_5_class), ]
      subtype_colorscale <- lapply(seq_along(present_levels), function(i) {
        list((i-1)/(length(present_levels)-1), lund_colors[present_levels[i]])
      })
      ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
        z = matrix(as.numeric(ann$subtype_5_class), nrow = 1),
        x = sorted_samples,
        y = "Subtype 5 Class",
        type = "heatmap",
        showscale = FALSE,
        colorscale = subtype_colorscale,
        hoverinfo = "text",
        text = paste("Subtype:", ann$subtype_5_class)
      )
      legends[[length(legends) + 1]] <- tags$div(
        tags$b("Subtype 5 Class:"),
        tags$ul(
          lapply(present_levels, function(lvl) {
            tags$li(
              tags$span(style = paste0("display:inline-block;width:15px;height:15px;background:", lund_colors[lvl], ";margin-right:5px;"), ""),
              lvl
            )
          })
        )
      )
    }
    
    # Add meta-gene annotation tracks
    if (length(mg_exprs) > 0) {
      for (mg_name in names(mg_exprs)) {
        mg_vec <- mg_exprs[[mg_name]]
        # Ensure meta-gene vector is named by sample IDs
        if (is.null(names(mg_vec)) || !all(sorted_samples %in% names(mg_vec))) {
          names(mg_vec) <- colnames(expr_mat)
        }
        mg_vec <- mg_vec[sorted_samples]
        ann_tracks[[length(ann_tracks) + 1]] <- plot_ly(
          z = matrix(mg_vec, nrow = 1),
          x = sorted_samples,
          y = mg_name,
          type = "heatmap",
          colorscale = custom_colorscale,
          zmin = -2,
          zmax = 2,
          showscale = FALSE,
          hoverinfo = "text",
          text = paste(mg_name, ":", signif(mg_vec, 3))
        )
      }
    }
    
    # Add user-selected annotation tracks
    ann_vars <- input$heatmap_ann_vars
    if (!is.null(ann_vars) && length(ann_vars) > 0) {
      for (i in seq_along(ann_vars)) {
        var <- ann_vars[i]
        vals <- as.character(ann[[var]])
        
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