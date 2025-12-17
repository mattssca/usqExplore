survival_server <- function(input, output, session, data_filtered) {
  
  # Keep track of previous gene choices
  prev_gene_choices <- reactiveVal(NULL)
  
  # Update gene choices with meta-genes
  observe({
    # This would need access to metagene_exprs from heatmap module
    # For now, just use expr_gene_names
    gene_choices <- expr_gene_names
    current_sel <- input$km_gene
    valid_sel <- current_sel[current_sel %in% gene_choices]
    if (!identical(prev_gene_choices(), gene_choices)) {
      updateSelectizeInput(session, "km_gene", choices = gene_choices, selected = if (length(valid_sel) > 0) valid_sel else NULL)
      prev_gene_choices(gene_choices)
    }
  })
  
  # KM Survival Plot
  output$km_plot <- renderPlot({
    gene_selected <- input$km_gene
    cut_method <- input$km_cut
    endpoint <- input$km_endpoint
    filter_subtypes <- input$km_filter_subtypes
    genes_all <- expr_gene_names
    
    req(gene_selected, cut_method, endpoint, filter_subtypes)
    
    df <- data_filtered()
    gene <- gene_selected
    
    # Filter samples by selected subtypes
    df <- df[df$subtype_5_class %in% filter_subtypes, , drop = FALSE]
    
    # Endpoint mapping (adjust column names as needed)
    endpoint_map <- list(
      "Biological Progression" = list(time_col = "progression_biological_time", event_col = "progression_biological_event"),
      "Clinical Progression" = list(time_col = "progression_clinical_time", event_col = "progression_clinical_event"),
      "Clinical Progression During Follow-up" = list(time_col = "progression_fu_time", event_col = "progression_fu_event"),
      "BCG Any" = list(time_col = "bcg_any_time", event_col = "bcg_any_event"),
      "BCG Adequate" = list(time_col = "bcg_adequate_time", event_col = "bcg_adequate_event"),
      "Recurrance" = list(time_col = "recurrence_time", event_col = "recurrence_event")
    )
    
    ep <- endpoint_map[[endpoint]]
    if (is.null(ep)) return(NULL)
    
    # Filter out samples with NA in the event column
    df <- df[!is.na(df[[ep$event_col]]), , drop = FALSE]
    
    valid_samples <- df$sample_id[df$sample_id %in% colnames(expr_mat)]
    
    # Get expression data - check if it's a gene in expr_mat
    if (gene %in% rownames(expr_mat)) {
      expr <- expr_mat[gene, valid_samples]
    } else {
      return(NULL)
    }
    
    df <- df[df$sample_id %in% valid_samples, , drop = FALSE]
    
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
    subtype <- droplevels(df$subtype_5_class[valid])
    
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
    if (length(filter_subtypes) == length(levels(df$subtype_5_class))) {
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
}