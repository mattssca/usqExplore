forest_server <- function(input, output, session, data_filtered) {
  
  # Forest plot results reactive
  forest_results <- reactive({
    req(input$forest_endpoint, input$forest_filter_subtypes)
    df <- data_filtered()
    df <- df[df$subtype_5_class %in% input$forest_filter_subtypes, , drop = FALSE]
    valid_samples <- df$sample_id[df$sample_id %in% colnames(expr_mat)]
    df <- df[df$sample_id %in% valid_samples, , drop = FALSE]
    
    endpoint_map <- list(
      "Biological Progression" = list(time_col = "progression_biological_time", event_col = "progression_biological_event"),
      "Clinical Progression" = list(time_col = "progression_clinical_time", event_col = "progression_clinical_event"),
      "Clinical Progression During Follow-up" = list(time_col = "progression_fu_time", event_col = "progression_fu_event"),
      "BCG Any" = list(time_col = "bcg_any_time", event_col = "bcg_any_event"),
      "BCG Adequate" = list(time_col = "bcg_adequate_time", event_col = "bcg_adequate_event")
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
  
  # Forest plot table
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
  
  # Download forest table
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
  
  # Forest plot
  output$forest_plot <- renderPlot({
    forest_df <- forest_results()
    req(forest_df)
    library(ggplot2)
    forest_df$color <- "black"
    forest_df$color[forest_df$P < 0.05] <- "red"
    forest_df$color[forest_df$Bonferroni_P < 0.05] <- "blue"
    all_subtypes <- length(input$forest_filter_subtypes) == length(levels(metadata$subtype_5_class))
    # Use sample count from the filtered data
    df <- data_filtered()
    df <- df[df$subtype_5_class %in% input$forest_filter_subtypes, , drop = FALSE]
    valid_samples <- df$sample_id[df$sample_id %in% colnames(expr_mat)]
    df <- df[df$sample_id %in% valid_samples, , drop = FALSE]
    
    endpoint_map <- list(
      "Biological Progression" = list(time_col = "progression_biological_time", event_col = "progression_biological_event"),
      "Clinical Progression" = list(time_col = "progression_clinical_time", event_col = "progression_clinical_event"),
      "Clinical Progression During Follow-up" = list(time_col = "progression_fu_time", event_col = "progression_fu_event"),
      "BCG Any" = list(time_col = "bcg_any_time", event_col = "bcg_any_event"),
      "BCG Adequate" = list(time_col = "bcg_adequate_time", event_col = "bcg_adequate_event"),
      "Recurrance" = list(time_col = "recurrence_time", event_col = "recurrence_event")
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