exploration_server <- function(input, output, session, data_filtered) {
  
  # Sample count display
  output$n_samples <- renderText({
    n <- nrow(data_filtered())
    paste("Samples after filtering:", n)
  })
  
  # Main exploration plot
  output$plot <- renderPlotly({
    df <- data_filtered()
    x <- input$xvar
    y <- if (input$yvar != "None") input$yvar else NULL
    color <- if (input$color != "None") input$color else NULL
    show_prop <- isTRUE(input$show_proportion)
    color_map <- if (!is.null(color) && color %in% c("subtype_5_class", "Prediction.7.Class")) lund_colors else NULL
    
    if (input$sort_x && !is.null(x) && !is.null(y)) {
      ord <- order(df[[y]], na.last = TRUE)
      df <- df[ord, , drop = FALSE]
      df[[x]] <- factor(df[[x]], levels = unique(df[[x]]))
    }
    
    hide_x_ticks <- identical(x, "sample_id")
    
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
}