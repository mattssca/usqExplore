correlation_server <- function(input, output, session, data_filtered) {
  
  # Correlation scatter plot logic
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
}