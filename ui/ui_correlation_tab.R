correlation_tab_ui <- function() {
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
  )
}