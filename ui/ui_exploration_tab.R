exploration_tab_ui <- function() {
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
                    prettySwitch("show_proportion", "Show proportions", value = FALSE)
             ),
             column(7,
                    div(class = "plot-container", plotlyOutput("plot", height = "565px")),
                    div(class = "after-plot-text", textOutput("n_samples"))
             )
           )
  )
}
