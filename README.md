<img src="man/figures/logo.png" align="right" height="280" alt="usqExplore logo" />

# usqExplore

Welcome to the UROSCANSEQ Explorer Shiny App!

## What does the app do?

UROSCANSEQ Explorer is an interactive Shiny application for exploring gene expression, clinical metadata, and survival analysis in urothelial cancer. The app provides:

- Dynamic heatmaps of gene expression and meta-genes
- Interactive filtering and annotation of samples by clinical and molecular features
- Survival analysis (Kaplan-Meier plots) by gene, meta-gene, or LundTax signature
- Forest plots of hazard ratios for LundTax signatures using Cox models
- Downloadable sample lists and tables

## How to Download and Run Locally

### Step-by-step guide

1. **Clone the Repository**
Clone the repository to your local machine:
```bash
git clone https://github.com/mattssca/usqUtils.git
```
2. **Download data**
   - Log in to LU Box and locate the associated data folder for this package `(/Git/data/usqExplore)`.
   - Download the contents of the LU Box data folder.
   - Add the folder (with its contents) to the cloned version of usqExplore's data fodler `usqExplore/data/`
   
3. **Install required R packages**
In the RStudio Console, run:
     ```r
     install.packages(c("shiny", "plotly", "bslib", "DT", "shinyWidgets", "RColorBrewer", "pheatmap", "later", "colourpicker", "survival", "survminer"))
     ```
6. **Run the app**
   - Set your working directory to the cloned repo.
   - In the RStudio Console, run:
     ```r
     shiny::runApp('app.R')
     ```
   - The app will open in your web browser

## Need help?

If you have any issues or questions, please open an issue on GitHub or contact the maintainer.

Enjoy exploring your data!
