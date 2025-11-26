<img src="man/figures/logo.png" align="right" height="280" alt="UroscanseqExplorer logo" />

# UROSCANSEQ Explorer
Welcome to the UROSCANSEQ Explorer Shiny App!

## What does the app do?

UROSCANSEQ Explorer is an interactive Shiny application for exploring gene expression, clinical metadata, and survival analysis in urothelial cancer. The app provides:
- Dynamic heatmaps of gene expression and meta-genes
- Interactive filtering and annotation of samples by clinical and molecular features
- Survival analysis (Kaplan-Meier plots) by gene, meta-gene, or LundTax signature
- Forest plots of hazard ratios for LundTax signatures using Cox models
- Downloadable sample lists and tables

## Data Used

The app uses the following data files (included in the repository):
- `uroscanseq_meta.Rdata`: Clinical and molecular metadata for all samples
- `expr_met_uroscanseq.Rdata`: Gene expression matrix
- `genes_to_plot_df.Rdata`: Predefined gene signatures
- `exp_genes.Rdata`: List of available genes
- `var_categories.Rdata`: Metadata variable categories
- Additional files for sample order and filtering

## How to Download and Run Locally

### Step-by-step guide

1. **Download the repository**
   - Click the green "Code" button on the GitHub page and select "Download ZIP"
   - Unzip the folder to your computer

2. **Open RStudio**
   - Launch RStudio on your computer

3. **Open the project folder**
   - Go to `File` > `Open Project...` or `Open Folder...`
   - Navigate to the unzipped `uroscanseq_exploration` folder and open it

4. **Install required R packages**
   - In the RStudio Console, run:
     ```r
     install.packages(c("shiny", "plotly", "bslib", "DT", "shinyWidgets", "RColorBrewer", "pheatmap", "later", "colourpicker", "survival", "survminer"))
     ```

5. **Run the app**
   - In the RStudio Console, run:
     ```r
     shiny::runApp('app.R')
     ```
   - The app will open in your web browser

## Need help?
If you have any issues or questions, please open an issue on GitHub or contact the maintainer.

Enjoy exploring your data!
