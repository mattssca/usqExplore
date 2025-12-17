library(shiny)
library(plotly)
library(bslib)
library(DT)
library(shinyWidgets)
library(RColorBrewer)
library(pheatmap)
library(later)
library(colourpicker)
library(survival)
library(survminer)

# Load required data
load("data/LU_box/metadata.Rdata")
load("data/LU_box/expr_mat.Rdata")
expr_mat = as.matrix(expr_mat)
load("data/sample_order.Rdata")
names(sample_order) <- colnames(expr_mat)[sample_order]
load("data/genes_to_plot_df.Rdata")
load("data/expr_gene_names.Rdata")
load("data/LU_box/var_categories.Rdata")
load("data/predefined_sample_lists.Rdata")

# Global variables
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
  list(0, "#4DF76F"),
  list(0.5, "black"),
  list(1, "#F74D4D")
)

forest_signature_names <- c(
  "Proliferation Score", "Mol. grade (WHO1999)", "Mol. grade (WHO2022)", "Progression Score", 
  "Immune 141_UP", "NK Cells", "T Cells", "CD8+ T Cells", "Cytotoxicity Score", "B Cells", 
  "Myeloid DCs", "Monocytic Lineage", "Macrophages", "M2 Macrophages", "Neutrophils", 
  "Stromal 141_UP", "Fibroblasts", "Endothelial Cells", "Smooth Muscle"
)

forest_metadata_names <- c(
  "proliferation_score", "molecular_grade_who_1999_score", "molecular_grade_who_2022_score", "progression_score", 
  "immune141_up", "nk_cells", "t_cells", "t_cells_cd8", "cytotoxicity_score", "b_cells", 
  "myeloid_dendritic_cells", "monocytic_lineage", "macrophages", "m2_macrophage", "neutrophils", 
  "stromal141_up", "fibroblasts", "endothelial_cells", "smooth_muscle"
)

# Bin Numeric Variables (internal)
int_bin_numeric_variables <- function(this_data = NULL, num_bins = NULL){
  numeric_columns = sapply(this_data, is.numeric)
  bin_column = function(column, num_bins){
    range_min = min(column, na.rm = TRUE)
    range_max = max(column, na.rm = TRUE)
    binned_column = cut(column, 
                        breaks = seq(range_min, range_max, length.out = num_bins + 1), 
                        labels = FALSE, 
                        include.lowest = TRUE)
    return(binned_column)
  }
  this_data[numeric_columns] <- lapply(this_data[numeric_columns], bin_column, num_bins = num_bins)
  return(this_data)
}