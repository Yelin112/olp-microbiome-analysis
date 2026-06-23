# =============================================================================
# Gene Expression Heatmap
# Input : ../fpkm-2.xlsx        (FPKM expression matrix)
#         ../group (1).xlsx     (sample → group annotation)
# Output: gene_expression_heatmap.pdf / .pptx
# =============================================================================

required_packages <- c("tidyverse", "pheatmap", "readxl", "officer", "rvg")
new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(new_packages)
}

library(tidyverse)
library(pheatmap)
library(readxl)
library(officer)
library(rvg)

# =============================================================================
# Paths
# =============================================================================
data_dir <- "E:\\打工人\\Zeng\\2023-3 OLP Microbiome\\分析\\help_others\\LHY\\Other"
output_dir <- file.path(data_dir, "heatmap")

# =============================================================================
# Color palettes
# =============================================================================
# Heatmap gradient: blue → white → red
heat_colors <- colorRampPalette(c("#20a0e5", "white", "#e73838"))(100)

# Group annotation colors (consistent with CIBERSORTx analysis)
group_colors <- c(
  "Control" = "#888888",
  "Ble" = "#55b7ec",
  "Are+Ble" = "#ee6161",
  "Are" = "#4fb984",
  "IA" = "#ed9548",
  "MA" = "#fdd456",
  "MP" = "#6358cb"
)

# PPTX helper (base R / grid graphics)
save_pptx_base <- function(plot_expr, file_path) {
  doc <- officer::read_pptx() |>
    officer::add_slide(layout = "Blank", master = "Office Theme") |>
    officer::ph_with(
      value = rvg::dml(code = plot_expr),
      location = officer::ph_location_fullsize()
    )
  print(doc, target = file_path)
  invisible(file_path)
}

# =============================================================================
# 1. Load data
# =============================================================================
expr_raw <- read_excel(file.path(data_dir, "fpkm-2.xlsx")) |>
  column_to_rownames("gene_id") |>
  mutate(across(everything(), as.numeric))

group_df <- read_excel(file.path(data_dir, "group (1).xlsx")) |>
  column_to_rownames("Sample")

# =============================================================================
# 2. Order samples by group, then Z-score normalize rows
# =============================================================================
group_order <- c("Control", "Ble", "Are+Ble", "Are", "IA", "MA", "MP")
sample_order <- rownames(group_df)[
  order(match(group_df$Group, group_order))
]

expr_ordered <- as.matrix(expr_raw[, sample_order])
mat_scaled <- t(scale(t(expr_ordered))) # row-wise Z-score

# Clip extreme values for display
mat_scaled <- pmax(pmin(mat_scaled, 2), -2)

# =============================================================================
# 3. Column annotation
# =============================================================================
ann_col <- group_df[sample_order, , drop = FALSE]
ann_colors <- list(Group = group_colors)

# =============================================================================
# 4. Shared pheatmap arguments (square cells: cellwidth = cellheight)
# =============================================================================
cell_size <- 22 # pt; same value ensures square cells

heatmap_args <- list(
  mat = mat_scaled,
  color = heat_colors,
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_method = "ward.D2",
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  cellwidth = cell_size,
  cellheight = cell_size,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9,
  border_color = "white",
  show_colnames = TRUE,
  show_rownames = TRUE,
  angle_col = 45,
  legend_breaks = c(-2, -1, 0, 1, 2),
  legend_labels = c("-2", "-1", "0", "1", "2"),
  main = "Gene Expression Heatmap (Z-score)",
  silent = TRUE
)

# =============================================================================
# 5. Draw and save
# =============================================================================
ph <- do.call(pheatmap, heatmap_args)

# --- PDF ---
pdf_path <- file.path(output_dir, "gene_expression_heatmap.pdf")
pdf(pdf_path, width = 14, height = 7)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()
cat("✓ PDF saved:", pdf_path, "\n")

# --- PPTX ---
pptx_path <- file.path(output_dir, "gene_expression_heatmap.pptx")
save_pptx_base(
  {
    ph_pptx <- do.call(pheatmap, heatmap_args)
    grid::grid.newpage()
    grid::grid.draw(ph_pptx$gtable)
  },
  pptx_path
)
cat("✓ PPTX saved:", pptx_path, "\n")

cat("\n✓ Done.\n")
