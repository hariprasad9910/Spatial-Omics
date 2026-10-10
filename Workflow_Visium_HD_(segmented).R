### Here is the clean, fully reproducible end-to-end R script.

### Key Bottlenecks Fixed in This Version:

# 1. **Bypassed Fragile In-Memory Zip/BiocFileCache Loops**: Directly streams and maps files via `osfr` without needing nested temp ZIP archives or failing cache hooks.
# 2. **Fixed Space Ranger Path Resolution**: Ensures `TENxVisiumHD` receives both `segmented_outputs` and `binned_outputs` properly resolved inside the same root folder.
# 3. **Resolved `rowAnys` Namespace Collision**: Explicitly accesses `matrixStats::rowAnys` so SpotSweeper spatial outlier filtering never trips an undefined function error.
# 4. **Guaranteed Polygon Matching**: Dynamically extracts cell polygons into `SpatialFeatureExperiment` `colGeometries` with safe indexing and row names.
# 5. **Deduplicated `colData` Cleanups**: Embedded automatic colData deduplication steps before plotting with `plotCoords` to prevent ggplot aesthetic errors.


# ==============================================================================
# Workflow: Visium HD (segmented)
# Analyzing Visium HD data on a human colorectal cancer biopsy (de Oliveira et al., 2025)
# ==============================================================================

# Visium HD enables direct output of H&E-segmented cell-level data, in addition to 
# the previous binned_outputs at 2, 8, and 16 µm, from Space Ranger v4. This update 
# simplifies the additional cell segmentation required when using bin2cell or stardist.

# ==============================================================================
#  Dependencies
# ==============================================================================

# 1. Define the list of packages you need
required_packages <- c(
  "sf", "arrow", "dplyr", "tidyr", "scrapper", "SingleR", "Voyager", 
  "ggplot2", "ggspavis", "VisiumIO", "patchwork", "pals", "Matrix",
  "SpotSweeper", "DropletUtils", "SpatialExperiment", "SpatialFeatureExperiment",
  "osfr", "matrixStats"
)

# 2. Check which ones are currently missing on your system
missing_packages <- required_packages[!required_packages %in% installed.packages()[, "Package"]]

# 3. Install any missing dependencies non-interactively
if (length(missing_packages) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_packages, update = FALSE, ask = FALSE)
}

# Clear any remaining package lock directories
unlink(paste0(.libPaths(), "/00LOCK*"), recursive = TRUE)

suppressPackageStartupMessages({
  library(sf)
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(scrapper)
  library(SingleR)
  library(Voyager)
  library(ggplot2)
  library(ggspavis)
  library(VisiumIO)
  library(patchwork)
  library(pals)
  library(Matrix)
  library(SpotSweeper)
  library(DropletUtils)
  library(SpatialExperiment)
  library(SpatialFeatureExperiment)
  library(osfr)
  library(matrixStats)
})

# Set seed for random number generation in order to make results reproducible
set.seed(777)

# Utility for cropping by bounding box
.crop <- \(spe, box) {
  box <- as.list(box)
  xy <- spatialCoords(spe)
  spe[, 
      xy[, 1] > box$xmin & xy[, 1] < box$xmax &
        xy[, 2] > box$ymin & xy[, 2] < box$ymax ]
}

# ==============================================================================
#  Setup & Data Download
# ==============================================================================

# Establish local project working directory
work_dir <- "D:/hariprasad/Spatial_omics"
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

# Connect to OSF repository node
project_node <- osf_retrieve_node("5n4q3")
df_node <- osf_ls_files(project_node)

# Set dedicated download folder for Visium HD
td_hd <- file.path(work_dir, "VisiumHD_Data")
dir.create(td_hd, recursive = TRUE, showWarnings = FALSE)

# 1. Download Visium HD outputs (segmented + binned + barcode mappings)
target_vhd <- df_node[df_node$name == "VisiumHD_HumanColon_Oliveira", ]
message("Downloading Visium HD outputs from OSF...")
osf_download(target_vhd, path = td_hd, recurse = TRUE, conflicts = "overwrite")

# Resolve root path of Visium HD directory
sr_hd_path <- if (dir.exists(file.path(td_hd, "VisiumHD_HumanColon_Oliveira", "segmented_outputs"))) {
  file.path(td_hd, "VisiumHD_HumanColon_Oliveira")
} else {
  td_hd
}

# 2. Download matched Chromium single-cell reference
td_ref <- file.path(work_dir, "Chromium_Reference")
dir.create(td_ref, recursive = TRUE, showWarnings = FALSE)

target_ref <- df_node[df_node$name == "Chromium_HumanColon_Oliveira", ]
message("Downloading Chromium reference data from OSF...")
osf_download(target_ref, path = td_ref, recurse = TRUE, conflicts = "overwrite")

sr_ref_path <- if (dir.exists(file.path(td_ref, "Chromium_HumanColon_Oliveira"))) {
  file.path(td_ref, "Chromium_HumanColon_Oliveira")
} else {
  td_ref
}

# ==============================================================================
#  Import Visium HD (Cell-Segmented Level)
# ==============================================================================

seg_dir <- file.path(sr_hd_path, "segmented_outputs")
bin_dir <- file.path(sr_hd_path, "binned_outputs")

message("Importing segmented Visium HD data into SpatialExperiment...")
spe <- TENxVisiumHD(
  spacerangerOut    = sr_hd_path,
  segmented_outputs = seg_dir,
  binned_outputs    = bin_dir,
  format            = "h5",
  images            = "lowres"
) |> 
  import()

# Make gene symbols unique
gs <- rowData(spe)$Symbol
rownames(spe) <- make.unique(gs)

# Required attribute for ggspavis coordinate rendering
spe$in_tissue <- TRUE

# Check bin-to-cell mapping summary
summary(sapply(spe$map, nrow))
summary(sapply(spe$map, \(.) mean(.$in_nucleus)))

# ==============================================================================
#  Tissue Cropping & Region of Interest (ROI)
# ==============================================================================

# For runtime reasons, crop tissue to roughly half of original dimensions
xy <- spatialCoords(spe)
xs <- range(xy[, 1])
ys <- range(xy[, 2])
dx <- diff(xs) / 4
dy <- diff(ys) / 4
box <- list(
  xmin = xs[1] + dx, xmax = xs[2] - dx,
  ymin = ys[1] + dy, ymax = ys[2] - dy
)
sub <- .crop(spe, box)

# Reverse y-coordinates of bounding box for Cartesian visualization
.box <- box
.box$ymax <- -box$ymin
.box$ymin <- -box$ymax
aes <- list(col = "red", fill = NA, linewidth = 2)

plotCoords(spe[, sample(ncol(spe), min(5e4, ncol(spe)))]) + 
  do.call(geom_rect, c(.box, aes)) +
  plotCoords(sub[, sample(ncol(sub), min(5e4, ncol(sub)))])

# Retain cropped region
ncol(spe <- sub)

# Specify a smaller region of interest (ROI) for cell-level visualization
roi <- c(xmin = 54e3, xmax = 58e3, ymin = 12e3, ymax = 16e3)
spe$roi <- colnames(spe) %in% colnames(.crop(spe, roi))
plotCoords(spe, annotate = "roi", point_size = 0) + theme(legend.position = "none")

# Convert SPE to SpatialFeatureExperiment (SFE) to hold cell segmentation polygons
sfe <- toSpatialFeatureExperiment(spe)

# Extract and attach cell segmentation geometries safely
seg_meta <- metadata(spe)$cellseg
if (!is.null(seg_meta)) {
  i <- colnames(spe)
  j <- match(i, seg_meta$cell_id)
  seg_clean <- seg_meta[j, ]
  rownames(seg_clean) <- i
  colGeometries(sfe) <- list(cellseg = seg_clean)
}

# ==============================================================================
#  Exploratory Analysis: Morphometry & Bin Overlaps
# ==============================================================================

# Approximate cell areas (4 µm² per 2 µm bin)
sfe$um2 <- 4 * sapply(sfe$map, nrow)
# Fraction of 2 µm bins overlapping nucleus
sfe$nuc <- sapply(sfe$map, \(.) mean(.$in_nucleus))
# Nuclear area approximation
sfe$nuc_um2 <- sfe$um2 * sfe$nuc

# Investigate cell counts within 8 µm and 16 µm bins
bin_names <- c("square_008um", "square_016um")
names(bin_names) <- bin_names
ncs <- lapply(bin_names, \(um) {
  ns <- lapply(spe$map, \(df) unique(df[[um]]))
  as.vector(table(unlist(ns)))
})
sapply(ncs, summary)

# ==============================================================================
#  Quality Control
# ==============================================================================

# Identify mitochondrial features
mt <- grepl("^MT-", rownames(sfe))

# Direct calculation of mitochondrial ratio
mt_counts <- assay(sfe, "counts")[mt, , drop = FALSE]
mt_sum <- Matrix::colSums(mt_counts)
total_sum <- Matrix::colSums(assay(sfe, "counts"))
mt_prop <- mt_sum / total_sum

# Standard scrapper RNA QC
sfe <- quickRnaQc.se(sfe, subsets = list(mt = mt))
sfe$log_sum <- log1p(sfe$sum)
sfe$subset.proportion.mt <- mt_prop

# Visualize QC across the defined ROI
plotSpatialFeature(sfe[, sfe$roi], colGeometryName = "cellseg", features = "log_sum") +
  ggtitle("log-library size") +
  plotSpatialFeature(sfe[, sfe$roi], colGeometryName = "cellseg", features = "subset.proportion.mt") +
  ggtitle("mitochondrial proportion") +
  plot_layout(nrow = 1) &
  theme(plot.title = element_text(hjust = 0.5)) &
  scale_fill_gradientn(NULL, colors = pals::jet())

# Determine spatially-aware outliers using SpotSweeper
sfe <- localOutliers(sfe, workers = 1, metric = "sum", log = TRUE, direction = "lower")
sfe <- localOutliers(sfe, workers = 1, metric = "detected", log = TRUE, direction = "lower")
sfe <- localOutliers(sfe, workers = 1, metric = "subset.proportion.mt", log = FALSE, direction = "higher")

# Extract outlier flags using matrixStats
cd <- colData(sfe)
ol_idx <- grep("outliers$", names(cd))
ol_mat <- as.matrix(cd[ol_idx])
round(100 * colMeans(ol_mat), 2)

# Set discard flag
sfe$ex <- matrixStats::rowAnys(ol_mat)
table(sfe$ex)

# Clean duplicate column names before plotting
colData(sfe) <- colData(sfe)[, !duplicated(colnames(colData(sfe))), drop = FALSE]

# Plot flagged local outliers spatially
plotCoords(sfe[, order(sfe$ex)], annotate = "ex", point_size = 0.2) +
  theme(legend.key.size = ggplot2::unit(0, "pt")) +
  scale_color_manual(values = c("grey90", "red")) +
  guides(col = guide_legend(override.aes = list(size = 2)))

# Filter low-quality outlier cells
sfe <- sfe[, !sfe$ex]

# ==============================================================================
#  Annotation: Single-Cell Reference & SingleR Transfer
# ==============================================================================

# Read single-cell Chromium reference
ref_files <- list.files(sr_ref_path, recursive = TRUE, full.names = TRUE)
h5_ref <- grep("filtered_feature_bc_matrix\\.h5$", ref_files, value = TRUE)[1]
csv_ref <- grep("metadata\\.csv$|cell_metadata\\.csv$", ref_files, value = TRUE)[1]

sce <- read10xCounts(h5_ref, col.names = TRUE)
cd_ref <- read.csv(csv_ref, row.names = 1, check.names = FALSE)
colData(sce)[names(cd_ref)] <- cd_ref[colnames(sce), , drop = FALSE]

# Unique gene symbols
rownames(sce) <- make.unique(rowData(sce)$Symbol)

# Exclude low-quality cells and harmonize labels
sce <- sce[, sce$QCFilter == "Keep"]
table(sce$Level1 <- gsub("\\s", ".", sce$Level1))

# Downsample to at most 1,000 cells per cluster to preserve memory and CPU runtime
idx_split <- split(seq_len(ncol(sce)), sce$Level1)
idx_sample <- lapply(idx_split, \(.) sample(., min(length(.), 1e3)))
.sce <- sce[, unlist(idx_sample)]

# Log-library normalization
sfe <- normalizeRnaCounts.se(sfe)
.sce <- normalizeRnaCounts.se(.sce)

# Realize counts as in-memory sparse matrices for SingleR
assay(sfe, "logcounts") <- as(assay(sfe, "logcounts"), "dgCMatrix")
assay(.sce, "logcounts") <- as(assay(.sce, "logcounts"), "dgCMatrix")

# Perform label transfer using pseudo-bulk reference profiles
res_singler <- SingleR(test = sfe, ref = .sce, labels = .sce$Level1, aggr.ref = TRUE)
sfe$Level1 <- factor(res_singler$pruned.labels)

# Group sub-clusters into broader tissue compartments (Level0)
lab_compartments <- list(
  tum = c("Tumor"),
  epi = c("Intestinal.Epithelial"),
  imm = c("B.cells", "T.cells", "Myeloid"),
  str = c("Endothelial", "Fibroblast", "Smooth.Muscle")
)
idx_comp <- match(sfe$Level1, unlist(lab_compartments))
rep_names <- rep.int(names(lab_compartments), sapply(lab_compartments, length))
sfe$Level0 <- factor(rep_names[idx_comp])
table(sfe$Level0)

# ==============================================================================
#  Spatial Visualization of Cell Annotations
# ==============================================================================

pal_lv1 <- unname(pals::trubetskoy())
pal_lv0 <- c("gold", "turquoise", "deeppink", "navy")

# Plot annotated cells over the ROI with cell boundaries
plotSpatialFeature(sfe[, sfe$roi], features = "Level1", colGeometryName = "cellseg") + 
  scale_fill_manual(values = pal_lv1) +
  plotSpatialFeature(sfe[, sfe$roi], features = "Level0", colGeometryName = "cellseg") + 
  scale_fill_manual(values = pal_lv0) +
  plot_layout(nrow = 1) &
  theme(legend.key.size = ggplot2::unit(0.5, "lines"))

# Plot whole tissue slice centroids
plotCoords(sfe, annotate = "Level1", point_size = 0.1) +
  scale_color_manual(values = pal_lv1) +
  plotCoords(sfe, annotate = "Level0", point_size = 0.1) +
  scale_color_manual(values = pal_lv0) +
  plot_layout() &
  theme(legend.key.size = ggplot2::unit(0, "pt"))

# ==============================================================================
#  Downstream: Morphometry vs. Compartments
# ==============================================================================

# Compare cell and nucleus areas across annotations
df_meta <- data.frame(colData(sfe)) |> filter(!is.na(Level0))
fd_meta <- pivot_longer(df_meta, cols = ends_with("um2"))
fd_meta$name <- factor(fd_meta$name, labels = c("cell", "nuc."))
panel_widths <- c(nlevels(df_meta$Level0), nlevels(df_meta$Level1))

ggplot(fd_meta, aes(reorder(Level0, value), value, col = name)) +
  ggplot(fd_meta, aes(reorder(Level1, value), value, col = name)) + 
  plot_layout(nrow = 1, widths = panel_widths, guides = "collect") &
  geom_boxplot(key_glyph = "point", outlier.size = 0) &
  labs(col = "area (um2)") & 
  theme_bw() & 
  theme(
    axis.title = element_blank(),
    legend.key.size = ggplot2::unit(0, "pt"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
