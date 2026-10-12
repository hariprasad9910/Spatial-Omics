# ==============================================================================
# Complete End-to-End Visium HD (8 µm & 16 µm) Analysis Pipeline
# Human Colorectal Cancer Biopsy (de Oliveira et al., 2025)
# ==============================================================================

# ==============================================================================
# 1. System Dependencies and Individual Package Diagnostics
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cran.rstudio.com")

# Isolate critical dependency installations for testing
critical_pkgs <- c("magick", "leidenAlg")
for (pkg in critical_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("=== Debugging installation of critical package:", pkg, "==="))
    tryCatch({
      install.packages(pkg, repos = "https://cran.rstudio.com", verbose = TRUE)
    }, error = function(e) {
      message(paste("FAIL:", pkg, "-", e$message))
    })
  }
}

# Install the rest of the spatial analysis packages
packages <- c(
  "SpatialExperiment", "SingleCellExperiment", "DropletUtils", "VisiumIO",
  "scrapper", "SpotSweeper", "Banksy", "spacexr", "ggspavis", "patchwork",
  "pheatmap", "pals", "osfr", "Matrix", "ggplot2", "Statial"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("=== Attempting to install:", pkg, "==="))
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE, verbose = FALSE)
    }, error = function(e) {
      message(paste("Error installing package:", pkg, "-", e$message))
    })
  } else {
    message(paste("Package already available:", pkg))
  }
}

suppressPackageStartupMessages({
  library(Banksy)
  library(dplyr)
  library(DropletUtils)
  library(ggplot2)
  library(ggspavis)
  library(magick)
  library(OSTA.data)
  library(patchwork)
  library(pheatmap)
  library(scrapper)
  library(sf)
  library(spacexr)
  library(SpatialExperiment)
  library(SpotSweeper)
  library(tidyr)
  library(VisiumIO)
  library(Matrix)
  library(pals)
  library(osfr)
  library(Statial)
})

# ==============================================================================
# 2. Download Datasets Directly from OSF (Bypassing OSTA.data_load bug)
# ==============================================================================

library(osfr)
library(VisiumIO)
library(SpatialExperiment)

# 1. Set destination folder on your working drive
dest_dir <- "D:/hariprasad/Spatial_omics/VisiumHD_HumanColon_Oliveira"
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Retrieve the OSF project node (5n4q3)
project <- osf_retrieve_node("5n4q3")
files <- osf_ls_files(project)

# 3. Locate and download the Visium HD folder recursively
target_folder <- files[files$name == "VisiumHD_HumanColon_Oliveira", ]

message("Downloading Visium HD dataset directly via osfr...")
osf_download(
  target_folder, 
  path = dest_dir, 
  recurse = TRUE, 
  conflicts = "overwrite"
)

# 4. Resolve the root directory path containing binned_outputs / segmented_outputs
sr_path <- if (dir.exists(file.path(dest_dir, "VisiumHD_HumanColon_Oliveira", "binned_outputs"))) {
  file.path(dest_dir, "VisiumHD_HumanColon_Oliveira")
} else {
  dest_dir
}

message("Data ready at: ", sr_path)

# Example: Import 8 µm bins
vhd8 <- TENxVisiumHD(
  spacerangerOut = sr_path,
  processing     = "filtered",
  format         = "h5",
  images         = "lowres",
  bin_size       = "008"
) |> import()

# use gene symbols as feature names
gs <- rowData(vhd8)$Symbol
rownames(vhd8) <- make.unique(gs)
vhd8

# These data come with bin-level annotations from deconvolution estimates by RCTD implemented in spacexr. 
# Specifically, two sets of labels are available that correspond to the most and second most frequent of 38 cell types per bin, respectively:

# retrieve annotations
# 1. Point to the file using your defined sr_path
gz <- "binned_outputs/square_008um/deconvolution.csv.gz"
df <- read.csv(file.path(sr_path, gz), row.names = 2)
head(df <- df[complete.cases(df), -1])

# keep only annotated bins
vhd8 <- vhd8[, rownames(df)]
names(df) <- c("DeconClass", "DeconLabel1", "DeconLabel2")
colData(vhd8)[names(df)] <- df


lab <- grep("Label", names(cd <- colData(vhd8)), value=TRUE)
pal <- hcl.colors(length(unique(unlist(cd[lab]))), "Spectral")
lapply(lab, \(.) {
  plotCoords(vhd8, annotate=., point_size=0.2, point_shape=15) + ggtitle(.)
}) |>
  wrap_plots(nrow=1, guides="collect") &
  scale_color_manual(NULL, values=pal) &
  theme(legend.key.size=unit(0, "lines"))

# Another resolution of Visium HD with 16 µm bins can be read in as follows:

# Example: Import 16 µm bins
vhd16 <- TENxVisiumHD(
  spacerangerOut = sr_path,
  processing     = "filtered",
  format         = "h5",
  images         = "lowres",
  bin_size       = "016"
) |> import()


# use symbols as feature names
gs <- rowData(vhd16)$Symbol
rownames(vhd16) <- make.unique(gs)
vhd16

# Remove bins overlaying empty tissue
vhd16$libsize <- colSums(counts(vhd16))

# We can visualize such bins by zooming into an empty tissue region. 
# On the left-hand side, one region shows a clear gap with low counts

p1 <- plotVisium(vhd16, 
                 annotate="libsize", point_shape=22,
                 zoom=TRUE, show_axes=TRUE, point_size=1.3) + 
  xlim(c(422, 442)) + ylim(c(318, 330)) + 
  ggtitle("16 µm - original")

vhd16 <- vhd16[, vhd16$libsize > 400]

p2 <- plotVisium(vhd16, 
                 annotate="libsize", point_shape=22,
                 zoom=TRUE, show_axes=TRUE, point_size=1.3) + 
  xlim(c(422, 442)) + ylim(c(318, 330)) + 
  ggtitle("16 µm - post > 400 UMI QC")

p1 | p2

# Clustering results generated by Space Ranger are provided for both resolutions

# Locate the clustering annotations using sr_path
csv <- list.files(sr_path, pattern = "clustering\\.csv\\.gz$", recursive = TRUE, full.names = TRUE)

# Sort files to match resolutions (8 µm then 16 µm)
csv_8um <- grep("008um", csv, value = TRUE)[1]
csv_16um <- grep("016um", csv, value = TRUE)[1]

df_8 <- read.csv(csv_8um, row.names = "barcode")
df_16 <- read.csv(csv_16um, row.names = "barcode")

colData(vhd8)$cluster <- factor(df_8[colnames(vhd8), "cluster"])
colData(vhd16)$cluster <- factor(df_16[colnames(vhd16), "cluster"])


# We can visualize the given clustering results at 8 µm and 16 µm resolutions. 
# Later, we will compare the clustering result at 16 µm with that from Banksy

(plotCoords(vhd8, 
            annotate="cluster", point_size=0.05, point_shape=15,
            pal=unname(pals::kelly())) + ggtitle("8 µm")) |
  (plotCoords(vhd16[, !is.na(vhd16$cluster)], 
              annotate="cluster", point_size=0.1, point_shape=15, 
              pal=unname(pals::kelly())) + ggtitle("16 µm")) 

# To keep runtimes low, our downstream analysis will be performed on a subset area of 16 µm bins 
# (1/64 of the data coverage area), which we read in here; there are no annotations available for this resolution in the public domain


.rng <- \(spe) {
  xy <- spatialCoords(spe)*scaleFactors(spe)
  xs <- range(xy[, 1]); ys <- nrow(imgRaster(spe))-range(xy[, 2])
  x1 <- xs[1] + 4*(xs[2]-xs[1])/8; x2 <- xs[2] - 3*(xs[2]-xs[1])/8
  y1 <- ys[1] + 3*(ys[2]-ys[1])/8; y2 <- ys[2] - 4*(ys[2]-ys[1])/8
  list(box=c(x1, x2, y1, y2), cov=c(xs, ys))
}
vhd8r <- .rng(vhd8)
# vhd16r <- .rng(vhd16) 
# use 8 µm bounding boxes to subset SPEs at both resolutions;
# it is similar to the 16 µm box range
.box <- \(roi, lty, lwd, col) {
  geom_rect(
    xmin=vhd8r[[roi]][1], xmax=vhd8r[[roi]][2], 
    ymin=vhd8r[[roi]][3], ymax=vhd8r[[roi]][4],
    col=col, fill=NA, linetype=lty, linewidth=lwd)
}
cov <- .box(roi="cov", lty=2, lwd=1/2, col="grey")
box <- .box(roi="box", lty=4, lwd=2/3, col="black")
# plotting
.lim <- \(spe) list(
  xlim(spe[["box"]][c(1, 2)]),
  ylim(spe[["box"]][c(4, 3)]))

plotVisium(vhd8, spots=FALSE, point_shape=22) + cov + box +
  ggtitle("grey: data coverage\n black: zoomed region") + 
  plotVisium(vhd8, point_size=0.8, zoom=TRUE, point_shape=22) + .lim(vhd8r) +
  ggtitle("8 µm bins in black box") + 
  plotVisium(vhd16, point_size=1.6, zoom=TRUE, point_shape=22) + .lim(vhd8r) +
  ggtitle("16 µm bins in black box") + 
  plot_layout(nrow=1, widths=c(1.5, 1, 1)) & facet_null() &
  theme(plot.title=element_text(hjust=0.5, vjust=0.5))


# helper for subsetting
.sub <- function(spe, rng, roi="box") {
  xs <- rng[[roi]][c(1, 2)]
  ys <- nrow(imgRaster(spe)) - rng[[roi]][c(3,4)]
  xy <- spatialCoords(spe)*scaleFactors(spe)
  spe[, xy[, 1] > xs[1] & xy[, 1] < xs[2] & 
        xy[, 2] > ys[1] & xy[, 2] < ys[2] ]
}

# Subset the 16 µm object to bins in the black box:

.vhd16 <- .sub(spe=vhd16, rng=vhd8r)
dim(.vhd16)

# ==============================================================================
# 4. Import Visium HD 16 µm and Run SpotSweeper Quality Control
# ==============================================================================
# From now on, we are operating on 2416 16 µm bins.


# calculate per-bin QC metrics
mt <- grepl("^MT-", rownames(.vhd16))
.vhd16 <- quickRnaQc.se(.vhd16, subsets=list(mt=mt))
# determine outliers based on 
# - low log-library size
# - few uniquely detected features
# - high mitochondrial count proportion
.vhd16 <- localOutliers(.vhd16, metric="sum", direction="lower", log=TRUE)
.vhd16 <- localOutliers(.vhd16, metric="detected", direction="lower", log=TRUE)
.vhd16 <- localOutliers(.vhd16, metric="subset.proportion.mt", direction="higher", log=FALSE)
.vhd16$discard <- 
  .vhd16$sum_outliers | 
  .vhd16$detected_outliers | 
  .vhd16$subset.proportion.mt_outliers
# tabulate number of bins retained 
# vs. removed by any criterion 
table(.vhd16$discard)

# Deduplicate colData columns if re-running
colData(.vhd16) <- colData(.vhd16)[, !duplicated(colnames(colData(.vhd16))), drop = FALSE]

# Ensure sum_log exists for visualization
if (!"sum_log" %in% colnames(colData(.vhd16))) {
  .vhd16$sum_log <- log1p(.vhd16$sum)
}

plotCoords(.vhd16, annotate="sum_log", point_size=0.2, point_shape=15) + 
  ggtitle("log-library size") +
  plotCoords(.vhd16, annotate="subset.proportion.mt", point_size=0.2, point_shape=15) +
  ggtitle("mitochondrial proportion") +
  plot_layout(nrow=1) & theme(
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines")) &
  scale_color_gradientn(colors=pals::jet())


# We can visualize the local outliers identified by SpotSweeper in space

plotCoords(.vhd16, point_shape=15, annotate="discard") + ggtitle("discard") +
  plotCoords(.vhd16, point_shape=15, annotate="sum_outliers") + ggtitle("low_lib_size") +
  plotCoords(.vhd16, point_shape=15, annotate="detected_outliers") + ggtitle("low_n_features") +
  plot_layout(nrow=1, guides="collect") & 
  theme(
    plot.title=element_text(hjust=0.5),
    legend.key.size=unit(0, "lines")) &
  guides(col=guide_legend(override.aes=list(size=3))) & 
  scale_color_manual("discard", values=c("lavender", "purple"))


# Lastly, we subset the Visium HD 16 µm object to bins that passed QC:

.vhd16 <- .vhd16[, !.vhd16$discard]
dim(.vhd16)


# ==============================================================================
# 5. Normalization, Banksy Spatially-Aware PCA & Leiden Clustering
# ==============================================================================

# First, we identify highly variable genes (HVGs) in the target area

.vhd16 <- normalizeRnaCounts.se(.vhd16)
assay(.vhd16, "logcounts") <- as(assay(.vhd16, "logcounts"), "dgCMatrix")

.vhd16 <- chooseRnaHvgs.se(
  .vhd16, top=3e3,
  more.var.args=list(use.min.width=TRUE))
hvg <- rownames(.vhd16)[rowData(.vhd16)$hvg]


# set seed for random number generation
# in order to make results reproducible
set.seed(112358)
# 'Banksy' parameter settings
k <- 8   # consider first order neighbors
l <- 0.2 # use little spatial information
a <- "logcounts"
xy <- c("array_row", "array_col")
# restrict to selected features
tmp <- .vhd16[hvg, ]
# compute spatially aware 'Banksy' PCs
tmp <- computeBanksy(tmp, assay_name=a, coord_names=xy, k_geom=k)
tmp <- runBanksyPCA(tmp, lambda=l, npcs=20)
reducedDim(.vhd16, "PCA") <- reducedDim(tmp)
## run UMAP (for visualization purposes only)
# .vhd16 <- runUmap.se(.vhd16, reddim.type="PCA")
# build cellular shared nearest-neighbor (SNN) graph
# and cluster using Leiden community detection algorithm
.vhd16 <- clusterGraph.se(
  .vhd16, num.neighbors=20,
  method="leiden", resolution=1.2,
  more.build.args=list(weight.scheme="jaccard"),
  more.cluster.args=list(leiden.objective="modularity", seed=112358),
  reddim.type="PCA", output.name="Banksy")
table(.vhd16$Banksy)


# We can score marker genes for each cluster

# marker gene scoring
mgs <- scoreMarkers.se(.vhd16, groups=.vhd16$Banksy)
# select for a few markers per cluster
top <- lapply(mgs, \(df) rownames(df)[df$cohens.d.min.rank <= 2])
top <- unique(unlist(top))
# average expression by clusters
pbs <- aggregateAcrossCells.se(
  .vhd16[top, ], factors=.vhd16$Banksy, assay.type="logcounts")
means <- t(t(assay(pbs, "sums")) / pbs$counts)


# The marker genes in each cluster can be visualized with a heatmap:

# visualize averages z-scaled across clusters
pheatmap(
  mat=t(means), scale="column", breaks=seq(-2, 2, length=101),
  cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)


# we can visualize the bin-wise expression of selected markers in space

gs <- c("MMP2", "PIGR", "IGHG1")
ps <- lapply(gs, \(.) plotCoords(.vhd16, annotate=., point_shape=15, 
                                 point_size=0.8, assay_name="logcounts"))
wrap_plots(ps, nrow=1) & theme(
  legend.key.width=unit(0.5, "lines"),
  legend.key.height=unit(1, "lines")) &
  scale_color_gradientn(colors=rev(hcl.colors(9, "Rocket")))

# ==============================================================================
# 6. Single-Cell Reference Preparation
# ==============================================================================

# We will now perform deconvolution on the 16 µm binned Visium HD data. 
# First, we retrieve matching (Chromium) scRNA-seq data, which includes 
# low- (Level1) and high-resolution (Level2) annotations of cells into 9 
# and 31 clusters, respectively. In order to ensure similar transcriptional 
# profiles across technologies, we filter the reference population to patient 2 only

library(DropletUtils)
library(SingleCellExperiment)
library(osfr)

# 1. Target directory specifically for Chromium reference
target_dir <- "D:/hariprasad/Spatial_omics/Chromium_Reference"
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Look ONLY for files inside paths containing 'Chromium'
h5_local <- list.files(
  c(target_dir, tempdir(), "C:/Users/acer/AppData/Local/Temp"),
  pattern = "filtered_feature_bc_matrix\\.h5$",
  recursive = TRUE, 
  full.names = TRUE
)
h5_chrom <- grep("Chromium", h5_local, value = TRUE, ignore.case = TRUE)

if (length(h5_chrom) > 0) {
  h5 <- h5_chrom[1]
  ref_dir <- dirname(h5)
  message("Using existing local Chromium reference from: ", ref_dir)
} else {
  message("Chromium reference not found locally. Downloading via osfr...")
  node <- osf_retrieve_node("5n4q3")
  files <- osf_ls_files(node)
  target_folder <- files[files$name == "Chromium_HumanColon_Oliveira", ]
  osf_download(target_folder, path = target_dir, recurse = TRUE, conflicts = "overwrite")
  
  ref_dir <- if (dir.exists(file.path(target_dir, "Chromium_HumanColon_Oliveira"))) {
    file.path(target_dir, "Chromium_HumanColon_Oliveira")
  } else {
    target_dir
  }
}

# 3. Locate the Chromium matrix and metadata inside ref_dir
fs <- list.files(ref_dir, recursive = TRUE, full.names = TRUE)
h5 <- grep("filtered_feature_bc_matrix\\.h5$", fs, value = TRUE)[1]
csv <- grep("metadata.*\\.csv$", fs, value = TRUE, ignore.case = TRUE)[1]

# Safety check
if (is.na(csv) || !file.exists(csv)) {
  stop("Metadata CSV not found in: ", ref_dir, "\nFiles present: ", paste(basename(fs), collapse = ", "))
}

message("Reading count matrix from: ", h5)
message("Reading cell metadata from: ", csv)

# 4. Read into SingleCellExperiment
sce <- read10xCounts(h5, col.names = TRUE)

# 5. Attach cell metadata
cd <- read.csv(csv, row.names = 1, check.names = FALSE)
colData(sce)[names(cd)] <- cd[colnames(sce), , drop = FALSE]

# 6. Gene symbols and patient filtering
rownames(sce) <- make.unique(rowData(sce)$Symbol)
sce <- sce[, sce$QCFilter == "Keep"]
sce <- sce[, grepl("P2", sce$Patient)]
sce$Level1 <- gsub("\\s", ".", sce$Level1)

# Inspect cleaned reference
sce
table(sce$Level1)

# We can merge the provided deconvolution annotation "DeconLabel1" 
# (most likely cell type for singlets, most dominant type for doublets) into a lower resolution using the single-cell reference annotations:
i <- match(vhd8$DeconLabel1, sce$Level2)
j <- match(vhd8$DeconLabel2, sce$Level2)
vhd8$.DeconLabel1 <- sce$Level1[i]
vhd8$.DeconLabel2 <- sce$Level1[j]
vhd8 <- vhd8[, !is.na(vhd8$.DeconLabel1)]
table(vhd8$.DeconLabel1)

# According to RCTD, approximately 75% of the 8 µm bins are singlets:
round(100*mean(vhd8$DeconClass == "singlet"), 2)

# Let us check the top-5 common doublet pairs according to the "doublet_certain" class:

dbl <- vhd8$DeconClass == "doublet_certain"
lab <- c(".DeconLabel1", ".DeconLabel2")
df <- data.frame(colData(vhd8)[dbl, lab])
# sort as to ignore order
df <- apply(df, 1, sort)
df <- do.call(rbind, df)
names(df) <- lab
# count unique pairs
ij <- paste(df[, 1], df[, 2], sep=";")
head(sort(table(ij), decreasing=TRUE), 5)


# We can perform the same subsetting procedure to 8 µm bins and compare the 
# difference in the number of bins in the zoomed (black box) area.

# subset zoomed (black box) area
.vhd8 <- .sub(spe=vhd8, rng=vhd8r)
# compare number of bins between resolutions
(n <- ncol(.vhd8))

(m <- ncol(.vhd16))

round(n/m, 2)

# ==============================================================================
# 7. RCTD Deconvolution on 16 µm Bins
# ==============================================================================
# We also assume there are at most two cell types in a 16 µm bin as in 8 µm, and thus we would expect a smaller proportion of singlets

library(spacexr)
library(SpatialExperiment)
library(SingleCellExperiment)

# 1. Create .sce by downsampling the loaded reference 'sce'
# (Crucial: prevents 'object .sce not found' error)
set.seed(112358)
cs <- split(seq_len(ncol(sce)), sce$Level1)
cs <- lapply(cs, \(.) sample(., min(length(.), 4e3)))
.sce <- sce[, unlist(cs)]

message("Downsampled reference created: ", ncol(.sce), " cells.")

# 2. Define createRctd helper function
createRctd <- function(spe, sce, cell_type_col = "Level1", max_cores = 1) {
  ref_counts <- assay(sce, "counts")
  ct <- factor(colData(sce)[[cell_type_col]])
  names(ct) <- colnames(sce)
  nUMI_ref <- colSums(ref_counts)
  
  ref <- Reference(counts = ref_counts, cell_types = ct, nUMI = nUMI_ref)
  
  sp_counts <- assay(spe, "counts")
  coords <- as.data.frame(spatialCoords(spe))
  colnames(coords) <- c("x", "y")
  rownames(coords) <- colnames(spe)
  
  valid <- colSums(sp_counts) > 0
  sp <- SpatialRNA(coords = coords[valid, ], counts = sp_counts[, valid])
  
  create.RCTD(spatialRNA = sp, reference = ref, max_cores = max_cores)
}

# 3. Define runRctd helper function
runRctd <- function(rctd, rctd_mode = "doublet", max_cores = 1) {
  rctd@config$max_cores <- max_cores
  spacexr::run.RCTD(rctd, doublet_mode = rctd_mode)
}

# 4. Build RCTD data object (using max_cores = 1)
rctd_data <- createRctd(.vhd16, .sce, cell_type_col = "Level1", max_cores = 1)

# Optional safety save to disk
saveRDS(rctd_data, "rctd_data_16um.rds")

# 5. Run RCTD deconvolution with max_cores = 1 (prevents Windows socket crash)
res <- runRctd(rctd_data, max_cores = 1, rctd_mode = "full")

# 6. Verify completion
table(res@results$results_df$spot_class, useNA = "ifany")

# Weights inferred by RCTD correspond to proportions of cell types, such that they sum to 1 for a given observation

# Extract raw weights matrix (spots x cell types)
ws_raw <- as.matrix(res@results$weights)

# Calculate sum per spot and normalize valid spots to 1
spot_sums <- rowSums(ws_raw)
valid_spots <- spot_sums > 0
ws_norm <- ws_raw
ws_norm[valid_spots, ] <- sweep(ws_raw[valid_spots, , drop = FALSE], 1, spot_sums[valid_spots], "/")

# 1. Ensure ws_df has spots as rows and cell types as columns
ws_df <- as.data.frame(ws_norm)

# 2. Subset .vhd16 to the deconvolved spots
overlap_spots <- intersect(colnames(.vhd16), rownames(ws_df))
.vhd16 <- .vhd16[, overlap_spots]
ws_df <- ws_df[overlap_spots, ]

# 3. Add each cell type proportion as a metadata column in colData
for (ct_name in colnames(ws_df)) {
  colData(.vhd16)[[ct_name]] <- ws_df[[ct_name]]
}

# Verify columns are attached
head(colData(.vhd16)[, colnames(ws_df)])

library(ggspavis)
library(patchwork)
library(ggplot2)

# Plot proportions for each cell type
cell_types <- colnames(ws_df)

ps <- lapply(cell_types, \(ct) {
  plotCoords(.vhd16, annotate = ct, point_size = 0.8, point_shape = 15) +
    ggtitle(ct)
})

wrap_plots(ps, nrow = 2, guides = "collect") &
  theme(
    legend.key.width  = unit(0.5, "lines"),
    legend.key.height = unit(1, "lines"),
    plot.title        = element_text(hjust = 0.5, size = 10)
  ) &
  scale_color_gradientn(colors = rev(hcl.colors(9, "Rocket")))


# Find the dominant cell type per bin
majority_idx <- apply(ws_df, 1, which.max)
dominant_type <- colnames(ws_df)[majority_idx]

# Flag rejected bins where total weight was 0
is_reject <- rowSums(ws_df) == 0
dominant_type[is_reject] <- NA

# Harmonize names with .vhd8 (replace '.' with space if .vhd8 uses spaces, or keep identical)
if (any(grepl(" ", levels(factor(.vhd8$.DeconLabel1))))) {
  dominant_type <- gsub("\\.", " ", dominant_type)
}

# Assign to colData as factor
.vhd16$.DeconLabel1 <- factor(dominant_type)

# Summary table of majority-vote calls
table(.vhd16$.DeconLabel1, useNA = "ifany")


# 4. Distribution of cell type frequency ratios (8 µm vs 16 µm)
# Align factors between resolutions to calculate exact ratios
common_labels <- intersect(names(table(.vhd8$.DeconLabel1)), names(table(.vhd16$.DeconLabel1)))

n8 <- table(.vhd8$.DeconLabel1)[common_labels]
n16 <- table(.vhd16$.DeconLabel1)[common_labels]

ratios <- as.numeric(n8 / n16)
df_ratios <- data.frame(
  cell_type = common_labels,
  ratio = ratios
)

# Plot density of the ratio distribution
ggplot(df_ratios, aes(x = ratio)) + 
  geom_density(fill = "grey90", color = "black") + 
  geom_rug() +
  xlab("Cell type frequency ratio (8 µm / 16 µm)") + 
  ggtitle("Distribution of cell type frequency\nratio between 8 vs. 16 µm bins") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5))


# ==============================================================================
# 8. Comparison of Deconvolution at 8 µm vs 16 µm
# ==============================================================================

plotVisium(.vhd8, 
           annotate=".DeconLabel1", zoom=TRUE, 
           point_size=0.8, point_shape=22) + 
  ggtitle("8 µm") + 
  plot_spacer() +
  plotVisium(.vhd16, 
             annotate=".DeconLabel1", zoom=TRUE, 
             point_size=1.6, point_shape=22) + 
  ggtitle("16 µm") +
  plot_layout(nrow=1, guides="collect", widths=c(1, 0.05, 1)) &
  guides(col=guide_legend(override.aes=list(size=2))) &
  scale_fill_manual(values=unname(pals::trubetskoy())) & 
  facet_null() & theme(
    plot.title=element_text(hjust=0.5),
    legend.key.size=unit(0, "lines")) 


# ==============================================================================
# 6 Concordance of results
# ==============================================================================

plotVisium(.vhd16, 
           annotate="cluster", zoom=TRUE, 
           point_shape=22, point_size=1.6, 
           pal=unname(pals::kelly())) + 
  plot_spacer() +
  plotVisium(.vhd16, 
             annotate="Banksy", zoom=TRUE, 
             point_shape=22, point_size=1.6, 
             pal=unname(pals::kelly())) + 
  plot_spacer() +
  plotVisium(.vhd16, 
             annotate=".DeconLabel1", zoom=TRUE, 
             point_shape=22, point_size=1.6, 
             pal=unname(pals::trubetskoy())) +
  plot_layout(nrow=1, widths=c(1, 0.05, 1, 0.05, 1)) & 
  facet_null() & theme(legend.key.size=unit(0, "lines"))


# We observe that Banksy delineates tissue borders very well, while deconvolution 
# provides direct biological insight into the underlying cluster

fq <- prop.table(table(.vhd16$Banksy, .vhd16$.DeconLabel1), 1)
pheatmap(fq, cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)


# From the above heatmap, clusters 5 and 7 are mostly immune cells; 
# cluster 9 maps to endothelial; clusters 1 and 13 correspond to intestinal epithelial 
# and smooth muscle cells, respectively; fibroblasts map to cluster 10; tumor spans clusters 3, 11 and 12


# ==============================================================================
# 7 Neighborhood analysis
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Statial non-interactively
BiocManager::install("Statial", update = FALSE, ask = FALSE)

# If not available on your current Bioconductor release repo, install via remotes:
if (!requireNamespace("Statial", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("SydneyBioX/Statial")
}


# we want to quantify the abundance between T cells, B cells, or Myeloid bins to Fibroblast or Tumor bins in 8 µm zoomed region

library(Statial)
library(SpatialExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(tidyr)
library(Matrix)

# Remove NA labels
.vhd8 <- .vhd8[, !is.na(.vhd8$.DeconLabel1)]

# Ensure space-separated names to match: "T cells", "B cells", "Smooth Muscle"
.vhd8$.DeconLabel1 <- gsub("\\.", " ", as.character(.vhd8$.DeconLabel1))
.vhd8$.DeconLabel1 <- factor(.vhd8$.DeconLabel1)

# Check label representation
table(.vhd8$.DeconLabel1)

# Let’s have a look at the abundance and distance results, also including deconvolution results:
# Extract coordinates and sync with colData
xy <- data.frame(spatialCoords(.vhd8))
colData(.vhd8)[names(xy)] <- xy

# Ensure sample_id is present
if (!"sample_id" %in% colnames(colData(.vhd8))) {
  .vhd8$sample_id <- "sample01"
}

# Run getAbundances on 1 core to prevent Windows socket crashes
.vhd8 <- getAbundances(
  .vhd8, 
  spatialCoords = names(xy),
  cellType      = ".DeconLabel1", 
  imageID       = "sample_id", 
  r             = 200, 
  nCores        = 1
)

# Extract abundances data frame
df <- as.data.frame(reducedDim(.vhd8, "abundances"))
df$CT <- .vhd8$.DeconLabel1
head(df[, 1:min(5, ncol(df))])

# In this region, we observe more fibroblast than tumor bins around each of the immune bins. 
# However, note that in this zoomed region, we have nearly 1.4 times more fibroblast bins compared to tumor bins. 
# Thus, we need to scale abundance values by the number of target cell type bins observed in the region.

# Subset to immune bins
source_cells <- df$CT %in% c("T cells", "B cells", "Myeloid")
target_cols  <- intersect(c("Tumor", "Fibroblast"), colnames(df))

fd <- pivot_longer(
  df[source_cells, c(target_cols, "CT")], 
  cols = all_of(target_cols),
  names_to = "target", 
  values_to = "n"
)

# Scale by the background counts of Tumor and Fibroblast in the ROI
n_tum <- sum(.vhd8$.DeconLabel1 == "Tumor", na.rm = TRUE)
n_fib <- sum(.vhd8$.DeconLabel1 == "Fibroblast", na.rm = TRUE)

fd$p <- ifelse(fd$target == "Tumor", fd$n / n_tum, fd$n / n_fib)

# Plot raw counts vs. normalized relative abundance
p_raw <- ggplot(fd, aes(x = CT, y = n, fill = target)) + 
  geom_boxplot(key_glyph = "point", outlier.size = 0.5) +
  labs(y = "# cells within 200px") +
  theme_classic() +
  theme(axis.title.x = element_blank())

p_norm <- ggplot(fd, aes(x = CT, y = p, fill = target)) + 
  geom_boxplot(key_glyph = "point", outlier.size = 0.5) +
  labs(y = "relative abundance") +
  theme_classic() +
  theme(axis.title.x = element_blank())

(p_raw + p_norm) +
  plot_layout(nrow = 1, guides = "collect") & 
  scale_fill_manual(values = c("Fibroblast" = "pink", "Tumor" = "yellow")) & 
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) &
  theme(legend.key.size = unit(0.5, "lines"))


# After correction, we see amplified differences in the relative abundance of fibroblast bins 
# compared to tumor bins around each immune cell type in the zoomed region.

# We can also investigate the “cell-cell-marker” relationship based on the minimum distance between 
# two cell types within a given radius. The matrix of minimum distance between all pairwise cell types is 
# stored in the reduced dimension, similar to how the abundance matrix was previously stored

# Convert to SingleCellExperiment to bypass Statial issue #16
tmp <- as(.vhd8, "SingleCellExperiment")

tmp <- getDistances(
  tmp, 
  cellType      = ".DeconLabel1", 
  imageID       = "sample_id",
  spatialCoords = names(xy), 
  maxDist       = 200, 
  nCores        = 1
)

reducedDim(.vhd8, "distances") <- reducedDim(tmp, "distances")

# Ensure counts and logcounts assays are realized in-memory as dgCMatrix
# This avoids slow HDF5 disk lookups ("TENxMatrixSeed" warnings) and completes instantly
assay(.vhd8, "counts") <- as(assay(.vhd8, "counts"), "dgCMatrix")

if (!"logcounts" %in% assayNames(.vhd8)) {
  .vhd8 <- scrapper::normalizeRnaCounts.se(.vhd8)
}
assay(.vhd8, "logcounts") <- as(assay(.vhd8, "logcounts"), "dgCMatrix")

# Here, we aim to determine whether the expression level of the tumor marker 
# changes in proximity to fibroblast. Using FTH1 as an example, we observe a 
# decrease in FTH1 expression in tumor bins as the distance to fibroblast increase. 
# This trend could be visualized using a scatter plot.

# Plot FTH1 expression as a function of distance to Fibroblasts
sample_name <- unique(.vhd8$sample_id)[1]

p_state <- plotStateChanges(
  cells         = .vhd8, 
  image         = sample_name, 
  from          = "Tumor", 
  to            = "Fibroblast",
  marker        = "FTH1", 
  cellType      = ".DeconLabel1", 
  imageID       = "sample_id",
  spatialCoords = names(xy)
)

(p_state$image + facet_null()) | p_state$scatter

# ==============================================================================
# 8 Conclusion
# ==============================================================================
# In this workflow, we demonstrated how to perform a standard analysis pipeline on a subset of 
# Visium HD data binned at 16 µm. Most methods and visualization techniques developed for 
# Visium can be applied to Visium HD. However, more sparsity is expected, as each data point 
# (e.g. 2, 8, or 16 µm) has smaller area compared to Visium (55 µm). 
# Each bin cannot be interpreted as a cell, but rather, a segmentation-free approach is applied here