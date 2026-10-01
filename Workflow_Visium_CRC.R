#Analyzing Visium data on a human colorectal cancer biopsy from de Oliveira et al. (2025)
#===============================================================================
# 1.Dependencies
#===============================================================================

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AUCell")

# Load the library again
library(AUCell)
library(BiocParallel)
library(DropletUtils)
library(ggspavis)
library(igraph)
library(jsonlite)
library(msigdbr)
library(OSTA.data)
library(patchwork)
library(pheatmap)
library(scater)
library(scrapper)
library(spacexr)
library(SpatialExperiment)
library(VisiumIO)
# specify whether/how to 
# perform parallelization
bp <- MulticoreParam(th <- 4)
# set seed for random number generation
# in order to make results reproducible
set.seed(194849)


#===============================================================================
# 2 Data import
#===============================================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STexampleData")

library(osfr)
library(VisiumIO)

# 1. Establish the clean local workspace folder structure
local_dir <- "Visium_HumanColon_Oliveira_Data"
spatial_dir <- file.path(local_dir, "spatial")
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Query the exact OSF project node for the OSTA datasets
project <- osf_retrieve_node("5n4q3")

# 3. Locate the 'Visium_HumanColon_Oliveira' dataset subfolder on OSF
osf_files <- osf_ls_files(project)
target_folder <- osf_files[osf_files$name == "Visium_HumanColon_Oliveira", ]

# 4. Access the contents of the 'outs' subfolder inside it
outs_folder <- osf_ls_files(target_folder)
outs_contents <- osf_ls_files(outs_folder[outs_folder$name == "outs", ])

# 5. Loop and download every file directly to the correct destination
message("Downloading files directly from OSF...")
for (i in seq_len(nrow(outs_contents))) {
  file_item <- outs_contents[i, ]
  fnm <- file_item$name
  
  # Send spatial image metadata files to the 'spatial' subdirectory
  spatial_extensions <- c(".png", ".json", ".csv", ".jpg")
  is_spatial <- any(sapply(spatial_extensions, function(ext) grepl(ext, fnm))) && !grepl("matrix", fnm)
  
  dest_path <- if (is_spatial) spatial_dir else local_dir
  osf_download(file_item, path = dest_path, conflicts = "overwrite")
}

# 6. Safety check: Rename the matrix file to standard naming if needed
all_base_files <- list.files(local_dir)
h5_match <- all_base_files[grep("\\.h5$", all_base_files)]
if (length(h5_match) > 0 && h5_match != "filtered_feature_bc_matrix.h5") {
  file.rename(from = file.path(local_dir, h5_match), 
              to = file.path(local_dir, "filtered_feature_bc_matrix.h5"))
}

message("Download complete! Constructing object...")

# 7. Construct and Import into the 'spe' Object
obj <- TENxVisium(
  spacerangerOut = local_dir, 
  format = "h5", 
  processing = "filtered", 
  images = "lowres"
)

spe <- import(obj)

# Output your built spatial matrix
print(spe)

#===============================================================================
# 3 Quality control
#===============================================================================
# 1. Install scrapper if you haven't already
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scrapper")


# 1. Install patchwork if it is completely missing from your system
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}

# 2. Make sure ggplot2 is loaded for theme/scale elements
library(ggplot2)
library(patchwork)

# 2. Load the library
library(scrapper)

# 2. Load the library
library(patchwork)

# 3. Use gene symbols as feature names
rownames(spe) <- make.unique(rowData(spe)$Symbol)

# 4. Add quality control metrics & determine outliers
# Note: scrapper automatically assigns the logical filter vector to spe$keep
sub <- list(mt = grep("^MT-", rownames(spe)))
spe <- quickRnaQc.se(spe, subsets = sub)

spe$discard <- !spe$keep
spe$log_sum <- log1p(spe$sum)
spe$mt_prop <- spe$subset.proportion.mt


# tabulate # & % of cells that'd be 
# discarded for different reasons
ths <- metadata(spe)$qc$thresholds
ols <- data.frame(
  low_sum=spe$sum < ths$sum,
  low_detected=spe$detected < ths$detected,
  high_mt_prop=spe$mt_prop > ths$subset.proportion,
  discard=spe$discard)
data.frame(
  check.names=FALSE,
  `#`=apply(ols, 2, sum), 
  `%`=round(100*apply(ols, 2, mean), 2))


# 1. Load the essential spatial visualization libraries
library(ggspavis)
library(ggplot2)
library(patchwork)

# 1. Clean up colData by keeping only the uniquely named columns
colData(spe) <- colData(spe)[, !duplicated(colnames(colData(spe))), drop = FALSE]

# 2. Assign the outlier metrics cleanly
colData(spe)[names(ols)] <- ols

# 3. Re-run your patchwork visualization panel
library(ggspavis)
library(ggplot2)
library(patchwork)

lapply(names(ols), \(.) {
  ggspavis::plotCoords(spe, annotate = .) + ggplot2::ggtitle(.)
}) |>
  patchwork::wrap_plots(nrow = 1, guides = "collect") &
  ggplot2::guides(col = guide_legend(override.aes = list(size = 3))) &
  ggplot2::scale_color_manual("discard", values = c("lavender", "purple")) &
  ggplot2::theme(plot.title = element_text(hjust = 0.5), legend.key.size = unit(0, "lines"))


#===============================================================================
# 4 Processing
#===============================================================================

# log-library size normalization
spe <- normalizeRnaCounts.se(spe)
# highly variable feature selection
spe <- chooseRnaHvgs.se(spe, top=2e3,
                        more.var.args=list(use.min.width=TRUE))
# principal component analysis
spe <- runPca.se(spe, features=rowData(spe)$hvg)

#===============================================================================
# 5 Clustering
#===============================================================================

# PCA-based shared nearest-neighbor (SNN) graph;
# cluster via Leiden community detection algorithm
spe <- clusterGraph.se(spe,
                       output.name="Leiden",
                       method="leiden", resolution=0.5,
                       more.build.args=list(weight.scheme="jaccard"))
table(spe$Leiden)

#===============================================================================
# 6 Deconvolution
#===============================================================================

library(osfr)
library(DropletUtils)
library(spacexr)
library(SpatialExperiment)

# 1. Download dataset via OSF Node ID (Workaround for OSTA.data_load cache bug)
dir.create(td <- tempfile())
project <- osf_retrieve_node("5n4q3")
osf_files <- osf_ls_files(project)
target_folder <- osf_files[osf_files[["name"]] == "Chromium_HumanColon_Oliveira", ]
osf_download(target_folder, path = td, recurse = TRUE, conflicts = "overwrite")
#==============================================================================
  # 2. Read into 'SingleCellExperiment'
  fs <- list.files(td, recursive = TRUE, full.names = TRUE)
h5 <- grep("h5$", fs, value = TRUE)

# Defensive check: verify the file structure safely
library(rhdf5)
h5_contents <- h5ls(h5)

# If it's a molecule info file, construct the matrix manually to bypass the error
if (any(grepl("molecule_info|barcode_idx", h5_contents$group))) {
  message("Detected Molecule Info H5 structure. Building SingleCellExperiment manually...")
  
  # Extract molecule info components safely
  mol_info <- DropletUtils::read10xMolInfo(h5)
  
  # Generate the count matrix from the raw gene and cell data indexes
  counts_mat <- DropletUtils::makeCountMatrix(
    gene = mol_info$data$gene, 
    cell = mol_info$data$cell, 
    value = mol_info$data$reads
  )
  
  # Assign clean dimnames based on parsed features
  dimnames(counts_mat) <- list(
    mol_info$genes[seq_len(nrow(counts_mat))], 
    mol_info$cell[seq_len(ncol(counts_mat))]
  )
  
  # Initialize the SingleCellExperiment container object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_mat))
  
} else {
  # Standard execution if the H5 formatting matches a standard matrix profile
  sce <- read10xCounts(h5, col.names = TRUE)
}

# Inspect output dims to ensure data loaded completely
dim(sce)

#==============================================================================
# 3. Add cell metadata
# FIX 1: Set recursive = TRUE to find files nested inside folders
# FIX 2: Change "csv\\$" to "csv$" to properly capture the extension anchor
fs_all <- list.files(td, recursive = TRUE, full.names = TRUE)
csv_file <- grep("csv$", fs_all, value = TRUE)

# Defensive Check: Verify the file path was extracted successfully before reading
if (length(csv_file) == 0) {
  stop("CSV metadata file not found. Check list.files(td, recursive=TRUE)")
}

# Read the file path smoothly
cd <- read.csv(csv_file, row.names = 1)
colData(sce)[names(cd)] <- cd[colnames(sce), ]


# 4. Use gene symbols as feature names
rownames(sce) <- make.unique(rowData(sce)[["Symbol"]])

# 5. Exclude cells deemed to be of low-quality
sce <- sce[, colData(sce)[["QCFilter"]] == "Keep"]

# 6. Tabulate subpopulations
table(colData(sce)[["Level1"]])

# ==================== DECONVOLUTION STEP ====================

# ==================== DECONVOLUTION STEP ====================

# 7. Prep reference data (Chromium); subset cells from same patient
.sce <- sce[, grepl("P2", colData(sce)[["Patient"]])]

# 8. Downsample to at most 2,000 cells per cluster
cs <- split(seq_len(ncol(.sce)), colData(.sce)[["Level1"]])
cs <- lapply(cs, function(x) sample(x, min(length(x), 2e3)))
.sce <- .sce[, unlist(cs)]

# 9. Run 'RCTD' deconvolution (GitHub Version Syntax)
th <- 4  # Adjust based on your available CPU cores

library(spacexr)

# -------------------------------------------------------------
# 9a-1. Construct the SpatialRNA object (Visium query target)
# -------------------------------------------------------------
spatial_counts <- assay(spe, "counts")
spatial_coords <- as.data.frame(spatialCoords(spe))
colnames(spatial_coords) <- c("x", "y")
rownames(spatial_coords) <- colnames(spe)

spatial_obj <- SpatialRNA(spatial_coords, spatial_counts)

# -------------------------------------------------------------
# 9a-2. Construct the Reference object (Single-cell Chromium)  <-- THIS WAS MISSING
# -------------------------------------------------------------
# Extract raw counts and make sure dimnames are intact
ref_counts <- assay(.sce, "counts")

# Extract cell types as a named factor
cell_types <- factor(colData(.sce)[["Level1"]])
names(cell_types) <- colnames(.sce)

# Calculate library sizes per cell (spacexr requires nUMI)
nUMI <- colSums(ref_counts)

# Build the spacexr Reference object
reference_obj <- Reference(
  counts = ref_counts, 
  cell_types = cell_types, 
  nUMI = nUMI
)

# -------------------------------------------------------------
# 9b. Run deconvolution
# -------------------------------------------------------------
rctd_data <- create.RCTD(spatial_obj, reference_obj, max_cores = th)
res_obj <- run.RCTD(rctd_data, doublet_mode = "full")


# 10. Scale RCTD weights so that proportions for each spot sum to 1

ws <- as.matrix(res_obj@results$weights)

# Normalize each spot across the 9 cell-type estimates
ws <- sweep(ws, 1, rowSums(ws), "/")

# Convert to data.frame
ws <- as.data.frame(ws)

# 11. Add proportion estimates as metadata

for (celltype in colnames(ws)) {
  colData(spe)[[celltype]] <- ws[[celltype]]
}

# 12. Assign each spot the cell type with the highest RCTD estimate

ids <- colnames(ws)[apply(ws, 1, which.max)]

colData(spe)[["RCTD"]] <- factor(ids)

# Cross-tabulate RCTD assignments with Leiden clusters
table(
  colData(spe)[["RCTD"]],
  colData(spe)[["Leiden"]]
)

# 13. Compartmentalize tissue domains

lab <- list(
  tum = "Tumor",
  epi = "Intestinal Epithelial",
  imm = c("B cells", "T cells", "Myeloid"),
  str = c("Endothelial", "Fibroblast", "Smooth Muscle")
)

idx <- match(
  colData(spe)[["RCTD"]],
  unlist(lab)
)

domain_names <- rep.int(
  names(lab),
  sapply(lab, length)
)

colData(spe)[["Domain"]] <- factor(domain_names[idx])

table(colData(spe)[["Domain"]])

#===============================================================================
# 7 Exploratory
#===============================================================================

#coloring by the proportion of a given cell type estimated to fall within a given spot:

lapply(names(ws), \(.) 
       plotCoords(spe, annotate=.)) |>
  wrap_plots(nrow=3) & theme(
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines")) &
  scale_color_gradientn(colors=pals::jet())

lapply(c("Leiden", "Domain", "RCTD"), 
       \(.) plotCoords(spe, annotate=.)) |>
  wrap_plots(nrow=1) &
  theme(legend.key.size=unit(0, "lines")) &
  scale_color_manual(values=unname(pals::trubetskoy()))

# characterize subpopulations from unsupervised clustering

cd <- data.frame(colData(spe))
df <- as.data.frame(with(cd, table(RCTD, Leiden)))
fd <- as.data.frame(with(cd, table(Domain, Leiden)))
ggplot(df, aes(Freq, RCTD, fill=Leiden)) + ggtitle("RCTD") +
  ggplot(fd, aes(Freq, Domain, fill=Leiden)) + ggtitle("Domain") +
  plot_layout(nrow=1, guides="collect") &
  labs(x="Proportion", y=NULL) &
  coord_cartesian(expand=FALSE) &
  geom_col(width=1, col="white", position="fill") &
  scale_fill_manual(values=unname(pals::trubetskoy())) &
  theme_minimal() & theme(aspect.ratio=1,
                          legend.key.size=unit(2/3, "lines"),
                          plot.title=element_text(hjust=0.5))

# inspect the key drivers of (expression) variability in terms of PCs
# add PCs as cell metadata
pcs <- reducedDim(spe, "PCA")
colnames(pcs) <- paste0("PC", seq(ncol(pcs)))
colData(spe)[colnames(pcs)] <- pcs

# visualize PCs 1-6 spatially
lapply(head(colnames(pcs), 6), 
       \(.) plotCoords(spe, annotate=.) +
         scale_color_gradientn(., colors=pals::jet())) |>
  wrap_plots(nrow=2) & theme(
    plot.title=element_blank(),
    legend.key.width=unit(0.5, "lines"),
    legend.key.height=unit(1, "lines"))

# clustering of low-quality spots seen earlier

lapply(c("detected", "log_sum", "mt_prop"), \(.)
       plotColData(spe, x=., y="Leiden", color_by="discard", point_size=0.1) +
         scale_x_discrete(limits=names(sort(by(spe[[.]], spe$Leiden, median))))) |>
  wrap_plots(nrow=1, guides="collect") &
  scale_color_manual("discard", values=c("lavender", "purple")) &
  guides(col=guide_legend(override.aes=list(alpha=1, size=3))) &
  theme_minimal() & theme(
    panel.grid.minor=element_blank(), 
    legend.key.size=unit(0, "lines"))

#===============================================================================
# 8 Signatures
#===============================================================================

# evaluate the expression of sets of genes


# retrieve hallmark gene sets from 'MSigDB'
db <- msigdbr(species="Homo sapiens", collection="H")
# get list of gene symbols, one element per set
gs <- split(db$ensembl_gene, db$gs_name)
# simplify set identifiers (drop prefix, use lower case)
names(gs) <- tolower(gsub("HALLMARK_", "", names(gs)))
# how many sets?
length(gs) 

# how many genes in each?
range(sapply(gs, length))

# Realize sparse gene expression matrix
mtx <- as(logcounts(spe), "dgCMatrix")

# Use Ensembl identifiers as feature names
rownames(mtx) <- rowData(spe)$ID

# Build per-spot gene rankings
rnk <- AUCell_buildRankings(
  mtx,
  BPPARAM = bp,
  plotStats = FALSE,
  verbose = FALSE
)

# Calculate AUC scores
auc <- AUCell_calcAUC(
  geneSets = gs,
  rankings = rnk,
  nCores = 1,
  verbose = FALSE
)

# Extract AUC results
res <- t(assay(auc))

# Make sure rows correspond to spots
res <- as.data.frame(res)

# Check dimensions and spot names
dim(res)
head(rownames(res))
head(colnames(spe))

# Add each gene-set AUC score to colData(spe)
for (set_name in colnames(res)) {
  colData(spe)[[set_name]] <- res[[set_name]]
}


# Calculate variance across spots
var <- colVars(as.matrix(res))

# Select the 8 gene sets with highest variance
top <- names(tail(sort(var), 8))

# Plot the top variable gene sets
lapply(top, \(x) {
  spe[[x]] <- as.numeric(scale(spe[[x]]))
  plotCoords(spe, annotate = x)
}) |>
  wrap_plots(nrow = 2, guides = "collect") &
  scale_color_gradientn(
    colors = pals::jet(),
    oob = scales::squish,
    limits = c(-2.5, 2.5)
  ) &
  theme(
    legend.key.width = unit(0.5, "lines"),
    legend.key.height = unit(1, "lines")
  )


for (. in c("Leiden", "RCTD")) {
  # aggregate AUC values by cluster
  ks <- spe[[.]]
  pb <- aggregateAcrossCells.se(auc[top, ], ks, assay.type="AUC")
  mu <- sweep(assay(pb, "sums"), 2, pb$counts, `/`)
  colnames(mu) <- levels(ks)
  # visualize as (cluster x set) heatmap
  pheatmap(
    mat=t(mu), scale="column", col=pals::coolwarm(), main=.,
    cellwidth=10, cellheight=10, treeheight_row=5, treeheight_col=5)
}


# correlate 'AUCell' signature scores with subpopulation
# proportion estimates from deconvolution with 'RCTD'
cm <- cor(as.matrix(ws), t(assay(auc[top, ])))


pheatmap(cm, 
         col=pals::coolwarm(),
         breaks=seq(-1, 1, length=25),
         cellwidth=10, cellheight=10, 
         treeheight_row=5, treeheight_col=5)




