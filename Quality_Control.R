# ==============================================================================
# STEP 1: LOAD NECESSARY BIOCONDUCTOR & CRAN PACKAGES
# ==============================================================================

library(STexampleData)      # Contains curated spatial transcriptomics datasets
library(ggspavis)           # Visualization tools specifically for spatial data
library(ggplot2)            # Standard data visualization engine
library(scater)             # Single-cell analysis and QC utilities
library(scrapper)           # High-throughput RNA quality control functions
library(SpotSweeper)        # Advanced local outlier and artifact detection
library(patchwork)          # Arranging and aligning multiple ggplot panels
library(SpatialExperiment)  # Core data structure for spatial data representation

# ==============================================================================
# STEP 2: LOAD RUNTIME DATASET (10x GENOMICS VISIUM DLPFC)
# ==============================================================================

# Import postmortem human dorsolateral prefrontal cortex dataset
spe <- Visium_humanDLPFC()

# Retain a baseline copy of the original object to save at the workflow's conclusion
spe_save <- spe

# ==============================================================================
# STEP 3: INITIAL FILTRATION & MITOCHONDRIAL METRIC COMPUTATION
# ==============================================================================

# Subset the experiment to retain only spots over the physical tissue section.
# Background spots capture cellular debris/ambient RNA and distort calculations.
spe <- spe[, spe$in_tissue == 1]

# Identify mitochondrial genes using regular expressions matching 'MT-' or 'mt-' prefixes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)

# Compute per-spot QC metrics (Library Size [sum], Unique Genes [detected], and Mito Proportion)
# We handle spots identically to individual cells at this operational phase.
spe <- quickRnaQc.se(spe, subsets = list(mito = is_mito))


# ==============================================================================
# STEP 4: GLOBAL THRESHOLD APPLICATION & PLOTTING
# ==============================================================================

# Define basic global thresholds based on initial exploratory visualizations:
# Spots with fewer than 600 UMIs, < 400 genes, or > 30% mitochondrial reads are flagged.
spe$qc_lib_size   <- spe$sum < 600
spe$qc_detected   <- spe$detected < 400
spe$qc_mito_prop  <- spe$subset.proportion.mito > 0.30

# Optional: Visualize global diagnostic metrics against total tissue cell count
p1 <- plotObsQC(spe, plot_type = "scatter", x_metric = "cell_count", y_metric = "sum", y_threshold = 600)
p2 <- plotObsQC(spe, plot_type = "scatter", x_metric = "cell_count", y_metric = "subset.proportion.mito", y_threshold = 0.30)
p1 | p2

# Generate spatial mapping plots to ensure global filters are not systematically 
# stripping out distinct biological anatomical features (e.g., specific cortical layers).
p_spatial_lib   <- plotObsQC(spe, plot_type = "spot", annotate = "qc_lib_size")
p_spatial_genes <- plotObsQC(spe, plot_type = "spot", annotate = "qc_detected")
p_spatial_mito  <- plotObsQC(spe, plot_type = "spot", annotate = "qc_mito_prop")
p_spatial_lib | p_spatial_genes | p_spatial_mito

# ==============================================================================
# STEP 5: NEIGHBORHOOD-BASED LOCAL OUTLIER DETECTION (via SpotSweeper)
# ==============================================================================

# Run local outlier assessment by evaluating spots against immediate surrounding tissue coordinates.
# For Visium's hexagonal layout, k = 36 computes 3 full rings of concentric neighbors.
# Note: For square grid arrangements (e.g., STOmics, Visium HD), use k = 48.
# Log transformations are applied to count parameters to force normal distributions.
spe <- localOutliers(spe, metric = "sum", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "detected", direction = "lower", log = TRUE)
spe <- localOutliers(spe, metric = "subset.proportion.mito", direction = "higher", log = FALSE)

# Generate comparative diagnostic visualizations matching log features against calculated local outliers
p_lib_log <- plotCoords(spe, annotate = "sum_log")
p_lib_out <- plotObsQC(spe, plot_type = "spot", in_tissue = "in_tissue", annotate = "sum_outliers", point_size = 0.2)
(p_lib_log / p_lib_out)

# ==============================================================================
# STEP 6: SPOT SELECTION ACCORDING TO CELL DENSITY AND HISTOLOGY
# ==============================================================================

# Spots with excessively high cell counts coupled with low gene detection are indicative
# of physical compression, tissue damage, or severe computational segmentation error.
# We establish an empirical cap of 10 cells per spot for standard cortical sections.
spe$qc_cell_count <- spe$cell_count > 10

# Crucial Architectural Rule: Do NOT blindly delete spots with a cell count of 0!
# In brain tissue, 0-cell spots map to neuropil (dense networks of axons and dendrites) 
# which contain distinct, biologically essential subcellular transcript profiles.

# ==============================================================================
# STEP 7: MULTISCALE LOCAL VARIANCE SCANNING AND HANGNAIL ARTIFACT REMOVAL
# ==============================================================================

# Load a dedicated verification dataset containing a confirmed hangnail artifact
data(DLPFC_artifact)
spe.hangnail <- DLPFC_artifact

# Flag spatial mitochondrial mapping profiles on the testing object
is_mito_hn   <- grepl("(^MT-)|(^mt-)", rowData(spe.hangnail)$gene_name)
spe.hangnail <- quickRnaQc.se(spe.hangnail, subsets = list(mito = is_mito_hn))

# Quantify local mitochondrial variance to map structural anomalies
spe.hangnail <- localVariance(spe.hangnail, n_neighbors = 36, 
                              metric = "subset.proportion.mito", name = "local_mito_variance_k36")

# Execute multiscale classification tool across 7 concentric rings (n_order = 7).
# Note: Ensure shape corresponds exactly to your spatial physical substrate configuration.
spe.hangnail <- findArtifacts(spe.hangnail, mito_percent = "expr_chrM_ratio", 
                              mito_sum = "expr_chrM", n_order = 7, 
                              shape = "hexagonal", name = "artifact")

# Perform hard filtration of the technical artifact layer
spe.hangnail <- spe.hangnail[, !spe.hangnail$artifact]

# ==============================================================================
# STEP 8: HARVEST DATASET AND EXPORT FOR DOWNSTREAM COMPUTE PIPELINES
# ==============================================================================

# Tune final mitochondrial cutoff filter threshold
spe$qc_mito <- spe$subset.proportion.mito > 0.28

# Combine distinct logical vectors across all outlier methodologies
spe$global_outliers <- spe$qc_lib_size | spe$qc_detected | spe$qc_mito | spe$qc_cell_count
spe$local_outliers  <- spe$sum_outliers | spe$detected_outliers | spe$subset.proportion.mito_outliers

# Construct final intersection vector for filtration processing
spe$discard <- spe$global_outliers | spe$local_outliers

# Slice out all poor quality/technical artifact coordinates from dataset matrix
spe <- spe[, !spe$discard]

# Drop clean unexpressed genes (rows containing entirely zero count occurrences across spots)
spe <- spe[rowSums(counts(spe)) > 0, ]

# Apply identical processing steps to baseline storage copy and write data to disk
ex <- rowSums(counts(spe_save)) != 0
spe_save <- spe_save[ex, ]
saveRDS(spe_save, "seq-spe_qc.rds")


