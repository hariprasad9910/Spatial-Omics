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
# STEP 3: CALCULATE QC METRICS
# ==============================================================================

# subset to keep only spots over tissue
spe <- spe[, spe$in_tissue == 1]
dim(spe)

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)

rowData(spe)$gene_name[is_mito]

# calculate per-spot QC metrics and store in colData
spe <- quickRnaQc.se(spe, subsets=list(mito=is_mito))

head(colData(spe))

# keep copy of object to save later
spe_save <- spe

# ==============================================================================
# STEP 4: GLOBAL OUTLIER DETECTION
# ==============================================================================

# To select thresholds for several QC metrics in our human DLPFC dataset: 
#(i) library size, (ii) number of expressed genes, (iii) proportion of mitochondrial reads, and 
#(iv) number of cells per spot.

par(mfrow=c(1, 4))
hist(spe$sum, xlab="sum", main="UMIs per spot")
hist(spe$detected, xlab="detected", main="Genes per spot")
hist(spe$subset.proportion.mito, xlab="proportion mito", main="Proportion mito UMIs")
hist(spe$cell_count, xlab="no. cells", main="No. cells per spot")

par(mfrow=c(1, 1))


# plot library size vs. number of cells per spot
p1 <- plotObsQC(spe, 
                plot_type="scatter", 
                x_metric="cell_count", 
                y_metric="sum",
                y_threshold=600) + 
  ggtitle("Library size vs. cells per spot")

# plot mito proportion vs. number of cells per spot
p2 <- plotObsQC(spe, 
                plot_type="scatter", 
                x_metric="cell_count", 
                y_metric="subset.proportion.mito",
                y_threshold=0.30) +
  ggtitle("Mito proportion vs. cells per spot")

p1 | p2

# select QC thresholds for library size, 
# detected features & mito. proportion
spe$qc_lib_size <- spe$sum < 600
spe$qc_detected <- spe$detected < 400
spe$qc_mito_prop <- spe$subset.proportion.mito > 0.30

# tabulate number of cells kept/flagged by each
qc <- grep("^qc", names(colData(spe)))
sapply


# check spatial pattern of discarded spots
p1 <- plotObsQC(spe, 
                plot_type="spot", annotate="qc_lib_size") + 
  ggtitle("Library size (< 600 UMI)")

p2 <- plotObsQC(spe, 
                plot_type="spot", annotate="qc_detected") + 
  ggtitle("Detected genes (< 400 genes)")

p3 <- plotObsQC(spe, 
                plot_type="spot", annotate="qc_mito_prop") + 
  ggtitle("Mito proportion (> 0.30)")

p1 | p2 | p3


# check spatial pattern of discarded spots if threshold is too high
spe$qc_lib_size_2000 <- spe$sum < 2000

# plot the spots flagged with the high threshold
p1 <- plotObsQC(spe, 
                plot_type="spot", annotate="qc_lib_size_2000") + 
  ggtitle("Library size (< 2000 UMI)")

# plot manually annotated reference layers
p2 <- plotCoords(spe, 
                 annotate="ground_truth", pal="libd_layer_colors") + 
  ggtitle("Manually annotated layers")

# plot library size by manual annotation
p3 <- plotColData(spe, 
                  x="ground_truth", y="sum", colour_by="ground_truth") + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ggtitle("Library size by layer") + xlab("")

p1 | p2 | p3

# library size and outliers
p1 <- plotObsQC(spe, 
                plot_type="violin", x_metric="sum", 
                annotate="qc_lib_size", point_size=0.5) + 
  xlab("Library size")

# detected genes and outliers
p2 <- plotObsQC(spe, 
                plot_type="violin", x_metric="detected", 
                annotate="qc_detected", point_size=0.5) + 
  xlab("Detected genes") 

# mito proportion and outliers
p3 <- plotObsQC(spe, 
                plot_type="violin", x_metric="subset.proportion.mito",
                annotate="qc_mito_prop", point_size=0.5) +
  xlab("Mito proportion")

p1 | p2 | p3

# ==============================================================================
# STEP 5: LOCAL OUTLIER DETECTION
# ==============================================================================

# detect local outliers based on library size, unique genes, mito. proportion
spe <- localOutliers(spe, metric="sum", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="detected", direction="lower", log=TRUE)
spe <- localOutliers(spe, metric="subset.proportion.mito", direction="higher", log=FALSE)

# spot plot of log-transformed library size
p1 <- plotCoords(spe, 
                 annotate="sum_log") + 
  ggtitle("log2(Library Size)")

p2 <- plotObsQC(spe, 
                plot_type="spot", in_tissue="in_tissue", 
                annotate="sum_outliers", point_size=0.2) + 
  ggtitle("Local Outliers (Library Size)")

# spot plot of log-transformed detected genes
p3 <- plotCoords(spe, 
                 annotate="detected_log") + 
  ggtitle("log2(Detected)")

p4 <- plotObsQC(spe, 
                plot_type="spot", in_tissue="in_tissue", 
                annotate="detected_outliers", point_size=0.2) + 
  ggtitle("Local Outliers (Detected)")

# spot plot of mitochondrial proportion
p5 <- plotCoords(spe, 
                 annotate="subset.proportion.mito") +
  ggtitle("Mito Proportion")

p6 <- plotObsQC(spe, 
                plot_type="spot", in_tissue="in_tissue", 
                annotate="subset.proportion.mito_outliers", point_size=0.2) +
  ggtitle("Local Outliers (Mito Prop)")

# plot using patchwork
(p1 / p2) | (p3 / p4) | (p5 / p6)


# z-transformed library size and outliers
p1 <- plotObsQC(spe, 
                plot_type="violin", x_metric="sum_z", 
                annotate="sum_outliers", point_size=0.5) + 
  xlab("sum_outliers")

# z-transformed detected genes and outliers
p2 <- plotObsQC(spe, 
                plot_type="violin", x_metric="detected_z", 
                annotate="detected_outliers", point_size=0.5) + 
  xlab("detected_outliers")

# z-transformed mito proportion and outliers
p3 <- plotObsQC(spe, 
                plot_type="violin", x_metric="subset.proportion.mito_z",
                annotate="subset.proportion.mito_outliers", point_size=0.5) +
  xlab("mito_outliers")

p1 | p2 | p3

# ==============================================================================
# STEP 6: REMOVE LOW QUALITY SPOTS 
# ==============================================================================

# select updated threshold for mito proportion
spe$qc_mito <- spe$subset.proportion.mito > 0.28
table(spe$qc_mito)

# combine global/local outliers
spe$global_outliers <- 
  spe$qc_lib_size | 
  spe$qc_detected | 
  spe$qc_mito
spe$local_outliers <- 
  spe$sum_outliers | 
  spe$detected_outliers | 
  spe$subset.proportion.mito_outliers

rbind( # tabulate kept/flagged cells
  global=table(spe$global_outliers),
  local=table(spe$local_outliers))


# check spatial pattern of combined set of discarded spots
plotObsQC(spe, plot_type="spot", annotate="global_outliers") +
  plotObsQC(spe, plot_type="spot", annotate="local_outliers")


# combine local and global outliers &
# remove combined set of low-quality spots
spe$discard <- 
  spe$global_outliers | 
  spe$local_outliers
spe <- spe[, !spe$discard]

# remove features with all 0 counts
spe <- spe[rowSums(counts(spe)) > 0, ]
dim(spe)
