# This workflow analyzes a 10x Genomics Visium dataset consisting of one sample (Visium capture area) 
# of postmortem human brain tissue from the dorsolateral prefrontal cortex (DLPFC) region, 
# originally described by Maynard et al. (2021).



library(SpatialExperiment)
library(STexampleData)
library(ggspavis)
library(patchwork)
library(scrapper)
library(pheatmap)

#Load sample 151673 from the DLPFC dataset.
spe <- Visium_humanDLPFC()
dim(spe)

#Plot data
plotCoords(spe)

# 1 Quality control 
# subset to keep only spots over tissue
spe <- spe[, spe$in_tissue == 1]
dim(spe)


# identify mitochondrial genes
nms <- rowData(spe)$gene_name
is_mito <- grepl("(^MT-)|(^mt-)", nms)
table(is_mito)

nms[is_mito]


# calculate per-spot QC metrics and store in colData
spe <- quickRnaQc.se(spe, subsets=list(mito=is_mito))


names(colData(spe))

#Select global filtering thresholds for the QC metrics by examining distributions using histograms.

par(mfrow=c(1, 4))
hist(spe$sum, xlab="sum", main="UMIs per spot")
hist(spe$detected, xlab="detected", main="Genes per spot")
hist(spe$subset.proportion.mito, xlab="proportion mito", main="Proportion mito UMIs")
hist(spe$cell_count, xlab="no. cells", main="No. cells per spot")

par(mfrow=c(1, 1))

# select global QC thresholds
spe$qc_lib_size <- spe$sum < 600
spe$qc_detected <- spe$detected < 400
spe$qc_mito <- spe$subset.proportion.mito > 0.28

# tabulate flagged cells
cd <- colData(spe)
qc <- grep("^qc", names(cd))
sapply(cd[qc], table)


#Plot the spatial distributions of the potentially identified low-quality spots

# plot spatial distributions of discarded spots
p1 <- plotObsQC(spe, 
                plot_type="spot", 
                annotate="qc_lib_size") + 
  ggtitle("Library size (< threshold)")
p2 <- plotObsQC(spe, 
                plot_type="spot", 
                annotate="qc_detected") +
  ggtitle("Detected genes (< threshold)")
p3 <- plotObsQC(spe, 
                plot_type="spot", 
                annotate="qc_mito") + 
  ggtitle("Mito proportion (> threshold)")

wrap_plots(p1, p2, p3, nrow=1, guides="collect") & labs(col="discard")

#Select spots to discard by combining the sets of identified low-quality spots according to each metric
# number of identified spots for each metric
ex <- cbind(spe$qc_lib_size, spe$qc_detected, spe$qc_mito)
apply(ex, 2, sum)

# combined set of identified spots
spe$discard <- rowSums(ex) > 0
table(spe$discard)

# Plot the spatial distribution of the combined set of identified low-quality spots to discard, 
# to again confirm that they do not correspond to any clearly biologically meaningful regions

# check spatial pattern of discarded spots
plotObsQC(spe, plot_type="spot", annotate="discard")

#Filter out the low-quality spots

# 2 Normalization
# calculate logcounts using library size factors
spe <- normalizeRnaCounts.se(spe)

summary(sf <- sizeFactors(spe))

hist(sf, breaks=20, main="Histogram of size factors")

assayNames(spe)

# Feature selection (HVGs)

# 3. Apply feature selection methods to identify a set of top highly variable genes (HVGs)

# remove mitochondrial genes
spe <- spe[!is_mito, ]
dim(spe)

# fit mean-variance relationship and select top HVGs
spe <- chooseRnaHvgs.se(
  spe,
  top=ceiling(0.1 * nrow(spe)),
  more.var.args=list(use.min.width=TRUE))

# select top HVGs
hvg <- rownames(spe)[rowData(spe)$hvg]
length(hvg)

# 4.Dimensionality reduction

# using 'scrapper' package
set.seed(123)
spe <- runPca.se(spe, features=hvg, number=50)

colnames(reducedDim(spe, "PCA")) <- paste0(
  "PC", seq_len(ncol(reducedDim(spe, "PCA"))))
spe <- runUmap.se(spe, reddim.type="PCA")
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# embeddings are matrices with
# rows = cells, columns = dims.
sapply(reducedDims(spe), dim)

# 5 Clustering

# graph-based clustering
set.seed(123)
spe <- clusterGraph.se(
  spe,
  num.neighbors=10,
  method="leiden",
  resolution=1,
  reddim.type="PCA",
  output.name="label")
table(spe$label)

# store cluster labels in column 'label' in colData
colLabels(spe) <- factor(spe$label)


# plot cluster labels & annotated reference labels in space
plotCoords(spe, annotate="label", pal="libd_layer_colors") +
plotCoords(spe, annotate="ground_truth", pal="libd_layer_colors")

# plot cluster labels in UMAP dimensions
plotDimRed(spe, plot_type="UMAP", annotate="label", pal="libd_layer_colors")

# 6 Marker genes

# using scrapper package
mgs <- scoreMarkers.se(spe, groups=spe$label)

top <- lapply(mgs, \(df) rownames(df)[df$cohens.d.min.rank <= 2])
length(top <- unique(unlist(top)))

pbs <- aggregateAcrossCells.se(
  spe[top, ], factors=spe$label, assay.type="logcounts")

means <- t(t(assay(pbs, "sums")) / pbs$counts)

# use gene symbols as feature names
mtx <- t(means)
colnames(mtx) <- rowData(pbs)$gene_name

# plot using pheatmap package
pheatmap(mat=mtx, scale="column")

