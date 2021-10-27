library(SingleCellExperiment)
library(scater)
library(Seurat)
library(ggplot2)

## load data, it will load two objects counts table: "count.raw.f" and meta data table: "cell.annot.f"
load("./rawData/nature_brain3008_raw_annot.RData")
# assign row names for meta data table
row.names(cell.annot.f)<-cell.annot.f$cell
# create SingleCellExperiment(sce) object
# about ingleCellExperiment(sce) object, further reading: https://bioconductor.org/packages/3.13/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html 
dat<-SingleCellExperiment(assays=list(counts = as.matrix(count.raw.f)),
                          colData=cell.annot.f)
## add per-cell metrics, add mitochondrial gene information
dat <- addPerCellQC(dat, subsets=list(Mito=grep("mt-", rownames(dat))))

# cell QC
# histogram for gene number 
hist(colData(dat)[,"detected"],breaks=100,
     xlab="gene number per cell",main="Histogram of numbere of genes per cell")

# histogram for reads counts
# disable scientific notation
options(scipen=999)
hist(colData(dat)[,"sum"],breaks=100,
     xlab="read counts per cell",main="Histogram of read counts per cell")

## detect outliers
## Sometimes it’s hard to come up with an obvious filtering cutoff. 
#In this case, adaptive threshold can help identify points that are more than the number  (e.g. 3) of median absolute deviations (MADs) away from the median in any of the variables we use for QC. 

# check library lower than 3 MADs (default)
qc.lib2 <- isOutlier(dat$sum, log=TRUE, type="lower")
attr(qc.lib2, "thresholds")

hist(colData(dat)[,"sum"],breaks=100,
     xlab="read counts per cell",main="Histogram of read counts per cell")
abline(v=155488,col="red")

# check gene numbers lower than 3 MADs
qc.nexprs2 <- isOutlier(dat$detected, log=TRUE, type="lower")
attr(qc.nexprs2, "thresholds")

hist(colData(dat)[,"detected"],breaks=100,
     xlab="gene number per cell",main="Histogram of numbere of genes per cell")
abline(v=2465.915,col="red")

qc.mito2 <- isOutlier(dat$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")

## scatter plot for total read counts per cell vs. proportion of mitochondrial genes
dat$qc.mito2<-qc.mito2
plotColData(dat, x = "sum", y="subsets_Mito_percent",colour_by="qc.mito2",add_legend=FALSE) 

## identifies low-quality cells based on frequently used QC metrics (from "isOutlier" command):
reasons <- quickPerCellQC(colData(dat), sub.fields=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))

## add discard into meta data 
dat$discard <- reasons$discard

# plot discard and remain cells
plotColData(dat, x="sum", y="detected", colour_by="discard")

## filtering features
#Let’s keep the genes which were detected (expression value > 1) in 2 or more cells. 
keep_feature <- nexprs(dat,byrow = TRUE,detection_limit = 1) >= 2
rowData(dat)$discard <- ! keep_feature
table(rowData(dat)$discard)

# After QC diagnosis, filtering out poor QC cells and genes
datPreQc <- dat
dat <- datPreQc[keep_feature,!datPreQc$discard]

# convert sce object to Seurat object for downstream analysis
dat.seurat<-as.Seurat(dat,counts = "counts",data="counts")
## normalization
dat.seurat <- NormalizeData(dat.seurat, normalization.method = "LogNormalize", scale.factor = 1000000)
## identify variable genes
dat.seurat <- FindVariableFeatures(dat.seurat, selection.method = "vst", nfeatures = 2000)
# scaling data
all.genes <- rownames(dat.seurat)
dat.seurat <- ScaleData(dat.seurat, features = all.genes)
# PCA reduction
dat.seurat <- RunPCA(dat.seurat, features = VariableFeatures(object = dat.seurat))
## select number of components
ElbowPlot(dat.seurat,ndims = 30)

## clustering
dat.seurat <- FindNeighbors(dat.seurat, dims = 1:20)
dat.seurat <- FindClusters(dat.seurat, resolution = 0.5)
dat.seurat <- RunUMAP(dat.seurat, dims = 1:20)
DimPlot(dat.seurat,label = T)&
  theme(
    text=element_text(size=8)
  )
# identify marker genes in each cluster
cluster.markers<-FindAllMarkers(dat.seurat,logfc.threshold = 0.5)

# sessionInfo
sessionInfo()