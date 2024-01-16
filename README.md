# AP-mice
#AP&amp; INK mice/ snRNA-seq
#the project has been done on Compute Canada/scRNA-seq of skin mice with physical activity
$ module load StdEnv/2020 r/4.2.2
 library(Seurat)
 library(SoupX)
 library(sctransform)
 library(harmony)
 library(scDblFinder)
 library(tidyverse)
 library(dplyr)
 library(ggplot2)
 
 #soupX_Ambient RNA correction
 toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED10_MPS12345537_E08_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED10_MPS12345537_E08_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED10_MPS12345537_E08_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster)
sc = autoEstCont(sc) 
WT_SED10_MPS12345537_E08_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP2_MPS12345537_C07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP2_MPS12345537_C07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP2_MPS12345537_C07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_AP2_MPS12345537_C07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP10_MPS12345537_D08_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP10_MPS12345537_D08_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP10_MPS12345537_D08_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_AP10_MPS12345537_D08_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP3_MPS12345537_A07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP3_MPS12345537_A07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP3_MPS12345537_A07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_AP3_MPS12345537_A07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP5_MPS12345537_F07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP5_MPS12345537_F07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP5_MPS12345537_F07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_AP5_MPS12345537_F07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP6_MPS12345537_H07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP6_MPS12345537_H07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP6_MPS12345537_H07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_AP6_MPS12345537_H07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP9_MPS12345537_B08_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP9_MPS12345537_B08_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-AP9_MPS12345537_B08_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_AP9_MPS12345537_B08_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED2_MPS12345537_D07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED2_MPS12345537_D07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED2_MPS12345537_D07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_SED2_MPS12345537_D07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED3_MPS12345537_B07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED3_MPS12345537_B07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED3_MPS12345537_B07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_SED3_MPS12345537_B07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED5_MPS12345537_G07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED5_MPS12345537_G07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED5_MPS12345537_G07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_SED5_MPS12345537_G07_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED6_MPS12345537_A08_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED6_MPS12345537_A08_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED6_MPS12345537_A08_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_SED6_MPS12345537_A08_9852 = adjustCounts(sc)

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED9_MPS12345537_C08_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED9_MPS12345537_C08_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/WT-SED9_MPS12345537_C08_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
WT_SED9_MPS12345537_C08_9852 = adjustCounts(sc)

#creating the object
APskin_samples <- CreateSeuratObject(counts = cbind(WT_SED10_MPS12345537_E08_9852 , WT_AP2_MPS12345537_C07_9852, WT_AP10_MPS12345537_D08_9852, WT_AP3_MPS12345537_A07_9852, WT_AP5_MPS12345537_F07_9852, WT_AP6_MPS12345537_H07_9852, WT_AP9_MPS12345537_B08_9852, WT_SED2_MPS12345537_D07_9852, WT_SED3_MPS12345537_B07_9852, WT_SED5_MPS12345537_G07_9852, WT_SED6_MPS12345537_A08_9852, WT_SED9_MPS12345537_C08_9852, project = "APmouseskin", min.cells = 5)

#specifying the different condition
APskin_samples@meta.data$sample <- c(rep("SED10", ncol(WT_SED10_MPS12345537_E08_9852)), rep("AP2", ncol(WT_AP2_MPS12345537_C07_9852)), rep("AP10", ncol(WT_AP10_MPS12345537_D08_9852)), rep("AP3", ncol(WT_AP3_MPS12345537_A07_9852)), rep("AP5", ncol(WT_AP5_MPS12345537_F07_9852)), rep("AP6", ncol(WT_AP6_MPS12345537_H07_9852)), rep("AP9", ncol(WT_AP9_MPS12345537_B08_9852)), rep("SED2", ncol(WT_SED2_MPS12345537_D07_9852)), rep("SED3", ncol(WT_SED3_MPS12345537_B07_9852)), rep("SED5", ncol(WT_SED5_MPS12345537_G07_9852)), rep("SED6", ncol(WT_SED6_MPS12345537_A08_9852)), rep("SED9", ncol(WT_SED9_MPS12345537_C08_9852)))

APskin_samples@meta.data$condition <- c(rep("SED", ncol(WT_SED10_MPS12345537_E08_9852)), rep("PhA", ncol(WT_AP2_MPS12345537_C07_9852)), rep("PhA", ncol(WT_AP10_MPS12345537_D08_9852)), rep("PhA", ncol(WT_AP3_MPS12345537_A07_9852)), rep("PhA", ncol(WT_AP5_MPS12345537_F07_9852)), rep("PhA", ncol(WT_AP6_MPS12345537_H07_9852)), rep("PhA", ncol(WT_AP9_MPS12345537_B08_9852)), rep("SED", ncol(WT_SED2_MPS12345537_D07_9852)), rep("SED", ncol(WT_SED3_MPS12345537_B07_9852)), rep("SED", ncol(WT_SED5_MPS12345537_G07_9852)), rep("SED", ncol(WT_SED6_MPS12345537_A08_9852)), rep("SED", ncol(WT_SED9_MPS12345537_C08_9852)))

#add a column to the metadata
APskin_samples[["percent.mt"]] <- PercentageFeatureSet(APskin_samples, pattern = "^mt-") 
APskin_samples <- PercentageFeatureSet(APskin_samples, pattern = "^mt-", col.name = "percent.mt")

#save the object
saveRDS(APskin_samples, file = "APskinMice.rds")

# QC steps
library(Seurat) 
library(sctransform) 
library(harmony) 
library(scDblFinder) 
library(tidyverse) 
library(dplyr) 
library(ggplot2)

APskin_samples <- readRDS("APskinMice.rds")

# Visualize QC metrics as a violin plot

VlnPlot(APskin_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(APskin_samples, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(APskin_samples, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

#changing and adding to the metadata to caracterize quality control metrics

#This is a common normalization method used to account for differences in sequencing depth between samples. The resulting values are often referred to as "counts per million" (CPM) or "reads per kilobase per million" (RPKM), and are used to compare gene expression levels between samples.
APskin_samples$log10GenesPerUMI <- log10(APskin_samples$nFeature_RNA) / log10(APskin_samples$nCount_RNA)

# Create metadata dataframe
metadata <- APskin_samples@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
                      
APskin_samples@meta.data <- metadata

saveRDS(APskin_samples, file = "APskinMice.rds")

# Visualize the number of cell counts per sample
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
   
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
   
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  	ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  	geom_boxplot() + 
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells vs NGenes")
   
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)
   
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=sample, x=percent.mt, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)
  
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)
#we can repeat QC visualization steps after filter step to see the difference that filtering dataset will make 

##Filter
##orig.ident: this often contains the sample identity if known, but will default to project as we had assigned it
#nCount_RNA: number of UMIs per cell(nUMI)
#nFeature_RNA: number of genes detected per cell(nGene)
# cell_level:Filter out low quality reads using selected thresholds - these will change with experiment
#logGENE/UMI :The number of genes per UMI for each cell is quite easy to calculate, and we will log10 transform the result for better comparison between samples
APskin_samplesF <- subset(x = APskin_samples, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (percent.mt < 20))

#gene_level: # Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = APskin_samplesF, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
APskin_samplesF <- CreateSeuratObject(filtered_counts, meta.data = APskin_samplesF@meta.data)

saveRDS(APskin_samplesF, file = "APskinMiceF.rds")

#normalising
#perform normalization & UMAP before finding doublets
#normalizing(in compute canada)

APskinMiceF <- SCTransform(APskinMiceF, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2")

APskinMiceF <- RunPCA(APskinMiceF, verbose = FALSE)

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = APskinMiceF, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = APskinMiceF, features = "PC_1", group.by = "sample", pt.size = .1)
p1 + p2

ElbowPlot(APskinMiceF)

options(repr.plot.height = 2.5, repr.plot.width = 6)
APskinMiceF <- APskinMiceF %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT")

harmony_embeddings <- Embeddings(APskinMiceF, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = APskinMiceF, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = APskinMiceF, features = "harmony_1", group.by = "sample", pt.size = .1)
p1 + p2

APskinMiceF <- APskinMiceF %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

saveRDS(APskinMiceF, file = "APskinMiceFn.rds")
#perform UMAP before doublet finding to have the original one for comparing to after doublet correction one

DimPlot(APskinMiceFn, reduction = "umap", label = TRUE, label.size = 5)


APskinMiceFn <- APskinMiceFn %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()
##doublet correction

#filter out doublets: scDblFinder-->integrated data (note: if you have very large differences in number of cells between samples the scores will not be directly comparable. We are working on improving this, but in the meantime it would be preferable to stratify similar samples and threshold the sets separately)

## load libraries
#BiocManager::install("plger/scDblFinder")
#BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

##finding doublets
DefaultAssay(APskinMiceFn) <- "RNA"
sce <- as.SingleCellExperiment(APskinMiceFn)
set.seed(2022)
sce <- scDblFinder(sce, clusters = "seurat_clusters", samples="sample")
sce$scDblFinder.score
sce$scDblFinder.class
APskinMiceFn$scDblFinder.score <- sce$scDblFinder.score
FeaturePlot(APskinMiceFnv, features = "scDblFinder.score")
table(sce$scDblFinder.class)
metadata(sce)$scDblFinder.stats

#add scDblFinder.class column from sce in soup_skin_samples object
APskinMiceFn@meta.data$isDoublet <- sce$scDblFinder.class

#Doublets visualization
DimPlot(APskinMiceFn, reduction = 'umap', group.by = "isDoublet")

#create a new object(fnd_soup_skin_samples) and delete the doublet
APskinMiceFn <- subset(APskinMiceFn, isDoublet == "doublet")
DimPlot(APskinMiceFn, reduction = 'umap', group.by = "isDoublet")
DimPlot(APskinMiceFn, reduction = 'umap', label = TRUE, label.size = 3)

#save new object
saveRDS(APskinMiceFn, file = "C:/Users/pegah/Desktop/skinProject/fromTheTop/soupCorrection/APskinMiceFnd.rds")

# after doublet correction we normalize again
APskinMiceFnd <- SCTransform(APskinMiceFn, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2")

APskinMiceFnd <- RunPCA(APskinMiceFnd, verbose = FALSE)

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = APskinMiceFnd, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = APskinMiceFnd, features = "PC_1", group.by = "sample", pt.size = .1)
p1 + p2

ElbowPlot(APskinMiceFnd)

options(repr.plot.height = 2.5, repr.plot.width = 6)
APskinMiceFnd <- APskinMiceFnd %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT")

harmony_embeddings <- Embeddings(APskinMiceFnd, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = APskinMiceFnd, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = APskinMiceFnd, features = "harmony_1", group.by = "sample", pt.size = .1)
p1 + p2

APskinMiceFnd <- APskinMiceFnd %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

saveRDS(APskinMiceFnd, file = "APskinMiceFndn.rds")



#to find a better resolution for clustering we try to find markers of different resolutions and compare.The best resolution is the one that clusters for resolution less or more than that do not show major different in annotation.
#look into database enrichr for anotation the cluster with specified resolution
#res0.04
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(enrichR)
library(SingleCellExperiment)
library(Matrix)

ClusterMarkers <- read.csv("APskinMiceFndn.rds")
head(ClusterMarkers)

setEnrichrSite("Enrichr") # Human and mouse genes

# list of all the databases

dbs <- listEnrichrDbs()
# this will list the possible libraries
dbs

# select libraries with cell types
db <- c('CellMarker_Augmented_2021',
        'Azimuth_Cell_Types_2021',
        'DSigDB',
        'Mouse_Gene_Atlas',
        'KEGG_2019_Mouse',
        'Reactome_2022',
        'Azimuth_2023',
        'GO_Biological_Process_2023',
        'GO_Cellular_Component_2023',
        'GeneSigDB',
        'GO_Molecular_Function_2023',
        'HDSigDB_Mouse_2021',
        'IDG_Drug_Targets_2022',
        'PanglaoDB_Augmented_2021',
        'Tabula_Muris',
        'RNAseq_Automatic_GEO_Signatures_Mouse_Down',
        'RNAseq_Automatic_GEO_Signatures_Mouse_Up',
        'TF_Perturbations_Followed_by_Expression')

#https://maayanlab.cloud/Enrichr/#libraries Aging_Perturbations_from_GEO_down, Aging_Perturbations_from_GEO_up

#Here is a small function to run easily on each cluster and find the cell type library predictions
checkCelltypes <- function(cluster_num = 0){
  clusterX <- ClusterMarkers %>% filter(cluster == cluster_num & avg_log2FC > 0.25)
  genes <- clusterX$gene
  # the cell type libraries
  # get the results for each library
  clusterX.cell <- enrichr(genes, databases = db)
  # visulize the results
  print(plotEnrich(clusterX.cell[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'CellMarker_Augmented_2021'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'Azimuth_Cell_Types_2021'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'DSigDB'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'Mouse_Gene_Atlas'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'KEGG_2019_Mouse'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'Reactome_2022'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'Azimuth_2023'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'GO_Biological_Process_2023'))
  print(plotEnrich(clusterX.cell[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'GO_Cellular_Component_2023'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'GeneSigDB'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'GO_Molecular_Function_2023'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'HDSigDB_Mouse_2021'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'IDG_Drug_Targets_2022'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'PanglaoDB_Augmented_2021'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'Tabula_Muris'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'RNAseq_Automatic_GEO_Signatures_Mouse_Down'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'RNAseq_Automatic_GEO_Signatures_Mouse_Up'))
  print(plotEnrich(clusterX.cell[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", title = 'TF_Perturbations_Followed_by_Expression'))
  
}

cluster0 <- checkCelltypes(cluster_num = 0)
cluster1 <- checkCelltypes(cluster_num = 1) 
cluster2 <- checkCelltypes(cluster_num = 2) 
#it has to be run for each cluster ....
cluster16 <- checkCelltypes(cluster_num = 16)



#chose the best res and assign cluster names/ in our dataset res=0.04 was the best
APskin_seurat0.4 <- readRDS("APskinMiceFndn.rds")
APskin_seurat0.4 <- APskin_seurat0.4 %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

#DimPlot(APskin_seurat0.4, reduction = "umap", label = TRUE, pt.size = 1.2)
## to remove krt+basal cells
#APskin_seurat0.4wokrt <- readRDS("APskinMiceFndn.rds")

APskin_seurat0.4wokrt <- subset(APskin_seurat0.4, subset = seurat_clusters %in% c('2','3','4','6','7','8','9','11','13','16','17','18','21'))
DimPlot(APskin_seurat0.4wokrt, reduction = "umap", label = TRUE, pt.size = 1.2)

APskin_seurat0.4wokrt <- APskin_seurat0.4wokrt %>%
  NormalizeData() %>%
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()
APskin_seurat0.4wokrt.markers <- FindAllMarkers(APskin_seurat0.4wokrt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(APskin_seurat0.4wokrt.markers, "APskin_seurat0.4wokrt.markers.csv")

DimPlot(APskin_seurat0.4wokrt, reduction = "umap", label = TRUE, label.size = 5)
saveRDS(APskin_seurat0.4wokrt, file = "APskin_seurat0.4wokrt.rds")

APskin_seurat0.4wokrt <- readRDS("APskin_seurat0.4wokrt.rds")

## to visualize top 5 markers
ClusterMarkers <- read.csv("APskin_seurat0.4wokrt.markers.csv")
head(ClusterMarkers)

library(dplyr)

top5 <- ClusterMarkers %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  labs(title = "wokrt_0.4")

genes_to_plot <- top5[[1]]$gene
genes_to_plot <- unique(genes_to_plot)
p <- DotPlot(APskin_seurat0.4wokrt, features = genes_to_plot)
p <- p + labs(title = top5$title)
# Rotate x-axis labels to be vertical
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
###assign names to clusters
new.cluster.ids <- c("Macrophage", "Fibroblast", "Endothelial","Endothelial", "Fibroblast", "Macrophage", "Fibroblast", "Endothelial", "T cell", "Macrophage", "Mast cell", "Macrophage", "Smooth muscle cell", "Fibroblast")
names(new.cluster.ids) <- levels(APskin_seurat0.4wokrt)
APskin_seurat0.4wokrt <- RenameIdents(APskin_seurat0.4wokrt, new.cluster.ids)
DimPlot(APskin_seurat0.4wokrt, reduction = "umap", label = TRUE, pt.size = 1.5, label.size = 3)
saveRDS(APskin_seurat0.4wokrt, file = "APskin_seurat0.4wokrtAnnotated.rds")


#######key gene identification / replace the name of the gene of interest in features = (") 
########feature plot for specific gene

#VlnPlot(APskin_seurat0.4wokrt, features = ("Krt14")) + NoLegend()
#FeaturePlot(APskin_seurat0.4wokrt, features = ("Krt14"))#, order = T)


ggsave(filename = "pwokrt_0.4.pdf", plot = p, width=30, height=10, scale=0.7)

## The next step will be sub-clustering clusters to clean and validate the clusters of interest/I brought the sub-clustering of lymphatic endothelial cells in this script. the process should repeat for each cluster.
##to remove some clusters(subclustering)
##to subset target clusters and remove others
#you can subset a cluster to subcluster to

APskin_seurat0.4wokrtc_lEC <- subset(APskin_seurat0.4wokrt, subset = seurat_clusters %in% c('2','7'))

APskin_seurat0.4wokrtc_lEC <- APskin_seurat0.4wokrtc_lEC %>% 
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, vst.flavor = "v2") %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% 
  RunHarmony("sample", plot_convergence = TRUE, assay.use = "SCT", reduction.save = "harmony2") %>% 
  RunUMAP(reduction = "harmony2", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony2", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(APskin_seurat0.4wokrtc_lEC, reduction = "umap", label = TRUE, label.size = 5) 
APskin_seurat0.4wokrtc_lEC.markers <- FindAllMarkers(object = APskin_seurat0.4wokrtc_lEC, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) 
write.csv(APskin_seurat0.4wokrtc_lEC.markers, "APskin_seurat0.4wokrtc_lEC.markers.csv") 
saveRDS(APskin_seurat0.4wokrtc_lEC, file = "APskin_seurat0.4wokrtc_lEC.rds")


APskin_seurat0.4wokrtc_lEC <- readRDS("APskin_seurat0.4wokrtc_lEC.rds")
#pseudo-bulk workflow -----------------
  # Acquiring necessary metrics for aggregation across cells in a sample
  # 1. counts matrix - sample level
  # counts aggregate to sample level
  seu.filtered <- readRDS("APskin_seurat0.4wokrtAnnotated.rds")
DimPlot(seu.filtered, reduction = "umap", label = TRUE, pt.size = 1.2)

View(seu.filtered@meta.data)

# Extract metadata to create new object

metadata <- seu.filtered@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$new.cluster.ids <- factor(seu.filtered@active.ident)
metadata$samples <- paste0(seu.filtered$condition, seu.filtered$sample)

View(metadata)
seu.filtered@meta.data <- metadata

DefaultAssay(seu.filtered)
View(seu.filtered)

cts <- AggregateExpression(seu.filtered, 
                           group.by = c("new.cluster.ids", "samples"),
                           assays = 'SCT',
                           slot = "counts",
                           return.seurat = FALSE)


cts <- cts$SCT

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split

splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame

cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

#gsub('.*_(.*)', '\\1', '2-EC_ntrs39Q')

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

head(cts.split.modified)
###to do DESeq2
#1.get the count metrix

counts_EC <- cts.split.modified$`Endothelial`

#2.generate sample level matadata

colData <- data.frame(samples = colnames(counts_EC), row.names = NULL)
#get condition column
colData <- colData %>%
  mutate(condition = c('AP','AP','AP','AP','AP','AP','ctrl','ctrl','ctrl','ctrl','ctrl','ctrl'))%>%
  column_to_rownames(var = 'samples')

colData

###Running DESeq2
#creating DESeq2 object

dds <- DESeqDataSetFromMatrix(countData = counts_EC,
                              colData = colData,
                              design = ~ condition)

#filter all the genes that have lower than 10 reads
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

#running DESeq

dds <- DESeq(dds)

#seeing the results
resultsNames(dds)#showing us the different groups that we can compare

res <- results(dds, name = "condition_ctrl_vs_AP")
res

# output tables
write.table(res, file="DESeqECmainAP.csv", quote=FALSE, sep=",")

#to visualize

library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)

saveRDS(res, 'APmainECddsRES.rds')
de <- readRDS('APmainECddsRES.rds')

#### Create a volcano plot in ggplot2 to show significantly expressed genes
# Convert DESeq2 results to a data frame
de <- data.frame(
  gene = rownames(de),
  log2FoldChange = de$log2FoldChange,
  padj = de$padj,
  pvalue = de$pvalue
)

# Display the first few rows of the resulting data frame
head(de)

#to visualize 

p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$pvalue < 0.05 & de$log2FoldChange > 0.6] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[ de$pvalue < 0.05 & de$log2FoldChange < -0.6] <- "DOWN"
# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text(data = subset(de, padj <= 0.05), aes(label = gene), vjust = 1, hjust = 1, size = 4)

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

###Pathway analysis
##################################    ##########################################

library(Seurat)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("enricher")
library(enrichR)
library(dplyr)
library(patchwork)
library(sctransform)
library(harmony)
library(tidyverse)
library(ggplot2)
install.packages("scDblFinder")
library(scDblFinder)
library(SingleCellExperiment)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ExperimentHub")
library(ExperimentHub)
library(DESeq2)
library(pheatmap)
library(SummarizedExperiment)
library(RColorBrewer)
library(msigdbr)
if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("DESeq2")
library(clusterProfiler)
library(fgsea)
library(data.table)
library(parallel)

##enrichment analysis
##load data
EC <- read.csv(file = "DESeqECmainAP.csv" , header = T, sep=",")

##Enrichment##
#get gene database
#all gene

all_gene_sets = msigdbr(species = "Mus musculus")
head(all_gene_sets)

#halmark
h_gene_sets = msigdbr(species = "mouse", category = "H")
head(h_gene_sets)

#to check if it is a data.frame
class(h_gene_sets)

#to define significant genes
# to look at padj

ggplot(EC, aes(x=padj)) +
  geom_histogram(bins = 10) +
  theme_classic() 

#to see if we have padj==0
table(EC$padj == 0)

#cutoff non-significant pvalue since padj significancy is not a proper treshold here 
signif <- EC %>%
  filter(pvalue <= 0.05)

table(signif$pvalue >= 0.05)

###to get matching columns between my dataset and databese

gene.name <- unique(signif$gene)
H.gene.symbol <- select(h_gene_sets, gs_name, gene_symbol)

#run enrichment
enrich.H<- enricher(gene = gene.name, TERM2GENE = H.gene.symbol)

##Extract results
class(enrich.H)

head(enrich.H@result)

class(enrich.H@result$GeneRatio)

##format the results to be able to plot
#separate ratios into 2 columns of data (to see how big the overlap between the gene set and halmark database)

enrich.H.df <- enrich.H@result %>% 
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  
  #Calculate k/K to see how many significant genes are in the gene set relative to haw big the geneset is
  mutate("k.K"=size.overlap.term/size.term) 



##visualize results##
enrich.H.df %>% 
  filter(pvalue <= 0.05) %>%
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Significant genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Mtb significant genes \nenriched in Hallmark gene sets")

################################https://github.com/hawn-lab/workshops_UW_Seattle/commit/b8fa2bf94eff5693d8ab2d79c8f3dc3423796f07
######GSEA########
#this method could not be used for more than two condition comparison, if your data has more than two groups or condition separate group of interest before

#load DEGenes data
DEG <- read.csv(file = "DESeqECmainAP.csv", sep=",")


#Get gene set database
H = msigdbr(species = "mouse", category = "H")
head(H)

class(H)

##or load your manual pathway database

#H <- read.csv(file = "H.gene.name2.csv" , header = T, sep=",")

#format gene set database
head(H)

H.gene.name.ls <- H %>% 
  select(gs_name, gene_symbol) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(gene_symbol))) %>% 
  deframe()


#GSEA does not need cutoff non-significant pvalue since padj significancy is not a proper treshold here 
##ranking the gene dataset
DEG <- sign(DEG$log2FoldChange)*(-log10(DEG$pvalue)) # we will use the signed p values from spatial DGE as ranking
names(DEG) <- DEG$gene # genes as names#
head(DEG)
DEG <- sort(DEG, decreasing = TRUE) # sort genes by ranking
plot(DEG)

#to check if we have infinite ranking because fgsea does not accept non-finit numbers
max(DEG)
min(DEG)

#table(signif$pvalue >= 0.05)

#Extract expression data
FC <- as.data.frame(DEG$log2FoldChange) %>% 
  as.data.frame(DEG$gene)%>% 
  #Move gene IDs from rownames to a column
  rownames_to_column("gene")

#format for gsea 
FC.vec <- FC$"signif$log2FoldChange"
names(FC.vec) <- FC$gene

#set score type 
#min(FC.vec)
#max(FC.vec)

#bc we have + & - LFC standardiz it 

#scoreType <- "std"

#Run GSEA#####
#the number of nperm=1000 can increase till 1milion depending on the pvalue that you have
#gsea.H <- fgseaSimple(pathways = H.gene.name.ls,
#                      stats = FC.vec,
#                      scoreType = scoreType,
#                      nperm=1000)





###Sahel fgsea##
gsea.H <- fgseaSimple(pathways = H.gene.name.ls,
                      stats = FC.vec,
                      minSize  = 1,
                      maxSize  = 1000,
                      nperm    = 1e7,
                      nproc    = ceiling( detectCores() / 2 ))


#### plot gsea results ####
class(gsea.H)

gsea.H %>% 
  filter(pval <= 0.05) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(pathway = gsub("HALLMARK_","", pathway),
         pathway = gsub("_"," ", pathway)) %>% 
  
  ggplot(aes(x=reorder(pathway, NES), #Reorder gene sets by NES values
             y=NES)) +
  geom_col() +
  theme_classic() +
  #Force equal max min
  lims(y=c(-3.2,3.2)) +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Normalized enrichment score (NES), pval <= 0.05",
       x="Gene set",
       title = "Hallmark GSEA \nDown in +Mtb <--         --> Up in +Mtb")

##output##
#to write a table
gsea.H <- apply(gsea.H,2,as.character)
write.table(gsea.H, file="GSEA_DESeqECmainAP.csv", sep=",", row.names = F)

# Save a single object to a file
#saveRDS(gsea.H, "GSEA_DESeqECmainAP.csv")


