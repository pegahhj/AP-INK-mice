# AP-INK-mice
#AP&amp; INK mice/ snRNA-seq
#the project has been done on Compute Canada/It is the integration of two datasets regarding skin microvascular of two types of mice
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
   

