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

toc = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/ATX71AV_MPS12345537_E07_9852/outs/filtered_feature_bc_matrix')) 
tod = Seurat::Read10X(file.path('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/ATX71AV_MPS12345537_E07_9852/outs/raw_feature_bc_matrix')) 
sc = SoupChannel(tod, toc) 
cluster = read.table('/lustre03/project/6003727/pegahhj/skinInkattacApInteg/cellrangerOut/ATX71AV_MPS12345537_E07_9852/outs/analysis/clustering/gene_expression_graphclust/clusters.csv', sep=',', header=T) 
sc = setClusters(sc,cluster$Cluster) 
sc = autoEstCont(sc) 
ATX71AV_MPS12345537_E07_9852 = adjustCounts(sc)

#creating the object
APskin_samples <- CreateSeuratObject(counts = cbind(WT_SED10_MPS12345537_E08_9852 , WT_AP2_MPS12345537_C07_9852, WT_AP10_MPS12345537_D08_9852, WT_AP3_MPS12345537_A07_9852, WT_AP5_MPS12345537_F07_9852, WT_AP6_MPS12345537_H07_9852, WT_AP9_MPS12345537_B08_9852, WT_SED2_MPS12345537_D07_9852, WT_SED3_MPS12345537_B07_9852, WT_SED5_MPS12345537_G07_9852, WT_SED6_MPS12345537_A08_9852, WT_SED9_MPS12345537_C08_9852, ATX71AV_MPS12345537_E07_9852), project = "APmouseskin", min.cells = 5)

#specifying different condition
APskin_samples@meta.data$sample <- c(rep("SED10", ncol(WT_SED10_MPS12345537_E08_9852)), rep("AP2", ncol(WT_AP2_MPS12345537_C07_9852)), rep("AP10", ncol(WT_AP10_MPS12345537_D08_9852)), rep("AP3", ncol(WT_AP3_MPS12345537_A07_9852)), rep("AP5", ncol(WT_AP5_MPS12345537_F07_9852)), rep("AP6", ncol(WT_AP6_MPS12345537_H07_9852)), rep("AP9", ncol(WT_AP9_MPS12345537_B08_9852)), rep("SED2", ncol(WT_SED2_MPS12345537_D07_9852)), rep("SED3", ncol(WT_SED3_MPS12345537_B07_9852)), rep("SED5", ncol(WT_SED5_MPS12345537_G07_9852)), rep("SED6", ncol(WT_SED6_MPS12345537_A08_9852)), rep("SED9", ncol(WT_SED9_MPS12345537_C08_9852)), rep("ATX", ncol(ATX71AV_MPS12345537_E07_9852)))

APskin_samples@meta.data$condition <- c(rep("SED", ncol(WT_SED10_MPS12345537_E08_9852)), rep("PhA", ncol(WT_AP2_MPS12345537_C07_9852)), rep("PhA", ncol(WT_AP10_MPS12345537_D08_9852)), rep("PhA", ncol(WT_AP3_MPS12345537_A07_9852)), rep("PhA", ncol(WT_AP5_MPS12345537_F07_9852)), rep("PhA", ncol(WT_AP6_MPS12345537_H07_9852)), rep("PhA", ncol(WT_AP9_MPS12345537_B08_9852)), rep("SED", ncol(WT_SED2_MPS12345537_D07_9852)), rep("SED", ncol(WT_SED3_MPS12345537_B07_9852)), rep("SED", ncol(WT_SED5_MPS12345537_G07_9852)), rep("SED", ncol(WT_SED6_MPS12345537_A08_9852)), rep("SED", ncol(WT_SED9_MPS12345537_C08_9852)), rep("ATX", ncol(ATX71AV_MPS12345537_E07_9852)))

#add a column to the metadata
APskin_samples[["percent.mt"]] <- PercentageFeatureSet(APskin_samples, pattern = "^mt-") 
APskin_samples <- PercentageFeatureSet(APskin_samples, pattern = "^mt-", col.name = "percent.mt")

#save the object
saveRDS(APskin_samples, file = "APskin_samples.rds")
