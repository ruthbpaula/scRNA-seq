## scRNA-seq processing of raw counts in Seurat
## import raw count matricies into Seurat
## QC, normalization
## Library integration via CCA to account for batch effects, mapping of similar cell types across libraries
## Doublet detection
## Path to fully processed Seurat object (aSyn and Control) with cell annotation metadata

#rm(list=ls(all=T))
setwd("/home/humble_local2_25t/ruth/data/scRNA-seq")
library(Seurat)
library(ggplot2)
library(gridExtra)

#library(future)
#plan(multiprocess, workers = 4)
#options(future.globals.maxSize = 8000*1024^2)
#options(future.globals.maxSize = 20000*1024^2)

## ------------ ##
# # import raw count matrix from cellranger output
# aSyn_b1 <- Read10X('/home/humble_local2_25t/ruth/data/26394_aSyn/filtered_feature_bc_matrix')
# aSyn_b2 <- Read10X('/home/humble_local2_25t/ruth/data/26403_aSyn/filtered_feature_bc_matrix')
# aSyn_b3 <- Read10X('/home/humble_local2_25t/ruth/data/26412_aSyn/filtered_feature_bc_matrix')
# elav_b1 <- Read10X('/home/humble_local2_25t/ruth/data/26392_elav_gal4_F/filtered_feature_bc_matrix')
# elav_b2 <- Read10X('/home/humble_local2_25t/ruth/data/26399_elav_gal4_F/filtered_feature_bc_matrix')
# elav_b3 <- Read10X('/home/humble_local2_25t/ruth/data/26408_elav_gal4_F/filtered_feature_bc_matrix')
# 
# ## ------------ ##
# #create seurat object
# aSyn_b1 <- CreateSeuratObject(counts = aSyn_b1, project = 'aSyn_b1')
# aSyn_b2 <- CreateSeuratObject(counts = aSyn_b2, project = 'aSyn_b2')
# aSyn_b3 <- CreateSeuratObject(counts = aSyn_b3, project = 'aSyn_b3')
# elav_b1 <- CreateSeuratObject(counts = elav_b1, project = 'elav_b1')
# elav_b2 <- CreateSeuratObject(counts = elav_b2, project = 'elav_b2')
# elav_b3 <- CreateSeuratObject(counts = elav_b3, project = 'elav_b3')
# 
# ## ------------ ##
# ### merging seurat objects into 1
# data <- merge(x = aSyn_b1, y=c(aSyn_b2, aSyn_b3, elav_b1, elav_b2, elav_b3), add.cell.ids= c('T1','T10','T20','E1','E10','E20'), project = 'aSynRW_timeSeries')
# 
# ## ------------ ##
# ## view number of cells captured per library
# cc <- as.data.frame(table(data@meta.data$orig.ident))
# colnames(cc) <- c('condition','count')
# cc$var <- NA
# cc$var[grep('b2', cc$condition)] <- 'b2'
# cc$var[grep('b3', cc$condition)] <- 'b3'
# cc$var[is.na(cc$var)] <- 'b1'
# cc$var <- as.factor(cc$var)
# cc$geno <- NA
# cc$geno[grep('elav', cc$condition)] <- 'elav'
# cc$geno[grep('aSyn', cc$condition)] <- 'aSyn'
# 
# p1 <- ggplot(cc, aes(x=var, y=count, group=geno, fill = geno)) +
#   geom_bar(position = position_dodge() , stat = 'identity', color = 'black') +
#   theme_minimal(base_size = 16) +
#   labs(x=NULL, y='raw cell counts') +
#   scale_fill_manual(name = 'genotype', breaks = c('elav','aSyn'), values = c('navy', 'orange'))

## ------------ ##
## quantify technical covariates, mitochondrial read percentage, QC filtering, visualization

## FOR POSITIVE CONTROL: run the same pipeline with TauWT instead of aSyn

data_raw <- readRDS("/home/humble_local_25t/timothyw/scRNA/multiModel/analysis/mergedSeurat_unprocessed_010822.rds")

# subsetting for aSyn and controls
Idents(object = data_raw) <- "orig.ident"
data <- subset(x = data_raw, idents = c("26392", "26399", "26408", "26395", "26404", "26413"))

data <- PercentageFeatureSet(data, pattern = '^mt:', col.name = 'percent.mt')
data$log_ncount <- log10(data$nCount_RNA)
qc1 <- VlnPlot(data, features = c("nFeature_RNA", "log_ncount", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
qc2 <- CombinePlots(plots = list(plot1, plot2))

pdf(file = 'mtDNA_hist_TauWT_Elav_F.pdf', height = 6, width = 6)
hist(data@meta.data$percent.mt, breaks = 100)
dev.off()
pdf(file='nCount_hist_TauWT_Elav_F.pdf', height =6, width = 6)
hist(data@meta.data$nCount_RNA, breaks = 1000)
dev.off()
pdf(file='nFeature_hist_TauWT_Elav_F.pdf', height = 6, width = 6)
hist(data@meta.data$nFeature_RNA, breaks = 100)
dev.off()

##
qc1a <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) + geom_hline(yintercept = 200, col = 'red', linetype = 'dashed') + geom_hline(yintercept = 3000, col = 'red', linetype = 'dashed') + NoLegend() + labs(x=NULL)
qc1b <- VlnPlot(data, features = "log_ncount",  pt.size = 0) + NoLegend()+ labs(x=NULL)
qc1c <- VlnPlot(data, features = "percent.mt",  pt.size = 0) + geom_hline(yintercept = 20, col = 'red', linetype = 'dashed') + NoLegend()+ labs(x=NULL)

saveRDS(qc1a, "qc1a_TauWT.rds")
saveRDS(qc1b, "qc1b_TauWT.rds")
saveRDS(qc1c, "qc1c_TauWT.rds")

## ------------ ##
## batch correction or data integration across libraries by Seurat CCA
## note: memory intensive workflow
#split object by condition
data.list <- SplitObject(data, split.by = 'orig.ident')

#sctransform regression based normalization per library, QC cutoff, mt dna regression
for (i in 1:length(data.list)) {
  data.list[[i]] <- subset(data.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
  data.list[[i]] <- SCTransform(data.list[[i]], variable.features.n = 5000, vars.to.regress = 'percent.mt')
}

#prep integration features per condition, top 5000 features
data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 5000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features)


# integrate cells across conditions using selected integration anchors
data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = data.features)

saveRDS(data.list, "data.list.rds")
rm(data.list)
gc()
saveRDS(data.features, "data.features.rds")
rm(data.features)
gc()

saveRDS(data.anchors, "data.anchors.rds")

data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT")


##### Log and size factor normalization of raw counts for gene expression analyses - non-integrated counts ######
# Save in 'RNA' assay
DefaultAssay(data.integrated) <- 'RNA'
data.integrated <- NormalizeData(data.integrated, normalization.method = "LogNormalize", scale.factor = 10000, assay = 'RNA')

saveRDS(data.integrated, file = "data.integrated_TauWT_processed.rds")

## ------------ ##
## PCA and unsupervised clustering of integrated dataset
## CCA 5k
## --- pca/umap/tsne --- ##

setwd("/home/humble_local2_25t/ruth/data/scRNA-seq")
library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)

#sub_aSyn <- readRDS("seuratObj_aSyn_elavCtrl_F.rds")

DefaultAssay(data.integrated) <- 'integrated'
data.integrated <- RunPCA(data.integrated, npcs = 100, verbose = T) 
data.integrated <- RunUMAP(data.integrated, dims = 1:100, n_threads = 8, n_sgd_threads = 'auto')
data.integrated <- RunTSNE(data.integrated, dims = 1:100, num_threads = 8, perplexity = 30, verbose = T)

## --- clustering --- ##
data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:100)
data.integrated <- FindClusters(data.integrated, resolution = 2)   # if not good, change PCA parameters
gc()

## visualization
p1 <- DimPlot(data.integrated, reduction = "umap", group.by = 'orig.ident', label = F, repel = T) 
p1.0 <- DimPlot(data.integrated, reduction = "umap", split.by = 'orig.ident', group.by = 'seurat_clusters', label = F, repel = F, ncol = 3) + NoLegend()
p1.1 <- DimPlot(data.integrated, reduction = "umap", group.by = 'seurat_clusters', label = T) + NoLegend()
ggsave(p1, filename = 'umap_cca_conditions_resolution2.pdf', width = 8, height = 6, device = 'pdf')
ggsave(p1.0, filename = 'umap_cca_splitConditions_resolution2.pdf', width = 22, height = 12, device = 'pdf')
ggsave(p1.1, filename = 'umap_cca_clusters_resolution2.pdf', device = 'pdf', width = 8, height = 6)

## Getting the number of clusters and cells
levels(data.integrated@meta.data$seurat_clusters) # 144 clusters for resolution = 2; 114 clusters for resolution = 0.8; 151 clusters for TauWT res2
nrow(as.data.frame(data.integrated@meta.data$seurat_clusters)) # 64943 rows or cells; 77595 for TauWT

saveRDS(data.integrated, "data_integrated_res2_TauWT_Elav_F.rds")

## ------------ ##
## quantify doublets in scRNA data using doubletfinder
## see tutorial: https://github.com/chris-mcginnis-ucsf/DoubletFinder
setwd("/home/humble_local2_25t/ruth/data/scRNA-seq")
library(Seurat)
library(DoubletFinder) # Installation: remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

data.integrated <- readRDS("data_integrated_res2_aSyn_Elav_F.rds")

# doublet rate based on linear increase - values provided by 10x Genomics manual; based on # cells recovered
doublet_rate <- function(x) { 
  y = 5.272e-04 + 7.589e-06*x
  return(y)
}  

#SCT as default for doublet estimation
DefaultAssay(data.integrated) <- 'SCT' 
byGT <- SplitObject(data.integrated, split.by = 'orig.ident')
conditions <- names(byGT)
meta.list <- as.list(c())
# DF for each condition / GEM
for (i in 1:length(conditions)) {
  #est pKa
  sweep.res.list_data <- paramSweep_v3(byGT[[i]], PCs = 1:40, sct = T) # , num.cores = 20
  sweep.stats_data <- summarizeSweep(sweep.res.list_data, GT = FALSE)
  bcmvn_data <- find.pK(sweep.stats_data)
  #pick pKa
  a <- as.numeric(as.character(bcmvn_data$pK[which(bcmvn_data$BCmetric == max(bcmvn_data$BCmetric))]))
  if (a > 0.12) {
    bcmvn_data <- find.pK(sweep.stats_data)[1:12,]
    a <- as.numeric(as.character(bcmvn_data$pK[which(bcmvn_data$BCmetric == max(bcmvn_data$BCmetric))]))
  }
  # estimate *homotypic* doublets. - using manual annotations will reduce over-estimation of doublets (nExp_poi.adj) - can't distingush homotypics
  # tentatively use assigned seurat clusters as cell 'annotation', assuming each cluster == a unique cell identity
  homotypic.prop <- modelHomotypic(byGT[[i]]@meta.data$seurat_clusters)   
  
  # estimate *total* doublet formation, based on cell loading density, based on 10x genomics metrics for loading density rate (linear function above)
  nExp_poi <- round(doublet_rate(ncol(byGT[[i]]))*nrow(byGT[[i]]@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder w/o adjustment
  ph1 <- doubletFinder_v3(byGT[[i]], PCs = 1:40, pN = 0.25, pK = a, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  #with adjustment for homotypic estimation (seurat cluster annotations)
  ph2 <- doubletFinder_v3(ph1, PCs = 1:40, pN = 0.25, pK = a, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_",a,'_',nExp_poi), sct = T)
  df <- ph2@meta.data
  df$DF_param <- paste0('pN_0.25_pK_',a,'_nExp_poi_',nExp_poi,'_nExp_poi_adj_', nExp_poi.adj)
  colnames(df)[9:11] <- c('pANN','DF_class','DF_class_adj')
  meta.list[[i]] <- df
  #names(meta.list)[i] <- conditions[i]
}

saveRDS(meta.list, "meta.list.rds")

setwd("/home/humble_local2_25t/ruth/data/scRNA-seq")
library(Seurat)

meta.list <- readRDS("meta.list.rds")
data.integrated <- readRDS("data_integrated_res2_aSyn_Elav_F.rds")
DefaultAssay(data.integrated) <- 'SCT'

library(data.table)
#df.final <- do.call(rbind, meta.list)
df.final <- as.data.frame(rbindlist(meta.list, fill=FALSE))

## -- ##
library(future.apply)
#plan(multiprocess, workers = 10)
#rownames(df.final) <- future_sapply(1:nrow(df.final), function(x) strsplit(rownames(df.final), '\\.')[[x]][2])
rownames(df.final) <- rownames(data.integrated@meta.data)
#merged_df.final <- merge(df.final, as.data.frame(data.integrated@meta.data), by = "row.names")   # careful if using this line -> may mess up with cell labels

#saveRDS(df.final, "df.final_new_label.rds")
#saveRDS(data.integrated@meta.data, "meta.data.rds")
#saveRDS(merged_df.final, "merged_df.final2.rds")
## -- ##

## ------------ ##
#save DF results to RDS object
#data.integrated@meta.data <- merged_df.final
data.integrated@meta.data <- df.final
sub_aSyn <- data.integrated

#sub_aSyn <- RenameCells(sub_aSyn, new.names = rownames(sub_aSyn@meta.data))   # careful if using this line -> may mess up with cell labels
Idents(object = sub_aSyn) <- "DF_class_adj"
sub_aSyn <- subset(x = sub_aSyn, idents = "Singlet")
#sub_aSyn <- SplitObject(sub_aSyn, split.by = 'DF_class_adj')   # alternative way to define idents
#sub_aSyn <- sub_aSyn[['Singlet']]   # alternative way to subset Seurat object

library(dplyr)
nrow(as.data.frame(sub_aSyn@meta.data$DF_class_adj))   # counting the number of singlets left - 59260 singlets -> number of doublets: 5863

saveRDS(sub_aSyn, 'sub_aSyn_res2.rds')
rm(data.integrated)
rm(byGT)
gc()

## the above code is for practice, demonstrating how to process raw scRNA gene counts after alignment using 10x Genomic's CellRanger
## ------------ ##
## ------------ ##
# for downstream analysis, e.g., DGE, cell abundance, deconvolution, use integrated dataset here: ~/home/humble_local_25t/timothyw/scRNA/multiModel/analysis/multiModel_Seurat_rPCA_refControl_Integrated_alignmentV1_010822.rds
# this integrated dataset contains multiple genotypes, make sure to subset samples that are alpha-synuclean and elav-Gal4 F control. 











