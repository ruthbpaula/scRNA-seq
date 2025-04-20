## scRNA gene lookup function - Adapted from Tim's code
## 2/28/23
## per genotype, lookup expression of individual gene: UMAP, neuron v glia violin, ranked boxplots, DEG list

setwd('/home/humble_local_25t/timothyw/scRNA/multiModel/analysis')
ruth_wd <- "/home/humble_local2_25t/ruth/data/scRNA-seq/new_Justin_res2_120PCs_labeled_data/re_generate/"

library(dplyr)
library(ggplot2)
library(splines)
library(DESeq2)
library(gridExtra)
library(Seurat)
library(matrixStats)
library(SeuratData)
library(SeuratDisk)
library(gprofiler2)
library(reshape2)
library(ComplexHeatmap)
library(seriation)
library(circlize)
library(cowplot)
library(openxlsx)

## required objects ----
# Seurat Obj (should consider re-writing fcn so don't have to load this entire rds)
mmFinal <- readRDS('mm_d10_SeurObj_optimized_res2_120pcs_260k.rds')
colnames(mmFinal@meta.data)
mmFinal$antonly <- NA
mmFinal$antonly[is.na(as.numeric(mmFinal$final_annotation))] <- mmFinal$final_annotation[is.na(as.numeric(mmFinal$final_annotation))]
Idents(mmFinal) <- 'antonly'
#write.table(mmFinal@reductions$umap@cell.embeddings, 'MM_scRNA_UMAP_embeddings.txt', sep = '\t', col.names = T, row.names = T, quote = F)

# non-0 prop rankings per gene for 3 genotypes
#gt.res <- readRDS('non-zero-prop.rds')
#rankings <- as.list(c())
#for (i in 1:length(gt.res)){
#  df <- gt.res[[i]]
#  cellranks <- data.frame(gene = rownames(df), rank = NA)
#  for (j in 1:nrow(df)) {
#    cellranks[j,'rank'] <- paste(colnames(df)[order(df[j,], decreasing = T)],collapse=',')
#  }
#  rankings[[i]] <- cellranks
#  names(rankings)[i] <- names(gt.res)[i]
#}
#saveRDS(rankings, 'non-zero-prop_rankings.rds')
rankings <- readRDS('non-zero-prop_rankings.rds')

# MAST DEG results
mast.batch <- read.delim('MAST_MM_DEGs_aws_withbatch.txt', sep = '\t', header= T, stringsAsFactors = F)
mast.sig <- mast.batch[mast.batch$padj < 0.05,]

# Loom object
mmFinal.loom <- Connect(filename = 'MM_res2_PC120_final_SeurObj-to-loom.loom', mode = 'r')

#cell classifications
categories <- read.delim('broad_cell_type.csv', sep =',', header = T, stringsAsFactors = F)

## GOIs ----
#goi <- c('Csp','Gba1a','CG8027','Ids','CG6201','LManII','Npc1a','Gba1b','beta-Man','CG6984')
# goi <- c('Ect3','Snmp1','CG15533','Lip4','CG10104','MFS3')#'Dsb'  # sent by Jinghan to Ruth on 05/15/2023
# goi_1 <- unique(read.csv(paste0(ruth_wd, "Jinghan_remaining_LSD_homologs_to_analyze_05232023_short.txt"), sep="\t")[,2])  # fly symbol
# goi_1 <- unique(read.csv(paste0(ruth_wd, "Jinghan_remaining_LSD_homologs_to_analyze_05232023_less_conserved_homologs.txt"), sep="\t")[,2])  # fly symbol
# goi_1 <- goi_1[!grepl("^CG", goi_1)]
# # goi_2 <- unique(read.csv(paste0(ruth_wd, "Jinghan_remaining_LSD_homologs_to_analyze_05232023_short.txt"), sep="\t")[,3])  # CG number
# goi_2 <- unique(read.csv(paste0(ruth_wd, "Jinghan_remaining_LSD_homologs_to_analyze_05232023_less_conserved_homologs.txt"), sep="\t")[,3])  # CG number
# table_1 <- cbind(as.data.frame(goi_1), as.data.frame(goi_1 %in% rownames(mmFinal)))
# goi_1 <- table_1[table_1[,2] %in% "TRUE",][,1]
# table_2 <- cbind(as.data.frame(goi_2), as.data.frame(goi_2 %in% rownames(mmFinal)))
# goi_2 <- table_2[table_2[,2] %in% "TRUE",][,1]
# goi <- c(goi_1, goi_2)
goi <- read.table("Jinghan_gene_graphs_to_re-generate.txt")[,1]
goi %in% rownames(mmFinal)  # all should be TRUE

# genelookup function ----
# make sure that all genes are present in dataset, check ID version / annotation symbol

topClusters <- function(gene, genotype) {
  gt <- rankings[[genotype]]
  pos <- which(gt$gene %in% gene)
  top <- sapply(1:length(pos), function(y) strsplit(gt[pos[y], 'rank'],',')[[1]][1:20]) 
  colnames(top) <- gt$gene[pos]
  return(top)
}

## Optional: select the genes here
gene <- 'Fuca'
genotype <- 'aSyn'

setwd('/home/humble_local2_25t/ruth/data/scRNA-seq/new_Justin_res2_120PCs_labeled_data/re_generate/')

lookup <- function(gene, genotype) {
  # paired UMAP gene expression
  p <- FeaturePlot(mmFinal, features = gene, 
                   cells = rownames(mmFinal@meta.data)[mmFinal$genotype%in%c('elav-gal4 F', genotype)], 
                   reduction = 'umap', 
                   split.by = 'genotype', 
                   cols = c('grey', 'red'), label = T, repel = T, raster = T, order = T, raster.dpi = c(1500,1500), slot = 'data')
  
  # Boxplot of top expression (based on prop non-0)
  gene.list <- topClusters(gene, genotype)
  gene <- colnames(gene.list)
  cell <- gene.list[,1]
  cell.positions <- mmFinal.loom[['col_attrs/final_annotation']][] %in% cell & mmFinal.loom[['col_attrs/genotype']][] %in% c(genotype, 'elav-gal4 F')
  gene.pos <- mmFinal.loom[['row_attrs/Gene']][] == gene
  counts <- as.data.frame(mmFinal.loom[['matrix']][cell.positions, gene.pos])
  counts$genotype <- mmFinal.loom[['col_attrs/genotype']][cell.positions]
  counts$cell <- mmFinal.loom[['col_attrs/final_annotation']][cell.positions]
  rownames(counts) <- mmFinal.loom[['col_attrs/CellID']][cell.positions]
  colnames(counts)[1] <- 'Counts'
  
  counts$cell <- factor(counts$cell, levels = c(aggregate(counts$Counts, by = list(counts$cell), FUN = median)[,1][order(aggregate(counts$Counts, by = list(counts$cell), FUN = median)[,2], decreasing = T)]))
  counts.non0 <- counts[counts$Counts>0,]
  counts.non0$genotype <- factor(counts.non0$genotype, levels = c('elav-gal4 F', genotype))
  
  p2 <- ggplot(counts.non0, aes(x=cell, y=Counts, fill=genotype)) + 
    geom_boxplot(outlier.size = 0)+
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.04, vjust = 1.03)) +
    labs(x=NULL, title = paste0(gene, ' top 20 clusters non-0 prop')) +
    scale_fill_manual(breaks = c('elav-gal4 F', genotype), values = c('grey70','red'))
  
  # Neuron vs glia mean expression (all cells)
  cell.all <- mmFinal.loom[['col_attrs/genotype']][] %in% c(genotype, 'elav-gal4 F')
  counts.all <- as.data.frame(mmFinal.loom[['matrix']][cell.all, gene.pos])
  counts.all$genotype <- mmFinal.loom[['col_attrs/genotype']][cell.all]
  counts.all$cell <- mmFinal.loom[['col_attrs/final_annotation']][cell.all]
  rownames(counts.all) <- mmFinal.loom[['col_attrs/CellID']][cell.all]
  colnames(counts.all)[1] <- 'Counts'
  counts.all <- left_join(counts.all, categories, by = c('cell' = 'Cluster'))
  counts.all$class <- 'Neuron'
  counts.all$class[counts.all$Type == 'Glia'] <- 'Glia'
  counts.all.non0 <- counts.all[counts.all$Counts>0,] ## dataframe of cells with non-0 gene expression of gene of interest
  counts.all.non0$condition <- paste(counts.all.non0$genotype, counts.all.non0$class, sep = '_')
  # Plotting  Violin plot ----
  p3 <- ggplot(counts.all.non0, aes(x=class, y=Counts, group = condition, fill = genotype)) + 
    geom_violin() +
    geom_boxplot(fill='white', width =0.3 , outlier.size = 0, outlier.color = NA, position = position_dodge(0.9)) +
    geom_point(aes(fill=genotype), position=position_jitterdodge(dodge.width=0.9, jitter.width = 0.1), size = 0.3, alpha = 0.6) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1.04, vjust = 1.03)) +
    labs(x=NULL, title = paste0(gene, ' ,neuron v glia, 0-count cells removed, \n%non-0 of total cells: ', 
                                round(100*(nrow(counts.all.non0)/nrow(counts.all)), 3),'%')) +
    scale_fill_manual(breaks = c('elav-gal4 F', genotype), values = c('grey70','red'))
  
  # DEG plot
  if (gene %in% mast.batch$gene[mast.batch$genotype == genotype]) {
    all.mast <- mast.batch[mast.batch$genotype == genotype & mast.batch$gene == gene,]
    lfc <- all.mast[,c('avg_log2FC','cluster')]
    pval <- all.mast[,c('padj','cluster')]
    pval$padj <- -log10(pval$padj)
    pval$padj[pval$padj < -log10(0.05)] <- NA
    #
    ord <- order(lfc$avg_log2FC, decreasing = T)
    lfc <- lfc[ord,]
    pval <- pval[ord,]
    rownames(lfc) <- lfc$cluster
    lfc$cluster <- NULL
    lfc <- as.matrix(lfc)
    rownames(pval) <- pval$cluster
    pval$cluster <- NULL
    pval <- as.matrix(pval)
    
    # circle heatmap 1
    # deg num
    col_fun1 = colorRamp2(c(0,1), c("white", "black"))
    col_fun2 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    circos.clear()
    c1 <- function() {
      circos.par(gap.after = 20 , start.degree=(120-20))
      circos.heatmap(mat = pval,
                     col = col_fun1, 
                     rownames.side = "outside" , cell.border = 'grey70', cluster = F , rownames.cex = 1)
      
      circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        if(CELL_META$sector.numeric.index == 1) { # the last sector
          cn = colnames(pval)
          n = length(cn)
          circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                      n - 1:n + 0.5, cn, 
                      cex = 1.2, adj = c(0, 0.5), facing = "inside")
        }
      }, bg.border = NA)
      
      # LFC
      circos.heatmap(mat = lfc,
                     col = col_fun2, 
                     cell.border =  'grey70')
      circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        if(CELL_META$sector.numeric.index == 1) { # the last sector
          cn = colnames(lfc)
          n = length(cn)
          circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                      n - 1:n + 0.5, cn, 
                      cex = 1.2, adj = c(0, 0.5), facing = "inside")
        }
      }, bg.border = NA)
      title(paste0(genotype, ' ', gene, ' differential expression'))
      lgd1 = Legend(title = "-log10(padj)", col_fun = col_fun1)
      lgd2 = Legend(title = 'log2FC',col_fun = col_fun2 )
      draw(lgd1, x = unit(0.7, "in"), y= unit(1, "in") )
      draw(lgd2, x = unit(1.4, "in"), y = unit(1, "in"))
    }
    circos.clear()
    pdf(paste0(ruth_wd, genotype,'_',gene,'_scRNA_expression.pdf'), width = 14, height = 22)
    # Move to a new page
    grid.newpage()
    # Create layout : nrow = 3, ncol = 2
    pushViewport(viewport(layout = grid.layout(ncol = 2,nrow = 4)))
    # A helper function to define a region on the layout
    define_region <- function(row, col){
      viewport(layout.pos.row = row, layout.pos.col = col)
    } 
    # Arrange the plots
    print(p[[2]], vp = define_region(1, 1))
    print(p[[1]], vp = define_region(1, 2))
    print(p3, vp = define_region(2, 1))
    print(p2, vp=define_region(3, 1:2))
    print(c1(), vp = define_region(4,1:2))
    dev.off()
  } else {
    pdf(paste0(ruth_wd, genotype,'_',gene,'_scRNA_expression.pdf'), width = 14, height = 20)
    # Move to a new page
    grid.newpage()
    # Create layout : nrow = 3, ncol = 2
    pushViewport(viewport(layout = grid.layout(ncol = 2,nrow = 3)))
    # A helper function to define a region on the layout
    define_region <- function(row, col){
      viewport(layout.pos.row = row, layout.pos.col = col)
    } 
    # Arrange the plots
    print(p[[2]], vp = define_region(1, 1))
    print(p[[1]], vp = define_region(1, 2))
    print(p3, vp = define_region(2, 1))
    print(p2, vp=define_region(3, 1:2))
    dev.off()
  }
}

###################################################################################################################################
# Optional - select specific genes
goi <- 'Fuca'

# loop through goi for aSyn
for (i in 1:length(goi)) {
  lookup(gene = goi[i], genotype = 'aSyn')
}

# OPTIONAL - loop through genotypes for 1 gene
levels(factor(mmFinal$genotype))
gts <- c('AB42', 'TauWT', 'aSyn')
for (i in 1:length(gts)) {
  lookup(gene = goi[10], genotype = gts[i])
}
