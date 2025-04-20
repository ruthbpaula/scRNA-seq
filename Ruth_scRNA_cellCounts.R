### Cell count analysis from scRNA-seq library

setwd("/home/humble_local2_25t/ruth/data/scRNA-seq/res2_120PCs_labeled_data")

library(DESeq2)
library(Seurat)
library(reshape2)
library(dplyr)

# import data
#data.integrated <- readRDS("testing_new_cluster/data_integrated_res2_TauWT_Elav_F.rds")  ## old clustering data
#data.integrated <- readRDS("seuratObj_aSyn_elavCtrl_F_label.rds")
#Idents(object = data.integrated) <- "orig.ident"
#sub_aSyn <- subset(x = data.integrated, idents = c("26392", "26399", "26408", "26395", "26404", "26413")) # c("26392", "26399", "26408", "26394", "26403", "26412"))

sub_aSyn <- readRDS("seuratObj_aSyn_elavCtrl_F_label.rds")
sub_aSyn@meta.data$genotype <- factor(sub_aSyn$orig.ident, levels = c("26392", "26394", "26399", "26403", "26408", "26412"), labels = c("elav-gal4 F","aSyn","elav-gal4 F","aSyn","elav-gal4 F","aSyn"))

sub_aSyn@meta.data$genotype <- gsub("elav-gal4 F",'elav', sub_aSyn@meta.data$genotype)

# stratify by condition / or library - libraries are indicated as factors in orig.ident column
conditions <- levels(as.factor(sub_aSyn$orig.ident))
count.list <- as.list(c())

for (i in 1:length(conditions)) {
  ph <- sub_aSyn@meta.data[sub_aSyn$orig.ident == conditions[i],]
  count.list[[i]] <- count_(ph, 'final_annotation')  # count conditions in column that contains cell annotations
  count.list[[i]]$freq <- 100*(count.list[[i]]$n/sum(count.list[[i]]$n)) #relative composition within 
  count.list[[i]]$condition <- conditions[i]
}

count.df <- do.call(rbind, count.list)
count.df$genotype <- as.factor(sapply(1:nrow(count.df), function(x) strsplit(count.df$condition, '_')[[x]][1]))
#count.df$age <- as.numeric(gsub('d','',sapply(1:nrow(count.df), function(x) strsplit(count.df$condition, '_')[[x]][2])))
levels(as.factor(count.df$final_annotation))
#barplot, individual conditions, grouped by genotype

# cell counts over time - line plots
count.df$final_annotation <- as.factor(count.df$final_annotation)
count.df$cell <- paste0(count.df$final_annotation, '_', count.df$genotype)

saveRDS(count.df, file = "count.df_res2_aSyn_label.rds")

# # # log2FC, standard error (timepoints as replicates) # # #


## DESEQ2 approach in analyzing count data
# normalization of count data (use deseq2 median of ratios, then run nb.binomial glm)
# median of ratios: account for total # cells captured per condition, and variable cell composition per condition

#count.df <- readRDS("count.df.rds")

count.df.rawCounts <- dcast(data = count.df, formula = final_annotation ~ condition,  value.var = 'n')
hist(count.df.rawCounts[,3], breaks = 100)

rownames(count.df.rawCounts) <- count.df.rawCounts$final_annotation
count.df.rawCounts$final_annotation <- NULL
count.df.rawCounts <- as.matrix(count.df.rawCounts)
count.df.rawCounts[is.na(count.df.rawCounts)] <- 0
as.integer(count.df.rawCounts)
count.df.rawCounts <- count.df.rawCounts[, c("26392", "26399", "26408", "26394", "26403", "26412")] # controls first, mutants after

info <- data.frame(sample = colnames(count.df.rawCounts), genotype = factor(c(rep('elav',3), rep('aSyn',3)), levels = c('elav','aSyn')))  # age = as.factor(c(rep(c(1,10,20),2)))

# reordering columns - see if results improve   -> did not change
#count.df.rawCounts <- count.df.rawCounts[, c("26392", "26399", "26408", "26394", "26403", "26412")]
#info <- data.frame(sample = colnames(count.df.rawCounts), genotype = factor(c(rep('elav',3), rep('aSyn',3)), levels = c('elav','aSyn')))  # age = as.factor(c(rep(c(1,10,20),2)))

# adding batch column into info   -> improved slightly
info$batch <- as.factor(c(1,2,3,1,2,3))

# creating DESeq object and running DESeq
dds <- DESeqDataSetFromMatrix(countData = count.df.rawCounts, colData = info, design = ~ genotype + batch)  # add age as a covariate as well when you have different age time points
#keep <- rowSums(counts(dds)) >= 10  # this sampling is only used for bulk RNA-seq data
#dds <- dds[keep,]
dds <- DESeq(dds)

# median-to-ratio normalized cell counts 
cell_norm_count <- DESeq2::counts(dds, normalized = T)
hist(cell_norm_count[,1]) #negative binomial dist 
write.table(cell_norm_count, file = 'integrated_cellCounts_DESeq2_norm_aSyn_Elav_F_res2_label.txt', sep = '\t', col.names = T, row.names = T, quote = F)

#Wald test for significance 
library(ashr)

res_geno <- results(object = dds, contrast = c('genotype','aSyn','elav'))
res_geno_shrink <- as.data.frame(lfcShrink(dds, res=res_geno, type = 'ashr', coef = 'genotype_aSyn_vs_elav'))
summary(res_geno)
write.table(res_geno_shrink, file = 'cell_count_DESeq2_lfcshrink_integrated_snn_res.2_res2_aSyn_label.txt', sep = '\t' , col.names = T, row.names = T, quote = F)
res_geno_shrink$cluster <- rownames(res_geno_shrink)
#res_geno_shrink$CImin <- res_geno_shrink$log2FoldChange-1.96*res_geno_shrink$lfcSE
#res_geno_shrink$CImax <- res_geno_shrink$log2FoldChange+1.96*res_geno_shrink$lfcSE

library(ggplot2)
setwd("/home/humble_local2_25t/ruth/data/scRNA-seq/res2_120PCs_labeled_data")

res_geno_shrink <- read.table('cell_count_DESeq2_lfcshrink_integrated_snn_res.2_res2_aSyn_label.txt')
res_geno_shrink$cluster <- rownames(res_geno_shrink)

ggplot(res_geno_shrink, aes(y=reorder(cluster, log2FoldChange), x=log2FoldChange, fill = -log10(padj))) +
  geom_bar(stat='identity') +
  labs(y='clusters', x = 'log2FC(aSyn/elav) of cell frequency') +
  theme_minimal(base_size = 20) +
  coord_flip() +
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'grey70', midpoint = 1.3, na.value = 'grey', guide = 'colourbar') +
  theme(axis.text.x = element_text(angle = 90)) +
  xlim(-0.1,0.1)
#geom_errorbar(aes(ymin = CImin, ymax = CImax), color = 'black')

ggsave("DESeq2_cellCounts_aSyn_Elav_F_res2_horiz_rescaled_label.png", width = 40, height = 7)
