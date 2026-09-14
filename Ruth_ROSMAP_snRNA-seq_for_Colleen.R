setwd("/mnt/belinda_local/ruth/data/snRNA-seq/ROSMAP")

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# For donor-aware statistical testing
library(lme4)
library(lmerTest)
library(emmeans)


### ROSMAP snRNA-seq graphs for Colleen
# I am trying to check whether genes within her list of AD genetic interactions screen are 
# 1) highly expressed, or 
# 2) enriched in expression compared to other cell types. 
# Also, I would like to know which cell types (top 1 or 2) each of these genes is most highly expressed in.

###
# 1. Defining the approach ----
###

# Checking what is inside the .rds objects
x <- readRDS("Astrocytes.rds")

class(x)
dim(x)

# if Seurat:
Assays(x)
DefaultAssay(x)
x

slotNames(x[["RNA"]])

head(x@meta.data)
colnames(x@meta.data)


# First, we need to define high expression
# Highly expressed: gene is in the top 10% of genes by expression within that cell type and is detected in ≥25% of nuclei of that cell type.
# Very highly expressed: top 5% and detected in ≥50% of nuclei.

# For enrichment, I would use donor-level pseudobulk, because your projid column gives us the biological replicate. This avoids treating 150,000 nuclei as 150,000 independent samples.

# Input:
#   Multiple cell-type-specific Seurat .rds files
#
# Required metadata:
#   projid
#   cell_type_high_resolution
#
# Outputs:
#   - absolute expression metrics
#   - % cells expressing each gene
#   - expression percentile within each cell type
#   - high-expression classification
#   - donor-level pseudobulk expression
#   - cell-type enrichment
#   - mixed-model enrichment statistics
#   - top 1 / top 2 cell types per gene
#   - plots


###
# 2. User settings ----
###

# Folder containing all downloaded RDS files
rds_dir <- "/mnt/belinda_local/ruth/data/snRNA-seq/ROSMAP"

# Output folder
out_dir <- "ROSMAP_gene_expression_analysis"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Files
rds_files <- list.files(
  rds_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

# Your gene list
genes <- c(
  "SLC12A9",
  "ECHDC3",
  "HS3ST5",
  "PSEN1",
  "PSEN2",
  "MADD",
  "DOC2A",
  "SNX32",
  "ADAM17",
  "WWOX",
  "MAF",
  "ATP8B3",
  "BIN1",
  "JAZF1",
  "CTSH",
  "NCK2",
  "PTK2B",
  "PICALM",
  "CCDC47",
  "STX6",
  "TSPAN14",
  "WDR81",
  "RIN3",
  "APP",
  "ICA1",
  "SNX1",
  "ANK3",
  "RTF2"
)

genes <- unique(genes)

# Metadata columns
donor_col <- "projid"
celltype_col <- "cell_type_high_resolution"


# Expression thresholds

# "High expression":
# top 10% of genes AND detected in >=25% of cells
high_percentile <- 90
high_detection <- 25

# "Very high expression":
# top 5% AND detected in >=50% of cells
very_high_percentile <- 95
very_high_detection <- 50


# Pseudobulk QC

# Require at least this many nuclei from a donor/cell-type
# before using that donor in enrichment analyses
min_cells_per_pseudobulk <- 20

# Require at least this many donors for a cell type
min_donors_per_celltype <- 10


###
# 3. Containers ----
###

pb_counts_list <- list()
pb_meta_list <- list()

cell_metrics_list <- list()

reference_features <- NULL


###
# 4. Process each RDS object ----
###

# OPTIONAL - First, inspect file splits (such as excitatory neurons set1, set2, set3...) (takes some time to run)
# lapply(rds_files, function(f) {
#   x <- readRDS(f)
#   
#   data.frame(
#     file = basename(f),
#     n_cells = ncol(x),
#     n_donors = length(unique(x$projid)),
#     n_celltypes = length(unique(x$cell_type_high_resolution)),
#     celltypes = paste(
#       unique(x$cell_type_high_resolution),
#       collapse = " | "
#     )
#   )
# }) %>%
#   bind_rows()


# Process the objects
for (f in rds_files) {
  
  cat("\nProcessing:", basename(f), "\n")
  
  obj <- readRDS(f)
  
  # Raw counts
  counts <- GetAssayData(
    obj,
    assay = "RNA",
    slot = "counts"
  )
  
  # Verify that all objects use the same feature set/order
  if (is.null(reference_features)) {
    
    reference_features <- rownames(counts)
    
  } else {
    
    if (!identical(reference_features, rownames(counts))) {
      
      stop(
        paste(
          "Feature names/order differ in:",
          basename(f)
        )
      )
    }
  }
  
  
  # Metadata
  
  meta <- obj@meta.data
  
  # Make absolutely sure metadata and matrix columns match
  meta <- meta[colnames(counts), , drop = FALSE]
  
  keep <- (
    !is.na(meta[[donor_col]]) &
      !is.na(meta[[celltype_col]]) &
      meta[[celltype_col]] != ""
  )
  
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  
  donor <- as.character(meta[[donor_col]])
  celltype <- factor(meta[[celltype_col]])
  
  
  cat(
    "  Cells:", ncol(counts),
    "\n  Donors:", length(unique(donor)),
    "\n  Cell types:", length(unique(celltype)),
    "\n"
  )
  
  print(table(celltype))
  
  ## 4A. Donor x cell-type ----
  #     pseudobulk counts
  
  pb_group <- factor(
    paste(
      donor,
      as.character(celltype),
      sep = "||"
    )
  )
  
  # Sparse design matrix:
  # cells x donor-celltype groups
  Z_pb <- Matrix::sparse.model.matrix(
    ~ 0 + pb_group
  )
  
  colnames(Z_pb) <- levels(pb_group)
  
  # genes x donor-celltype
  pb_counts <- counts %*% Z_pb
  
  colnames(pb_counts) <- levels(pb_group)
  
  # Pseudobulk metadata
  
  first_cell <- match(
    levels(pb_group),
    as.character(pb_group)
  )
  
  n_cells_pb <- table(pb_group)
  
  pb_meta <- data.frame(
    pb_id = levels(pb_group),
    projid = donor[first_cell],
    cell_type = as.character(celltype[first_cell]),
    n_cells = as.integer(
      n_cells_pb[levels(pb_group)]
    ),
    stringsAsFactors = FALSE
  )
  
  pb_counts_list[[length(pb_counts_list) + 1]] <- pb_counts
  pb_meta_list[[length(pb_meta_list) + 1]] <- pb_meta
  
  
  ## 4B. Cell-level descriptive metrics for target genes ----
  
  genes_present_here <- intersect(
    genes,
    rownames(counts)
  )
  
  if (length(genes_present_here) > 0) {
    
    target_counts <- counts[
      genes_present_here,
      ,
      drop = FALSE
    ]
    
    # Total UMI count of each nucleus
    library_size <- Matrix::colSums(counts)
    
    # Common normalization across ALL RDS objects
    #
    # log(1 + counts / total_counts * 10,000)
    
    norm_target <- target_counts %*%
      Matrix::Diagonal(
        x = 10000 / pmax(library_size, 1)
      )
    
    # Sparse matrix: zeros remain zero
    norm_target@x <- log1p(norm_target@x)
    
    
    # Aggregate by high-resolution cell type
    
    # ct_levels <- levels(celltype)
    # 
    # Z_ct <- Matrix::sparse.model.matrix(
    #   ~ 0 + celltype
    # )
    # 
    # colnames(Z_ct) <- ct_levels
    
    ct_levels <- levels(celltype)
    
    if (length(ct_levels) == 1) {
      
      # All cells belong to one cell type
      Z_ct <- Matrix::Matrix(
        1,
        nrow = length(celltype),
        ncol = 1,
        sparse = TRUE
      )
      
      colnames(Z_ct) <- ct_levels
      
    } else {
      
      Z_ct <- Matrix::sparse.model.matrix(
        ~ 0 + celltype
      )
      
      colnames(Z_ct) <- ct_levels
    }
    
    # Number of cells expressing each target gene
    detected <- (target_counts > 0) %*% Z_ct
    
    # Sum of normalized expression
    norm_sum <- norm_target %*% Z_ct
    
    n_cells_ct <- as.numeric(
      table(celltype)[ct_levels]
    )
    
    # Save long-format summaries
    
    for (j in seq_along(ct_levels)) {
      
      cell_metrics_list[[length(cell_metrics_list) + 1]] <- data.frame(
        
        gene = genes_present_here,
        
        cell_type = ct_levels[j],
        
        n_cells = n_cells_ct[j],
        
        n_detected = as.numeric(
          detected[, j]
        ),
        
        sum_logCP10K = as.numeric(
          norm_sum[, j]
        )
      )
    }
  }
  
  # Free memory
  
  rm(
    obj,
    counts,
    meta,
    pb_counts,
    Z_pb
  )
  
  gc()
}


# Processing: Astrocytes.rds 
# Cells: 149558 
# Donors: 427 
# Cell types: 3 
# celltype
# Ast CHI3L1  Ast DPP10   Ast GRM3 
# 24252      32779      92527 
# 
# Processing: Excitatory_neurons_set1.rds 
# Cells: 296936 
# Donors: 427 
# Cell types: 1 
# celltype
# Exc L2-3 CBLN2 LINC02306 
# 296936 
# 
# Processing: Excitatory_neurons_set2.rds 
# Cells: 421529 
# Donors: 425 
# Cell types: 4 
# celltype
# Exc L3-4 RORB CUX2    Exc L3-5 RORB PLCH1   Exc L4-5 RORB GABRG1 Exc L4-5 RORB IL1RAPL2 
# 184784                  37949                  79361                 119435 
# 
# Processing: Excitatory_neurons_set3.rds 
# Cells: 324765 
# Donors: 426 
# Cell types: 9 
# celltype
# Exc L5 ET Exc L5-6 RORB LINC02196        Exc L5/6 IT Car3             Exc L5/6 NP               Exc L6 CT      Exc L6 THEMIS NFIA                 Exc L6b                Exc NRGN           Exc RELN CHD7 
# 3454                   22343                   18371                   17247                   23073                   66676                   25055                   45859                  102687 
# 
# Processing: Immune_cells.rds 
# Cells: 83889 
# Donors: 426 
# Cell types: 5 
# celltype
# CAMs  Mic MKI67 Mic P2RY12   Mic TPT1    T cells 
# 2167        866      73061       5261       2534 
# 
# Processing: Inhibitory_neurons.rds 
# Cells: 329699 
# Donors: 423 
# Cell types: 25 
# celltype
# Inh ALCAM TRPM3              Inh CUX2 MSR1           Inh ENOX2 SPHKAP          Inh FBN2 EPB41L4A              Inh GPC5 RIT2            Inh L1 PAX6 CA4         Inh L1-2 PAX6 SCGN        Inh L1-6 LAMP5 CA13          Inh L3-5 SST MAFB 
# 10897                      24885                      12396                       6769                       4788                       5070                       1617                      15060                      23294 
# Inh L5-6 PVALB STON2            Inh L5-6 SST TH             Inh L6 SST NPY   Inh LAMP5 NRG1 (Rosehip)             Inh LAMP5 RELN          Inh PTPRK FAM19A1 Inh PVALB CA8 (Chandelier)             Inh PVALB HTR4            Inh PVALB SULF1 
# 5286                       4502                       1404                      26065                       8813                      10438                      14644                      42032                      21231 
# Inh RYR3 TSHZ2             Inh SGCD PDE3A             Inh SORCS1 TTN             Inh VIP ABI3BP             Inh VIP CLSTN2             Inh VIP THSD7B              Inh VIP TSHZ2 
# 24230                       4030                       9610                      14448                      17631                       8211                      12348 
# 
# Processing: Oligodendrocytes.rds 
# Cells: 645142 
# Donors: 427 
# Cell types: 1 
# celltype
# Oli 
# 645142 
# 
# Processing: OPCs.rds 
# Cells: 90502 
# Donors: 427 
# Cell types: 1 
# celltype
# OPC 
# 90502 
# 
# Processing: Vasculature_cells.rds 
# Cells: 17974 
# Donors: 423 
# Cell types: 5 
# celltype
# End  Fib FLRT2 Fib SLC4A4        Per        SMC 
# 6514       3728        819       5308       1605 



###
# 5. Combine pseudobulk data ----
###

pb_all <- do.call(
  cbind,
  pb_counts_list
)

pb_meta <- bind_rows(
  pb_meta_list
)

stopifnot(
  identical(
    colnames(pb_all),
    pb_meta$pb_id
  )
)

if (anyDuplicated(pb_meta$pb_id)) {
  
  stop(
    "Duplicate donor/cell-type pseudobulk IDs found."
  )
}


saveRDS(pb_counts_list, "ROSMAP_gene_expression_analysis/pb_counts_list.rds")
saveRDS(pb_meta_list, "ROSMAP_gene_expression_analysis/pb_meta_list.rds")
saveRDS(pb_all, "ROSMAP_gene_expression_analysis/pb_all.rds")
saveRDS(pb_meta, "ROSMAP_gene_expression_analysis/pb_meta.rds")


###
# 6. Which requested genes exist? ----
###

genes_present <- intersect(
  genes,
  rownames(pb_all)
)

genes_missing <- setdiff(
  genes,
  rownames(pb_all)
)

cat(
  "\nGenes found:\n",
  paste(genes_present, collapse = ", "),
  "\n"
)

if (length(genes_missing) > 0) {
  
  cat(
    "\nGenes NOT found:\n",
    paste(genes_missing, collapse = ", "),
    "\n"
  )
}


###
# 7. Cell-level expression summary ----
###

cell_metrics <- bind_rows(
  cell_metrics_list
) %>%
  
  group_by(gene, cell_type) %>%
  
  summarise(
    
    n_cells = sum(n_cells),
    
    n_detected = sum(n_detected),
    
    sum_logCP10K = sum(sum_logCP10K),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    pct_expr =
      100 * n_detected / n_cells,
    
    mean_logCP10K =
      sum_logCP10K / n_cells
  )

write.table(cell_metrics, "ROSMAP_gene_expression_analysis/Colleen_genes_cell_metrics.txt", sep="\t", row.names = F, col.names = T)


###
# 8. Transcriptome-wide expression percentile ----
#
# This is what lets us define whether a target gene is
# "highly expressed" relative to all ~33,000 genes.
###

celltypes_all <- sort(
  unique(pb_meta$cell_type)
)

ct_factor <- factor(
  pb_meta$cell_type,
  levels = celltypes_all
)

Z_ct_all <- Matrix::sparse.model.matrix(
  ~ 0 + ct_factor
)

colnames(Z_ct_all) <- celltypes_all


# Sum all donor pseudobulks belonging to each cell type
# genes x cell-types
ct_counts <- pb_all %*% Z_ct_all


# Total counts in each cell type
ct_library_size <- Matrix::colSums(
  ct_counts
)

# Convert to CPM

ct_cpm <- as.matrix(ct_counts)

ct_cpm <- sweep(
  ct_cpm,
  2,
  ct_library_size,
  "/"
)

ct_cpm <- ct_cpm * 1e6

# Percentile of EVERY gene within each cell type

expression_percentile <- apply(
  
  ct_cpm,
  
  2,
  
  function(v) {
    
    100 *
      rank(
        v,
        ties.method = "average"
      ) /
      length(v)
  }
)

# Extract only genes of interest

target_cpm <- ct_cpm[
  genes_present,
  ,
  drop = FALSE
]

target_percentile <- expression_percentile[
  genes_present,
  ,
  drop = FALSE
]


absolute_expression <- data.frame(
  
  gene = rep(
    rownames(target_cpm),
    times = ncol(target_cpm)
  ),
  
  cell_type = rep(
    colnames(target_cpm),
    each = nrow(target_cpm)
  ),
  
  pooled_CPM = as.vector(
    target_cpm
  ),
  
  expression_percentile = as.vector(
    target_percentile
  )
)


###
# 9. Add donor counts ----
###

celltype_donor_counts <- pb_meta %>%
  
  group_by(cell_type) %>%
  
  summarise(
    total_donors =
      n_distinct(projid),
    
    .groups = "drop"
  )


###
# 10. Absolute-expression classification ----
###

expression_summary <- cell_metrics %>%
  
  left_join(
    absolute_expression,
    by = c(
      "gene",
      "cell_type"
    )
  ) %>%
  
  left_join(
    celltype_donor_counts,
    by = "cell_type"
  ) %>%
  
  mutate(
    
    expression_class = case_when(
      
      expression_percentile >=
        very_high_percentile &
        
        pct_expr >=
        very_high_detection
      ~ "Very high",
      
      expression_percentile >=
        high_percentile &
        
        pct_expr >=
        high_detection
      ~ "High",
      
      TRUE ~ "Not high"
    )
  )

write.table(expression_summary, "ROSMAP_gene_expression_analysis/Colleen_genes_expression_summary.txt", sep="\t", row.names = F, col.names = T)


###
# 11. Prepare donor-level pseudobulk data ----
###

# Minimum cell count per donor/cell-type
keep_pb <- (
  pb_meta$n_cells >=
    min_cells_per_pseudobulk
)

# Require enough donors per cell type

valid_celltypes <- pb_meta[
  keep_pb,
] %>%
  
  count(
    cell_type,
    name = "n_donors"
  ) %>%
  
  filter(
    n_donors >=
      min_donors_per_celltype
  ) %>%
  
  pull(cell_type)


keep_pb <- (
  keep_pb &
    pb_meta$cell_type %in%
    valid_celltypes
)


pb_meta_qc <- pb_meta[
  keep_pb,
  ,
  drop = FALSE
]

pb_qc <- pb_all[
  ,
  keep_pb,
  drop = FALSE
]


cat(
  "\nCell types passing pseudobulk QC:\n"
)

print(
  table(pb_meta_qc$cell_type)
)

write.table(pb_meta_qc, "ROSMAP_gene_expression_analysis/All_genes_pb_meta_qc.txt", sep="\t", row.names = F, col.names = T)


###
# 12. Donor-level CPM ----
#
# IMPORTANT:
# library size comes from ALL genes,
# not just genes of interest.
###

pb_library_size <- Matrix::colSums(
  pb_qc
)

target_pb_counts <- pb_qc[
  genes_present,
  ,
  drop = FALSE
]

target_pb_cpm <- target_pb_counts %*%
  Matrix::Diagonal(
    x = 1e6 /
      pmax(pb_library_size, 1)
  )

# Convert small target-gene matrix to ordinary matrix

pb_cpm_matrix <- t(
  as.matrix(
    target_pb_cpm
  )
)

rownames(pb_cpm_matrix) <-
  colnames(target_pb_counts)

# Make long table

pb_expression <- as.data.frame(
  pb_cpm_matrix
) %>%
  
  rownames_to_column(
    "pb_id"
  ) %>%
  
  left_join(
    pb_meta_qc,
    by = "pb_id"
  ) %>%
  
  pivot_longer(
    
    cols =
      all_of(genes_present),
    
    names_to = "gene",
    
    values_to = "CPM"
  ) %>%
  
  mutate(
    
    log2CPM1 =
      log2(CPM + 1)
  )


###
# 13. Donor-balanced expression per cell type ----
###

donor_summary <- pb_expression %>%
  
  group_by(
    gene,
    cell_type
  ) %>%
  
  summarise(
    
    n_donors =
      n_distinct(projid),
    
    mean_donor_CPM =
      mean(CPM),
    
    median_donor_CPM =
      median(CPM),
    
    mean_log2CPM1 =
      mean(log2CPM1),
    
    median_log2CPM1 =
      median(log2CPM1),
    
    .groups = "drop"
  )


###
# 14. Descriptive cell-type enrichment ----
#
# Compares each cell-type mean with the mean expression
# across the OTHER cell types.
#
# Equal weighting of cell types.
###

donor_summary <- donor_summary %>%
  
  group_by(gene) %>%
  
  mutate(
    
    other_celltype_mean_CPM =
      
      (
        sum(mean_donor_CPM) -
          mean_donor_CPM
      ) /
      (n() - 1),
    
    fold_enrichment =
      
      (
        mean_donor_CPM + 0.1
      ) /
      (
        other_celltype_mean_CPM + 0.1
      ),
    
    log2_fold_enrichment =
      log2(fold_enrichment),
    
    specificity_z =
      
      ifelse(
        
        sd(mean_log2CPM1) > 0,
        
        (
          mean_log2CPM1 -
            mean(mean_log2CPM1)
        ) /
          sd(mean_log2CPM1),
        
        0
      )
  ) %>%
  
  ungroup()

write.table(donor_summary, "ROSMAP_gene_expression_analysis/Colleen_genes_donor_summary.txt", sep="\t", row.names = F, col.names = T)


###
# 15. Donor-aware statistical enrichment test ----
#
# Mixed model:
#
# expression ~ cell type + (1 | donor)
#
# For every gene, each cell type is compared against
# the equally weighted average of all OTHER cell types.
###

mixed_model_results <- list()


for (g in genes_present) {
  
  cat(
    "\nMixed model:",
    g,
    "\n"
  )
  
  dd <- pb_expression %>%
    filter(gene == g)
  
  
  # Only cell types represented in this gene's data
  ct_levels <- sort(
    unique(dd$cell_type)
  )
  
  if (length(ct_levels) < 2) {
    next
  }
  
  
  dd$cell_type <- factor(
    dd$cell_type,
    levels = ct_levels
  )
  
  # Random intercept for donor
  
  fit <- lmerTest::lmer(
    
    log2CPM1 ~
      0 + cell_type +
      (1 | projid),
    
    data = dd,
    
    REML = FALSE
  )
  
  # Estimated marginal means
  
  emm <- emmeans::emmeans(
    fit,
    ~ cell_type
  )
  
  # Build:
  #
  # Cell type i
  #      versus
  # mean(other cell types)
  
  K <- length(ct_levels)
  
  contrast_list <- list()
  
  for (i in seq_len(K)) {
    
    weights <- rep(
      -1 / (K - 1),
      K
    )
    
    weights[i] <- 1
    
    contrast_list[[paste0(ct_levels[i]," vs rest")]] <- weights
  }
  
  
  contrast_result <- emmeans::contrast(
    
    emm,
    
    method = contrast_list,
    
    adjust = "none"
  )
  
  
  tmp <- as.data.frame(
    contrast_result
  )
  
  tmp$gene <- g
  
  tmp$cell_type <- sub(
    " vs rest$",
    "",
    tmp$contrast
  )
  
  
  mixed_model_results[[length(mixed_model_results) + 1]] <- tmp
}


# What is happening is that emmeans sees >3,000 pseudobulk observations and decides not to calculate 
# Satterthwaite/Kenward-Roger degrees of freedom because those calculations can become expensive. 
# Instead, it effectively uses asymptotic inference (infinite df / normal approximation).
# 
# With ~12,000 donor × cell-type observations, that approximation is very reasonable. 
# The difference between a large-df t distribution and a normal distribution becomes negligible.


# Combine statistical results

mixed_model_results <- bind_rows(
  mixed_model_results
)

# Global FDR across all gene x cell-type tests

mixed_model_results <- mixed_model_results %>%
  
  mutate(
    FDR = p.adjust(
      p.value,
      method = "BH"
    )
  ) %>%
  
  dplyr::select(
    gene,
    cell_type,
    enrichment_estimate = estimate,
    enrichment_SE = SE,
    any_of(c("df", "t.ratio", "z.ratio")),
    p.value,
    FDR
  )


###
# 16. Build final master table ----
###

final_results <- expression_summary %>%
  
  left_join(
    
    donor_summary,
    
    by = c(
      "gene",
      "cell_type"
    )
  ) %>%
  
  left_join(
    
    mixed_model_results,
    
    by = c(
      "gene",
      "cell_type"
    )
  )

write.table(final_results, "ROSMAP_gene_expression_analysis/Colleen_genes_final_results.txt", sep="\t", row.names = F, col.names = T)


###
# 17. Top 1 / Top 2 cell type for each gene ----
#
# Ranking based on donor-balanced mean log2(CPM + 1).
###

top2_long <- donor_summary %>%
  
  group_by(gene) %>%
  
  arrange(
    desc(mean_log2CPM1),
    .by_group = TRUE
  ) %>%
  
  slice_head(
    n = 2
  ) %>%
  
  mutate(
    rank = row_number()
  ) %>%
  
  ungroup()

write.table(top2_long, "ROSMAP_gene_expression_analysis/Colleen_genes_top2_long.txt", sep="\t", row.names = F, col.names = T)


# Wide, easy-to-read table

top2_table <- top2_long %>%
  
  dplyr::select(
    gene,
    rank,
    cell_type,
    mean_donor_CPM,
    mean_log2CPM1,
    log2_fold_enrichment
  ) %>%
  
  tidyr::pivot_wider(
    names_from = rank,
    
    values_from = c(
      cell_type,
      mean_donor_CPM,
      mean_log2CPM1,
      log2_fold_enrichment
    ),
    
    names_glue = "top{rank}_{.value}"
  )

write.table(top2_table, "ROSMAP_gene_expression_analysis/Colleen_genes_top2_table.txt", sep="\t", row.names = F, col.names = T)


###
# 18. Save tables ----
###

# write.csv(
#   
#   final_results,
#   
#   file.path(
#     out_dir,
#     "gene_celltype_full_results.csv"
#   ),
#   
#   row.names = FALSE
# )
# 
# 
# write.csv(
#   
#   top2_table,
#   
#   file.path(
#     out_dir,
#     "gene_top2_celltypes.csv"
#   ),
#   
#   row.names = FALSE
# )


write.csv(
  
  mixed_model_results,
  
  file.path(
    out_dir,
    "Colleen_genes_gene_celltype_enrichment_statistics.csv"
  ),
  
  row.names = FALSE
)


###
# 19. Visualization 1: ----
#     classic single-cell dot plot
#
# Size = % nuclei expressing gene
# Color = mean normalized expression
###

# Assign each high-resolution cell type to a broad class

celltype_info <- final_results %>%
  
  dplyr::distinct(cell_type) %>%
  
  dplyr::mutate(
    
    major_class = dplyr::case_when(
      
      grepl("^Mic", cell_type) ~ "Microglia",
      grepl("^Oli", cell_type) ~ "Oligodendrocytes",
      
      grepl("^Exc", cell_type) ~ "Excitatory neurons",
      grepl("^Inh", cell_type) ~ "Inhibitory neurons",
      
      grepl("^CAM", cell_type) ~ "CAMs",
      grepl("^End", cell_type) ~ "Endothelial",
      
      grepl("^Ast", cell_type) ~ "Astrocytes",
      grepl("^Fib", cell_type) ~ "Fibroblasts",
      
      grepl("^Per", cell_type) ~ "Pericytes",
      grepl("^SMC", cell_type) ~ "SMC",
      
      TRUE ~ "Other"
    )
  )


# Specify biological order of the broad classes

major_class_order <- c(
  "Microglia",
  "Oligodendrocytes",
  "Excitatory neurons",
  "Inhibitory neurons",
  "CAMs",
  "Endothelial",
  "Astrocytes",
  "Fibroblasts",
  "Pericytes",
  "SMC",
  "Other"
)


# Calculate mean expression for sorting WITHIN each class

celltype_order_df <- final_results %>%
  
  dplyr::group_by(cell_type) %>%
  
  dplyr::summarise(
    
    overall = mean(
      mean_logCP10K,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::left_join(
    celltype_info,
    by = "cell_type"
  ) %>%
  
  dplyr::mutate(
    
    major_class = factor(
      major_class,
      levels = major_class_order
    )
  ) %>%
  
  # broad biological group first;
  # then high-expression subtypes first within each group
  dplyr::arrange(
    major_class,
    dplyr::desc(overall)
  )


# Top-to-bottom ordering

celltype_order <- celltype_order_df$cell_type

# Factor levels actually used by the y-axis
y_levels <- rev(celltype_order)

# Match each displayed cell type to its major class
y_info <- data.frame(
  cell_type = y_levels,
  y = seq_along(y_levels)
) %>%
  dplyr::left_join(
    celltype_info %>% dplyr::select(cell_type, major_class),
    by = "cell_type"
  )

# Find boundaries where the major cell class changes
group_boundaries <- which(
  y_info$major_class[-1] !=
    y_info$major_class[-nrow(y_info)]
) + 0.5


# Plot itself

p_dot <- ggplot(
  final_results,
  aes(
    x = gene,
    y = factor(
      cell_type,
      levels = rev(celltype_order)
    )
  )
) +
  
  # Lines separating major cell-type groups
  geom_hline(
    yintercept = group_boundaries,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_point(
    aes(
      size = pct_expr,
      color = mean_logCP10K
    )
  ) +
  
  scale_size_continuous(
    range = c(1, 10)
  ) +
  
  scale_color_viridis_c(
    option = "plasma",
    direction = 1
  ) +
  
  labs(
    x = NULL,
    y = "Cell type",
    size = "% nuclei\nexpressing",
    color = "Mean\nlog(CP10K + 1)",
    title = "ROSMAP snRNA-seq expression"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    ),
    panel.grid.major = element_line(
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank()
  )

p_dot


ggsave(
  
  file.path(
    out_dir,
    "01_expression_dotplot.pdf"
  ),
  
  p_dot,
  
  width = 12,
  height = 9
)


###
# 20. Visualization 2: ----
#     cell-type enrichment heatmap
#
# Positive = enriched
# Negative = depleted
###

p_enrich <- ggplot(
  
  final_results,
  
  aes(
    x = gene,
    y = factor(
      cell_type,
      levels = rev(celltype_order)
    ),
    fill = log2_fold_enrichment
  )
  
) +
  
  geom_tile() +
  
  # Separate major cell-type groups
  geom_hline(
    yintercept = group_boundaries,
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  
  labs(
    x = NULL,
    y = "Cell type",
    fill = "log2 fold\nenrichment",
    title = "Cell-type expression enrichment"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    ),
    panel.grid = element_blank()
  )

p_enrich

ggsave(
  
  file.path(
    out_dir,
    "02_celltype_enrichment_heatmap.pdf"
  ),
  
  p_enrich,
  
  width = 12,
  height = 9
)


###
# 21. Visualization 3: ----
#     expression percentile heatmap
#
# Makes the "highly expressed" definition intuitive
###

p_percentile <- ggplot(
  
  final_results,
  
  aes(
    x = gene,
    y = factor(
      cell_type,
      levels = rev(celltype_order)
    ),
    fill = expression_percentile
  )
  
) +
  
  geom_tile() +
  
  # Separate major cell-type groups
  geom_hline(
    yintercept = group_boundaries,
    color = "black",
    linewidth = 0.5
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    direction = 1,
    limits = c(0, 100)
  ) +
  
  labs(
    x = NULL,
    y = "Cell type",
    fill = "Expression\npercentile",
    title = "Expression percentile within each cell type"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    ),
    panel.grid = element_blank()
  )

p_percentile

ggsave(
  
  file.path(
    out_dir,
    "03_expression_percentile_heatmap.pdf"
  ),
  
  p_percentile,
  
  width = 12,
  height = 9
)


###
# 22. Visualization 4: ----
#     top 2 cell types per gene
###

p_top2 <- ggplot(
  
  top2_long,
  
  aes(
    x = mean_log2CPM1,
    y = reorder(
      cell_type,
      mean_log2CPM1
    )
  )
  
) +
  
  geom_col() +
  
  facet_wrap(
    ~ gene,
    scales = "free"
  ) +
  
  labs(
    
    x =
      "Mean donor log2(CPM + 1)",
    
    y = NULL,
    
    title =
      "Top two expressing cell types for each gene"
  ) +
  
  theme_bw()

p_top2

ggsave(
  
  file.path(
    out_dir,
    "04_top2_celltypes.pdf"
  ),
  
  p_top2,
  
  width = 12,
  height = 8
)


###
# 23. Visualization 5: ----
#     High-expression diagnostic plot
#
# Upper-right region corresponds to:
# high transcriptome percentile + broadly detected
###

p_high <- ggplot(
  
  final_results,
  
  aes(
    x = expression_percentile,
    y = pct_expr
  )
  
) +
  
  geom_point(
    aes(
      shape = expression_class
    ),
    size = 3
  ) +
  
  geom_vline(
    xintercept = high_percentile,
    linetype = 2
  ) +
  
  geom_hline(
    yintercept = high_detection,
    linetype = 2
  ) +
  
  facet_wrap(
    ~ gene
  ) +
  
  labs(
    
    x =
      "Expression percentile within cell type",
    
    y =
      "% nuclei expressing gene",
    
    shape =
      "Expression class",
    
    title =
      "Classification of highly expressed genes"
  ) +
  
  theme_bw()

p_high

ggsave(
  
  file.path(
    out_dir,
    "05_high_expression_diagnostic.pdf"
  ),
  
  p_high,
  
  width = 12,
  height = 9
)



cat(
  "\nAnalysis complete.\nResults saved to:\n",
  normalizePath(out_dir),
  "\n"
)
