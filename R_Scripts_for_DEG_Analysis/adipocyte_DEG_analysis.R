# This script is to process the adipocyte dataset associated with the publication Harms et al. Cell Reports 2020

# The paper can be found at https://doi.org/10.1016/j.celrep.2019.03.026 and has the PMID:30943403


# 1. Workspace set up -----------------------------------------------------

library(dplyr)
library(DESeq2)
library(biomaRt)

datafiles <- dir('GSE115020_RAW/')

counts_matrix <- as.data.frame(read.table(paste0('GSE115020_RAW/',datafiles[1]), header = TRUE))

for (file in datafiles) {
x <- as.data.frame(read.table(paste0('GSE115020_RAW/',file), header = TRUE))
counts_matrix <- merge(counts_matrix, x)
}

# Removing the Raw counts from the column names
sample_names <- colnames(counts_matrix)[3:ncol(counts_matrix)]
colnames(counts_matrix)[3:ncol(counts_matrix)] <- lapply(sample_names, function(x) gsub(pattern='Raw_count_Hs_',replacement = '',x=x))

# Setting up the metadata which tells which condition is what
sample_names <- colnames(counts_matrix)[3:ncol(counts_matrix)]
conditions = unlist(lapply(sample_names, function(x) strsplit(x,split = '_rep')[[1]][[1]]))
metaData = data.frame(Condition = conditions, row.names = sample_names)


# Starting DEG Analysis ---------------------------------------------------

proc_data <- counts_matrix
row.names(proc_data) =proc_data[,'id']
proc_data <- proc_data[, !(names(proc_data) %in% c('id','symbol'))]

# Dropping any genes which were extremely low counts (<10) across all samples
proc_data <- proc_data[rowSums(proc_data) > 10,]

dds <- DESeqDataSetFromMatrix(proc_data, metaData, design = ~Condition)
dds <- DESeq(dds)

# Simple check to make sure there is some reproduction of the original paper's result
vsd <- vst(dds, blind = T)
plotPCA(vsd, intgroup = 'Condition')


# Doing some targeted DEG analysis ------------------------------------------

# To limit some of the data that will be analyzed just staying at the 2 week time points and some of the naive cell states

# Differentiated pre-adipocytes versus adipocytes cultured in a floating culture
res <- results(dds, contrast = c('Condition','float_2_wk','Diff_PA_2_wk'), alpha = .05)
res <- data.frame(res)

suffix = '_Float_v_Diff'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res <- res

# FLoat versus transwells at 2 weeks
res <- results(dds, contrast = c('Condition','float_2_wk','TW_2_wk'), alpha = .05)
res <- data.frame(res)

suffix = '_Float_v_TW'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# FLoat versus explant at 2 weeks
res <- results(dds, contrast = c('Condition','float_2_wk','explant_2_wk'), alpha = .05)
res <- data.frame(res)

suffix = '_Float_v_Explant'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# FLoat versus fresh tissue
res <- results(dds, contrast = c('Condition','float_2_wk','whole_tissue_fresh'), alpha = .05)
res <- data.frame(res)

suffix = '_Float_v_Fresh'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# FLoat versus packed cells
res <- results(dds, contrast = c('Condition','float_2_wk','packed_cells_D0'), alpha = .05)
res <- data.frame(res)

suffix = '_Float_v_Packed'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# Packed versus Differentiated pre-adipocytes
res <- results(dds, contrast = c('Condition','packed_cells_D0', 'Diff_PA_2_wk'), alpha = .05)
res <- data.frame(res)

suffix = '_Packed_v_Diff'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# Whole fresh versus explant
res <- results(dds, contrast = c('Condition','explant_2_wk','whole_tissue_fresh'), alpha = .05)
res <- data.frame(res)

suffix = '_Explant_v_Whole'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# Differentiated versus non-differentiated
res <- results(dds, contrast = c('Condition','Diff_PA_2_wk','Diff_PA_not_differentiated'), alpha = 1e-13)
res <- data.frame(res)

suffix = '_Diff_v_Undiff'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')


# Exporting the results ---------------------------------------------------

write.table(full_deg_res, 'Adipocyte_Culture_Methods_DEG_Results.csv',sep = ',',row.names = FALSE)

