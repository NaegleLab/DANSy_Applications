# This script is to process the adipocyte dataset associated with the publication Harms et al. Cell Reports 2020

# The paper can be found at https://doi.org/10.1016/j.celrep.2019.03.026 and has the PMID:30943403


# 1. Workspace set up -----------------------------------------------------

library(dplyr)
library(DESeq2)
library(biomaRt)
library(GEOquery)

counts_matrix <-read.table('GSE129290_raw_counts_GRCh38.p13_NCBI.tsv.gz', sep = '\t', header=TRUE)
full_metadata = getGEO('GSE129290')
phenotype_data = full_metadata$GSE129290_series_matrix.txt.gz@phenoData@data
metaData <- phenotype_data['title']
metaData['conditions'] = lapply(phenotype_data['title'], function(x) substr(x, 7,nchar(x)-2))
metaData = metaData['conditions']
# Performing DEG analysis -------------------------------------------------

# For this dataset I will be splitting up the counts matrices into 2 separate experiments. This is due to the original paper using and conducting the sequencing separately. Further, the iPSCs used came from likely different batches.
proc_data <- counts_matrix[,1:13]
row.names(proc_data)  =proc_data[,'GeneID']
proc_data = proc_data[,-1]

# Prefiltering to reduce runtime
proc_data = proc_data[rowSums(proc_data) > 10,]
metaData_temp = metaData %>% dplyr::slice(1:12)
dds <- DESeqDataSetFromMatrix(proc_data, metaData_temp, design = ~conditions)
dds <- DESeq(dds)

# Checking some quality assurance and which comparisons may not be worth doing
vsd <- vst(dds, blind = T)
plotPCA(vsd, intgroup = 'conditions')

# QC again on just a handful of the comparisons to ensure no weird issues
res <- results(dds, contrast = c('conditions','2.4dyn','0dyn_Static'), alpha = .05)


# Data frame set up -------------------------------------------------------
res <- data.frame(res)

suffix = '_2-4dyn_v_0dyn'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res <- res

res <- results(dds, contrast = c('conditions','0.5dyn','0dyn_Static'), alpha = .05)
res <- data.frame(res)
suffix = '_0-5dyn_v_0dyn'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','0.01dyn','0dyn_Static'), alpha = .05)
res <- data.frame(res)
suffix = '_0-01dyn_v_0dyn'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','0.5dyn','0.01dyn'), alpha = .05)
res <- data.frame(res)
suffix = '_0-5dyn_v_0-01dyn'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','2.4dyn','0.01dyn'), alpha = .05)
res <- data.frame(res)
suffix = '_2-4dyn_v_0-01dyn'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','2.4dyn','0.5dyn'), alpha = .05)
res <- data.frame(res)
suffix = '_2-4dyn_v_0-5dyn'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# Setting up export but first setting up ID conversions -------------------

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror='useast')
gene_ids <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name','entrezgene_id'),
                  filters = 'entrezgene_id',
                  values = full_deg_res['gene_id'],
                  mart = ensembl)

# Merging the results with the gene conversion matrix
data_to_export = merge(gene_ids, full_deg_res, by.x = 'entrezgene_id', by.y = 'gene_id')
write.table(data_to_export, file= 'iBMEC_Flow_DEG_Results.csv', sep = ',')


# Same as above but with the co-cultures ----------------------------------

proc_data <- counts_matrix[,c(1,14:ncol(counts_matrix))]
row.names(proc_data)  =proc_data[,'GeneID']
proc_data = proc_data[,-1]

# Prefiltering to reduce runtime
proc_data = proc_data[rowSums(proc_data) > 10,]
metaData_temp = metaData %>% dplyr::slice(13:nrow(metaData))
dds <- DESeqDataSetFromMatrix(proc_data, metaData_temp, design = ~conditions)
dds <- DESeq(dds)

# Checking some quality assurance and which comparisons may not be worth doing
vsd <- vst(dds, blind = T)
plotPCA(vsd, intgroup = 'conditions')

# QC again on just a handful of the comparisons to ensure no weird issues
res <- results(dds, contrast = c('conditions','cocultured_AP','no_coculture'), alpha = .05)

# Data frame set up -------------------------------------------------------
res <- data.frame(res)

suffix = '_AP_v_no'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res <- res

res <- results(dds, contrast = c('conditions','cocultured_iNeural','no_coculture'), alpha = .05)
res <- data.frame(res)
suffix = '_iNeural_v_no'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','cocultured_AP','cocultured_iNeural'), alpha = .05)
res <- data.frame(res)
suffix = '_AP_v_iNeural'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

gene_ids <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name','entrezgene_id'),
                  filters = 'entrezgene_id',
                  values = full_deg_res['gene_id'],
                  mart = ensembl)

# Merging the results with the gene conversion matrix
data_to_export = merge(gene_ids, full_deg_res, by.x = 'entrezgene_id', by.y = 'gene_id')
write.table(data_to_export, file= 'iBMEC_Coculture_DEG_Results.csv', sep = ',')
