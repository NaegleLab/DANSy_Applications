# This script is to process the adipocyte dataset associated with the publication Harms et al. Cell Reports 2020

# The paper can be found at https://doi.org/10.1016/j.celrep.2019.03.026 and has the PMID:30943403


# 1. Workspace set up -----------------------------------------------------

library(dplyr)
library(DESeq2)
library(biomaRt)
library(GEOquery)

counts_matrix <-read.table('GSE180568_raw_counts_GRCh38.p13_NCBI.tsv.gz', sep = '\t', header=TRUE)

# Full metadata that is available for these given that the NCBI version is automatically run and only maintians th sample ID
full_metadata = getGEO('GSE180568')
phenotype_data = full_metadata$GSE180568_series_matrix.txt.gz@phenoData@data
metaData <- phenotype_data['title']
metaData['conditions'] = lapply(phenotype_data['title'], function(x) substr(x, 1,nchar(x)-3))
metaData = metaData['conditions']

# Performing DEG analysis -------------------------------------------------
proc_data <- counts_matrix
row.names(proc_data)  =proc_data[,'GeneID']
proc_data = proc_data[,-1]

# Prefiltering to reduce runtime
proc_data = proc_data[rowSums(proc_data) > 10,]

dds <- DESeqDataSetFromMatrix(proc_data, metaData, design = ~conditions)
dds <- DESeq(dds)

# Simple check to make sure there is some reproduction of the original paper's result
vsd <- vst(dds, blind = T)
plotPCA(vsd, intgroup = 'conditions')

res <- results(dds, contrast = c('conditions','MN_4.11','MN_A375'), alpha = .05,lfcThreshold = 1)
res <- data.frame(res)

suffix = '_4-11_v_A375'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res <- res

res <- results(dds, contrast = c('conditions','MN_2.18_0','MN_A375'), alpha = .05)
res <- data.frame(res)
suffix = '_2-18_v_A375'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_4.11','MN_2.18_0'), alpha = .05)
res <- data.frame(res)
suffix = '_4-11_v_2-18'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT34_nodrug','MN_A375'), alpha = .05)
res <- data.frame(res)

suffix = '_CRT34_v_A375'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT35_nodrug','MN_A375'), alpha = .05)
summary(res)
res <- data.frame(res)

suffix = '_CRT35_v_A375'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT35_nodrug','MN_CRT34_nodrug'), alpha = .05)
summary(res)
res <- data.frame(res)

suffix = '_CRT35_v_CRT34'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT35_nodrug','MN_CRT35_drug'), alpha = .05)
summary(res)
res <- data.frame(res)

suffix = '_CRT35_v_CRT35drug'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT34_nodrug','MN_CRT34_drug'), alpha = .05)
summary(res)
res <- data.frame(res)

suffix = '_CRT34_v_CRT34drug'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT35_nodrug','MN_2.18_0'), alpha = .05)
summary(res)
res <- data.frame(res)

suffix = '_CRT35_v_2-18'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT35_nodrug','MN_4.11'), alpha = .05)
summary(res)
res <- data.frame(res)


suffix = '_CRT35_v_4-11'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT34_nodrug','MN_2.18_0'), alpha = .05)
res <- data.frame(res)
suffix = '_CRT34_v_2-18'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','MN_CRT34_nodrug','MN_4.11'), alpha = .05)
summary(res)
res <- data.frame(res)

suffix = '_CRT34_v_4-11'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

# Setting up export but first setting up ID conversions -------------------

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ids <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name','entrezgene_id'),
                  filters = 'entrezgene_id',
                  values = full_deg_res['gene_id'],
                  mart = ensembl)

# Merging the results with the gene conversion matrix
data_to_export = merge(gene_ids, full_deg_res, by.x = 'entrezgene_id', by.y = 'gene_id')
write.table(data_to_export, file= 'SOX10_invasion_DEG_Results.csv', sep = ',')
                       