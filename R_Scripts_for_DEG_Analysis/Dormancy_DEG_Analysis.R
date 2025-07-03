# This script is to process the adipocyte dataset associated with the publication Harms et al. Cell Reports 2020

# The paper can be found at https://doi.org/10.1016/j.celrep.2019.03.026 and has the PMID:30943403


# 1. Workspace set up -----------------------------------------------------

library(dplyr)
library(DESeq2)
library(biomaRt)
library(GEOquery)

counts_matrix <-read.table('GSE182890_raw_counts_GRCh38.p13_NCBI.tsv.gz', sep = '\t', header=TRUE)

# Full metadata that is available for these given that the NCBI version is automatically run and only maintians th sample ID
full_metadata = getGEO('GSE182890')
phenotype_data = full_metadata$GSE182890_series_matrix.txt.gz@phenoData@data
metaData <- phenotype_data['title']
metaData['conditions'] = unlist(lapply(strsplit(phenotype_data[['title']],'\\['), function(x) substr(x[2], 1,nchar(x[2])-2)))
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

res <- results(dds, contrast = c('conditions','DHEP','THEP'), alpha = .05)
res <- data.frame(res)

suffix = '_DHEP_v_THEP'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res <- res

res <- results(dds, contrast = c('conditions','DHEsC','THEP'), alpha = .05)
res <- data.frame(res)
suffix = '_DHEPsC_v_THEP'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','DHEsC','DHEP'), alpha = .05)
res <- data.frame(res)
suffix = '_DHEPsC_v_DHEP'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','DHEPsD','DHEP'), alpha = .05)
res <- data.frame(res)
suffix = '_DHEPsD_v_DHEP'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','DHEPsD','DHEsC'), alpha = .05)
res <- data.frame(res)
suffix = '_DHEPsD_v_DHEPsC'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res = merge(full_deg_res, res, by='gene_id')

res <- results(dds, contrast = c('conditions','DHEPsD','THEP'), alpha = .05)
res <- data.frame(res)
suffix = '_DHEPsD_v_THEP'
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
write.table(data_to_export, file= 'Dormancy_DDR_DEG_Results.csv', sep = ',')
