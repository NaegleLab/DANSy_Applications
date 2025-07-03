# This script is to process the adipocyte dataset associated with the publication Harms et al. Cell Reports 2020

# The paper can be found at https://doi.org/10.1016/j.celrep.2019.03.026 and has the PMID:30943403


# 1. Workspace set up -----------------------------------------------------

library(dplyr)
library(DESeq2)
library(biomaRt)
library(GEOquery)

counts_matrix <-read.table('GSE127887_raw_counts_GRCh38.p13_NCBI.tsv.gz', sep = '\t', header=TRUE)
full_metadata = getGEO('GSE127887')
phenotype_data = full_metadata$GSE127887_series_matrix.txt.gz@phenoData@data
metaData <- phenotype_data['title']
metaData['conditions'] = unlist(lapply(strsplit(phenotype_data[['title']],'-'), function(x) x[2]))
metaData = metaData['conditions']
# Performing DEG analysis -------------------------------------------------
proc_data <- counts_matrix
row.names(proc_data)  =proc_data[,'GeneID']
proc_data = proc_data[,-1]

# Prefiltering to reduce runtime
proc_data = proc_data[rowSums(proc_data) > 10,]

dds <- DESeqDataSetFromMatrix(proc_data, metaData, design = ~conditions)
dds <- DESeq(dds)

# Checking some quality assurance and which comparisons may not be worth doing
vsd <- vst(dds, blind = T)
plotPCA(vsd, intgroup = 'conditions')

# QC again on just a handful of the comparisons to ensure no weird issues
res <- results(dds, contrast = c('conditions','soft','stiff'), alpha = .05)
resLFC <- lfcShrink(dds, coef="conditions_stiff_vs_soft", type="apeglm")
plotMA(resLFC,ylim=c(-5,5), alpha=0.05)


# Data frame set up -------------------------------------------------------
res <- data.frame(res)

suffix = '_soft_stiff'
res <-  res %>%
  dplyr::select(log2FoldChange, stat,padj) 
colnames(res) <- paste0(colnames(res), suffix)
res['gene_id'] = rownames(res)
full_deg_res <- res

# Setting up export but first setting up ID conversions -------------------

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror='useast')
gene_ids <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol','external_gene_name','entrezgene_id'),
                  filters = 'entrezgene_id',
                  values = full_deg_res['gene_id'],
                  mart = ensembl)

# Merging the results with the gene conversion matrix
data_to_export = merge(gene_ids, full_deg_res, by.x = 'entrezgene_id', by.y = 'gene_id')
write.table(data_to_export, file= 'Mechanical_Conditioning_DEG_Results.csv', sep = ',')
