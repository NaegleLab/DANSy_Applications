# R Scripts for generating the differential gene expression analysis results.

Each dataset that was analyzed has a corresponding R script which generates a csv file containing that matches what is provided in the [DEG_data](DEG_data) folder. For each of dataset, the raw counts are retrieved from the NCBI GEO repository using the accession number provided in the Table below. We retrieved the raw RNA-sequencing counts matrix using the NCBI-generated version and these were used for the analysis. 

The differential gene expression analysis was performed using the DESeq2 pipeline. For each dataset, in the R script a principle component analysis (PCA) was performed to visualize how variable the transcriptome of each condition was globally.

We have provided a renv environment to help reproduce the complete R project.

Below is the table that contains GEO accession information, links to the study, and links to obtain the raw reads.

| Dataset Name | Original Study | PMID | GEO accession | Link to data download |
| - | - | - | - | - |
| Adipocyte Culturing Methods | [Harms et al. Cell Reports 2019](https://doi.org/10.1016/j.celrep.2019.03.026) | [30943403](https://pubmed.ncbi.nlm.nih.gov/30943403/) | [GSE115020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115020) | Download the GSE115020_RAW.tar file
| iBMEC Tissue Chips (Both flow and co-culture datasets) | [Vatine et al. Cell Stem Cell 2019](https://doi.org/10.1016/j.stem.2019.05.011) | [31173718](https://pubmed.ncbi.nlm.nih.gov/31173718/) | [GSE129290](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129290) | [Download the NCBI-generated Series raw counts matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129290)
| Breast Cancer Cell Dormancy and Type III Collagen | [Di Martino et al. Nature Cancer 2022](https://doi.org/10.1038/s43018-021-00291-9) | [35121989](https://pubmed.ncbi.nlm.nih.gov/35121989/) | [GSE182890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182890) | [Download the NCBI-generated Series raw counts matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE182890)
| Mechanical Conditioning of Breast Cancer Cells | [Watson et al. Cell Reports 2021](https://doi.org/10.1016/j.celrep.2021.109293) | [34192535](https://pubmed.ncbi.nlm.nih.gov/34192535/) | [GSE127887](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127887) |[Download the NCBI-generated Series raw counts matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE127887)
| SOX10-deficient melanoma cells | [Capparelli et al. Nature Communications 2022](https://doi.org/10.1038/s41467-022-28801-y) | [35296667](https://pubmed.ncbi.nlm.nih.gov/35296667/) | [GSE180568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180568) |[Download the NCBI-generated Series raw counts matrix](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE180568) 