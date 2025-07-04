# DANSy_Applications


This is our analysis that applies the linguistic technique n-gram analysis with network theory to protein domain architectures, to represent the proteome as an abstracts the functional connections between proteins to describe either proteome-wide (base DANSy) or phenotype-specific changes from differential expression results (deDANSy). 

### DANSy Overview
![Overview of the general workflow](Figures/N%20gram%20network%20workflow.png)

### deDANSy Overview
![Overview of the deDANSy workflow](Figures/deDANSy%20Overview.png)

How to cite: Please cite our [bioRxiv paper](https://doi.org/10.1101/2024.12.04.626803), which contains further details on the methods and their applications. Here, we provide the code that produces the results and figures from the manuscript.

## Getting started

For this work, we recommend creating a local copy of this repository by inputting the following into your terminal:

    git clone https://github.com/NaegleLab/DANSy_Applications

Create a virtual environment containing all the dependencies for the analysis using the following code in a terminal.

    conda create env -f dansy_apps.yml

Activate the environment using `conda activate dansy_apps` for specific scripts or select the dansy_apps kernel for the Jupyter notebooks.

## Proteome Reference File

Several of the analysis provided here rely on reference files generated by [CoDIAC](https://github.com/NaegleLab/CoDIAC). We have provided the reference files for the analysis conducted in our manuscript in the [data folder](data/Current_Human_Proteome/). The copy of the reference file used for the analysis was generated on January 7th, 2025, and will be the default file used for analysis. 

If you wish to generate the most up to date reference file to use for analysis, you will to take the following steps. First download the SwissProt ID list from [Gencode](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.metadata.SwissProt.gz) and place in the main directory of your local copy of this repo. Then, go to the [whole_proteome_reference.py](scripts/whole_proteome_reference.py) file and change the reference file suffix variable to the current date. Finally, run the following code in a terminal to establish the environment that includes CoDIAC, which will query UniProt and InterPro for the domain architectures. (Note: This can take up to 2 hours after a fresh install as it will also establish a biomart sqlite database for use in other notebooks.)

    conda create env -f dansy_codiac.yml
    conda activate dansy_codiac
    python scripts/whole_proteome_reference.py
    conda deactivate

If you are going to use the new build, you will then have to tell each of the notebooks to use the current build of the reference file by changing the `generateCompleteProteome.import_reference_file()` command to `generateCompleteProteome.import_reference_file(new_reference_file_suffix)`.

## Additional Datasets

For more specific analysis, the following datasets are recommended to be downloaded from their source.

| Analysis | Dataset Source/Code |
| -------- | ------------------- |
| Fusion Gene Analysis|ChimerSeq Excel File from [ChimerDB](https://www.kobic.re.kr/chimerdb/download)|
| TCGA Clinical Data | [cBioPortal](https://www.cbioportal.org/)
| PTM Systems| Provided in the [Multispecies Reference Files](data/Current_Multispecies_Files/) and see the example [CoDIAC script](scripts/species_phospho_reference_files.py)|
| Differential Gene Expression Results | Provided in the [DEG_data folder](DEG_data) and see the [R Scripts README file](R_Scripts_for_DEG_Analysis/README.md) | 

## What is provided and how to run the code.

We have provided several Jupyter notebooks that serve as examples of either DANSy or deDANSy analysis. Below are short summaries of their applications and the results, which are discussed in our [manuscript](https://doi.org/10.1101/2024.12.04.626803 ).

For custom uses of DANSy or deDANSy, please visit our [DANSy repository](https://github.com/NaegleLab/DANSy), which provides more general use cases and information on how to get started with the DANSy package for new datasets.

### Complete human proteome analysis

Analyzing the domain architectures across the human proteome and broadly characterizing the resulting network and information encoded by different versions related to n-gram length. Here, we find limiting n-grams to 10-grams in length will recapitulate network characteristics of the complete proteome and that specific domains such as the protein kinase, zinc finger C2H2, and EGF-like domains are n-grams that frequently lie along the shortest path and are connected to the most other n-gram nodes in the network.

Associated notebooks:

- [human_model_comp.ipynb](human_model_comp.ipynb)
- [human_proteome.ipynb](human_proteome.ipynb)

### PTM System Analysis

Focused on reversible post-translational modification systems (e.g. phosphorylation, methylation, acetylation) that operate under a reader-writer-eraser paradigm. Characterizing broad properties of how individual components combine in domain architectures and identifying general grammatical rules where eraser domains do not require additional reader domains to modify their activity. Meanwhile, reader domains will frequently associate as "adverbs" to modify the activity of writer domains. Further, the catalytic domains that act as "verbs" must be on separate proteins.

Associated notebooks:

- [PTM_systems_analysis.ipynb](PTM_systems_analysis.ipynb)

### PTM System Evolution

Additional analysis on the phosphorylation systems that compares the network characteristics of domains associated with phosphotyrosine (pTyr) and phosphoserine/threonine (pSer/Thr) systems during the evolutionary period from yeast to humans. The pTyr system is evolutionarily younger and rapidly expanded during the transition to metazoans. Our analysis suggests that during this transition, species were sampling several configurations of the network before converging to a similar set of grammatical rules observed in other PTM systems.

Associated notebooks:

- [ngram_species_comparison.ipynb](ngram_species_comparison.ipynb)

### Fusion Gene Analysis

Here, we analyze how fusion genes may provide an avenue to explore new grammatical structures of domain architectures. Using network topology, we studied fusion genes reported in TCGA and find that the domain architectures of gene fusions still adhere to the same grammatical rules that were established for domain architectures within the complete proteome. We further explore kinase fusions and find that the diversity of domain combinations it is involved in both within the natural proteome and fusions reflect its involvement across biological functions by connecting distinct n-gram families.

Associated notebooks:

- [fusion_genes.ipynb](fusion_genes.ipynb)
- [fusion_gene_survival.ipynb](fusion_gene_survival.ipynb)

### Applications to Differential Gene Expression Analysis

In these notebooks, we provide an example case of how differential expression DANSy (deDANSy) describes the coordinated behavior of gene expression/protein interaction networks. We focus on the gene expression changes induced by SOX10-deficiency through either a CRISPR-mediated or a therapeutic-induced mechanism. We demonstrate these two perturbations, despite showing similar phenotypes, achieve these through different domain language subnetworks. The RNA-sequencing results come from [this study](https://rdcu.be/eujbz), which showed SOX10-deficient melanoma cells exhibit a dormant and invasive phenotype and sequencing results were retrieved from the NCBI Gene Omnibus repository with accession number [GSE180568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180568).

Additionally, we provide other example deDANSy results from other publicly available datasets. The differential gene expression results (csv files) are provided in the [DEG_data](DEG_data) folder, and the R scripts that use DESeq2 to obtain those results are provided in the [R_Scripts folder](R_Scripts_for_DEG_Analysis). In the Jupyter notebooks, we collect and plot the data for each dataset. We also provide additional scripts that describe considerations taken during the development of deDANSy these are in [deDANSy_development_considerations](deDANSy_development_considerations).

Through this work, we show the network separation and syntax enrichment analysis provided by deDANSy characterizes the molecular basis of cell phenotypes, and identifies novel information about the effectors which is distinct from current gene set enrichment analysis (GSEA) approaches.

Associated notebooks:

- [deDANSy_All_Example.ipynb](deDANSy_All_Example.ipynb)
- [sox10_deDANSy_Full_Analysis.ipynb](sox10_deDANSy_Full_Analysis.ipynb)

