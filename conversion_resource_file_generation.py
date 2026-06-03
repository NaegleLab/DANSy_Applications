import os
import pandas as pd
import numpy as np
import urllib.request
import uniprotIDMapping as unimapping
from pybiomart import Dataset

# This script downloads several pieces of data that speed up several parts of the analysis necessary for converting between gene names and the proteome.
data_directory = './'
uniprot_file_name = 'uniprot_id_gene_names.tsv'
ensembl_conv_file = 'ENSEMBL_Gene_Conversion.csv'
exon_filename = 'Gene_exon_information.csv'

# Chromosome names to help limiting any data downloaded from ENSEMBL via pybiomart
chroms = [str(x) for x in np.arange(1,23)] + ['X','Y','MT']

# A couple of helper functions
def unpack_uniprot_results(res):
    '''
    Unpacks the results that come from a ID mapping job submission to the UniProt API with the input query and the corresponding UniProt ID as individual outputs.

    Parameters
    ----------
    - res: dict
        The inner results dictionary that is returned when an ID mapping job is completed from UniProt

    Returns
    -------
    - query: str
        The original query input (e.g. the ENSEMBL ID)
    - uniprot_id: str
        The mapped UniProt ID for the input
    '''
    query = res['from']
    uniprot_id = res['to']['primaryAccession']
    return query, uniprot_id


def fill_missing_value(row, mapped_ids):
    '''
    This fills in missing values within a row of the pandas DataFrame returned from pybiomart that contains all the gene information based on mapped values returned from the UniProt ID mapping results. 

    Parameters
    ----------
        - row: pandas DataFrame row
            The row of the pandas DataFrame the method is to be applied to
        - mapped_ids: dict
            Dictionary of key-value pairs for ENSEMBL IDs to UniProt IDs. Should only contain missing conversions.

    Returns
    -------
        - new_val: str
            The UniProt ID that corresponds to the ENSEMBL ID of that row. Note: this will be the existing value if the ENSEMBL ID had a corrsponding UniProt ID already.
    '''
    idOI = row['Gene stable ID']
    if idOI in mapped_ids:
        new_val = mapped_ids[idOI]
    else:
        new_val = row['UniProtKB/Swiss-Prot ID']
    
    return new_val

def main():
    # First grabbing the UniProt IDs that are currently present within the human proteome. The url for the api come from UniProt after selecting only reviewed entries.
    uniprotfetch_exists = os.path.exists(''.join([data_directory, uniprot_file_name]))
    if not uniprotfetch_exists:
        url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&fields=accession%2Cid%2Cgene_primary&format=tsv&query=%28%28reviewed%3Atrue%29%29+AND+%28model_organism%3A9606%29'
        x = urllib.request.urlretrieve(url, filename = uniprot_file_name)

    # Now get the pybiomart dataset
    # pybiomart database for gene name conversions
    dataset = Dataset(host = 'http://useast.ensembl.org', name='hsapiens_gene_ensembl')

    gene_ID_conv = dataset.query(attributes=['ensembl_gene_id','external_gene_name','external_synonym','uniprotswissprot','transcript_is_canonical','chromosome_name','start_position','end_position'], 
                                filters={'transcript_is_canonical':True, 'chromosome_name':chroms})

    full_id_conv = gene_ID_conv.copy()
    full_id_conv.dropna(axis=0,subset=['Gene name','Gene Synonym'], inplace=True, how='all')

    # Unfortunately, pybiomart/ENSEMBL does not always map the corresponding gene ID to a UniProt ID despite UniProt mapping to the ENSEMBL. As a result, we have to double check the missing ones to ensure we have the full proteome.
    missing_uniprots = full_id_conv[full_id_conv['UniProtKB/Swiss-Prot ID'].isna()]['Gene stable ID'].drop_duplicates()

    # Need to submit a job to do the ID conversion with UniProt using a python toolset they provided in the API documentation (saved as uniprotIDMapping.py).
    job_id = unimapping.submit_id_mapping(
        from_db="Ensembl", to_db="UniProtKB-Swiss-Prot", ids=list(missing_uniprots)
    )

    running = True
    while running:
        link = unimapping.get_id_mapping_results_link(job_id)
        results = unimapping.get_id_mapping_results_search(link)
        running = unimapping.check_id_mapping_results_ready(job_id)

    success_mapped = results['results']
    id_mapping = {}
    for i in success_mapped:
        k,v = unpack_uniprot_results(i)
        id_mapping[k] = v

    full_human_uniprot = pd.read_csv('uniprot_id_gene_names.tsv', sep='\t')
    missing_genes = full_id_conv[full_id_conv['Gene stable ID'].isin(results['failedIds'])]['Gene name'].unique()

    mapped_missing = full_human_uniprot[full_human_uniprot['Gene Names (primary)'].isin(missing_genes)]
    # Now mapping them using the gene names and their corresponding ENSEMBL IDs
    temp = full_id_conv[(full_id_conv['Gene name'].isin(mapped_missing['Gene Names (primary)'])) & (full_id_conv['Gene stable ID'].isin(results['failedIds']))].filter(['Gene stable ID','Gene name']).drop_duplicates()
    mapped_missing = mapped_missing.merge(temp, right_on = 'Gene name',left_on = 'Gene Names (primary)')
    additional_ids = dict(zip(mapped_missing['Gene stable ID'],mapped_missing['Entry']))
    id_mapping = {**id_mapping, **additional_ids}


    full_id_conv['UniProtKB/Swiss-Prot ID'] = full_id_conv.apply(fill_missing_value, args=[id_mapping], axis = 1)

    full_id_conv.to_csv(ensembl_conv_file, index=False)

    # Now this is for the Exon information
    full_exon_information = dataset.query(attributes=['ensembl_gene_id','exon_chrom_start','exon_chrom_end','cds_start','strand','genomic_coding_start','genomic_coding_end','cds_end'],filters = {'transcript_is_canonical':True, 'chromosome_name':chroms})
    full_exon_information.to_csv(exon_filename,index=False)

if __name__ == '__main__':
    main()