import os
import pandas as pd
import numpy as np
from pybiomart import Dataset

CHROMOSOMES = [str(x) for x in np.arange(1,23)] + ['X','Y','MT']

def retrieve_cancer_fusion_data(fusion_db,cancertype):
    '''
    Retrieves cancer patient fusion gene data for patients with a specified cancer type.

    Parameters:
    -----------
        - fusion_db: pandas DataFrame
            - Dataframe containing all potential patients of interest
        - cancertype: list
            - list of cancer types to be extracted

    Returns:
    --------
        - cancer_fusion_data: pandas DataFrame
            - A filtered copy of the inputted dataframe.
    
    '''
    cancer_fusion_data = fusion_db.copy()
    cancer_fusion_data = cancer_fusion_data[cancer_fusion_data['Highly_Reliable_Seq'] == 'Seq+']
    cancer_fusion_data = cancer_fusion_data[cancer_fusion_data['Cancertype'].isin(cancertype)]
    cancer_fusion_data.drop_duplicates(inplace=True)

    return cancer_fusion_data

def extract_fusion_ensembl(cancer_fusion_data, gene_ID_conversion):
    '''
    Retrieves the ensembl IDs of genes for each gene fusion of interest.

    Parameters:
    -----------
        - cancer_fusion_data: pandas DataFrame
            - Dataframe containing the fusion gene names as retrieved from ChimerDB 
        - gene_ID_conversion: pandas DataFrame
            - Dataframe containing the gene name, synonym and Ensembl IDs

    Returns:
    --------
        -ensembl_dict: dict
            - dict containing the conversion between the gene name in ChimerDB and its Ensembl ID
    '''
    
    # Combining the genes of interest
    fusion_genes_OI = set(cancer_fusion_data.H_gene).union(cancer_fusion_data.T_gene)
    
    conv_OI = gene_ID_conversion[gene_ID_conversion['Gene name'].isin(fusion_genes_OI)].filter(['Gene stable ID','Gene name']).drop_duplicates()
    temp_conv = dict(zip(conv_OI['Gene name'],conv_OI['Gene stable ID']))

    # Now let's find how many did not get converted
    missing_genes = fusion_genes_OI.difference(temp_conv.keys())
    missing_genes = [x.upper() for x in missing_genes if isinstance(x,str)]
    alt_conv = gene_ID_conversion[gene_ID_conversion['Gene Synonym'].isin(missing_genes)].filter(['Gene stable ID','Gene name', 'Gene Synonym']).drop_duplicates()
    mis_temp = dict(zip(alt_conv['Gene Synonym'],alt_conv['Gene stable ID']))
    #unmappable = set(missing_genes).difference(mis_temp.keys())
    ensembl_id_conv = {**temp_conv, **mis_temp}

    return ensembl_id_conv

def verify_fusion_gene(fusion_info:dict, ensembl_conversion:dict, complete_gene_info:pd.DataFrame):
    '''
    This takes individual genes within the gene fusion and verifies that they can be converted into the correct ensembl ID.
    '''

    fusion_gene = fusion_info['Gene']
    chromloc = fusion_info['Chromosome']
    position = fusion_info['Position']

    # First check if the gene is in the baseline ensembl conversion dictionary
    check_loc = False
    upper_flag = False
    if fusion_gene in ensembl_conversion.keys():
        check_loc = True
    else:
        # Check if it needs to be uppercased for a proper check
        if fusion_gene.upper() in ensembl_conversion.keys():
            check_loc = True
            fusion_gene = fusion_gene.upper()
            fusion_info['Gene'] = fusion_gene
            upper_flag = True
    
    # Now double check that the fusion gene chromosome is an acceptable chromosome and breakpoint is within the bounds of the gene.
    verified = False
    ensembl_output = None
    if check_loc:
        ensembl_output = ensembl_conversion[fusion_gene]
        acceptable_chrom = chromloc in CHROMOSOMES
        if acceptable_chrom:
            gene_index = complete_gene_info[complete_gene_info['Gene stable ID'] == ensembl_conversion[fusion_gene]].index
            gene_information = complete_gene_info.iloc[gene_index[0]]
            correct_chrom = gene_information['Chromosome/scaffold name'] == chromloc
            within_bounds = (position >= gene_information['Gene start (bp)']-500) & (position <= gene_information['Gene end (bp)']+500)
            if correct_chrom and within_bounds:
                verified = True
                
            elif not correct_chrom:
                # If the chromosome is wrong need to make sure we get the correct gene instead
                potential_info = complete_gene_info[complete_gene_info['Gene name'] == fusion_gene].filter(['Gene stable ID', 'Gene name','Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)']).drop_duplicates()
                if len(potential_info) > 1:
                    overlap = potential_info[potential_info['Chromosome/scaffold name'] == chromloc]
                    new_ensembl = overlap['Gene stable ID'].values[0]
                else:
                    potential_info = complete_gene_info[complete_gene_info['Gene Synonym'] == fusion_gene].filter(['Gene stable ID', 'Gene name','Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)']).drop_duplicates()
                    if len(potential_info) > 1:
                        overlap = potential_info[potential_info['Chromosome/scaffold name'] == chromloc]
                        new_ensembl = overlap['Gene stable ID'].values[0]
                    else:
                        new_ensembl = None
                if new_ensembl is not None:
                    ensembl_output = new_ensembl
                    verified = True
            
            else:
                pass

    return verified, upper_flag,ensembl_output

def retrieve_ensembl_exon_data(filename, dataset:Dataset = None):
    '''
    Retrieves the ensembl ID information. It will query Ensembl multiple times if necessary due to a large number of IDs whose information needs to be retrieved.

    Parameters:
    -----------
        - dataset: pybiomart Dataset object
            - pybiomart dataset that connects to Ensembl
        - fusion_ensembl: list
            - list of ensembl IDs whose information is to be retrieved

    Returns:
    --------
        - gene_pos_info: pandas DataFrame
            - dataframe containing all exon position, coding sequence, strand, and exon length.
    '''

    # First check if the file exists that should contain all the exon information.
    if os.path.exists(filename):
        exon_information = pd.read_csv(filename)
    else:
        # If the file does not exist check if a pybiomart Dataset has been provided otherwise make a new one.
        if dataset == None:
            dataset = Dataset(host = 'http://useast.ensembl.org', name='hsapiens_gene_ensembl')
    
        exon_information = dataset.query(attributes=['ensembl_gene_id','exon_chrom_start','exon_chrom_end','cds_start','strand','transcript_length','cds_end'],
                                      filters = {'chromosome_name':CHROMOSOMES, 'transcript_is_canonical':True})
        exon_information.to_csv(filename, index=False)

    return exon_information