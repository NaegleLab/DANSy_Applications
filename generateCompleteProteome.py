from os import listdir
import os.path as path
import ngramUtilities
import pandas as pd
import numpy as np
import json

def generate_complete_adjacency(generate_files = False, readable_flag = False, filenames = None, ref_version = None, ref_dir = None):
    '''
    Base function that conducts the full n-gram analysis pipeline on the complete canonical human proteome and generates the adjacency matrix and the n-grams that were subsumed due to redundant information.

    Parameters:
    -----------
        - generate_files: bool
            Flag to determine if csv files should be generated in the project directory
        - readable_flag: bool
            Flag to determine if the adjacency matrix should use the InterPro IDs or the domain names
        - filename: list
            Strings for the name of files to be generated.
        - ref_version: str
            Version of reference files to use indicated by the suffix of the file name
        -ref_dir: str
            Directory for the csv file
    
    Returns:
    --------
        - full_adj: pandas DataFrame
            dataframe that contains the adjacency matrix of the complete proteome. Self-loops have been removed
        - removed_ngrams: list
            List of the n-grams that had redundant information that were subsumed during analysis.
    '''

    if filenames == None:
        adj_file = 'Current_Complete_Proteome_Adjacency.csv'
        ngram_file = 'Current_Complete_Proteome_Removed_Ngrams.csv'
    else:
        adj_file = filenames[0]
        ngram_file = filenames[1]

    ref_df, interpro_dict = import_proteome_files(ref_version, ref_dir)
    
    all_domains = [x for x in interpro_dict.keys()]
    adj_df, _, _, removed_ngrams, _ = ngramUtilities.full_ngram_analysis(ref_df, all_domains, min_arch=1, readable_flag=readable_flag, max_ngram=303,max_node_len=303)

    full_adj = adj_df.copy()
    for node in full_adj.columns:
        full_adj.loc[node,node]=0
   
    if generate_files:
        pd.DataFrame.to_csv(full_adj,adj_file)
        np.savetxt(ngram_file, removed_ngrams, delimiter = ',', fmt = '%s')
    return full_adj, removed_ngrams

def import_proteome_files(reference_file_suffix = None, ref_file_dir = None):
    '''
    Imports the files that are used for the generation of the reference dataframe of the complete canonical proteome.

    Note: Need to adjust this so it looks in only one folder from here on out.

    Parameters:
    -----------
        - reference_file_version: str
          String of the suffix of the reference files to be used  

    Returns:
    --------
        - ref_df: pandas DataFrame
            Dataframe containing the InterPro, UniProt, and PDB information of individual proteins as retrieved via CoDIAC
        - interpro_dict: dict
            dictionary containing the InterPro IDs and domain names for conversion purposes
    
    '''
    all_refs = []
    if ref_file_dir == None:
        d = 'data/Current_Human_Proteome/'
    else:
        d = ref_file_dir
    
    if reference_file_suffix == None:
        reference_file_suffix = 'Jan7_2025.csv' # The current version of the reference files
    
    ref_files = listdir(d)
    for fileName in ref_files:
        if fileName.endswith(reference_file_suffix):
            all_refs.append(d+fileName)

    ref_df = ngramUtilities.import_reference_file(all_refs)
    ref_df, interpro_dict = ngramUtilities.add_Interpro_ID_architecture(ref_df)

    # If the full proteome reference files already had Interpro IDs this can lead to an empty dictionary
    if len(interpro_dict) == 0:
        interpro_dict = ngramUtilities.generate_interpro_conversion(ref_df)

    return ref_df, interpro_dict

def import_adjacency(adj_file = None, adj_dir = 'data/'):
    '''
    Imports the adjacency matrix from either a designated file or by generating a new file. 

    Parameters:
    -----------
        adj_file: str
            File name of the adjacency file
        adj_dir: str (optional)
            Directory Name of where the adjacency file is located. By default it will place files in the data folder.
    
    Returns:
    --------
        full_adj: pandas DataFrame
            Full adjacency matrix of the n-gram analysis of the human proteome.

    '''
    if adj_file == None:
        adj_file  = 'Current_Complete_Proteome_Adjacency.csv'
    
    full_path = adj_dir+adj_file
    if path.exists(full_path):
        full_adj = pd.read_csv(full_path, header=0,index_col=0)
    else:
        print('Generating a new Adjacency Matrix but not a csv file.')
        full_adj, _ = generate_complete_adjacency()

    return full_adj

def import_removed_ngrams(rm_ngrams_file = None, rm_dir = '\data'):
    '''
    Imports the list of n-grams that were removed during n-gram analysis of the complete proteome from either the file or via the n-gram analysis.

    Parameters:
    -----------
        rm_ngrams_file: str (optional)
            File name if not provided will generate a File with the title: Current_Complete_Proteoeme_Removed_Ngrams.csv
        rm_dir: str (optional)
            Directory Name of where the file is located. By default it will place files in the data folder.
    
    Returns:
    --------
        removed_ngrams: list
            List of subsumed n-grams with redundant information

    '''
    if rm_ngrams_file == None:
        rm_ngrams_file  = 'Current_Complete_Proteome_Removed_Ngrams.csv'

    full_path = rm_dir + rm_ngrams_file
    if path.exists(full_path):
        temp = pd.read_csv(full_path, header = None)
        removed_ngrams = temp[0].tolist()

    else:
        _, removed_ngrams = generate_complete_adjacency()

    return removed_ngrams

def import_adjacency_from_json(adj_json = "Current_Adjacency.json", adj_dir = 'data/'):
    '''
    Imports the adjacency matrix from the designated json file. The JSON file should have keys correspond to the index of the final adjacency dataframe.
    
    Parameters:
    -----------
        adj_json: str
            File name of the adjacency json file
        adj_dir: str (optional)
            Directory Name of where the adjacency file is located. By default it will look for files in the data folder.
    
    Returns:
    --------
        adj: pandas DataFrame
            Full adjacency matrix of the n-gram analysis of the human proteome.
    '''

    full_path = adj_dir + adj_json
    if path.exists(full_path):
        js_adj = pd.read_json(full_path, orient = 'index')
        # The JSON file helps reduce the memory footprint of the adjacency matrix by omitting 0s so have to fill them in. Also ensuring all inputted data is an integer.
        js_adj.fillna(0, inplace=True)
        adj = js_adj.astype('int64')
    else:
        raise FileExistsError('The JSON Adjacency File does not exist.')
    
    return adj


def generate_adjacency_json(adj_file = None, adj_json = "Current_Adjacency.json", adj_dir = 'data/'):
    '''
    Imports the adjacency matrix from the designated json file. The JSON file should have keys correspond to the index of the final adjacency dataframe.
    
    Parameters:
    -----------
        adj_file: str
            File name for the csv file corresponding to the adjacency matrix.
        adj_json: str
            File name of the adjacency json file
        adj_dir: str (optional)
            Directory Name of where the adjacency file is located. By default it will look for files in the data folder.
    
    Returns:
    --------
        adj: pandas DataFrame
            Full adjacency matrix of the n-gram analysis of the human proteome.
    '''

    # Setting a default name if not provided.
    if adj_file == None:
        adj_file  = 'Current_Complete_Proteome_Adjacency.csv'

    base_adj = import_adjacency(adj_file=adj_file, adj_dir=adj_dir)
    adj_dict = base_adj.to_dict(orient='index')

    # Running through each entry and removing the 0s except the self-referencing ones to preserve adjacency matrix shape and reduce memory footprint.
    for ngram in adj_dict.keys():
    
        keys_2_rm = base_adj.index[base_adj.loc[ngram] == 0].tolist()
        for k in keys_2_rm:
            if k != ngram:
                del adj_dict[ngram][k]
    
    with open(adj_json,'w') as output:
        json.dump(adj_dict, output)
    
    return None