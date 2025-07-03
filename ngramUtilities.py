
from collections import Counter
import pandas as pd
import numpy as np
import networkx as nx


def reduce_ngram_dict(ngram_dict, num_arch =1):
    """
    Filters the n-gram dictionary to n-grams that meet the specified threshold.
    
    Parameters:
    -----------
        - ngram_dict: dict 
            key-value pairs of domain architecture ngram and a list of protein identifiers containing the ngram
        - num_arch: int
            threshold for number of instances of specific architecture to retain

    Returns:
    --------
        - ngram_dict_sub: dict 
            a reduced version of the inputted dictionary
    """
    
    ngram_dict_sub = {ngram:v for ngram,v in ngram_dict.items() if len(v) >= num_arch}
    
    return ngram_dict_sub

    
def return_ngram_adjacency(ngram_dict, interpro_conversion_dict, readable = 1):
    """
    Creates an adjacency matrix of all unique ngrams. If one is a subset of another, it will 
    put an entry in the matrix, which is the number of family members. This is sorted by length of the ngram 
    (i.e. 1-gram, 2-gram, 3-gram, etc.)

    Parameters:
    -----------
        - ngram_dict: dict
            n-grams produced from a dataframe and processed to remove redundant architectures
        - interpro_conversion_dict: dict
            key-value pairs of interpro names and IDs
        - readable: bool (Optional)
            flag to determine whether to pass a human readable adjacency matrix

    Returns:
    --------
        - ngram_adj_df: Pandas dataframe
            The adjacency matrix where values correspond to the number of instances the (sub)member occurs.
    """
    
    #first setup the unique and sorted architectures
    ngram_numdom_dict = {ngram:len(prot_list) for ngram,prot_list in ngram_dict.items()}
    
    # Setting up the sorted order of ngrams
    ngram_sorted = []
    tuples_sorted = sorted(ngram_numdom_dict.items(), key=lambda x:x[1], reverse=True)
    for tuple in tuples_sorted:
        ngram_sorted.append(tuple[0])
        
    # Setting up the adjacency matrix in a pandas DataFrame (may want to consider changing the base data structure later)
    ngram_adj_df = pd.DataFrame(index=ngram_sorted, columns=ngram_sorted,dtype=int)

    # To speed things up a little using the length of the ngram to determine for loop order
    ngrams = list(ngram_dict.keys())
    ngrams = sorted(ngrams, key = lambda x:len(x.split('|')), reverse=False)
    ngram_check = set(ngrams)

    for outer_ngram in ngrams:
        ngram_adj_df.loc[outer_ngram, outer_ngram] = ngram_numdom_dict[outer_ngram]

        # Don't check itself and remove from future loops
        ngram_check = ngram_check.difference([outer_ngram])
        for inner_ngram in ngram_check:
            
            if outer_ngram in inner_ngram:
                ngram_adj_df.loc[outer_ngram, inner_ngram] = ngram_numdom_dict[inner_ngram]

    # Convert the Interpro ID-based naming to the domain name to make it human readable in exported files
    if readable:
        ngram_sorted_names = []
        for arch in ngram_sorted:
            dom_arch = arch.split('|')
            dom_arch_convert = []
            for dom_id in dom_arch:
                dom_arch_convert.append(interpro_conversion_dict[dom_id])
            dom_arch_convert = '|'.join(dom_arch_convert)
            ngram_sorted_names.append(dom_arch_convert)
        
        # Creating dictionary for column and index conversion
        ngram_id_dict = dict(zip(ngram_sorted, ngram_sorted_names))
        ngram_adj_df.rename(columns=ngram_id_dict, index=ngram_id_dict, inplace=True)
    return ngram_adj_df

def concatenate_ngrams(ngram_dict, max_length = None):
    """ 
    Merge n-grams together in the dictionary if their architecture is a subarchitecture of another and they are fully overlapping
    (i.e. this will remove the architectures that are redundant as they only capture a single family of proteins that increased
    domain complexity, but for which there is no other overlap)

    Parameters:
    -----------
        - ngram_dict: dict 
            N-grams produced from a dataframe and processed to remove redundant architectures
        - max_length: int, optional
            Threshold of the maximum length of n-grams that are extracted and concatenated

    Returns:
    --------
        - ngram_dict_sub_concat: dict 
            Contains the nonredundant n-grams and their protein members
        - vals_to_delete_total: list
            Key values that were removed from the original dict
    """
    
    # Go through an initial time to get the first obvious ones
    ngram_dict_sub_concat, vals_to_delete = concatenate_ngrams_inner(ngram_dict.copy(), max_length)
    vals_to_delete_total = vals_to_delete

    # Repeat until no more n-grams can be collapsed.
    while vals_to_delete:
        ngram_dict_sub_concat, vals_to_delete = concatenate_ngrams_inner(ngram_dict_sub_concat.copy(), max_length)
        vals_to_delete_total+=vals_to_delete
        
    return ngram_dict_sub_concat, vals_to_delete_total

def concatenate_ngrams_inner(ngram_dict, max_length = None):
    """ 
    Inner function of the concatenating function that merges n-grams together in the dictionary if their architecture is a subarchitecture of another and they are fully overlapping.

    Parameters:
    -----------
        - ngram_dict: dict 
            N-grams produced from a dataframe and processed to remove redundant architectures
        - max_length: int (Optional)
            Threshold of the maximum length of n-grams that are extracted and concatenated

    Returns:
    --------
        - ngram_dict_sub_concat: dict 
            Contains the nonredundant n-grams and their protein members
        - ngrams_to_remove: list
            Key values that were removed from the original dict
    """
    if max_length == None:
        max_length = 303 # Within the human proteome this is the maximum length
    ngrams = list(ngram_dict.keys())
    ngrams_to_remove = []

    # Start off with longer n-grams so that they are run through quickly to then be concatenated earlier and reduce loop iterations. 
    ngrams = sorted(ngrams, key = lambda x:len(x.split('|')), reverse=True)

    ngrams_to_remove = set()
    for ngram in ngrams:
        if ngram in ngrams_to_remove: # In case an n-gram gets removed early just skip the process (In testing this has not been encountered.)
            continue
        else:
            # Inner loop only checks against ngrams that haven't been removed already
            ngrams_check = set(ngrams).difference(ngrams_to_remove) 
            for node in ngrams_check:
                if (ngram in node) and (ngram!=node): #it's a string subset
                    if set(ngram_dict[ngram]).difference(ngram_dict[node]) == set(): #and has matching protein lists
                        ngram_length = len(ngram.split('|'))
                        node_length = len(node.split('|'))
                        if ngram_length <= max_length and node_length <= max_length:
                            # Keep only the longer n-gram that encodes more information
                            ngrams_to_remove.add(ngram)              
                        else:
                            ngrams_to_remove.add(node)
    
    ngrams_to_remove = list(ngrams_to_remove)
    for ngram in ngrams_to_remove:
        del ngram_dict[ngram]
    return ngram_dict, ngrams_to_remove     

def return_all_n_grams_from_key(domain_architecture_str, domains, max_ngram = None):
    """
    Given a domain architecture find all n-grams that are possible from starting from the domain(s) of interest.

    Parameters:
    -----------
        - domain_architecture_str : str
            The domain architecture to be assessed
        - domains : str
            The domain(s) of interest that n-grams must contain
        - max_ngram: int (Optional)
            the maximum length an n-gram can be extracted. If no maximum is set then will determine it based on the length of the full architecture.

    Returns:
    --------
        - ngram_list : list
            all possible n-grams extracted from the complete architecture
    """
    
    # Ensure it is a list
    if isinstance(domains, str):
        domains = [domains]

    domain_list = domain_architecture_str.split('|')
    arch_length = len(domain_list)

    if max_ngram == None:
        max_length = arch_length
    else:
        max_length = max_ngram

    ngram_list = []
    
    for domain_of_interest in domains:
        found_ind = []
        for i in range(0, arch_length):
            if domain_list[i] == domain_of_interest:
                found_ind.append(i)
        
        #start by adding the individual times a domain of interest was found
    
        for target_ind in found_ind:
            i = target_ind-1 
            n_term_max = target_ind - max_length + 1 
            #Towards the N-term side
            while i >= np.max([0, n_term_max]):
                ngram = '|'.join(domain_list[i:target_ind+1])
                ngram_list.append(ngram)
                i-=1
                
            # Towards C-term Side
            i=target_ind+2
            
            while i <= np.min([arch_length, max_length]):
                
                ngram = '|'.join(domain_list[target_ind:i])
                ngram_list.append(ngram)
                i+=1

            # Expand towards both sides
            i = target_ind-1
            j = target_ind+2
            while i >= 0:
                while j <= np.min([arch_length, max_length]):
                    ngram = '|'.join(domain_list[i:j])
                    ngram_list.append(ngram)
                    j+=1
                i-=1

    #after uniquifying the ngram list, in case there are doubles, add back in the unique domains if more than one. 
        ngram_list = list(set(ngram_list))
        
        for ind in found_ind:
            ngram_list.append(domain_list[ind])
    return ngram_list

def get_ngrams_from_df(df, interpro_id, arch_num=1, max_ngram = None):
    """
    Given a dataframe containing the Interpro and Uniprot reference values return a dict of n-grams and their relative number limited to those with a specified number of instances.

    Parameters:
    -----------
        df: Pandas dataframe
            Contains all architectures and UniProt and InterPro IDs
        interpro_id: str
            The InterPro ID that seeds the n-grams. Note this has to be the ID not the name to prevent mismatching due to substrings within domain names
        arch_num: int (Optional)
            The minimum number of instances an architecture occurs to be kept in the final dict
        max_ngram: int (Optional)
            The maximum length an n-gram can be
    
    Returns:
    --------
        ngram_dict_sub: dict 
            Key value pairs of the n-grams and the UniProt IDs associated with them but limited to those with the minimum number of instances
        ngram_dict: dict
            Same as ngram_dict_sub but all extracted ngrams
    """

    ngram_dict = {}
    for _, row in df.iterrows():
        domain_arch_str = row['Interpro Domain Architecture IDs']
        ngram_list = return_all_n_grams_from_key(domain_arch_str, interpro_id, max_ngram)
        uniprot_id = row['UniProt ID']
        for ngram in ngram_list:
            if ngram not in ngram_dict:
                
                ngram_dict[ngram] = []
            ngram_dict[ngram].append(uniprot_id)
    ngram_dict_sub = reduce_ngram_dict(ngram_dict,arch_num)

    return ngram_dict_sub, ngram_dict


def get_all_ngrams_from_df(df, max_ngram = None):
    """
    Given a dataframe containing the Interpro architectures and IDs along with Uniprot IDs extracts all possible n-grams and the adjacency matrix.

    Parameters:
    -----------
        df: Pandas dataframe
            Reference dataframe containing all values of interest
        max_ngram: int (Optional)
            maximum length of an n-gram that can be found
    
    Returns:
    --------
        ngram_adj_df: Pandas dataframe
            The adjacency matrix where entries represent the number of instances a (sub)member occurs
        df: Pandas DataFrame
            A modified reference dataframe containing Interpro ID architectures
        ngram_dict_convert: dict
            Dictionary that contains the UniProt IDs for each n-gram
        removed_ngrams: list
            The n-grams that were collapsed into other n-grams due to complete overlap
    """

    # Find all domains within the dataframe
    domains = [arch.split('|') for arch in df['Interpro Domain Architecture IDs']]
    domain_list = Counter(sum(domains,[])).keys()

    ngram_adj_df, df, ngram_dict_convert, removed_ngrams, _ = full_ngram_analysis(df,domain_list,min_arch=1, max_ngram=max_ngram)
    return ngram_adj_df, df, ngram_dict_convert, removed_ngrams

def add_Interpro_ID_architecture(orig_df):
    """
    Add in a new column to a reference dataframe containing the domain architecture as the Interpro IDs rather than domain names to prevent substring mismatching due to similarities in domain names.

    Inputs:
        - df: Pandas dataframe
            Reference data frame containing the starting information
    
    Outputs:
        - df: Pandas dataframe
            Amended dataframe containing the new architecture string as the new column
        - interpro_conversion: dict 
            Key-value pairs of the Interpro IDs to the Interpro domain name. Note: If the ID architecture column exists already this will be an empty dictionary.
    """
    if 'Interpro Domain Architecture IDs' not in orig_df.columns:
        df = orig_df.copy()
        interpro_conversion = {}
        id_arch_list = []
        
        for _, row in df.iterrows():

            # Extract the string containing the Interpro Domain information
            interproStr = row['Interpro Domains']
            if interproStr != '':
                interproStr = interproStr.split(';')
                interpro_Domain_IDs = []
                
                for domain in range(0, len(interproStr)):
                    
                    # Extract the domain name and ID in the order they appear in the string
                    splitInterpro = interproStr[domain].split(':')
                    domain_name = splitInterpro[0]
                    domain_interpro_ID = splitInterpro[1]
                    interpro_Domain_IDs.append(domain_interpro_ID)
                    # Check if the domain is present in the ID to name conversion dictionary and add otherwise
                    if domain_interpro_ID not in interpro_conversion.keys():
                        interpro_conversion[domain_interpro_ID] = domain_name
            else:
                interpro_Domain_IDs = ''
            # Rejoin the domain architecture into a single string
            interpro_Domain_IDs_str = '|'.join(interpro_Domain_IDs)
            id_arch_list.append(interpro_Domain_IDs_str)

        # Add new column to the dataframe
        df['Interpro Domain Architecture IDs'] = id_arch_list
    else:
        df = orig_df.copy()
        interpro_conversion = {}
    return df, interpro_conversion

def generate_interpro_conversion(df):
    '''
    Generates the dict that contains the conversion from interpro ids to a human legible name.

    Parameters:
    -----------
        - df: pandas DataFrame
            reference dataframe that contains all the InterPro domain information
    
    Returns:
    --------
        - interpro_conversion: dict
            dict containing the InterPro IDs as keys with the human legible name as values.

    '''

    interpro_conversion = {}
    for _, row in df.iterrows():

            # Extract the string containing the Interpro Domain information
            interproStr = row['Interpro Domains']
            if interproStr != '':
                interproStr = interproStr.split(';')
                
                for domain in range(0, len(interproStr)):
                    
                    # Extract the domain name and ID in the order they appear in the string
                    splitInterpro = interproStr[domain].split(':')
                    domain_name = splitInterpro[0]
                    domain_interpro_ID = splitInterpro[1]
                    # Check if the domain is present in the ID to name conversion dictionary and add otherwise
                    if domain_interpro_ID not in interpro_conversion.keys():
                        interpro_conversion[domain_interpro_ID] = domain_name

    return interpro_conversion

def full_ngram_analysis(df, Interpro_ID,min_arch = 1, max_node_len = None, max_ngram = None, readable_flag = False, concat_flag = True, verbose = True):
    """
    Conduct n-gram analysis on a protein domain dataframe whose n-grams must include specific seeding domains.

    Parameters
    ----------
        df: Pandas DataFrame
            dataframe containing the reference file(s) data
        Interpro_ID: str | list
            Interpro IDs that seed the n-gram extraction
        max_ngram: int (Optional, Recommended)
            Maximum length of n-grams to be extracted. Can be omitted if the max_node_length is provided.
        min_arch: int (Optional)
            Minimum occurrences for an architecture to be retained (Default = 2)
        max_node_length: int (Optional)
            Maximum length of n-grams that can be concatenated into
        readable_flag: bool (Optional)
            flag as to whether the adjacency matrix should be human legible    
        concat_flag: bool (Optional)
            Flag as whether the n-gram collapsing is to be conducted.
        verbose: bool (Optional)
            Flag as to whether to print statements that indicate process of the analysis.

    Returns
    -------
        ngram_adj_df: Pandas DataFrame
            The adjacency matrix of n-grams
        df: Pandas DataFrame
            A modified reference dataframe containing Interpro ID architectures
        ngram_dict_convert: dict
            Dictionary that contains the UniProt IDs for each n-gram
        removed_ngrams: dict
            Dictionary containing the n-grams that were concatenated into other n-grams due to complete overlap
        interpro_conversion_dict: dict
            Dictionary that contains all the InterPro IDs and the corresponding names for conversion as needed

    """

    # Checking the parameters

    # If parameters were both provided comparing the two to ensure there are no weird issues
    if max_ngram != None and max_node_len != None:
        if max_ngram != max_node_len:
            if max_ngram < max_node_len:
                max_node_len = max_ngram
    
    # If only the max_ngram is provided will set the maximum node length to be equivalent
    elif max_ngram != None and max_node_len == None:
        max_node_len = max_ngram
    
    # If neither are provided then erroring out for the time being.
    else:
        raise ValueError('Please provide a maximum length of n-gram to be extracted.')

    # Getting n-gram information
    if verbose:
        print('Starting to fetch n-grams.')
    interpro_conversion_dict = generate_interpro_conversion(df)
    ngram_dict_sub, ngram_dict_comp = get_ngrams_from_df(df,Interpro_ID, min_arch, max_ngram)
    if verbose:
        print('Finished getting all n-grams')

    # Starting to collapse the n-grams to non-redundant n-grams
    if concat_flag:
        ngram_dict_concat, vals_to_delete = concatenate_ngrams(ngram_dict_sub,max_node_len)
    else:
        ngram_dict_concat = ngram_dict_comp
        vals_to_delete = []
    
    # Generating the adjacency matrix
    if verbose:
        print('Starting to generate adjacency')
    ngram_adj_df = return_ngram_adjacency(ngram_dict_concat, interpro_conversion_dict, readable_flag)
    if verbose:
        print('Finished building adjacency.')

    # Returning a human legible version of the results
    if readable_flag:
        # Pass a human-readable dictionary of the identifying
        ngram_dict_convert = {}
        for ngram in ngram_dict_concat.keys():
            gram = ngram.split('|')
            gram_convert = []
            for k in gram:
                k_con = interpro_conversion_dict[k]
                gram_convert.append(k_con)
            gram_convert_str = '|'.join(gram_convert)
            ngram_dict_convert[gram_convert_str] = ngram_dict_concat[ngram]

        # Passing a human-readable dictionary of the removed ngrams
        removed_ngrams = []
        for ngram in vals_to_delete:
            gram = ngram.split('|')
            gram_convert = []
            for k in gram:
                k_con = interpro_conversion_dict[k]
                gram_convert.append(k_con)
            gram_convert_str = '|'.join(gram_convert)
            removed_ngrams.append(gram_convert_str)
    else:
        ngram_dict_convert = ngram_dict_concat
        removed_ngrams = vals_to_delete

    ngram_adj_df.fillna(0,inplace=True)
    return ngram_adj_df, df, ngram_dict_convert, removed_ngrams, interpro_conversion_dict

def import_reference_file(reference_file):
    """
    Read and process the reference file(s) containing the interpro and uniprot information. This will also provide a seperate dataframe that provides background information (e.g. names/counts of domains, n-gram length) of the imported data.

    Paremeters:
    ----------- 
        -reference_file: str | list
            Path(s) to the reference files

    Returns:
    --------
        - df: Pandas dataframe
            Processed dataframe of the imported data
    """

    # Checking if the reference_files are a list or not
    if isinstance(reference_file, str):
        df = pd.read_csv(reference_file)
    elif isinstance(reference_file, list):
        df = pd.DataFrame()
        for ref_file in reference_file:
            temp_df = pd.read_csv(ref_file)
            df = pd.concat([temp_df, df])
        df.drop_duplicates(inplace=True, ignore_index=True, subset='UniProt ID')   

    # Just a bit of cleanup of the strings and adding in the Interpro ID based architecture
    df.fillna('',inplace=True,)

    # Adding in the Interpro Domain Architecture IDs here
    df,_ = add_Interpro_ID_architecture(df)

    return df

def get_ngram_cutoff(ngram_counts, thres = 0.8):
    """
    Get a cutoff for the maximum length of an n-gram to use in downstream analysis. Values are computed by the relative fraction of the cumulative sum described by the n-gram length. 

    NOTE: Plan is to deprecate this function as it is no longer required for downstream analysis and is redundant with other functions that are available.

    Parameters:
    -----------
        - ngram_counts: dict 
            The instances an n-gram occurs
        - thres: int
            The relative fraction of proteins to be considered

    Returns:
    --------
        - n_cutoff: int
            N-gram length cutoff
    """
    for i in range(1, np.max(list(ngram_counts.keys()))):
        if i not in ngram_counts.keys():
            ngram_counts[i] = 0
    ngram_counts = dict(sorted(ngram_counts.items(), key=lambda x:x[0], reverse=False))

    ngram_cumsum = np.array(list(ngram_counts.items())).cumsum(axis=0)[:,1]
    ngram_cumsum_norm = ngram_cumsum/np.max(ngram_cumsum)
    n_cutoff = np.where(ngram_cumsum_norm > thres)[0].tolist()
    n_cutoff = n_cutoff[0]+1

    # Returning the normalized cumulative sum to provide information about where the cutoff came from.
    ngram_cumsum_norm = dict(zip(ngram_counts.keys(), ngram_cumsum_norm))

    return n_cutoff, ngram_cumsum_norm