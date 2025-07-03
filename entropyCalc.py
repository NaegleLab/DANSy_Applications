import pandas as pd
import numpy as np
from collections import Counter

def get_token_list(corpus_list, delim = '|'):
    """
    Extracts the individual tokens from a corpus given a delimiter seperating the tokens

    Inputs:
        - corpus_list: list containing all the words/structures containing the tokens of interest
        - delim: str of the delimiter that separates the tokens of interest

    Outputs:
        - token_list: dict of all the identified tokens and their counts
    
    """
    split_corpus = [ngram.split(delim) for ngram in corpus_list]
    token_list = Counter(sum(split_corpus,[]))

    return token_list

def calc_F0(corpus_list, delim = '|'):
    """
    Calculates the entropy based on the number of units of interest. This refers to the 0-gram entropy that does not account for the relative frequency of the smallest unit (e.g. letters, domains, etc)

    Inputs:
        - corpus_list: list containing all the words/structures containing the tokens of interest
        - delim: str of the delimiter that separates the tokens of interest

    Outputs:
        - f0: float of the F0 entropy based on the tokens of interest rounded to 3 decimals
    """
    
    token_list = get_token_list(corpus_list, delim)
    f0 = np.round(np.log2(len(token_list)), 3)

    return f0

def calc_F1(corpus_list, delim = '|'):
    """
    Calculates the 1-gram entropy which is based on the relative frequency of the individual tokens

    Inputs:
        - corpus_list: list containing all the words/structures containing the tokens of interest
        - delim: str of the delimiter that separates the tokens of interest

    Outputs:
        - f1: float of the F1 entropy based on the token frequency rounded to 3 decimals
    
    """
    token_list = get_token_list(corpus_list, delim)
    
    # Initializing a token frequency dict
    token_freq = dict.fromkeys(token_list, 0)
    tot_toks = sum(token_list.values())

    for token in token_freq:
        token_freq[token] = token_list[token]/tot_toks

    int_step = [-p*np.log2(p) for p in token_freq.values()]
    f1 = round(sum(int_step),3)

    return f1

def calc_FN_from_adjacency(adj_df, corpus_list,n, removed_ngrams, delim = '|', db_check = 0):
    """
    Calculates the bigram entropy from an adjacency matrix which is calculated according the conditional probability associated with each bigram given the starting token and the probability of each bigram within the corpus. These probabilites are estimated by the relative frequency of the bigrams (and tokens) in the corpus 
    
    The inputted adjacency matrix is not required to be limited to bigrams only as the function will limit the analysis to only bigrams. Missing tokens from the adjacency matrix will also be filled in as necessary from the corpus list that was used to generate the adjacency matrix.

    Inputs:
        - adj_df: dataframe that is the adjacency matrix of n-grams of interest within the corpus list
        - corpus_list: list containing all the words/structures composed of the tokens of interest
        - n: int of n-gram length of interest
        - removed_ngrams: dict containing the n-grams that were removed due to complete overlap with larger n-grams. 
        - delim: str of the delimiter that separates the tokens of interest
        - db_check: bool to check whether a cumulative set of frequencies should be determined or not

    Outputs:
        - fn: float of the FN entropy based on the token frequency rounded to 3 decimals
    
    """

    # The adjacency matrix is based on the count of larger architectures that contain it and not the number of instances each architecture occurs.
    arch_OI_list = adj_df.columns.tolist()
    archs_OI = set(arch_OI_list).union(removed_ngrams)

    # Ensuring that the unigrams are accounted for in the architectures of interest
    unigrams = get_token_list(corpus_list,delim)
    archs_OI = archs_OI.union(unigrams)
    # Ensuring only the n-grams of interest are passed on
    arch_2_rm = []
    for ngram in archs_OI:
        arch_len = len(ngram.split(delim))

        # If a cumulative n-gram check is desired
        if db_check:
            if arch_len > n:
                arch_2_rm.append(ngram)
        else:
            if arch_len != n:
                arch_2_rm.append(ngram)
    
    # Removing n-grams not associated with the length of interest
    archs_OI = archs_OI.difference(arch_2_rm)

    # Retreiving the ngram frequencies to use for entropy calcualtions.
    ngram_freq = build_freq_dict(archs_OI, corpus_list,delim)
    fn = 0
    for ngram, f_dets in ngram_freq.items():
        p = f_dets['Count Freq.']
        pbi = f_dets['Conditional Freq.']
        if p == 0:
            fn += 0
        else:
            fn += -p*np.log2(pbi)

    return fn

def get_ngram_instance_count(ngram, corpus_list, delim):
    """
    Count the number of instances an n-gram occurs. Instance in this case refers to the actual times the n-gram occurs even within larger structures not just the number of structures it occurs in.

    Inputs:
        - ngram: str of the n-gram of interest
        - corpus_list: list containing the corpus to be analyzed for the ngram
        - delim: str of the delimiter separating individual tokens

    Outputs:
        - ngram_count: int of the number of instances
    """

    arch_length = len(ngram.split(delim))
    ngram_count = 0
    for raw_arch in corpus_list:
        if ngram in raw_arch:
            split_arch = raw_arch.split(delim)
            tot_ngrams = len(split_arch)-arch_length+1
            for x in range(0, tot_ngrams):
                sub_check = delim.join(split_arch[x:x+arch_length])
                if ngram == sub_check:
                    ngram_count += 1        

    return ngram_count

def build_freq_dict_from_adj(adj_df, corpus_list, removed_ngrams, delim = '|'):
    """
    Function that builds a dictionary containing the 3 different frequency values for each n-gram of interest within an adjacency matrix. The frequencies that are returned are:
        - adjacency frequency: relative frequency of the n-gram within the n-grams of interest only
        - n-length frequency: relative frequency of the n-gram amongst only same length n-grams in the adjacency matrix
        - conditional frequency: the conditional frequency of each n-gram given the preceding (n-1)-grams. For unigrams this will refer to the frequency of the unigram within the entire corpus (i.e. not limited to unigrams within the adjacency matrix only)

    Inputs:
        - adj_df: dataframe that is the adjacency matrix of n-grams of interest within the corpus list
        - corpus_list: list containing all the words/structures composed of the tokens of interest
        - n_max: int of the maximum n-gram length
        - removed_ngrams: dict containing the n-grams that were removed due to complete overlap with larger n-grams. 
        - delim: str of the delimiter that separates the tokens of interest
        
    Outputs:
        freq_dict: dict containing frequency values mapped to each n-gram
    """
    ngrams_OI = adj_df.columns.tolist()
    ngrams_OI = set(ngrams_OI).union(removed_ngrams)
    freq_dict = build_freq_dict(ngrams_OI, corpus_list, delim)

    return freq_dict

def build_freq_dict(ngrams_OI,corpus_list,delim = '|'):
    """
    Function that builds a dictionary containing the 3 different frequency values for each n-gram of interest within an adjacency matrix. The frequencies that are returned are:
        - adjacency frequency: relative frequency of the n-gram within the n-grams of interest only
        - n-length frequency: relative frequency of the n-gram amongst only same length n-grams in the adjacency matrix
        - conditional frequency: the conditional frequency of each n-gram given the preceding (n-1)-grams. For unigrams this will refer to the frequency of the unigram within the entire corpus (i.e. not limited to unigrams within the adjacency matrix only)

    Inputs:
        - adj_df: dataframe that is the adjacency matrix of n-grams of interest within the corpus list
        - corpus_list: list containing all the words/structures composed of the tokens of interest
        - n_max: int of the maximum n-gram length
        - removed_ngrams: dict containing the n-grams that were removed due to complete overlap with larger n-grams. 
        - delim: str of the delimiter that separates the tokens of interest
        
    Outputs:
        freq_dict: dict containing frequency values mapped to each n-gram
    """

    nlength_dict = {}
    ngram_count = {}

    # Extract individual n-gram counts and categorize based on n-gram length too
    for ngram in ngrams_OI:
        ngram_length = len(ngram.split(delim))
        ngram_count[ngram] = get_ngram_instance_count(ngram, corpus_list,delim)
        
        if ngram_length not in nlength_dict:
            nlength_dict[ngram_length] = {}
        nlength_dict[ngram_length][ngram] = ngram_count[ngram]
    
    # Get unigram counts in the total corpus
    token_counts = get_token_list(corpus_list, delim)
    #token_tot = sum(token_counts.values())

    # quick check on unigrams given that some may have character overlaps between them but are distinct from each other.
    for ngram in nlength_dict[1]:
        if token_counts[ngram] != ngram_count[ngram]:
            ngram_count[ngram] = token_counts[ngram]
            nlength_dict[ngram_length][ngram] = token_counts[ngram]
            print('Note: count for the token %s had to readjusted'%ngram)

    # Get n-gram length based total counts values
    nlength_count = {}
    for n in nlength_dict:
        nlength_count[n] = sum(nlength_dict[n].values())

    # Getting the total number of n-grams (limited up to the length of interest) in the corpus
    tot_ngram_count = sum(ngram_count.values())

    freq_dict = {}
    for ngram in ngrams_OI:
        ngram_length = len(ngram.split(delim))
        n_p = ngram_count[ngram]/tot_ngram_count
        
        # Conditional Frequency Results
        if ngram_length == 1:
            bi = ''
            pbi = token_counts[ngram]/tot_ngram_count
            #pbi = token_counts[ngram]/token_tot
            bi_count = None
            
        else:
            bi = ngram.split('|')[:-1]
            bi = '|'.join(bi) # Note this is mostly useless for bigrams
            try:
                if bi in ngram_count:
                    pbi = ngram_count[ngram]/ngram_count[bi]
                    bi_count = ngram_count[bi]
                else:
                    bi_count = get_ngram_instance_count(bi,corpus_list,delim)
                    ngram_count[bi] = bi_count
                    pbi = ngram_count[ngram]/ngram_count[bi]
            except:
                print('There is an issue with n-gram %s the n-1 gram count is troublesome:\n %s Count: %s'%(ngram,bi, ngram_count[bi]))
                pbi = 0
                bi_count = 0
                

        freq_dict[ngram]= {'Conditional Freq.':pbi,'Count':ngram_count[ngram],'n-1 gram':bi, 'n-1 gram Count':bi_count, 'Count Freq.':n_p, 'Total Length Count':nlength_count[ngram_length], 'Total N-gram Count':tot_ngram_count}

    return freq_dict

def ngram_freq_diversity(ngrams_OI,corpus_list,delim = '|'):
    """
    Function that builds a dictionary containing the 3 different frequency values for each n-gram of interest within an adjacency matrix. The frequencies that are returned are:
        - adjacency frequency: relative frequency of the n-gram within the n-grams of interest only
        - n-length frequency: relative frequency of the n-gram amongst only same length n-grams in the adjacency matrix
        - conditional frequency: the conditional frequency of each n-gram given the preceding (n-1)-grams. For unigrams this will refer to the frequency of the unigram within the entire corpus (i.e. not limited to unigrams within the adjacency matrix only)

    Inputs:
        - adj_df: dataframe that is the adjacency matrix of n-grams of interest within the corpus list
        - corpus_list: list containing all the words/structures composed of the tokens of interest
        - n_max: int of the maximum n-gram length
        - removed_ngrams: dict containing the n-grams that were removed due to complete overlap with larger n-grams. 
        - delim: str of the delimiter that separates the tokens of interest
        
    Outputs:
        freq_dict: dict containing frequency values mapped to each n-gram
    """

    nlength_dict = {}
    ngram_count = {}

    # Extract individual n-gram counts and categorize based on n-gram length too
    for ngram in ngrams_OI:
        ngram_length = len(ngram.split(delim))
        ngram_count[ngram] = get_ngram_instance_count(ngram, corpus_list,delim)
        
        if ngram_length not in nlength_dict:
            nlength_dict[ngram_length] = {}
        nlength_dict[ngram_length][ngram] = ngram_count[ngram]
    
    # Get unigram counts in the total corpus
    token_counts = get_token_list(corpus_list, delim)

    # quick check on unigrams given that some may have character overlaps between them but are distinct from each other.
    for ngram in nlength_dict[1]:
        if token_counts[ngram] != ngram_count[ngram]:
            ngram_count[ngram] = token_counts[ngram]
            nlength_dict[ngram_length][ngram] = token_counts[ngram]
            print('Note: count for the token %s had to readjusted'%ngram)

    freq_dict = {}
    for ngram in ngrams_OI:
        ngram_length = len(ngram.split(delim))

        # Conditional Frequency Results
        if ngram_length > 1:
            # Get the final domain
            fin_dom = ngram.split('|')[-1]
            
            # Getting just the n-grams with the same length
            ngrams_to_check = nlength_dict[ngram_length]

            # Getting the n-grams whose final domain is the same as that of the currently evaluated n-gram
            common_end = {}
            for n1_gram in ngrams_to_check:
                temp = n1_gram.split('|')
                if temp[-1] == fin_dom:
                    common_end[n1_gram] = ngram_count[n1_gram]

            try:
                tot_possible_ngram_end = sum([v for v in common_end.values()])
                div_freq = ngram_count[ngram]/tot_possible_ngram_end
            except:
                print('There is an issue with n-gram %s the n-1 gram count is troublesome'%(ngram))
                div_freq = 0
                

            freq_dict[ngram]= {'Diversity Freq.':div_freq,'Count':ngram_count[ngram],'n-1 grams':common_end}

    return freq_dict

def crossEntropy(p_dict, q_dict):
    '''Calculates the cross entropy based of the two provided frequency distributions.
    
    Parameters:
    -----------
        - p_dict: dict
            A dict of frequencies of the model being compared
        - q_dict: dict
            A dict of frequencies of the baseline model.
            
    Returns:
    --------
        ce : int
            Cross entropy of the two frequency distributions
    '''

    # Checking that all n-grams in the model to be determined are found in the baseline otherwise raising an exception error
    memb_check = set(p_dict.keys()).difference(q_dict.keys())
    if len(memb_check) > 0:
        raise Exception('Cross entropy cannot be determined as an n-gram is missing from the baseline set.')

    # Now calculating the cross entropy based on the conditional probabilities.
    ce = 0
    for ngram in p_dict:
        p = p_dict[ngram]['Count Freq.']
        qbi = q_dict[ngram]['Conditional Freq.']
        if p != 0:
            ce += -p*np.log2(qbi)

    return ce

def relativeEntropy(p_dict, q_dict):
    '''
    Calculates the relative entropy or the Kullback-Leibler divergence of the two provided frequency distributions.
    
    Parameters:
    -----------
        - p_dict: dict
            A dict of frequencies of the model being compared
        - q_dict: dict
            A dict of frequencies of the baseline model.
            
    Returns:
    --------
        kl : int
            Relative entropy of the two frequency distributions
    '''

    memb_check = set(p_dict.keys()).difference(q_dict.keys())
    if len(memb_check) > 0:
        raise Exception('Cross entropy cannot be determined as an n-gram is missing from the baseline set.')

    # Now calculating the cross entropy based on the conditional probabilities.
    kl = 0
    for ngram in p_dict:
        p = p_dict[ngram]['Count Freq.']
        pbi = p_dict[ngram]['Conditional Freq.']
        qbi = q_dict[ngram]['Conditional Freq.']
        if p != 0:
            kl += p*np.log2(pbi/qbi)

    return kl    

def add_nc_term_ngram(ref_df, readable = True):
    '''
    Function that adds an N-terminal or C-terminal unigram to the n-gram of an entire column of the reference dataframe.

    Parameters:
    -----------
        - ref_df: pandas DataFrame
            DataFrame that contains the reference data including the Interpro Domain Architecture
        - readable: bool (optional)
            Bool to determine whether it should be added to the human-readable or ID version of the domain architecture.
    '''

    if readable:
        colOI = 'Interpro Domain Architecture'
    else:
        colOI = 'Interpro Domain Architecture IDs'
    
    # Check that the column is within the DataFrame
    if colOI not in ref_df.columns:
        raise Exception('The reference DataFrame is missing the column of interest.')
    
    # Now appending the N and C term designations which will be the same regardless of whether using the legible or not version of the domain architecture.
    ref_df[colOI] = ref_df[colOI].apply(lambda x: f'<Nterm|{x}|Cterm>')

    return ref_df

def calc_ngram_simpson_diversity(domain, freq_dict, adj, mode):
    '''
    Calculates the Simpson Diversity Index for the specified domain using the inputted frequency dict and adjacency matrix. N-grams that will be used for the diversity index have to either start or end with the domain of interest. This is designated with the mode parameter where 'pro' designates all n-grams of interest that start, while 'retro' corresponds to those that end with the domain. 
    
    This uses the small dataset with replacement definition of the Simpson Index:
        l = \Sigma_i^R n_i*(n_i - 1)/(N*(N-1)) 
            Where R is the richness (i.e. total domain partners), n_i is the count of the bigram, N is the total count of bigrams involving the domain of interest

    And returns the diversity index as defined as:
        D = 1-l
    
    Parameters:
    -----------
        domain: str
            Domain of interest to calculate the diversity index of
        freq_dict: dict
            Dict of all frequencies within the corpus of n-grams in the dataset
        adj: pandas DataFrame
            DataFrame containing the adjacency matrix of n-grams
        mode: str
            The mode for the n-gram children to be searched for must be either 'pro' or 'retro'
    Returns:
    --------
        D: float
            Simpson Diversity Index of the domain of interest         
    '''

    # Verifying that the domain is a unigram as this is currently not implemented to take in longer n-grams
    
    # Verifying that the mode is provided.
    if mode not in ['retro', 'pro']:
        raise ValueError('Improper mode designation.')
    
    # Based on the mode retrieving the n-gram children that will be used in the diversity calculations based on if it is 1 domain longer and either ends or starts with the domain of interest.
    child_cand = adj.columns[adj.loc[domain] > 0].tolist()
    ngramOI_len = len(domain.split('|'))
    if mode == 'retro':
        ngram_children = [ngram for ngram in child_cand if (len(ngram.split('|')) == ngramOI_len + 1) and '|'.join(ngram.split('|')[-ngramOI_len:]) == domain]
    else:
        ngram_children = [ngram for ngram in child_cand if (len(ngram.split('|')) == ngramOI_len + 1) and '|'.join(ngram.split('|')[:ngramOI_len]) == domain]
    
    # Retrieving the individual and total counts
    child_freq = {k:v['Count'] for k,v in freq_dict.items() if k in ngram_children}
    N = sum([v for v in child_freq.values()])
    
    # Calculating the Diversity Index
    if N > 1:
        l = calc_simpson_index(child_freq.values(), N)
        D = 1 - l
    else: # If a domain only has one n-gram associated with it a diversity index cannot be calculated.
        D = np.nan

    return D

def calc_simpson_index(n_list, N):
    '''
    Calcutes the Simpson Index given a list of individual species numbers and the total count.
    
    This uses the small dataset with replacement definition of the Simpson Index:
        l = \Sigma_i^R n_i*(n_i - 1)/(N*(N-1)) 
            Where R is the richness (i.e. total domain partners), n_i is the count of the bigram, N is the total count of bigrams involving the domain of interest

    Parameters:
    -----------
        - n_list: list
            List of species counts
        - N: int
            Total number of species of interest

    Returns:
        - l: float
            The Simpson Index
    '''

    l = 0
    for n in n_list:
        l += (n*(n-1))
    
    l = l/(N*(N-1))

    return l