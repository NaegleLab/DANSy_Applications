import networkAnalysisUtilities as nAU
import pandas as pd
import networkx as nx
import numpy as np
import ngramUtilities
from pybiomart import Dataset
import math
from tqdm import tqdm

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
    cancer_fusion_data = cancer_fusion_data[cancer_fusion_data['Frame'] == 'In-Frame']
    cancer_fusion_data = cancer_fusion_data[cancer_fusion_data['Genome_Build_Version'] == 'hg38']
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
            - Dataframe containing the gene name, synonym and Ensembl IDs to allow for quick conversion. Ideally also contains versions of the IDs to allow for accurate conversions due to gene names/synonyms referring to different genes

    Returns:
    --------
        - fusion_gene_ensembl: list
            - list of ensembl IDs to be queried for further analysis
        -ensembl_dict: dict
            - dict containing the conversion between the gene name in ChimerDB and its Ensembl ID
    '''
    
    # Combining the genes of interest to get IDs for
    fusion_genes_OI = cancer_fusion_data['H_gene'].tolist() + cancer_fusion_data['T_gene'].tolist()
    fusion_genes_OI = list(set(fusion_genes_OI))
    
    # Initializing variables
    fusion_gene_ensembl = []
    ensembl_dict = {}
    ensembl_version = {}
    ensembl_primary = {}

    # Getting the conversions of interest
    for _, row in gene_ID_conversion.iterrows():
        g = row['Gene name']
        syn = row['Gene Synonym']
        version = row['Version (gene)']
        if g in fusion_genes_OI:
            
            # Due to some genes having multiple Ensembl IDs need to check if the conversion exists already and if so if it is the most recent version (as a shortcut). This shouldn't create new issues
            if g not in ensembl_dict:
                ensembl_dict[g] = row['Gene stable ID']
                ensembl_version[g] = version
                ensembl_primary[g] = 1 # Designating this as the primary name of the gene
            else:
                if (version > ensembl_version[g]) and ensembl_primary[g]:
                    ensembl_dict[g] = row['Gene stable ID']
                    ensembl_version[g] = version
                elif not ensembl_primary[g]: # For instances when a gene synonym matches the actual name of another gene (i.e. GARS1 has a synonym of SMAD1 which is not what is typically known as SMAD1)
                    ensembl_dict[g] = row['Gene stable ID']
                    ensembl_version[g] = version
                    ensembl_primary[g] = 1
            

            # Need to double-check if the gene synonym is also in as some patients will have the synonym instaed of the primary gene name
            if syn in fusion_genes_OI:
                if syn not in ensembl_dict:
                    ensembl_dict[syn] = row['Gene stable ID']
                    ensembl_version[syn] = version
                    ensembl_primary[syn] = 0 #Designating this as the synonym being found so that if a synonym matches an actual gene name it gets overwritten.
                else:
                    if (version > ensembl_version[syn]) and not ensembl_primary[syn]: # Since synonym will have the same verison
                        ensembl_dict[syn] = row['Gene stable ID']
                        ensembl_version[syn] = version

        elif syn in fusion_genes_OI:
            if syn in fusion_genes_OI:
                if syn not in ensembl_dict:
                    ensembl_dict[syn] = row['Gene stable ID']
                    ensembl_version[syn] = version
                    ensembl_primary[syn] = 0
                elif not ensembl_primary[syn]: # Do not touch the gene if it was an original gene name
                    if version > ensembl_version[syn]: # Since synonym will have the same verison
                        ensembl_dict[syn] = row['Gene stable ID']
                        ensembl_version[syn] = version

        # Final check if a strange issue where capitalizing all the letters matches the genes of interest    
        elif g.upper() in fusion_genes_OI:
            if g.upper not in ensembl_dict:
                ensembl_dict[g.upper()] = row['Gene stable ID']
                ensembl_version[g.upper()] = version
                ensembl_primary[g.upper()] = 1
            else:
                if version > ensembl_version[g.upper()]:
                    ensembl_dict[g.upper()] = row['Gene stable ID']
                    ensembl_version[g.upper()] = version
            
    fusion_gene_ensembl = list(set([v for v in ensembl_dict.values()]))

    return fusion_gene_ensembl, ensembl_dict

def get_pt_fusion(cancer_fusion_data, ensembl_dict):
    '''
    Gets the patient and their fusion information. If one of the fusion partners lacks an Ensembl ID then the patient is removed from the analysis pool.

    Parameters:
    -----------
        - canceer_fusion_data: pandas DataFrame
            - Dataframe containing all patient information
        - ensembl_dict: dict
            - Conversion dictionary between gene names and their corresponding Ensembl ID

    Returns:
    --------
        - pt_check: dict
            - dict containing the patient information including fusion gene locations and breakpoints along with very basic patient information.
    '''
    pt_check = {}
    for idx, row in cancer_fusion_data.iterrows():
        five_gene = {'Gene':row['H_gene'],
                    'Chromosome':row['H_chr'][3:],
                    'Position':row['H_position'],
                    'Strand':row['H_strand'],
                    'Fusion pos':5
                    }
        three_gene = {'Gene':row['T_gene'],
                    'Chromosome':row['T_chr'][3:],
                    'Position':row['T_position'],
                    'Strand':row['T_strand'],
                    'Fusion pos':3
                    }
        # Making sure the gene is included in the ensembl conversion otherwise removing the pt from analysis
        if (five_gene['Gene'] in ensembl_dict.keys()) and (three_gene['Gene'] in ensembl_dict.keys()):
            pos_check = {'Genes':{ensembl_dict[five_gene['Gene']]:five_gene,
                    ensembl_dict[three_gene['Gene']]:three_gene},
                    'Fusion Protein':row['Fusion_pair'],
                    'Genome Build':row['Genome_Build_Version'],
                    'Cancer Type':row['Cancertype'],
                    'Pt ID':row['BarcodeID']}
            pt_check[idx] = pos_check
        else:
            print('Patient %s has been removed from further analysis due to missing gene information.'%idx)

    return pt_check

def get_pt_exons(pt_check, exon_info):
    '''
    Retrieves the exon information of a patients gene fusions and whether it was included in the fusion gene and/or truncated. Can take in older genome build coordinates and converts to the current build.

    Parameters:
    -----------
        - pt_check: dict
            - dictionary containing the patient information
        - exon_info: pandas DataFrame
            - Dataframe containing the exon locations and protein coding sequence positions

    Returns:
    --------
        - pt_exon_check: dict
            - Dictionary with keys for each patient that contains a nested dictionary that includes pandas Dataframes of boolean values whether an exon was included within the fusion and/or truncated.
        printed statement of how many patients had misannotated strand information within the ChimerDB
    '''
    pt_exon_check = {}
    # For some patients they use hte old build of the human genome, while the pybiomart database uses the current build so creating a liftover object to convert coordinates.
    mismatch_flag = 0
    mismatch_cnt = 0
    for pt in pt_check:
        pt_OI = pt_check[pt]
        pt_exon_check[pt] = {}
        
        # Running through each member of the fusion gene and getting their information
        for transcript in pt_OI['Genes']:
            gene_OI = pt_OI['Genes'][transcript]
            gene_info = exon_info.get_group(transcript)
            
            # These should not vary for any exon
            gene_strand = gene_info['Strand'].tolist()[0]

            # The CDS start is a convenient way regardless of strand to order the exons
            gene_info = gene_info.sort_values(by = 'CDS start')
            
            # Retrieving the breakpoint position
            bp_position = gene_OI['Position']
            
            pt_check[pt]['Genes'][transcript]['hg38 Breakpoint Position'] = bp_position

            # Figuring out which exons are after the breakpoint
            gene_info['Exon Inclusion'] = gene_info['Exon region start (bp)'] >= bp_position
            gene_info['Exon Truncation'] = gene_info['Exon region end (bp)'] > bp_position

            # Correcting exon inclusion if it is the 5' gene on the + strand or the 3' gene on the  - strand
            if (gene_OI['Fusion pos'] == 5 and gene_strand == 1) or (gene_OI['Fusion pos'] == 3 and gene_strand == -1):
                gene_info['Exon Inclusion'] = ~gene_info['Exon Inclusion']
                gene_info['Exon Truncation'] = ~gene_info['Exon Truncation']
            
            # Excluding the exons that do not correspond to the CDS region
            gene_info['CDS NaN Check'] = np.isnan(gene_info['CDS start'])
            gene_info['CDS Exons'] = np.where(gene_info['CDS NaN Check'],False, gene_info['Exon Inclusion'])
            gene_info['CDS Truncation'] = np.where(gene_info['CDS NaN Check'],False, gene_info['Exon Truncation'])

            # There are some patient fusions where I notice a mismatch or missing strand annotations for the gene. For these I am taking the canonical strand orientation but printing a statement that there was an issue.
            if (gene_strand == 1 and gene_OI['Strand'] == '-') or (gene_strand == -1 and gene_OI['Strand'] == '+'):
                mismatch_flag = 1
                mismatch_cnt += 1

            # There are some exons which do not contain a CDS that are demarcated using NaN in the CDS start column
            bp_check = gene_info['CDS Exons']
            bp_trun_check = gene_info['CDS Truncation']
            pt_exon_check[pt][transcript]={'Inclusion':bp_check,'Truncation':bp_trun_check}   
    
    if mismatch_flag:
        print('There were %2.0f patient(s) had a mismatch of strand information'%mismatch_cnt)

    return pt_exon_check

def add_pt_peptide_info(pt_check, pt_exon_check,exon_info):
    '''
    Adds fusion peptide information based on the exon information provided and generates confidence scores (0-3) in the peptide annotations if there were issues with the exon information. 
    
    A score of 0 will indicate that there was some missing information or an overall issue that will lead to patient removal from downstream analysis. Non-integer lengths and truncation events will lead to decreases in confidence as this method uses solely position information to reduce information retrieved from Ensembl.

    Parameters:
    -----------
        - pt_check: dict
            - dictionary of patient information including gene fusion position information
        - pt_exon_check: dict
            - dictionary containing the exon information of genes of interest
        - exon_info: pandas DataFrame
            - dataframe containing the exon positional information
    
    Returns:
    --------
        - valid_pts: dict
            - modified pt_check with the same keys but with additional values regarding the amino acid position and confidence scores for each peptide.
    '''
    
    pts_2_remove =[]
    valid_pts = pt_check.copy()

    for pt in pt_exon_check:

        cur_pt = pt_exon_check[pt]
        for gene in cur_pt:
            pt_OI = valid_pts[pt]['Genes']
            exon_check = cur_pt[gene]['Inclusion']
            trunc_check = cur_pt[gene]['Truncation']
            gene_info = exon_info.get_group(gene).copy()
            gene_info = gene_info.sort_values(by = 'CDS start')
            
            if any(exon_check):
                if pt_OI[gene]['Fusion pos'] == 5:
                    exons_kept = exon_check[exon_check]
                    fusion_exons = gene_info.filter(exons_kept.index, axis=0)
                    seq_len = fusion_exons['CDS end'].tolist()[-1]
                    
                else:
                    
                    exons_kept = exon_check[exon_check]
                    fusion_exons = gene_info.filter(exons_kept.index, axis = 0)
                    cds_start = fusion_exons['CDS start'].tolist()
                    cds_end = fusion_exons['CDS end'].tolist()
                    seq_len = cds_end[-1] - cds_start[0] +1 - 3# The plus one is to account for length arithmetic while the -3 is for the stop codon
            else:
                seq_len = 0
                

            if any(exon_check != trunc_check):
                trunc_event = [exon_check != trunc_check]
                trunc_exon = gene_info.filter(exon_check.index[trunc_event[0]], axis = 0)
                trunc_diff = trunc_exon['Exon region end (bp)'].tolist()[0] - pt_OI[gene]['hg38 Breakpoint Position']
               
            else:
                trunc_diff = 0
               
            if seq_len > 0:
                aa_len = (seq_len - trunc_diff)/3
            else:
                aa_len = 0

            # Putting the amino acid position in for ease in converting to the domains
            if pt_OI[gene]['Fusion pos'] == 5:
                aa_pos = math.ceil(aa_len)
                
            elif seq_len > 0:
                aa_pos = math.ceil(cds_end[-1]/3 - aa_len - 1)
            else:
                aa_pos = 0
                print('Due to a likely alternative splicing product of gene %s patient %s will be excluded.'%(gene, pt))
                pts_2_remove.append(pt)
                
            
            valid_pts[pt]['Genes'][gene]['Fusion peptide length'] = aa_len
            valid_pts[pt]['Genes'][gene]['Fusion peptide aa position'] = aa_pos

    for pt in pts_2_remove:
        del valid_pts[pt]
    
    return valid_pts

def add_pt_domain_arch(pt_dict, fusion_genes_df, ensembl_2_uniprot):
    '''
    Adds the theoretical domain architecture is based on where the breakpoints within the genes are located for each provided patient.

    Parameters:
    -----------
        - pt_dict: dict
            - dictionary containing the patient and fusion gene information
        - fusion_genes_df: pandas DataFrame
            - dataframe containing the UniProt IDs and Interpro Domain architecture of proteins of interest
        - ensembl_2_uniprot: dict
            - dictionary containing conversions between the ENSEMBL ID to the UniProt ID
        
    Returns:
    --------
        - pt_check: dict
            - modified input dictionary with a new field containing the Interpro Domain Architecture information.
    '''
    # Going through each of the fusion proteins that were annotated and finding what the theoretical domain architecture is
    pts_2_remove = []
    pt_check = pt_dict.copy()
    for pt in pt_check:
        pt_OI = pt_check[pt]
        fusion_arch = []
        legible_arch = []
        arch_info = []
        for gene in pt_OI['Genes']:
            fusion_info = pt_OI['Genes'][gene]
            
            # Ensuring the gene has a UniProt ID associated with it
            if gene in ensembl_2_uniprot:
                if len(ensembl_2_uniprot[gene]) == 1: 
                    prot_info = fusion_genes_df[fusion_genes_df['UniProt ID'] == ensembl_2_uniprot[gene][0]]
                    prot_fusion_contrib = []
                    prot_fusion_contrib_ids = []
                    
                    # Checking that the protein was fetched during reference file generation and then flagging the patient to be removed.
                    if prot_info.empty:
                        domains = []
                        print('Due to missing domain annotations for %s patient %s will be excluded from further analysis.'%(gene, pt))
                        pts_2_remove.append(pt)
                    else:
                        pt_check[pt]['Genes'][gene]['UniProt ID'] = prot_info['UniProt ID'].tolist()[0]
                        domains = prot_info['Interpro Domains'].str.split(';').tolist()[0]
                        pt_OI['Genes'][gene]['Interpro Domain Architecture'] = prot_info['Interpro Domain Architecture'].tolist()[0]
                        pt_OI['Genes'][gene]['Interpro Domain Architecture IDs'] = prot_info['Interpro Domain Architecture IDs'].tolist()[0]
                        
                    
                    if domains and domains != ['']:
                        for d in domains:
                            dom_info = d.split(':')
                            if fusion_info['Fusion pos'] == 5:
                                if int(dom_info[3]) <= fusion_info['Fusion peptide aa position']:
                                    fusion_arch.append(dom_info[1])
                                    legible_arch.append(dom_info[0])
                                    arch_info.append(d)

                                    # Adding in the fusion protein contribution:
                                    prot_fusion_contrib.append(dom_info[0])
                                    prot_fusion_contrib_ids.append(dom_info[1])
                            else:
                                if int(dom_info[2]) >= fusion_info['Fusion peptide aa position']:
                                    fusion_arch.append(dom_info[1])
                                    legible_arch.append(dom_info[0])
                                    arch_info.append(d)
                                    
                                    # Adding in the fusion protein contribution:
                                    prot_fusion_contrib.append(dom_info[0])
                                    prot_fusion_contrib_ids.append(dom_info[1])
                    pt_OI['Genes'][gene]['Fusion Protein Contribution'] = '|'.join(prot_fusion_contrib)
                    pt_OI['Genes'][gene]['Fusion Protein Contribution IDs'] = '|'.join(prot_fusion_contrib_ids)
                else:
                    print('Gene %s has multiple UniProt ID due to alternative splicing and thus patient %s will be exlcuded from further analysis.'%(gene, pt))
                    pts_2_remove.append(pt) 
            else:
                print('Gene %s was not mapped to a UniProt ID and thus patient %s will be exlcuded from further analysis.'%(gene, pt))
                pts_2_remove.append(pt)

        pt_check[pt]['Interpro Domains'] = ';'.join(arch_info)
        pt_check[pt]['Interpro Domain Architecture'] = '|'.join(legible_arch)
        pt_check[pt]['Interpro Domain Architecture IDs'] = '|'.join(fusion_arch)
    
    for pt in set(pts_2_remove):
        del pt_check[pt]

    return pt_check

def add_pt_ngrams(pt_check,interpro_dict):
    '''
    Adds all possible n-grams generated by the fusion protein for each patient.

    Parameters:
    -----------
        - pt_check: dict
            - dictionary containing all patient information
        - interpro_dict: dict
            - dictionary containing the key-value pairs of domain name-InterPro ID
    
    Returns:
    --------
        - pt_check: dict
            - modified version of of the input dictionary with a new key-value pair of n-grams
    
    '''
    all_domains = [x for x in interpro_dict.keys()]

    for fusion_dets in pt_check.values():
        dom_arch = fusion_dets['Interpro Domain Architecture IDs']
        ngrams = ngramUtilities.return_all_n_grams_from_key(dom_arch, all_domains)
        ngram_legible = []
        # Converting into a human legible version
        for n in ngrams:
            gram = n.split('|')
            gram_convert = []
            for k in gram:
                k_con = interpro_dict[k]
                gram_convert.append(k_con)
            gram_convert_str = '|'.join(gram_convert)
            ngram_legible.append(gram_convert_str)
        fusion_dets['N-grams'] = ngrams
        fusion_dets['N-grams Legible'] = ngram_legible
    
    return pt_check


def calc_fusion_network_changes(pt_check,adj_df, removed_ngrams,orig_network_vals):
    '''
    Determines and calculates simple descriptors of network changes within caused by the presence of n-grams associated with fusion proteins. 
    
    There are two categories of network changes that are determined: gross topology changes and localized soft clustering changes. If the patient's fusion protein domain architecture already exists within the proteome it will not be assessed.
    
    Parameters:
    -----------
        - pt_check: dict
            - dictionary containing all patient information
        - adj_df: pandas DataFrame
            - dataframe that makes up the adjacency matrix of the complete (or at least reference) proteome of interst
        - removed_ngrams: list
            - list of n-grams that were subsumed by others within the adjacency matrix
        - orig_network_vals: list
            - list containing the reference network values of the number of connected components, isolates, and articulation points (in order)
        - k: int (Optional)
            - number of iterations to be performed for the community detection
        - seed: int (Optional)
            - random number generator seed

    Returns:
    --------
        - gross_changes: pandas DataFrame
            - dataframe containing the gross topological changes associated with the fusion proteins including the new nodes and edges. Negative numbers correspond to a gain for the non-node/edge measurements
        - comm_info: dict
            - dictionary containing the information for each patient of communities detected across all iterations. (Note this will change and is very memory intensive)
        printed statement of how many domain architectures were skipped.
    '''
    cc_orig = orig_network_vals[0]
    isol_orig = orig_network_vals[1]
    artic_orig = orig_network_vals[2]
    

    gross_changes = pd.DataFrame(columns = ['Connected Components','Isolates','Articulation Points','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges'],index=pt_check.keys())
    skip_cnt = 0 # Keeping count of how many are skipped for the entire dataset
    for pt in tqdm(pt_check):
        pt_adj = adj_df.copy()
        # For each patient there will be a new adjacency matrix that the n-grams from the fusion protein architecture will be added to

        # Prior to doing any checks on the n-gram will check if the domain architecture already exists. If so then no need to run through the entire analysis
        pt_arch = pt_check[pt]['Interpro Domain Architecture IDs']
        
        if pt_arch in pt_adj.index.tolist():
            gross_changes.loc[pt,:] = 0
            skip_cnt += 1
        
        else:
            
            # Checking if any of the n-grams have been removed previously
            fusion_ngrams = pt_check[pt]['N-grams']
            ngram_to_return = list(set(fusion_ngrams).intersection(removed_ngrams))
            new_ngrams = set(fusion_ngrams).difference(pt_adj.index.tolist()).difference(ngram_to_return)
            ngrams_to_remove = []
            
            # Checking for any n-grams that can be subsumed by a longer n-gram and then removing it
            for gram in ngram_to_return:
                for inner_gram in ngram_to_return:
                    if gram != inner_gram and gram in inner_gram:
                        ngrams_to_remove.append(gram)
                  
            if ngrams_to_remove:
                ngrams_to_remove = set(ngrams_to_remove)
                ngram_to_return = [x for x in ngram_to_return if x not in ngrams_to_remove]
            
            # Saving the easy gross topological changes 
            gross_changes.loc[pt,'Reintroduced N-grams'] = len(ngram_to_return)
            gross_changes.loc[pt,'New N-gram Count'] = len(new_ngrams)
            
            # Checking if any of the new n-grams can be subsumed by a longer new n-gram that was also added
            ngrams_to_remove = []
            for gram in new_ngrams:
                for inner_gram in new_ngrams:
                    if gram != inner_gram and gram in inner_gram:
                        ngrams_to_remove.append(gram)

            if ngrams_to_remove:
                ngrams_to_remove = set(ngrams_to_remove)
                new_ngrams = set([x for x in new_ngrams if x not in ngrams_to_remove])

            # Due to ease adding in the reintroduced n-grams but keeping track of the actual new ones
            new_ngrams_orig = new_ngrams.copy()
            new_ngrams = new_ngrams.union(ngram_to_return)
            new_edge = 0
            if new_ngrams:
                # Adding in all the new n-grams with default values of zero
                new_df = pd.DataFrame(columns=list(new_ngrams), index=list(new_ngrams)).fillna(0)
                pt_adj = pd.concat([pt_adj, new_df]).fillna(0)
                
                # Going through and adding connections in the adjacency matrix
                for n in new_ngrams:
                    new_connection_check = []
                    for j in pt_adj.columns.tolist():
                        if n in j and n != j:
                            pt_adj.loc[n,j] += 1
                            new_connection_check.append(j in new_ngrams_orig)
                            new_edge += 1
                        elif j in n and n!= j:
                            pt_adj.loc[j,n] += 1
                            new_connection_check.append(j in new_ngrams_orig)
                            new_edge += 1
                    
                    # Check if there was a completely subsumed n-gram that slipped through (Note: need to add what happens if this actually happens)
                    if all(new_connection_check):
                        print('There is an issue with the adjacency matrix for pt %s'%pt)
                        print('A subsumed n-gram made it through %s'%n)
                        chk = pt_adj.loc[n] > 0
                        
                        print(pt_adj.filter(pt_adj.index[chk]).columns.tolist())

            # Outputting the changes
            gross_changes.loc[pt,'New Nodes'] = len(new_ngrams)
            gross_changes.loc[pt,'New Edges'] = new_edge

            # Start of the network-based analysis and not simple node-edge counting
            if new_ngrams:
                G_pt = nx.from_pandas_adjacency(pt_adj)
                G_pt.remove_edges_from(nx.selfloop_edges(G_pt))
                # Get a couple of key parameters that can potentially change of the gross topology of the network
                cc_pt = nx.number_connected_components(G_pt)
                isol_pt = nx.number_of_isolates(G_pt)
                artic_pt = len(list(nx.articulation_points(G_pt)))
            
                gross_changes.loc[pt,'Connected Components'] = cc_orig - cc_pt
                gross_changes.loc[pt,'Isolates'] = isol_orig - isol_pt
                gross_changes.loc[pt,'Articulation Points'] = artic_orig - artic_pt
                pt_check[pt]['New N-grams Added'] = new_ngrams

            else:
                
                gross_changes.loc[pt,'Connected Components'] = 0
                gross_changes.loc[pt,'Isolates'] = 0
                gross_changes.loc[pt,'Articulation Points'] = 0
                
    print('There were %s architecture(s) skipped due to presence in the original adjacency matrix.'%skip_cnt)
    return gross_changes



def summarize_fusion_gross_changes(pt_check, uni_archs, cancertype, gross_changes, old_summary = {}):

    # Copying the input data frame if desired.
    fusion_categories = old_summary.copy()

    # Creating the initial categories for each
    cats = ['No Change','Interconnectedness','Connect Components','Empty Domain']
    
    for c in cancertype:
        fusion_categories[c] = {k:{'Pt Count':0, 'Pt Members':[], 'Unique Architectures':set()} for k in cats}
    
    print('Simple Gross Network Summary Generation')
    for idx, row in gross_changes.iterrows():
        pt_arch = pt_check[idx]['Interpro Domain Architecture']
        if sum(row) == 0:
            if pt_check[idx]['Interpro Domains'] == '':
                cat = 'Empty Domain'
            else:
                cat = 'No Change'
        elif row['Connected Components'] == 0 and row['Articulation Points'] <= 0:
            cat = 'Interconnectedness'
        elif row['Connected Components'] > 0 and row['Articulation Points'] < 0:
            cat = 'Connect Components'
        elif row['Connected Components'] == 0 and row['Articulation Points'] > 0:
            cat = 'Interconnectedness'
        elif row['Connected Components'] > 0 and row['Articulation Points'] == 0:
            cat = 'Connect Components'
        else:
            print('pt %s has wrecked havoc on the network likely due to using a domain name and not the Interpro ID'%idx)
        
        # Going through each patient associated with the n-gram of interest and accurately putting in each category
        for pt in uni_archs[pt_arch]['Patients']:
            pt_cancer = pt_check[pt]['Cancer Type']
            fusion_categories[pt_cancer][cat]['Pt Count'] += 1
            fusion_categories[pt_cancer][cat]['Pt Members'].append(pt)
        
        # Keeping track of unique n-gram domain architectures associated with each category
        if cat not in ['Empty Domain', 'No Change']:
            fusion_categories[pt_cancer][cat]['Unique Architectures'].add(pt_arch)               

    # adding in the empty domain patients that were omitted from the network calculation
    if '' in uni_archs:
        for pt in uni_archs['']['Patients']:
            pt_cancer = pt_check[pt]['Cancer Type']
            fusion_categories[pt_cancer]['Empty Domain']['Pt Count'] += 1
            fusion_categories[pt_cancer]['Empty Domain']['Pt Members'].append(pt)

    return fusion_categories

def retrieve_ensembl_exon_data(dataset, fusion_ensembl):
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
    if len(fusion_ensembl) <= 200:
        gene_pos_info = dataset.query(attributes=['ensembl_gene_id','exon_chrom_start','exon_chrom_end','cds_start','strand','transcript_length','cds_end'],filters = {'link_ensembl_gene_id':fusion_ensembl,'transcript_is_canonical':True})
    else:
        print('Need to split the ensembl list')
        gene_pos_info = pd.DataFrame()
        k = math.ceil(len(fusion_ensembl)/200)
        print('Will query ensembl %s times'%k)
        for x in range(0,k):
            temp_ensembl = fusion_ensembl[x*200:min([(x+1)*200,len(fusion_ensembl)])]
            temp_pos = dataset.query(attributes=['ensembl_gene_id','exon_chrom_start','exon_chrom_end','cds_start','strand','transcript_length','cds_end'],filters = {'link_ensembl_gene_id':temp_ensembl,'transcript_is_canonical':True})
            gene_pos_info = pd.concat([gene_pos_info, temp_pos], ignore_index=True)
        gene_pos_info.drop_duplicates(inplace=True)
        gene_pos_info.reset_index(drop=True)
    
    return gene_pos_info

def perform_fusion_analysis(fusion_db,cancertype,conv_df,dataset,ref_data,fusion_categories = {}, soft_clust_flag = 0):
    '''
    The main function that goes through each step of the fusion analysis for patients of a specific cancer type. Note for multiple cancer types it is recommended to limit only to 3 due to the ensembl querying during the analysis.

    Parameters:
    -----------
        - fusion_db: pandas DataFrame
            - dataframe containing all information regarding the fusion genes
        - cancertype: list
            - list containing all cancer types of interest
        - conv_df: pandas DataFrame
            - Dataframe containing the gene name, synonym, UniProt and Ensembl IDs to allow for quick conversion. Ideally also contains versions of the Ensembl IDs to allow for accurate conversions due to gene names/synonyms referring to different genes
        - dataset: pybiomart Dataset object
            - pybiomart Dataset object that ensembl information is derived from
        - ref_data: list
            - list containing the reference proteome data including (in order) the reference file DataFrame, adjacency matrix, removed n-gram list, interpro conversion dictionary, and gross topological information about the n-gram network
    
    Returns:
    --------
        - fusion_categories_ammended: dict
            - dictionary of the different categories a fusion gene falls under within the network changes according to gross topological changes
        - pt_check: dict
            - dictionary of all patient information (Note: this is planned to be removed in the next iteration of this function once some additional functions are generated.)
        - pt_comm: dict
            - dictionary containing all community detection results for each patient (Note: this is planned to be removed in the next iteration of this function once some additional functions are generated.)
 
    '''
    
    # Unpacking the reference data that is generated
    ref_df = ref_data[0]
    adj_df = ref_data[1]
    removed_ngrams = ref_data[2]
    interpro_dict = ref_data[3]
    proteome_gross_topo = ref_data[4]

    # Generating the background distributions of soft clusters for the reference data
    print('Generating reference network')
    G = nx.from_pandas_adjacency(adj_df)
    seeds = nAU.create_iteration_seeds(k = 100, seed = 882) #Will change these later for user decisions
    G = nAU.simplify_network(G)
    if soft_clust_flag:
        print('Performing soft clustering on the reference network')
        comm_check = nAU.generate_community_reference(G, seeds[0:50])
        orig_comm_info = nAU.extract_comm_iter_info(comm_check,50)
        comm_check_2 = nAU.generate_community_reference(G, seeds[50:100])
        orig_comm_info_2 = nAU.extract_comm_iter_info(comm_check_2,50)
        
        # Soft cluster results for the reference
        comp_freq = nAU.get_node_soft_clust_res(orig_comm_info,50)
        comp_freq_2 = nAU.get_node_soft_clust_res(orig_comm_info_2,50)
        t = nAU.get_freq_diff_noise_thres(comp_freq, comp_freq_2)
    else:
        print('No soft clustering being performed.')
    # Getting the basic inputs of the fusions
    fusion_data = retrieve_cancer_fusion_data(fusion_db,cancertype)


    print('Retrieving Ensembl IDs of interest')
    fusion_ensembl, ensembl_dict = extract_fusion_ensembl(fusion_data,conv_df)

    # Getting the basic exon information
    print('Adding exon and peptide information')
    exon_info = retrieve_ensembl_exon_data(dataset, fusion_ensembl)
    exon_info = exon_info.groupby('Gene stable ID')
    pt_check = get_pt_fusion(fusion_data, ensembl_dict)
    pt_exon_check = get_pt_exons(pt_check,exon_info)
    pt_check = add_pt_peptide_info(pt_check, pt_exon_check, exon_info)
    
    # Get UniProt IDs
    uni_ensembl_conv = {}
    uniprots_OI = set()
    for _, row in conv_df.iterrows():
        g = row['Gene stable ID']
        uni_id = row['UniProtKB/Swiss-Prot ID']
        if g in fusion_ensembl:
            uniprots_OI.add(uni_id)
            if str(uni_id) != 'nan':
                uni_ensembl_conv.setdefault(g, set()).add(uni_id)
   
    uni_ensembl_conv = {k:list(v) for k,v in uni_ensembl_conv.items()}

    uniprots_OI = [x for x in uniprots_OI if str(x) != 'nan']
    fusion_genes_df = ref_df.copy()
    fusion_genes_df = fusion_genes_df[fusion_genes_df['UniProt ID'].isin(uniprots_OI)]
    fusion_genes_df = fusion_genes_df.drop_duplicates(ignore_index=True)
    pt_removal = [pt for pt,details in pt_check.items() if 'Fusion Problem' in details.keys()]
    if pt_removal:
        for pt in pt_removal:
            print('Removing pt %s from downstream analysis'%pt)
            del pt_check[pt]
    pt_check = add_pt_domain_arch(pt_check,fusion_genes_df, uni_ensembl_conv)
    pt_check = add_pt_ngrams(pt_check,interpro_dict)

    # Here determining how many unique architectures there are and thus how many times new networks have to be created
    uni_archs = {}

    for pt in pt_check:
        pt_arch = pt_check[pt]['Interpro Domain Architecture'] # Will need to eventually change this to the IDs version
        if pt_arch in uni_archs:
            uni_archs[pt_arch]['Patients'].append(pt)
            uni_archs[pt_arch]['Count'] += 1
        else:
            uni_archs[pt_arch] = {'Patients':[pt],'Count':1}

    # From the unique architectures grabbing representative patients to limit the number of networks needing to be created
    rep_pts = [v['Patients'][0] for k,v in uni_archs.items() if k != '']
    rep_pt_check = pt_check.copy()
    pt_2_drop = [pt for pt in rep_pt_check if pt not in rep_pts]
    for pt in pt_2_drop:
        del rep_pt_check[pt]
    print('Reduced the network change calculation to %s unique architectures from %s.'%(len(rep_pt_check), len(pt_check)))
    
    print('Determining changes in the fusion network')
    gross_changes, pt_comm = calc_fusion_network_changes(rep_pt_check,adj_df,removed_ngrams,proteome_gross_topo, soft_cluster_flag= soft_clust_flag)
    fusion_categories_ammended = summarize_fusion_gross_changes(pt_check, uni_archs, cancertype, gross_changes, fusion_categories)
    
    # This is the soft clustering analysis
    if soft_clust_flag:
        diff_alt = {k:{} for k in cancertype}
  
        for pt in pt_comm:
            comm_for_check = pt_comm[pt]
            fusion_arch = pt_check[pt]['Interpro Domain Architecture']
            arch_count = uni_archs[fusion_arch]['Count']
            node_freq_pt = nAU.get_node_soft_clust_res(comm_for_check, 50)
            temp_comp = nAU.calc_soft_clust_membership_diff(node_freq_pt, comp_freq)
            diff_comms = nAU.compare_soft_cluster_freq(comp_freq, node_freq_pt, thres=t)
            for c in diff_comms:
                for inner_pt in uni_archs[fusion_arch]['Patients']:
                    pt_cancer = pt_check[inner_pt]['Cancer Type']
                    if c not in diff_alt[pt_cancer]:
                        diff_alt[pt_cancer][c] = 1
                    else:    
                        diff_alt[pt_cancer][c] += 1
            
            '''# Try to grab the fusion gene information
            if fusion_arch in temp_comp:    
                fusion_comm = temp_comp[fusion_arch]
                for c in fusion_comm:
                    if c not in diff_alt.index:
                        diff_alt.loc[c] = 0
                    diff_alt.loc[c, 'Fusion Count'] += arch_count'''
    else:
        diff_alt = {}
    return fusion_categories_ammended, pt_check, pt_comm, diff_alt

