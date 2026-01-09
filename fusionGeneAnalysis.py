import networkAnalysisUtilities as nAU
import pandas as pd
import networkx as nx
import numpy as np
import ngramUtilities
from pybiomart import Dataset
import math
from tqdm import tqdm
import liftover
import os
import dansy

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
    cancer_fusion_data = cancer_fusion_data[~cancer_fusion_data['Frame'].isin(['Out-of-Frame','Out-of-frame'])]
    #cancer_fusion_data = cancer_fusion_data[cancer_fusion_data['Genome_Build_Version'] == 'hg38']
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


def get_pt_fusion(cancer_fusion_data, ensembl_dict, conversion_df):
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
    removed_count = 0
    changed_convert = 0
    
    # Setting up a converter for the human genome
    converter = liftover.get_lifter('hg19', 'hg38', one_based=True)

    

    for idx, row in cancer_fusion_data.iterrows():

        genome = row['Genome_Build_Version']

        five_gene = {'Gene':row['H_gene'],
                        'Chromosome':row['H_chr'][3:],
                        'ChimerDB_Position':row['H_position'],
                        'Strand':row['H_strand'],
                        'Fusion pos':5
                        }
        three_gene = {'Gene':row['T_gene'],
                        'Chromosome':row['T_chr'][3:],
                        'ChimerDB_Position':row['T_position'],
                        'Strand':row['T_strand'],
                        'Fusion pos':3
                        }

        #   Performing a liftover if the genome is hg19
        if genome == 'hg19':
            for gene_OI in [five_gene, three_gene]:
                bp_position = gene_OI['ChimerDB_Position']
                liftover_issue = False
                chromosome = gene_OI['Chromosome']
                hg38info = converter.query(chromosome,bp_position)[0]
                if hg38info[0][3:] == chromosome:
                    bp_position = hg38info[1]
                else:
                    liftover_issue = True
                gene_OI['Position'] = bp_position
                gene_OI['liftover_issue'] = liftover_issue
        else:
            for gene_OI in [five_gene, three_gene]:
                gene_OI['Position'] = gene_OI['ChimerDB_Position']
                gene_OI['liftover_issue'] = False
        
        
        # Double check the conversion for uppercase instances of the gene name:
        # Gene checks to ensure can be used in downstream analysis
        # can_convert = (five_gene['Gene'] in ensembl_dict.keys()) and (three_gene['Gene'] in ensembl_dict.keys())
        # ck = False
        # if not can_convert:
        #     if five_gene['Gene'] not in ensembl_dict.keys():
        #         if five_gene['Gene'].upper() in ensembl_dict.keys():
        #             five_gene['Gene'] = five_gene['Gene'].upper()
        #             ck = True
        #     if three_gene['Gene'] not in ensembl_dict.keys():
        #         if three_gene['Gene'].upper() in ensembl_dict.keys():
        #             three_gene['Gene'] = three_gene['Gene'].upper()
        #             ck = True
        #     if ck:
        #         changed_convert += 1
        # # Gene checks to ensure can be used in downstream analysis
        # can_convert = (five_gene['Gene'] in ensembl_dict.keys()) and (three_gene['Gene'] in ensembl_dict.keys())

        # chrom_acceptable = (five_gene['Chromosome'] in CHROMOSOMES) and (three_gene['Chromosome'] in CHROMOSOMES)
        five_verified,fiveup,five_ensembl = verify_fusion_gene(five_gene,ensembl_dict,conversion_df)
        three_verified, threeup, three_ensembl = verify_fusion_gene(three_gene,ensembl_dict,conversion_df)

        # Making sure the gene is included in the ensembl conversion otherwise removing the pt from analysis
        if five_verified and three_verified:
            pos_check = {'Genes':{five_ensembl:five_gene,
                    three_ensembl:three_gene},
                    'Fusion Protein':row['Fusion_pair'],
                    'Genome Build':genome,
                    'Cancer Type':row['Cancertype'],
                    'Pt ID':row['BarcodeID']}
            pt_check[idx] = pos_check
        else:
            removed_count += 1
            
            #print('Patient %s has been removed from further analysis due to missing gene information.'%idx)

        if fiveup or threeup:
            changed_convert += 1

    print(f'There were {removed_count} fusion(s) dropped from further analysis due to missing gene information')
    print(f'There were {changed_convert} fusions that the gene name had to be changed to upppercase')
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
        genome = pt_OI['Genome Build']
        # Running through each member of the fusion gene and getting their information
        for transcript in pt_OI['Genes']:
            gene_OI = pt_OI['Genes'][transcript]
            gene_info = exon_info.get_group(transcript)
            
            # These should not vary for any exon
            gene_strand = gene_info['Strand'].tolist()[0]

            # The CDS start is a convenient way regardless of strand to order the exons
            gene_info = gene_info.sort_values(by = 'CDS start')
            
            liftover_issue = gene_OI['liftover_issue']
                    
            if liftover_issue:
                continue
            bp_position = pt_check[pt]['Genes'][transcript]['Position']

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
                exons_kept = exon_check[exon_check]
                fusion_exons = gene_info.filter(exons_kept.index, axis=0)
                cds_start = fusion_exons['CDS start'].tolist()
                cds_end = fusion_exons['CDS end'].tolist()
                
                if pt_OI[gene]['Fusion pos'] == 5:
                    seq_len = fusion_exons['CDS end'].tolist()[-1]    
                else:
                    seq_len = cds_end[-1] - cds_start[0] +1 - 3# The plus one is to account for length arithmetic while the -3 is for the stop codon
            else:
                seq_len = 0
                

            if any(exon_check != trunc_check):
                trunc_event = [exon_check != trunc_check]
                trunc_exon = gene_info.filter(exon_check.index[trunc_event[0]], axis = 0)
                trunc_diff = trunc_exon['Exon region end (bp)'].tolist()[0] - pt_OI[gene]['Position']
               
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
                # This represents the cases where the UTRs are spliced together for the 3' end so the amino acid position of the fusion is the CDS end
                cds_end = gene_info['CDS end'].tolist()
                aa_pos = max(cds_end)
                pts_2_remove.append(pt)
                
            
            valid_pts[pt]['Genes'][gene]['Fusion peptide length'] = aa_len
            valid_pts[pt]['Genes'][gene]['Fusion peptide aa position'] = aa_pos

    print(f'There were {len(pts_2_remove)} fusions with a UTR region as part of the fusion.')
    #for pt in pts_2_remove:
    #    del valid_pts[pt]

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
    missing_domain = 0
    for pt in pt_check:
        pt_OI = pt_check[pt]
        fusion_arch = []
        legible_arch = []
        arch_info = []
        for gene in pt_OI['Genes']:
            fusion_info = pt_OI['Genes'][gene]
            
            # Ensuring the gene has a UniProt ID associated with it
            if gene in ensembl_2_uniprot:
                #if len(ensembl_2_uniprot[gene]) == 1: 
                    prot_info = fusion_genes_df[fusion_genes_df['UniProt ID'] == ensembl_2_uniprot[gene]]
                    prot_fusion_contrib = []
                    prot_fusion_contrib_ids = []
                    
                    # Checking that the protein was fetched during reference file generation and then flagging the patient to be removed.
                    if prot_info.empty:
                        domains = []
                        #print('Due to missing domain annotations for %s patient %s will be excluded from further analysis.'%(gene, pt))
                        missing_domain += 1
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
                #else:
                #    print('Gene %s has multiple UniProt ID due to alternative splicing and thus patient %s will be exlcuded from further analysis.'%(gene, pt))
                #    pts_2_remove.append(pt) 
            else:
                #print('Gene %s was not mapped to a UniProt ID and thus patient %s will be exlcuded from further analysis.'%(gene, pt))
                pts_2_remove.append(pt)

        pt_check[pt]['Interpro Domains'] = ';'.join(arch_info)
        pt_check[pt]['Interpro Domain Architecture'] = '|'.join(legible_arch)
        pt_check[pt]['Interpro Domain Architecture IDs'] = '|'.join(fusion_arch)
    
    for pt in set(pts_2_remove):
        del pt_check[pt]

    print(f'Due to missing uniprot information {len(pts_2_remove)} fusion(s) were removed')
    print(f'There were {missing_domain} fusions that were removed due to domain annotations.')
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
    

    gross_changes = pd.DataFrame(columns = ['Connected Components','Articulation Points','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges'],index=pt_check.keys())
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
                #gross_changes.loc[pt,'Isolates'] = isol_orig - isol_pt
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
    cats = ['No Change','Interconnectedness','Connect Components','Empty Domain', 'Reintroduction']
    
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
        elif row['Connected Components'] == 0:
            cat = 'Interconnectedness'
        elif row['Connected Components'] > 0:
            cat = 'Connect Components'
        elif row['Connected Components'] == -1:
            cat = 'Reintroduction'
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

def perform_fusion_analysis(fusion_db,cancertype,conv_df,ref_data,fusion_categories = {}, version = 1, dansy_obj = None):
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
    G = nAU.simplify_network(G)
    # Getting the basic inputs of the fusions
    fusion_data = retrieve_cancer_fusion_data(fusion_db,cancertype)


    print('Retrieving Ensembl IDs of interest')
    ensembl_dict = extract_fusion_ensembl(fusion_data,conv_df)

    # Getting the basic exon information
    print('Adding exon and peptide information')
    exon_info = retrieve_ensembl_exon_data(filename='Gene_exon_information.csv')
    exon_info = exon_info.groupby('Gene stable ID')
    pt_check = get_pt_fusion(fusion_data, ensembl_dict,conv_df)
    pt_exon_check = get_pt_exons(pt_check,exon_info)
    pt_check = add_pt_peptide_info(pt_check, pt_exon_check, exon_info)
    
    # Get UniProt IDs
    uniprot_temp = conv_df.filter(['Gene stable ID', 'UniProtKB/Swiss-Prot ID']).drop_duplicates().dropna(axis=0, subset=['UniProtKB/Swiss-Prot ID'])
    uniprot_ensembl_conv = dict(zip(uniprot_temp['Gene stable ID'], uniprot_temp['UniProtKB/Swiss-Prot ID']))
    #uniprots_OI = set()
    #for _, row in conv_df.iterrows():
    #    g = row['Gene stable ID']
   #     uni_id = row['UniProtKB/Swiss-Prot ID']
  #      if g in [v for v in ensembl_dict.values()]:
 #           uniprots_OI.add(uni_id)
#            if str(uni_id) != 'nan':
#                uni_ensembl_conv.setdefault(g, set()).add(uni_id)
   
    #uni_ensembl_conv = {k:list(v) for k,v in uni_ensembl_conv.items()}

    uniprots_OI = [x for x in uniprot_ensembl_conv.values() if str(x) != 'nan']
    fusion_genes_df = ref_df.copy()
    fusion_genes_df = fusion_genes_df[fusion_genes_df['UniProt ID'].isin(uniprots_OI)]
    fusion_genes_df = fusion_genes_df.drop_duplicates(ignore_index=True)
    pt_removal = [pt for pt,details in pt_check.items() if 'Fusion Problem' in details.keys()]
    if pt_removal:
        for pt in pt_removal:
            print('Removing pt %s from downstream analysis'%pt)
            del pt_check[pt]
    pt_check = add_pt_domain_arch(pt_check,fusion_genes_df, uniprot_ensembl_conv)
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
    if version == 1:
        gross_changes = calc_fusion_network_changes(rep_pt_check,adj_df,removed_ngrams,proteome_gross_topo)
    elif version == 2 and dansy_obj is not None:
        gross_changes = calc_fusion_network_changes_v2(rep_pt_check, dansy_obj)
    fusion_categories_ammended = summarize_fusion_gross_changes(pt_check, uni_archs, cancertype, gross_changes, fusion_categories)
    
    return fusion_categories_ammended, pt_check, gross_changes
    
def calc_fusion_network_changes_v2(pt_check,proteome_dansy:dansy.dansy):
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

    Returns:
    --------
        - gross_changes: pandas DataFrame
            - dataframe containing the gross topological changes associated with the fusion proteins including the new nodes and edges. Negative numbers correspond to a gain for the non-node/edge measurements
        - comm_info: dict
            - dictionary containing the information for each patient of communities detected across all iterations. (Note this will change and is very memory intensive)
        printed statement of how many domain architectures were skipped.
    '''

    # from the dansy object getting the adjacency matrix and some key gross topology information

    removed_ngrams = proteome_dansy.collapsed_ngrams
    existing_ngrams = proteome_dansy.ngrams
    #isol_orig = orig_network_vals[1]
    #artic_orig = orig_network_vals[2]

    gross_changes = pd.DataFrame(columns = ['Connected Components','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges'],index=pt_check.keys())
    skip_cnt = 0 # Keeping count of how many are skipped for the entire dataset

    for pt in tqdm(pt_check):

        # Prior to doing any checks on the n-gram will check if the domain architecture already exists. If so then no need to run through the entire analysis
        pt_arch = pt_check[pt]['Interpro Domain Architecture IDs']
        
        if pt_arch in existing_ngrams:
            gross_changes.loc[pt,:] = 0
            skip_cnt += 1
        
        else:
            
            # Checking if any of the n-grams have been removed previously
            fusion_ngrams = pt_check[pt]['N-grams']
            ngram_to_return = list(set(fusion_ngrams).intersection(removed_ngrams))
            new_ngrams = set(fusion_ngrams).difference(existing_ngrams).difference(ngram_to_return)
            ngrams_to_remove = []
            
            # Checking for any n-grams that can be subsumed by a longer n-gram and then removing it
            for gram in ngram_to_return:
                for inner_gram in ngram_to_return:
                    if gram != inner_gram and gram in inner_gram:
                        ngrams_to_remove.append(gram)
                  
            if ngrams_to_remove:
                ngrams_to_remove = set(ngrams_to_remove)
                ngram_to_return = list(set(ngram_to_return).difference(ngrams_to_remove))
                #ngram_to_return = [x for x in ngram_to_return if x not in ngrams_to_remove]
            
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
                
                # N-grams to check focusing only on the connected components that the novel n-grams are within
                original_ngram_check = set()
                cc_impacted =0
                for i in nx.connected_components(proteome_dansy.G):
                    if len(set(fusion_ngrams).intersection(i)) > 0:
                        original_ngram_check.update(i)
                        cc_impacted += 1

                
                # Adding in all the new n-grams with default values of zero
                #new_df = pd.DataFrame(columns=list(new_ngrams), index=list(new_ngrams)).fillna(0)
                #pt_adj = pd.concat([pt_adj, new_df]).fillna(0)
                
                # Going through and adding connections in the adjacency matrix
                #for n in new_ngrams:
                #    new_connection_check = []
                #    for j in original_ngram_check:
                #        if n in j and n != j:
                #            new_connection_check.append(j in new_ngrams_orig)
                #            new_edge += 1
                #        elif j in n and n!= j:
                #            new_connection_check.append(j in new_ngrams_orig)
                #            new_edge += 1
                    
                    # Check if there was a completely subsumed n-gram that slipped through (Note: need to add what happens if this actually happens)
                 #   if all(new_connection_check):
                 #       print('There is an issue with the adjacency matrix for pt %s'%pt)
                  #      print('A subsumed n-gram made it through %s'%n)
                        #chk = pt_adj.loc[n] > 0
                        
                        #print(pt_adj.filter(pt_adj.index[chk]).columns.tolist())
                gross_changes.loc[pt,'Connected Components'] = cc_impacted - 1
                                    #gross_changes.loc[pt,'Isolates'] = isol_orig - isol_pt
                                    #gross_changes.loc[pt,'Articulation Points'] = artic_orig - artic_pt
                pt_check[pt]['New N-grams Added'] = new_ngrams
                
            else:
                
                gross_changes.loc[pt,'Connected Components'] = 0
                #gross_changes.loc[pt,'Isolates'] = 0
                #gross_changes.loc[pt,'Articulation Points'] = 0
        
            # Outputting the changes
                gross_changes.loc[pt,'New Nodes'] = len(new_ngrams)
                gross_changes.loc[pt,'New Edges'] = new_edge
                
    print('There were %s architecture(s) skipped due to presence in the original adjacency matrix.'%skip_cnt)


    return gross_changes