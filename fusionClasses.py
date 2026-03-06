import pandas as pd
import numpy as np
import liftover
import math
import dansy
import networkx as nx
import os
from typing import Literal
from tqdm import tqdm

CHROMOSOMES = [str(x) for x in np.arange(1,23)] + ['X','Y','MT']
LIFTER = liftover.get_lifter('hg19', 'hg38', one_based=True)
EXON_INFO = pd.read_csv('Gene_exon_information.csv')
EXON_INFO = EXON_INFO.groupby('Gene stable ID')
GENE_CONVERSION = pd.read_csv('ENSEMBL_Gene_Conversion.csv')
DANSY_REFERENCE_DF = dansy.import_proteome_files(ref_file_dir='data/Current_Human_Proteome',
                                                 ref_file_suffix='2026_0108.csv')
FUSION_POS_CONV = {5:'h', 3:'t'}
NT_BUFFER = 500
AA_BUFFER = 2

# Create a gene conversion dict to convert between ENSEMBL IDs and UniProt IDs
ENSEMBL2UNIPROT = GENE_CONVERSION.filter(['Gene stable ID', 'UniProtKB/Swiss-Prot ID']).drop_duplicates()
ENSEMBL2UNIPROT = dict(zip(ENSEMBL2UNIPROT['Gene stable ID'], ENSEMBL2UNIPROT['UniProtKB/Swiss-Prot ID']))


class fusionGene():
    '''
    Class to contain all relevant fusion gene information and predict a potential domain architecture for the fusion.

    Parameters
    ----------
    - name: str
        Name of the fusion pair (e.g. BCR-ABL1)
    - h_gene, t_gene: str
        The ENSEMBL ID of the head (or 5')/tail (or 3') gene
    - h_chr, t_chr: str
        The chromosome of the head/tail gene in UCSC annotation (i.e. chr1, chrX, etc.)
    - h_pos, t_pos: int
        The breakpoint position for the head/tail gene
    - h_strand, t_strand: str (Optional, Recommended)
        Which strand a gene is located. Can be omitted as it will be pulled from the ENSEMBL database
    - genome_build: str (Default: hg38)
        The genome build the positions have been derived from
    - id: str (Optional, Recommended)
        An ID for the gene fusion and breakpoint combination.
    - silent: bool
            Whether initial errors in fusion gene formation should result in an error or only log message.
    - add_info: dict
        Key-value pairs of additional information for the fusion that should be added
    '''

    def __init__(self,name:str,h_gene, t_gene, h_chr,h_pos, t_chr, t_pos, h_strand = None,t_strand = None, genome_build:Literal['hg38','hg19'] = 'hg38', id=None, silent = False, add_info = None):
        self.name = name
        self.h_gene = h_gene
        self.t_gene = t_gene
        self.errors = []
        c2r = True

        # Mapping the ENSEMBL genes to their UniProt IDs
        if all([h_gene in ENSEMBL2UNIPROT, t_gene in ENSEMBL2UNIPROT]):
            self.h_uniprot = ENSEMBL2UNIPROT[h_gene]
            self.t_uniprot = ENSEMBL2UNIPROT[t_gene]
        elif silent:
            self.errors.append('0:Unable to map ENSEMBL to UniProt.')
            c2r = False
        else:
            raise ValueError('Cannot map one of the genes to UniProt. Check that an ENSEMBL ID was provided.')

        # Genomic Information
        self.h_chr = h_chr
        self.t_chr = t_chr
        self.h_pos = h_pos
        self.t_pos = t_pos
        self.h_strand = h_strand
        self.t_strand = t_strand
        self.id = id


        if genome_build == 'hg19':
            self.original_bps = {'H_chr':h_chr, 'T_chr':t_chr,'H_pos':h_pos, 'T_pos':t_pos}
            self.h_chr, self.h_pos = LIFTER.query(h_chr, h_pos)[0][0:2]
            self.t_chr, self.t_pos = LIFTER.query(t_chr, t_pos)[0][0:2]

        if c2r:
            self.get_fusion_exons()
            self.add_peptide_information()
            self.create_domain_arch()
            
        
        if self.errors:
            pass
        else:
            self.get_domain_ngrams()
        
        if add_info is not None:
            self.add_info_attribs = []
            for k,v in add_info.items():
                if ' ' in k:
                    k = k.replace(' ','_')
                self.__setattr__(k, v)
                self.add_info_attribs.append(k)

    def get_fusion_exons(self):
        
        for fusion_pos, var_prefix in FUSION_POS_CONV.items():
            
            gene = self.__getattribute__(var_prefix+'_gene')
            fusion_strand = self.__getattribute__(var_prefix+'_strand')
            bp_position = self.__getattribute__(var_prefix+'_pos')
            
            # Get established genomic and exon information from the genome
            gene_info = EXON_INFO.get_group(gene)
            gene_info = gene_info.sort_values(by = 'CDS start')
            gene_strand = gene_info['Strand'].tolist()[0]
            
            
            # Double-check if there is a mismatch in strand annotations and setting a flag for logging purposes
            if 'strand_mismatch' in self.__dict__:
                if self.strand_mismatch:
                    pass
                else:
                    self.strand_mismatch = fusion_strand != gene_strand        
            else:
                self.strand_mismatch = fusion_strand != gene_strand

            # Figuring out which exons are after the breakpoint and if the breakpoint is within/before the exon end
            gene_info['Exon Inclusion'] = gene_info['Exon region start (bp)'] <= bp_position
            gene_info['Exon Truncation'] = gene_info['Exon region end (bp)'] > bp_position
                
            # Creating a correction flag since the above can have a small issue later on when we take the strand into account and if the fusion position is 5' or 3'
            trunc_correction = any((gene_info['Exon Inclusion'] == gene_info['Exon Truncation']) & gene_info['Exon Inclusion'])
            if trunc_correction:
                trunc_check = gene_info[(gene_info['Exon Inclusion'] == gene_info['Exon Truncation']) & gene_info['Exon Inclusion']].index
                
            # Correcting exon inclusion if it is the 5' gene on the + strand or the 3' gene on the  - strand
            if (fusion_pos == 5 and gene_strand == -1) or (fusion_pos == 3 and gene_strand == 1): #
                gene_info['Exon Inclusion'] = ~gene_info['Exon Inclusion']
                gene_info['Exon Truncation'] = ~gene_info['Exon Truncation']
            
            # For specific cases when the tail gene has the breakpoint at the start of an exon that would otherwise be removed with the correction above
            if all([trunc_correction, fusion_pos == 3, gene_strand == 1]):
                gene_info.loc[trunc_check,'Exon Inclusion'] = True
                gene_info.loc[trunc_check,'Exon Truncation'] = True
            
            if all([trunc_correction, fusion_pos == 5, gene_strand == -1]):
                gene_info.loc[trunc_check,'Exon Inclusion'] = True
                gene_info.loc[trunc_check,'Exon Truncation'] = True

            # Excluding the exons that do not correspond to the CDS region
            gene_info['CDS NaN Check'] = np.isnan(gene_info['CDS start'])
            gene_info['CDS Exons'] = np.where(gene_info['CDS NaN Check'],False, gene_info['Exon Inclusion'])
            gene_info['CDS Truncation'] = np.where(gene_info['CDS NaN Check'],False, gene_info['Exon Truncation'])

            self.__setattr__(var_prefix+'_exon_inclusion',gene_info['CDS Exons'])
            self.__setattr__(var_prefix+'_exon_truncation',gene_info['CDS Truncation'])
            
    def add_peptide_information(self):

        for fusion_pos, var_prefix in FUSION_POS_CONV.items():
            exon_check = self.__getattribute__(var_prefix+'_exon_inclusion')
            trunc_check = self.__getattribute__(var_prefix+'_exon_truncation')
            gene = self.__getattribute__(var_prefix+'_gene')
            gene_info = EXON_INFO.get_group(gene).copy().sort_values(by = 'Genomic coding start')
            rel_pos = self.__getattribute__(var_prefix+'_pos')
            
            # If any exons are included in the fusion (very rare for none to be but it happens with UTR regions only)
            if any(exon_check):
                exons_kept = exon_check[exon_check]
                fusion_exons = gene_info.filter(exons_kept.index, axis=0)
                cds_start = fusion_exons['CDS start'].tolist()
                cds_end = fusion_exons['CDS end'].tolist()
                seq_len = cds_end[-1] - cds_start[0] +1 # The plus one is to account for 1-indexing
            else:
                seq_len = 0

            # Now check if there was a truncation event in one of the included exons
            if any((exon_check == trunc_check) & exon_check):
                trunc_event = [exon_check == trunc_check]
                trunc_exon = gene_info.filter(exon_check.index[trunc_event[0]], axis = 0)
                rel_cds_start = trunc_exon['Genomic coding start'].tolist()[0]
                rel_cds_end = trunc_exon['Genomic coding end'].tolist()[0]
                
                # Default value of 0 because this also captures if the "truncated" exon is just the very first/last bp of the exon.
                trunc_diff = 0
                # Based on if the gene is in the 5' or 3' position we truncate ever so slightly differently relative to either end (5') or the start (3')
                if rel_cds_start < rel_pos and rel_cds_end > rel_pos:
                    if fusion_pos == 5:
                        trunc_diff =  rel_cds_end - rel_pos
                    elif fusion_pos ==3:
                        trunc_diff =  rel_pos - rel_cds_start           
            else:
                trunc_diff = 0
        
            # Now converting the nucleotide sequence length to an amino acid length with correction for 1-indexing
            if seq_len > 0:
                aa_len = (seq_len - trunc_diff)/3 - 1
            else:
                aa_len = 0

            # Putting the amino acid position in for ease in converting to the domains
            if fusion_pos == 5:
                aa_pos = math.ceil(aa_len)
            elif seq_len > 0:
                aa_pos = math.ceil(cds_end[-1]/3 - aa_len - 1)
            else:
                # This represents the cases where the UTRs are spliced together for the 3' end so the amino acid position of the fusion is the CDS end
                cds_end = gene_info['CDS end'].tolist()
                aa_pos = max(cds_end)
            
            self.__setattr__(var_prefix+'_aa_len',aa_len)
            self.__setattr__(var_prefix+'_aa_pos',aa_pos)

    def create_domain_arch(self):
    
        # Going through each of the fusion proteins that were annotated and finding what the theoretical domain architecture is
        fusion_arch = []
        brk_flag = False
        for fusion_pos, var_prefix in FUSION_POS_CONV.items():
            protOI = self.__getattribute__(var_prefix+'_uniprot')
            aa_pos = self.__getattribute__(var_prefix+'_aa_pos')
            prot_info = DANSY_REFERENCE_DF[DANSY_REFERENCE_DF['UniProt ID'] == protOI]
            prot_fusion_contrib_ids = []
                    
            # Check that the protein is in the reference file otherwise create a flag and break out of the loop
            if prot_info.empty:
                domains = []
                self.errors.append('1:Missing domain annotations')
                brk_flag = True
                prot_fusion_contrib_ids = None
                continue
            else:
                domains = prot_info['Interpro Domains'].str.split(';').tolist()[0]
            
            if domains and domains != ['']:
                for d in domains:
                    dom_info = d.split(':')
                    if fusion_pos == 5:
                        if int(dom_info[3]) <= aa_pos +AA_BUFFER:
                            fusion_arch.append(dom_info[1])
                            # Adding in the fusion protein contribution:
                            prot_fusion_contrib_ids.append(dom_info[1])
                    else:
                        if int(dom_info[2]) >= aa_pos - AA_BUFFER:
                            fusion_arch.append(dom_info[1])
                            prot_fusion_contrib_ids.append(dom_info[1])

            prot_fusion_contrib_ids = '|'.join(prot_fusion_contrib_ids)
            self.__setattr__(var_prefix+'_contribution', prot_fusion_contrib_ids)

        if brk_flag:
            self.domain_architecture = None
        else:
            self.domain_architecture =  '|'.join(fusion_arch)

    def get_domain_ngrams(self):
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
        all_domains = list(set(self.domain_architecture.split('|')))
        self.ngrams = dansy.ngramUtilities.return_all_n_grams_from_key(self.domain_architecture, all_domains)
    
    def _set_impact(self,cat):
        '''
        Internal function that should only be called from within the fusionCollection class.
        '''

        self.dansy_impact = cat

    def to_dict(self):
        '''
        Creates a dictionary version of the fusion information.
        '''
        
        attributes_4_output = ['id','name','h_gene', 't_gene','h_uniprot','t_uniprot','h_chr','t_chr','h_pos','t_pos','h_strand','t_strand', 'domain_architecture','dansy_impact'] + self.add_info_attribs
        return {k:v for k,v in self.__dict__.items() if k in attributes_4_output}

class fusionCollection():
    '''
    This class represents a collection of fusion genes that have had their domain architecture predicted and then are analyzed using DANSy networks to determine broad changes in the domain n-gram network.
    
    '''

    def __init__(self, x, mapper, add_cols=None, silent=False):
        if isinstance(x, pd.DataFrame):
            self.from_df(x,mapper=mapper,add_cols = add_cols, silent=silent)
        else:
            raise TypeError('Currently unsupported type has been provided.')

    def from_df(self,df:pd.DataFrame, mapper:dict,add_cols:list = None, silent = False):
        '''
        Builds the collection from a pandas DataFrame. 

        Parameters
        ----------
        - df: pandas DataFrame
            DataFrame containing all relevant information
        - mapper: dict
            Key-value pairs that map column names to parameters necessary to build fusionGene objects. Acceptable keys are: h_gene, t_gene, h_/t_chr, h_/t_pos, h_/t_strand, name, genome_build, id.
        - silent: bool
            Whether initial errors in fusion gene formation should result in an error or only log message.
        '''

        required_keys = ['name','h_gene','t_gene','h_chr','t_chr','h_pos', 't_pos']
        if set(required_keys).difference(mapper.keys()):
            raise ValueError('At least one required input for fusion gene domain architecture is missing.')
        
        optional_keys = ['h_strand','t_strand','genome_build', 'id']
        all_keys = required_keys+optional_keys
        param_mapper = {k:v for k,v in mapper.items() if k in all_keys}

        # If an id column has been provided double-checking that it is unique
        create_unique_ids = True
        formlen = math.ceil(np.log10(len(df)))
        if 'id' in param_mapper:
            duplicated_ids = df[param_mapper['id']].duplicated().any()
            if duplicated_ids:
                print('Duplicated IDs detected in the provided ID column will create a new set of IDs')
            else:
                create_unique_ids = False

        # Now create the collection and create a dict for each unique fusion name
        self.fusion_list = {}
        self.fusion_name_id = {k:[] for k in df[param_mapper['name']].unique()}
        pflag = len(df) > 1000
        if pflag:
            pbar = tqdm(total=len(df))

        for idx, row in df.iterrows():
            if create_unique_ids:
                temp = f'{idx}'.zfill(formlen)
                id = 'DANSY_'+temp
            else:
                id = row[param_mapper['id']]
            
            inputs = {k:row[v] for k,v in param_mapper.items()}
            if add_cols is not None:
                add_info = {k:row[k] for k in add_cols}
            
            f = fusionGene(**inputs, id = id,add_info=add_info, silent=silent)
            self.fusion_list[id] = f
            self.fusion_name_id[f.name].append(id)

            if pflag:
                pbar.update()
        
        if pflag:
            pbar.close()

    def get_fusion(self, f):
        '''
        Retrieves the fusion with the provided ID.
        '''

        return self.fusion_list[f]  
    
    def return_fusion_ids(self, f):
        '''
        Returns all the ids that a specific fusion name has.
        '''
        return self.fusion_name_id[f]
        
    def get_dansy_impacts(self, dansy_obj, debug_flag = False):

        # First go through and get the unique domain architectures wihtin the collection and get lists of fusions for each
        uni_archs = {}
        for f in self.fusion_list.values():
            if f.domain_architecture not in uni_archs:
                uni_archs[f.domain_architecture] = [f.id]
            else:
                uni_archs[f.domain_architecture].append(f.id)

        # Now go through and check the impact of each individual fusion.
        gc_dict = {}
        cnt = 0
        pflag = len(uni_archs.keys()) > 1000
        if pflag:
            pbar = tqdm(total=len(uni_archs.keys()))
        for arch, reps in uni_archs.items():
            if arch is not None:
                gc,skipped = self.calc_network_change(arch, dansy_obj, reps[0], debug_flag)
                gc_dict[arch] = gc

                if skipped:
                    cnt += 1
            if pflag:
                pbar.update()

        if pflag:
            pbar.close()

        print(f'There were {cnt} fusion domain architectures previously found in the proteome.')
        self.gc_dict = gc_dict

    def calc_network_change(self, arch,proteome_dansy:dansy.dansy,rep:str, debug_flag = False):
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
        num_cats = ['Connected Components','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges']
        if debug_flag:
            gross_changes = pd.Series(index = ['Connected Components','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges','New n-grams', 'Reintroduced N-grams List'])
        else:
            gross_changes = pd.Series(index = ['Connected Components','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges'],dtype=int).fillna(0)
        skipped = False # Keeping count of how many are skipped for the entire dataset

            
        if arch in existing_ngrams or arch == '':
            gross_changes[num_cats] = 0
            skipped = True
        
        else:
            
            # Checking if any of the n-grams have been removed previously
            f = self.get_fusion(rep)
            fusion_ngrams = f.ngrams

            # Here making sure we only focus on n-grams up to the maximum length of the dansy network
            fusion_ngrams = [f for f in fusion_ngrams if len(f.split('|')) <= proteome_dansy.n]

            ngram_to_return = list(set(fusion_ngrams).intersection(removed_ngrams))
            new_ngrams = set(fusion_ngrams).difference(existing_ngrams).difference(ngram_to_return)
            ngrams_to_remove = []
            
            # Checking for any n-grams that can be subsumed by a longer n-gram and then removing it
            for gram in ngram_to_return:
                for inner_gram in ngram_to_return:
                    if gram != inner_gram and gram in inner_gram and len(gram.split('|')) > 1:
                        ngrams_to_remove.append(gram)

            if ngrams_to_remove:
                ngrams_to_remove = set(ngrams_to_remove)
                ngram_to_return = list(set(ngram_to_return).difference(ngrams_to_remove))
                #ngram_to_return = [x for x in ngram_to_return if x not in ngrams_to_remove]
            
            # Saving the easy gross topological changes 
            gross_changes['Reintroduced N-grams'] = len(ngram_to_return)
            gross_changes['New N-gram Count'] = len(new_ngrams)
            
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
            
            if new_ngrams_orig: # If there were any new n-grams going through and checking all the fusion n-grams to determine which connected components were impacted.
                
                # N-grams to check focusing only on the connected components that the novel n-grams are within
                cc_impacted = 0
                for i in nx.connected_components(proteome_dansy.G):
                    overlap = set(fusion_ngrams).intersection(i)
                    if len(overlap) > 0:
                        
                        cc_impacted += 1
                    else:
                        # Check if any of the reintroduced n-grams are part of longer ones in the connected component
                        found_flag = False
                        for ngram in ngram_to_return:
                            if any(ngram in j for j in i):
                                found_flag = True
                        if found_flag:
                            cc_impacted +=1 

                cc_val = cc_impacted - 1
            
            elif ngram_to_return: # As a special case there can be only reintroduced n-grams so to designate that using a value that should not be identified with -1.
                cc_val = -1

            else:
                cc_val = 0
               
            # Outputting the changes
            gross_changes['Connected Components'] = cc_val
            gross_changes['New Nodes'] = len(new_ngrams)
            gross_changes['New Edges'] = new_edge

            if debug_flag:
                gross_changes['New n-grams'] = list(new_ngrams)
                gross_changes['Reintroduced N-grams List'] = list(ngram_to_return)

        return gross_changes, skipped
    
    def categorize_dansy_impact(self):

        self.fusion_arch_impacts = {}
        for arch, gc in self.gc_dict.items():
            if arch == '':
                cat = 'No Annotation'
            elif arch is None:
                cat = None
            else:
                if gc[['Connected Components','New N-gram Count','Reintroduced N-grams', 'New Nodes','New Edges']].sum() == 0:
                    cat = 'No Change'
                elif gc['Connected Components'] == 0:
                    cat = 'Reinforcement'
                elif gc['Connected Components'] > 0:
                    cat = 'Bridges'
                elif gc['Reintroduced N-grams'] > 0 and gc['Connected Components'] == -1:
                    cat = 'Reintroduction'
                else:
                    cat = 'Double-check'
            
            
            self.fusion_arch_impacts[arch] = cat

    def propagate_fusion_impacts(self):

        for f in self.fusion_list.values():
            if f.domain_architecture is not None:
                f._set_impact(self.fusion_arch_impacts[f.domain_architecture])
            else:
                f._set_impact(None)
                

    def perform_dansy_analysis(self, dansy_obj, debug_flag = False):
        self.get_dansy_impacts(dansy_obj=dansy_obj,debug_flag=debug_flag)
        self.categorize_dansy_impact()
        self.propagate_fusion_impacts()
        
    def summarize_dansy_results(self):

        x = []
        for f in self.fusion_list.values():
            x.append(f.to_dict())
        
        return pd.DataFrame(x)
        

            
