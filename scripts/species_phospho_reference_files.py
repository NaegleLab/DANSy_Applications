from CoDIAC import InterPro, UniProt
import CoDIAC
import pandas as pd
import requests

from pybiomart import Dataset, Server

servers = Server(host='http://www.ensembl.org')
mart = servers['ENSEMBL_MART_ENSEMBL'] 

# pTyr System
py_list = ['IPR000980','IPR000242','IPR006020','IPR020635']

#pST System
pst_list = ['IPR001245','IPR006186','IPR001932','IPR004274','IPR001202','IPR023410','IPR000253','IPR001132','IPR001357','IPR002713']

fetch_data = 'May12_2025'

species_names = {'Homo sapiens':'H_sapiens',
                 'Mus musculus':'M_musculus',
                 'Rattus norvegicus':'R_norvegicus',
                 'Oryctolagus cuniculus':'O_cuniculus',
                 'Gallus gallus':'G_gallus',
                 'Xenopus tropicalis':'X_tropicalis',
                 'Danio rerio':'D_rerio',
                 'Carassius auratus':'C_auratus',
                 'Ciona intestinalis':'C_intestinalis',
                 'Strongylocentrotus purpuratus':'S_purpuratus',
                 'Drosophila melanogaster':'D_melanogaster',
                 'Caenorhabditis elegans':'C_elegans', 
                 'Nematostella vectensis':'N_vectensis',
                 'Trichoplax adhaerens':'T_adhaerens',
                 'Monosiga brevicollis':'M_brevicollis',
                 'Capsaspora owczarzaki (strain ATCC 30864)':'C_owczarzaki',
                 'Sphaeroforma arctica JP610':'S_arctica',
                 'Acanthamoeba castellanii (strain ATCC 30010 / Neff)':'A_castellanii',
                 'Dictyostelium discoideum':'D_discoideum',
                 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)':'S_cerevisiae'}

special_case = {'Capsaspora owczarzaki (strain ATCC 30864)':'Capsaspora owczarzaki',
                'Acanthamoeba castellanii (strain ATCC 30010 / Neff)':'Acanthamoeba castellanii',
                'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)':'Saccharomyces cerevisiae'}

reviewed_species = ['Homo sapiens','Mus musculus','Rattus norvegicus',]

# Now fetching the and generating the files
for id_list, ptm_sys in zip([py_list,pst_list],['pTyr','pSerThr']):
    for species,spec_name in species_names.items():
        print('Starting fetching of IDs for %s'%species)
        comp_uniprots = set()
        for Interpro_ID in id_list:
            rev_flag =  species in reviewed_species
            if species not in special_case:
                uniprot_IDs, uniprot_dict = CoDIAC.InterPro.fetch_uniprotids(Interpro_ID, REVIEWED=rev_flag, species=species)
            else:
                spec_search = special_case[species]
                uniprot_IDs, uniprot_dict = CoDIAC.InterPro.fetch_uniprotids(Interpro_ID, REVIEWED=rev_flag, species=spec_search)
            
            
            # If a species is not one that has well reviewed UniProt entries go through and start cleaning up the IDs
            if not rev_flag:
                
                uniprot_IDs = [v for v in uniprot_dict[species].values()]
                uniprot_IDs = sum(uniprot_IDs,[])


                uniprot_cleaned = {}
                ids_to_remove = []

                for ID in uniprot_IDs:
                    url = f"https://rest.uniprot.org/uniprotkb/{ID}?fields=accession,gene_primary,protein_existence,xref_ensembl,xref_refseq,length,xref_geneid"
                    get_url = requests.get(url)
                    response = get_url.json()

                    # First grabbing a few key values
                    review = response['entryType']
                    prot_exist = int(response['proteinExistence'][0]) #First character is always the number
                    refseq = False
                    ncbi_geneID = False
                    ensembl = []

                    # For all UniProt IDs an ENSEMBL, RefSeq or GeneID entry must exist as those are the entries associated with actual genes while others will be redundant entries corresponding to mRNA molecules that may be partial CDS

                    if response['uniProtKBCrossReferences']:
                    
                        for i in range(len(response['uniProtKBCrossReferences'])):
                            if response['uniProtKBCrossReferences'][i]['database'] == 'RefSeq':
                                refseq = True
                            
                            if response['uniProtKBCrossReferences'][i]['database'] == 'Ensembl':
                                if response['uniProtKBCrossReferences'][i]['properties'][1]['value']:
                                    ensembl = response['uniProtKBCrossReferences'][i]['properties'][1]['value']
                            
                            if response['uniProtKBCrossReferences'][i]['database'] == 'GeneID':
                                ncbi_geneID = True

                    if refseq or ncbi_geneID or ensembl or review == 'UniProtKB reviewed (Swiss-Prot)':

                        # If there was not an ensembl ID then changing it to the gene name
                        if not ensembl:
                            if 'genes' in response.keys():
                                try:
                                    ensembl = response['genes'][0]['geneName']['value']
                                except:
                                    ensembl = ID
                            else:
                                ensembl = ID
                            
                        entry_len = response['sequence']['length']
                        temp_dict = {'review':review,'existence':prot_exist, 'refseq status':refseq, 'protein length':entry_len, 'ID':ID}

                        # Setting up a comparison if the ensembl ID does not already exist
                        if ensembl in uniprot_cleaned:
                            # first checking review status which takes ultimate precedence
                            if uniprot_cleaned[ensembl]['review'] != 'UniProtKB reviewed (Swiss-Prot)':
                                if review == 'UniProtKB reviewed (Swiss-Prot)':
                                    ids_to_remove.append(uniprot_cleaned[ensembl]['ID'])
                                    uniprot_cleaned[ensembl] = temp_dict.copy()
                                elif prot_exist < uniprot_cleaned[ensembl]['existence']:
                                    ids_to_remove.append(uniprot_cleaned[ensembl]['ID'])
                                    uniprot_cleaned[ensembl] = temp_dict.copy()
                                elif refseq and not uniprot_cleaned[ensembl]['refseq status']:
                                    ids_to_remove.append(uniprot_cleaned[ensembl]['ID'])
                                    uniprot_cleaned[ensembl] = temp_dict.copy()
                                elif entry_len > uniprot_cleaned[ensembl]['protein length']:
                                    ids_to_remove.append(uniprot_cleaned[ensembl]['ID'])
                                    uniprot_cleaned[ensembl] = temp_dict.copy()
                                else:
                                    ids_to_remove.append(ID)
                            else:
                                ids_to_remove.append(ID)
                                        
                        else:
                            uniprot_cleaned[ensembl] = temp_dict.copy()
                    else:
                        ids_to_remove.append(ID)

                if len(ids_to_remove) > 0:
                    print('A total of %d will be removed from the original ID list length of %d'%(len(ids_to_remove),len(uniprot_IDs)))

                # Removing the redundant uniprot IDs
                nonredundant_IDs = [id for id in uniprot_IDs if id not in ids_to_remove]
                comp_uniprots.update(nonredundant_IDs)
            else:
                comp_uniprots.update(uniprot_IDs)
            
        print('Generating the Reference File')
        reference_File = 'data/Current_Multispecies_Files/'+spec_name+'Full_'+ptm_sys+'_System_Reference_File_'+fetch_data+'.csv'
        _ = CoDIAC.UniProt.makeRefFile(comp_uniprots, reference_File)