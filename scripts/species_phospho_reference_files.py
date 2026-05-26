from CoDIAC import InterPro, UniProt
import CoDIAC
import pandas as pd
import requests
import gzip
from Bio import SeqIO

from pybiomart import Dataset, Server

servers = Server(host='http://www.ensembl.org')
mart = servers['ENSEMBL_MART_ENSEMBL'] 

# pTyr System
py_list = ['IPR000980','IPR000242','IPR006020','IPR020635']

#pST System
pst_list = ['IPR001245','IPR006186','IPR001932','IPR004274','IPR001202','IPR023410','IPR000253','IPR001132','IPR001357','IPR002713', 'IPR000719','IPR000340']

fetch_data = '2026_0324'

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

species_information = {'Homo sapiens':['H_sapiens','UP000005640',9606],
                 'Mus musculus':['M_musculus','UP000000589',10090],
                 'Rattus norvegicus':['R_norvegicus','UP000002494',10116],
                 'Oryctolagus cuniculus':['O_cuniculus','UP000001811',9986],
                 'Gallus gallus':['G_gallus','UP000000539',9031],
                 'Xenopus tropicalis':['X_tropicalis','UP000008143',8364],
                 'Danio rerio':['D_rerio','UP000000437',7955],
                 'Carassius auratus':['C_auratus','UP000515129',7957],
                 'Ciona intestinalis':['C_intestinalis','UP000008144',7719],
                 'Strongylocentrotus purpuratus':['S_purpuratus','UP000007110',7668],
                 'Drosophila melanogaster':['D_melanogaster','UP000000803',7227],
                 'Caenorhabditis elegans':['C_elegans', 'UP000001940', 6239],
                 'Nematostella vectensis':['N_vectensis','UP000001593',45351],
                 'Trichoplax adhaerens':['T_adhaerens','UP000009022',10228],
                 'Monosiga brevicollis':['M_brevicollis','UP000001357',81824],
                 'Capsaspora owczarzaki (strain ATCC 30864)':['C_owczarzaki','UP000008743',595528],
                 'Sphaeroforma arctica JP610':['S_arctica','UP000054560',667725],
                 'Acanthamoeba castellanii (strain ATCC 30010 / Neff)':['A_castellanii','UP000011083',1257118],
                 'Dictyostelium discoideum':['D_discoideum','UP000002195',44689],
                 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)':['S_cerevisiae','UP000002311',559292]}

reviewed_species = ['Homo sapiens','Mus musculus','Rattus norvegicus',]

# Now fetching the and generating the files
for id_list, ptm_sys in zip([py_list,pst_list],['pTyr','pSerThr']):
    for species,spec_name in species_names.items():
        print('Starting fetching of IDs for %s'%species)
        
        # First let's load the species reference proteome as we will use that to reduce the uniprot ID list to the ones with highest confidence
        full_proteome = set()
        proteome_id = species_information[species][1]
        taxon_id = species_information[species][2]
        with gzip.open(f'data/reference_proteomes/{proteome_id}_{taxon_id}.fasta.gz','rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                full_proteome.add(record.id.split('|')[1])
        
        # Now getting the complete set of UniProt IDs
        comp_uniprots = set()
        
        for Interpro_ID in id_list:
            
            rev_flag =  species in reviewed_species
            if species not in special_case:
                uniprot_IDs, uniprot_dict = CoDIAC.InterPro.fetch_uniprotids(Interpro_ID, REVIEWED=rev_flag, species=species)
            else:
                spec_search = special_case[species]
                uniprot_IDs, uniprot_dict = CoDIAC.InterPro.fetch_uniprotids(Interpro_ID, REVIEWED=rev_flag, species=spec_search)
        
            uniprot_IDs = [v for v in uniprot_dict[species].values()]
            uniprot_IDs = sum(uniprot_IDs,[])

            comp_uniprots.update(uniprot_IDs)
            ids_removed = comp_uniprots.difference(full_proteome)

            if len(ids_removed) > 0:
                print('A total of %d will be removed from the original ID list length of %d'%(len(ids_removed),len(uniprot_IDs)))

            # Removing the redundant uniprot IDs
            comp_uniprots = comp_uniprots.intersection(full_proteome)

            
        print('Generating the Reference File')
        reference_File = 'data/Current_Multispecies_Files/'+spec_name+'Full_'+ptm_sys+'_System_Reference_File_'+fetch_data+'.csv'
        _ = CoDIAC.UniProt.makeRefFile(comp_uniprots, reference_File)