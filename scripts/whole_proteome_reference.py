from CoDIAC import UniProt,InterPro
import CoDIAC
import pandas as pd
from pybiomart import Dataset, Server
import os

# Change these according to current build of each reference file as desired.
file_suffix = 'Reference_File_MMDD_YYYY.csv' # Change this to the current date
gencode = pd.read_table('gencode.v47.metadata.SwissProt',header=None) # Ensure this is the gencode version file that has been downloaded

# Change these if you want slightly different behavior
single_file = True # Change if you want to spread the proteome reference files across multiple files
method = 'gencode' # We recommend gencode due to stability in cross references but can also use pybiomart

# If pybiomart is to used these are the commands that will retrieve UniProt IDs
servers = Server(host='http://www.ensembl.org')
mart = servers['ENSEMBL_MART_ENSEMBL'] 
dataset = Dataset(host = 'http://useast.ensembl.org', name='hsapiens_gene_ensembl')
gene_cnv = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name','uniprotswissprot'])
gene_cnv.dropna(axis=0, inplace=True,subset=['UniProtKB/Swiss-Prot ID'])


# Current Interpro File Data Folder
data_folder = 'data/Current_Human_Proteome/'

# Ensuring if there are existing reference files that they are not overwritten and finding which IDs have not been fetched.
current_uniprots = set()
full_dir = os.listdir(data_folder)
for file in full_dir:
    if file.endswith(file_suffix):
        temp_df = pd.read_csv(data_folder+file)
        temp_uniprot = temp_df['UniProt ID'].tolist()
        current_uniprots.update(temp_uniprot)

# There are two potential ways to retrieve the UniProt IDs one is using the pybiomart package and the other is using data downloaded from GenCode
if method == 'pybiomart':
    uniprot_ids = gene_cnv['UniProtKB/Swiss-Prot ID'].tolist()
elif method == 'gencode':
    uniprot_ids = gencode[1].tolist()
else:
    raise ValueError()

uniprot_ids = list(set(uniprot_ids).difference(current_uniprots))
print('Out of %s there are a total of %s left'%(len(set(uniprot_ids))+len(current_uniprots), len(set(uniprot_ids))))

# Double-check in case there was a disruption during reference file generation.
cnt = 0
hs_prot_dir = os.listdir(data_folder)
for file in hs_prot_dir:
    if file.endswith(file_suffix):
        file_str = file.split('_')
        cnt_chck = int(file_str[1])
        if cnt_chck > cnt:
            cnt = cnt_chck

# Second double-check for complete overlap between two reference files as there was a re-run with the same uniprot ID list
rm_files = []
for file in hs_prot_dir :
    if file.endswith(file_suffix) and file not in rm_files:
        file_content_1 = pd.read_csv(data_folder+file)
        id_list1 = file_content_1['UniProt ID'].tolist()
        for file_2 in hs_prot_dir:
            if file_2.endswith(file_suffix):
                if file_2 != file and file_2 not in rm_files:
                    file_content_2 = pd.read_csv(data_folder+file_2)
                    id_list2 = file_content_2['UniProt ID'].tolist()
                    if set(id_list1).difference(id_list2) == set():
                        print('Found a match for %s'%file)
                        rm_files.append(file_2)
if rm_files:
    cnt += 1

if single_file:
    # Final set of IDs that are missing
    ref_file_name = data_folder+'Proteome_'+str(cnt)+'_'+file_suffix
    uniprot_df = CoDIAC.UniProt.makeRefFile(uniprot_ids, ref_file_name)
else:
    # Note for everything below the 300 IDs is once UniProt starts to throw out exceptions.
    currently_retrieved = set()
    for x in range(0,round(len(uniprot_ids)/300)):
        a = uniprot_ids[x*300:300*(x+1)]
        ref_file_name = data_folder+'Proteome_'+str(cnt)+'_'+file_suffix
        uniprot_df = CoDIAC.UniProt.makeRefFile(a, ref_file_name)
        cnt += 1
        currently_retrieved.update(a)

    # Final set of IDs that are missing
    a = set(uniprot_ids).difference(currently_retrieved)
    ref_file_name = data_folder+'Proteome_'+str(cnt)+'_'+file_suffix
    uniprot_df = CoDIAC.UniProt.makeRefFile(a, ref_file_name)
