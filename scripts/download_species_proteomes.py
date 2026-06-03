import subprocess

# This is fetching the reference proteomes for the species of interest in the DANSy PTM analysis
fetch_date = '2026_0324'

# This is the uniprot ftp server that we are getting the canonical, reference proteomes from which have 1 UniProt ID per gene
ftp_server = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/'
local_dir = 'data/reference_proteomes/'

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

for species,info in species_information.items():
    spec_name = info[0]
    ref_proteome = info[1]
    taxon_id = info[2]
    species_portion = f'{ref_proteome}/{ref_proteome}_{taxon_id}.fasta.gz'
    try:
        full_url = ftp_server+species_portion
        wget_args = ['wget', full_url, '-P',local_dir,'-t','10','-w','30']
        subprocess.run(wget_args,check=True)
        print(f'Successfully fetched the reference proteome for {spec_name}')
    except subprocess.CalledProcessError as e:
        print(f'Issue with {spec_name} the error was: {e.output}')
