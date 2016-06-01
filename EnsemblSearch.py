# Trying to use purely Ensembl Compara database to search young duplicated genes
# inspired by MBE paper "Accelerated Evolution after Gene Duplication:A Time-Dependent Process Affecting Just One Copy"
# Xiang Ji
# xji3@ncsu.edu

import csv

# http://www.ensembl.info/blog/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
# http://useast.ensembl.org/biomart/martview/
# Details
# Database : Ensembl Genes 84
# Dataset  : Homo sapiens genes (CRCh38.p5)
# Filter   : Homologue filters, only paralogous human genes
# File name: Ensembl84_Human_Paralogs.txt


if __name__ == '__main__':
    Human_paralog_file = '/Users/xji3/Downloads/Ensembl84_Human_Paralogs.txt'

    visited_genes = set()
    gene_families = list()
    
    with open(Human_paralog_file, 'rb') as f:
        paralog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')

        last_Ensembl_id = ''
        current_family = list()
        for row in reader:
            combined_gene_name = '_'.join(row[0], row[4])
            if row[-2] == '1':
                if row[0] in visited_genes:
                    continue
                
                if row[0] == last_Ensembl_id:
                    current_family.append(row[5])
                else:
                    gene_families.append(set(current_family))
                    for gene in current_family:
                        visited_genes.add(gene)
                    current_family = [row[0], row[5]]
                    last_Ensembl_id = row[0]
        
        
    # Now subtract all two paralog families
    two_paralog_families = [family for family in gene_families if len(family) == 2]


#####################################################################################################
#####################################################################################################

    # Now build ortholog mapping
    Human_primate_ortholog_file = '/Users/xji3/Downloads/Ensembl84_Human_Primate_Orthologs.txt'

    ortholog_mapping = dict()

    species_list = ['Chimpanzee', 'Gibbon', 'Gorilla', 'Macaque', 'Olive baboon', 'Orangutan']
    entry_number = [1, 4, 7, 10, 13, 17]

    with open(Human_primate_ortholog_file, 'rb') as f:
        ortholog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')

        for row in reader:
            ensembl_id = row[0]
            if ensembl_id in visited_genes:
                ortholog_mapping[ensembl_id] = {species_list[i]:row[entry_number[i]] for i in range(len(species_list))}


    Human_mouse_ortholog_file = '/Users/xji3/Downloads/Ensembl84_Human_Mouse_Orthologs.txt'

    with open( Human_mouse_ortholog_file, 'rb') as f:
        ortholog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')

        for row in reader:
            ensembl_id = row[0]
            if ensembl_id in ortholog_mapping:
                ortholog_mapping[ensembl_id]['Mouse'] = row[2]


#####################################################################################################
#####################################################################################################

    # Now select families that have only one paralog in mouse
    Mouse_outgroup_two_paralog_families = [family for family in two_paralog_families if sum([ortholog_mapping[i].has_key('Mouse') for i in family]) == 1]
