# Trying to use purely Ensembl Compara database to search young duplicated genes
# inspired by MBE paper "Accelerated Evolution after Gene Duplication:A Time-Dependent Process Affecting Just One Copy"
# Xiang Ji
# xji3@ncsu.edu

import csv, os, sys, requests, numpy

# http://www.ensembl.info/blog/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
# http://useast.ensembl.org/biomart/martview/
# Details
# Database : Ensembl Genes 84
# Dataset  : Homo sapiens genes (CRCh38.p5)
# Filter   : Homologue filters, only paralogous human genes
# File name: Ensembl84_Human_Paralogs.txt


# about query cds using Ensembl id
# http://rest.ensembl.org/documentation/info/sequence_id

def getENSsequence(ENS_id, seq_type = 'cds', use_ccds = False):
    server = "http://rest.ensembl.org"
    if use_ccds:
        #https://rest.ensembl.org/documentation/info/xref_id
        ext_ccds = '/xrefs/id/' + ENS_id + '?external_db=CCDS;all_levels=1'
        r = requests.get(server+ext_ccds, headers={ "Content-Type" : "application/json"})
        ccds_id = str(r.json()[0][u'display_id'])

        ext = '/sequence/id/' + ccds_id + '?object_type=transcript;db_type=otherfeatures;type=cds;species=human'
    else:
        ext = "/sequence/id/" + ENS_id + '?multiple_sequences=1;type=' + seq_type
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    return str(r.text)

def outputFamilies(families, output_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list):
    for family in families:
        paralog1, paralog2 = list(family)
        for outgroup in outgroup_species_list:
            if len(ortholog_mapping[paralog2][outgroup]) > 1:
                trans = paralog1
                paralog1 = paralog2
                paralog2 = trans
        family_folder = output_folder + '_'.join(family) + '/'
        if not os.path.isdir(family_folder):
            os.mkdir(family_folder)

        with open(family_folder + '_'.join(family) + '.fa', 'w+') as f:
            for outgroup in outgroup_species_list:                    
                f.write('>' + outgroup.replace(' ', '_') + paralog1 + '\n')
                ENS_id = ortholog_mapping[paralog1][outgroup]
                seq = getENSsequence(ENS_id)
                seq = seq.split('\n')[0]                
                f.write(seq + '\n')
            for paralog in [paralog1, paralog2]:
                f.write('>Human' + paralog + '\n')
                try:
                    seq = getENSsequence(paralog, 'CDS', True)
                except:
                    seq = getENSsequence(paralog, 'cds', False)
                if paralog == 'ENSG00000163602':
                    seq = 'ATG' + seq[3:]

                seq = seq.split('\n')[0]
                f.write(seq + '\n')
                for species in two_paralog_species_list:
                    f.write('>' + species.replace(' ', '_') + paralog + '\n')
                    ENS_id = ortholog_mapping[paralog][species]
                    
                    seq = getENSsequence(ENS_id)
                    seq = seq.split('\n')[0]
                    f.write(seq + '\n')


if __name__ == '__main__':
    home_path = '/Users/xji3/'
    #home_path = '/Users/Xiang/'
    Human_paralog_file = home_path + 'Downloads/Ensembl84_Human_Paralogs.txt'

    visited_genes = set()
    gene_families = list()
    
    with open(Human_paralog_file, 'rb') as f:
        paralog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')

        last_Ensembl_id = ''
        current_family = list()
        for row in reader:
            combined_gene_name = '_'.join([row[0], row[4]])
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

    # Now load gene ids that locate in segmental duplication regions
    sg_dup_gene_ids = numpy.loadtxt(home_path + 'GitFolders/PrimateMultigeneFamilySearch/Segmental_Dup_Gene_ids.txt', dtype = str)
    sg_dup_gene_ids = [i[1:-1] for i in sg_dup_gene_ids]


#####################################################################################################
#####################################################################################################

    # Now build ortholog mapping
    Human_primate_ortholog_file = home_path + 'Downloads/Ensembl84_Human_Primate_Orthologs.txt'

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


    Human_mouse_ortholog_file = home_path + 'Downloads/Ensembl84_Human_Mouse_Orthologs.txt'

    with open( Human_mouse_ortholog_file, 'rb') as f:
        ortholog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')

        for row in reader:
            ensembl_id = row[0]
            if ensembl_id in ortholog_mapping:
                ortholog_mapping[ensembl_id]['Mouse'] = row[2]

    # now add mouse key to all genes
    for ensembl_id in ortholog_mapping:
        if not ortholog_mapping[ensembl_id].has_key('Mouse'):
            ortholog_mapping[ensembl_id]['Mouse'] = ''

    # Dog, Cow ortholog mapping
    Human_cow_dog_ortholog_file = home_path + 'Downloads/Ensembl84_Human_Cow_Dog_Orthologs.txt'
    with open( Human_cow_dog_ortholog_file, 'rb') as f:
        ortholog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')

        for row in reader:
            ensembl_id = row[0]
            if ensembl_id in ortholog_mapping:
                ortholog_mapping[ensembl_id]['Cow'] = row[2]
                ortholog_mapping[ensembl_id]['Dog'] = row[5]
    for ensembl_id in ortholog_mapping:
        if not ortholog_mapping[ensembl_id].has_key('Cow'):
            ortholog_mapping[ensembl_id]['Cow'] = ''
        if not ortholog_mapping[ensembl_id].has_key('Dog'):
            ortholog_mapping[ensembl_id]['Dog'] = ''


    # Bushbaby, Marmoset, Mouse Lemur, Shrew ortholog mapping
    species_list = ['Bushbaby', 'Marmoset', 'Mouse Lemur', 'Shrew']
    entry_number = [2, 5, 8, 11]
    Human_4primate_ortholog_file = home_path + 'Downloads/Ensembl84_Human_4primate_Orthologs.txt'
    with open(Human_4primate_ortholog_file, 'rb') as f:
        ortholog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')
        
        for row in reader:
            ensembl_id = row[0]
            if ensembl_id in visited_genes:
                for i in range(len(species_list)):
                    ortholog_mapping[ensembl_id][species_list[i]] = row[entry_number[i]]

    # Tarsier, Tree Shrew, Vervet ortholog mapping
    species_list = ['Tarsier', 'Tree Shrew', 'Vervet']
    entry_number = [2, 5, 8]
    Human_3primate_ortholog_file = home_path + 'Downloads/Ensembl84_Human_3primate_Orthologs.txt'
    with open(Human_3primate_ortholog_file, 'rb') as f:
        ortholog_labels = f.readline()[:-1].replace('"', '').replace('[0 low, 1 high]', '').split(',')
        reader = csv.reader(f, delimiter = ',', quotechar = '"')
        
        for row in reader:
            ensembl_id = row[0]
            if ensembl_id in visited_genes:
                for i in range(len(species_list)):
                    ortholog_mapping[ensembl_id][species_list[i]] = row[entry_number[i]]

                
#####################################################################################################
#####################################################################################################

    # Now select families that have only one paralog in mouse
    Mouse_outgroup_two_paralog_families = [family for family in two_paralog_families if sum([ortholog_mapping[i]['Mouse'] != '' for i in family]) == 1]


    # Now restric Chimpanzee, Gorilla, Orangutan, Gibbon to have both copies
    Post_Mouse_families = [family for family in Mouse_outgroup_two_paralog_families
                            if sum([ortholog_mapping[i]['Chimpanzee'] != '' for i in family]) == 2
                            and sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                            and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                            and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 2]
    # exlude one family of small RNA
    Post_Mouse_families.remove(set(['ENSG00000238578', 'ENSG00000238597']))


    # Decided to remove Bushbaby
    Post_Lemur_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Vervet'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Tarsier'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Mouse Lemur'] != '' for i in family]) == 1
                           ]    
    print 'Post Lemur families: ', len(Post_Lemur_families), len([family for family in Post_Lemur_families if sum([ens in sg_dup_gene_ids for ens in family])])


    Post_Tarsier_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Vervet'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Tarsier'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Mouse Lemur'] != '' for i in family]) == 1
                           ]
    same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Mouse Lemur'] != '' for i in family])]
    Post_Tarsier_families = [family for family in Post_Tarsier_families if family in same_ortholog]
    Post_Tarsier_families.remove(set(['ENSG00000124788', 'ENSG00000224470'])) # Gibbon seq % 3 != 0
    print 'Post Tarsier families: ', len(Post_Tarsier_families), len([family for family in Post_Tarsier_families if sum([ens in sg_dup_gene_ids for ens in family])])
    
    Post_NWM_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Vervet'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Tarsier'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Mouse Lemur'] != '' for i in family]) == 1
                           ]
    same_ortholog = [family for family in Post_NWM_families
                     if sum([ortholog_mapping[i]['Tarsier'] != ''
                             and ortholog_mapping[i]['Mouse Lemur'] != ''
                             and ortholog_mapping[i]['Marmoset'] != ''
                             for i in family])]
    Post_NWM_families = [family for family in Post_NWM_families if family in same_ortholog]
    print 'Post NWM families: ', len(Post_NWM_families), len([family for family in Post_NWM_families if sum([ens in sg_dup_gene_ids for ens in family])])

    Post_OWM_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Vervet'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Tarsier'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Mouse Lemur'] != '' for i in family]) == 1
                           ]
    print 'Post OWM families: ', len(Post_OWM_families), len([family for family in Post_OWM_families if sum([ens in sg_dup_gene_ids for ens in family])])

    Post_Gibbon_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Vervet'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Tarsier'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Mouse Lemur'] != '' for i in family]) == 1
                           ]
    print 'Post Gibbon families: ', len(Post_Gibbon_families), len([family for family in Post_Gibbon_families if sum([ens in sg_dup_gene_ids for ens in family])])


# Note: EDN_ECP is in this family ['ENSG00000129538', 'ENSG00000169385', 'ENSG00000169397', 'ENSG00000165799', 'ENSG00000258818', 'ENSG00000169413', 'ENSG00000214274', 'ENSG00000173431']

#####################################################################################################
#####################################################################################################


    Ensembl_GF = './Ensembl_GeneFamilies/'
    if not os.path.isdir(Ensembl_GF):
        os.mkdir(Ensembl_GF)

    Post_Mouse_folder = Ensembl_GF + 'Post_Mouse/'
    if not os.path.isdir(Post_Mouse_folder):
        os.mkdir(Post_Mouse_folder)

    Post_Lemur_folder = Ensembl_GF + 'Post_Lemur/'
    if not os.path.isdir(Post_Lemur_folder):
        os.mkdir(Post_Lemur_folder)
        
    Post_Tarsier_folder = Ensembl_GF + 'Post_Tarsier/'
    if not os.path.isdir(Post_Tarsier_folder):
        os.mkdir(Post_Tarsier_folder)

    Post_NWM_folder = Ensembl_GF + 'Post_NWM/'
    if not os.path.isdir(Post_NWM_folder):
        os.mkdir(Post_NWM_folder)

    Post_OWM_folder = Ensembl_GF + 'Post_OWM/'
    if not os.path.isdir(Post_OWM_folder):
        os.mkdir(Post_OWM_folder)

    # output Post_Lemur list
    with open(Post_Lemur_folder + 'Post_Lemur_List.txt', 'w+') as f:
        for family in Post_Lemur_families:
            f.write('_'.join(family) + '\n')
            
    # output Post_Tarsier list
    with open(Post_Tarsier_folder + 'Post_Tarsier_List.txt', 'w+') as f:
        for family in Post_Tarsier_families:
            f.write('_'.join(family) + '\n')

    # output Post_NWM list
    with open(Post_NWM_folder + 'Post_NWM_List.txt', 'w+') as f:
        for family in Post_NWM_families:
            f.write('_'.join(family) + '\n')

    # output Post_OWM list
    with open(Post_OWM_folder + 'Post_OWM_List.txt', 'w+') as f:
        for family in Post_OWM_families:
            f.write('_'.join(family) + '\n')

    # Post Lemur
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Macaque',
                                'Olive baboon', 'Orangutan', 'Marmoset', 'Tarsier']
    outgroup_species_list = ['Mouse Lemur']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_Lemur_families, Post_Lemur_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    


    # Post Tarier
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Macaque',
                                'Olive baboon', 'Orangutan', 'Marmoset']
    outgroup_species_list = ['Tarsier', 'Mouse Lemur']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_Tarsier_families, Post_Tarsier_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    


    # Post NWM
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Macaque',
                                'Olive baboon', 'Orangutan']
    outgroup_species_list = ['Tarsier', 'Bushbaby', 'Marmoset']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_NWM_families, Post_NWM_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    

    # Post OWM
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Macaque']
    outgroup_species_list = ['Tarsier', 'Bushbaby', 'Olive baboon', 'Orangutan', 'Marmoset']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_OWM_families, Post_OWM_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    


##        
##    # output list
##    with open(Post_Mouse_folder + 'Post_Mouse_List.txt', 'w+') as f:
##        for family in Post_Mouse_families:
##            f.write('_'.join(family) + '\n')
##
##    two_paralog_species_list = ['Chimpanzee', 'Gibbon', 'Gorilla', 'Macaque', 'Olive baboon', 'Orangutan']
##    for family in Post_Mouse_families:
##        paralog1, paralog2 = list(family)
##        if len(ortholog_mapping[paralog2]['Mouse']) > 1:
##            trans = paralog1
##            paralog1 = paralog2
##            paralog2 = trans
##        family_folder = Post_Mouse_folder + '_'.join(family) + '/'
##        if not os.path.isdir(family_folder):
##            os.mkdir(family_folder)
##
##        with open(family_folder + '_'.join(family) + '.fa', 'w+') as f:
##            f.write('>Mouse' + paralog1 + '\n')
##            ENS_id = ortholog_mapping[paralog1]['Mouse']
##            seq = getENSsequence(ENS_id)
##            seq = seq.split('\n')[0]                
##            f.write(seq + '\n')
##            for paralog in [paralog1, paralog2]:
##                f.write('>Human' + paralog + '\n')
##                f.write(getENSsequence(paralog, 'CDS', True) + '\n')
##                for species in two_paralog_species_list:
##                    f.write('>' + species.replace(' ', '_') + paralog + '\n')
##                    ENS_id = ortholog_mapping[paralog][species]
##                    
##                    seq = getENSsequence(ENS_id)
##                    seq = seq.split('\n')[0]
##                    f.write(seq + '\n')

    
        
    
    
