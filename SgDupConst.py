# Construct Human paralog gene pairs in segmental duplication region
# Xiang Ji
# xji3@ncsu.edu

# almost duplicated EnsemblSearch.py
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
    Human_gene_pair_file = home_path + 'GitFolders/PrimateMultigeneFamilySearch/Seg_Dup_Gene_Pairs.txt'
    Human_sg_pairs = numpy.loadtxt(Human_gene_pair_file, dtype = str)

    # Now subtract all two paralog families
    two_paralog_families = [pair.replace('"', '').split('_') for pair in Human_sg_pairs]
    visited_genes = list()
    for pair in two_paralog_families:
        visited_genes.append(pair[0])
        visited_genes.append(pair[1])
    visited_genes = list(set(visited_genes))

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


    # Decided to remove Bushbaby, Mouse_Lemur, Tarsier, Vervet, Marmoset
    Post_NWM_families = [family for family in two_paralog_families
                         if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 1
                         ]
##    same_ortholog = [family for family in Post_NWM_families
##                     if sum([ortholog_mapping[i]['Tarsier'] != ''
##                             and ortholog_mapping[i]['Mouse Lemur'] != ''
##                             and ortholog_mapping[i]['Marmoset'] != ''
##                             for i in family])]
##    Post_NWM_families = [family for family in Post_NWM_families if family in same_ortholog]
    print 'Post NWM families: ', len(Post_NWM_families), len([family for family in Post_NWM_families if sum([ens in sg_dup_gene_ids for ens in family])])

    Post_OWM_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 1
                           ]
    same_ortholog = [family for family in Post_OWM_families
                     if sum([ortholog_mapping[i]['Macaque'] != ''
                             and ortholog_mapping[i]['Olive baboon'] != ''
                             for i in family])]
    Post_OWM_families = [family for family in Post_OWM_families if family in same_ortholog]    

    print 'Post OWM families: ', len(Post_OWM_families), len([family for family in Post_OWM_families if sum([ens in sg_dup_gene_ids for ens in family])])

    Post_Gibbon_families = [family for family in two_paralog_families
                           if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                           and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 1
                           and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 1
                           ]
    same_ortholog = [family for family in Post_Gibbon_families
                     if sum([ortholog_mapping[i]['Macaque'] != ''
                             and ortholog_mapping[i]['Olive baboon'] != ''
                             and ortholog_mapping[i]['Gibbon'] != ''
                             for i in family])]
    Post_Gibbon_families = [family for family in Post_Gibbon_families if family in same_ortholog]
    print 'Post Gibbon families: ', len(Post_Gibbon_families), len([family for family in Post_Gibbon_families if sum([ens in sg_dup_gene_ids for ens in family])])

    Post_Oran_families = [family for family in two_paralog_families
                          if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 1
                          and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                          and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 1
                          and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 1
                          ]
    same_ortholog = [family for family in Post_Oran_families
                     if sum([ortholog_mapping[i]['Macaque'] != ''
                             and ortholog_mapping[i]['Olive baboon'] != ''
                             and ortholog_mapping[i]['Gibbon'] != ''
                             and ortholog_mapping[i]['Orangutan'] != ''
                             for i in family])]
    Post_Oran_families = [family for family in Post_Oran_families if family in same_ortholog]
    print 'Post Orangutan families: ', len(Post_Oran_families), len([family for family in Post_Oran_families if sum([ens in sg_dup_gene_ids for ens in family])])


# Note: EDN_ECP is in this family ['ENSG00000129538', 'ENSG00000169385', 'ENSG00000169397', 'ENSG00000165799', 'ENSG00000258818', 'ENSG00000169413', 'ENSG00000214274', 'ENSG00000173431']

#####################################################################################################
#####################################################################################################

    # Construct datasets that include outgroup species

    Sg_GF = './Sg_GeneFamilies/'
    if not os.path.isdir(Sg_GF):
        os.mkdir(Sg_GF)

    Ensembl_GF = Sg_GF + 'withOutgroup/'
    if not os.path.isdir(Ensembl_GF):
        os.mkdir(Ensembl_GF)

    Post_NWM_folder = Ensembl_GF + 'Post_NWM/'
    if not os.path.isdir(Post_NWM_folder):
        os.mkdir(Post_NWM_folder)

    Post_OWM_folder = Ensembl_GF + 'Post_OWM/'
    if not os.path.isdir(Post_OWM_folder):
        os.mkdir(Post_OWM_folder)

    Post_Gibbon_folder = Ensembl_GF + 'Post_Gibbon/'
    if not os.path.isdir(Post_Gibbon_folder):
        os.mkdir(Post_Gibbon_folder)

    Post_Oran_folder = Ensembl_GF + 'Post_Oran/'
    if not os.path.isdir(Post_Oran_folder):
        os.mkdir(Post_Oran_folder)

    # output Post_NWM list
    with open(Post_NWM_folder + 'Post_NWM_List.txt', 'w+') as f:
        for family in Post_NWM_families:
            f.write('_'.join(family) + '\n')

    # output Post_OWM list
    with open(Post_OWM_folder + 'Post_OWM_List.txt', 'w+') as f:
        for family in Post_OWM_families:
            f.write('_'.join(family) + '\n')

    # output Post_Gibbon list
    with open(Post_Gibbon_folder + 'Post_Gibbon_List.txt', 'w+') as f:
        for family in Post_Gibbon_families:
            f.write('_'.join(family) + '\n')

    # output Post_Oran list
    with open(Post_Oran_folder + 'Post_Oran_List.txt', 'w+') as f:
        for family in Post_Oran_families:
            f.write('_'.join(family) + '\n')


    # Post NWM
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Macaque',
                                'Olive baboon', 'Orangutan']
    outgroup_species_list = ['Marmoset']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_NWM_families, Post_NWM_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    

    # Post OWM
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Orangutan']
    outgroup_species_list = ['Olive baboon', 'Macaque']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_OWM_families, Post_OWM_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    

    # Post Gibbon
    two_paralog_species_list = ['Gorilla', 'Orangutan']
    outgroup_species_list = ['Gibbon', 'Olive baboon', 'Macaque']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_Gibbon_families, Post_Gibbon_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    


    # Post Orangutan
    two_paralog_species_list = [ 'Gorilla']
    outgroup_species_list = ['Gibbon', 'Orangutan', 'Olive baboon', 'Macaque']
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Post_Oran_families, Post_Oran_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)
    

#####################################################################################################
#####################################################################################################

    # Construct datasets that don't have any outgroup species
    Before_Marmoset_families = [family for family in two_paralog_families
                         if sum([ortholog_mapping[i]['Gibbon'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Gorilla'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Orangutan'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Macaque'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Olive baboon'] != '' for i in family]) == 2
                         and sum([ortholog_mapping[i]['Marmoset'] != '' for i in family]) == 2
                         ]
    Before_Marmoset_families.remove(['ENSG00000207866', 'ENSG00000207867'])  # mRNA
    Before_Marmoset_families.remove(['ENSG00000207866', 'ENSG00000207868'])  # mRNA
    Before_Marmoset_families.remove(['ENSG00000207867', 'ENSG00000207868'])  # mRNA

    Before_Marmoset_families.remove(['ENSG00000251796', 'ENSG00000251925'])  # snRNA
    Before_Marmoset_families.remove(['ENSG00000251796', 'ENSG00000253042'])  # snRNA
    Before_Marmoset_families.remove(['ENSG00000251925', 'ENSG00000253042'])  # snRNA

    Before_Marmoset_families.remove(['ENSG00000104818', 'ENSG00000267335'])
    Before_Marmoset_families.remove(['ENSG00000104818', 'ENSG00000267631'])
    Before_Marmoset_families.remove(['ENSG00000104827', 'ENSG00000267335'])
    Before_Marmoset_families.remove(['ENSG00000189052', 'ENSG00000267335'])
    Before_Marmoset_families.remove(['ENSG00000196337', 'ENSG00000267335'])
    Before_Marmoset_families.remove(['ENSG00000205777', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000215269', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000215274', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000216649', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000224659', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000227488', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000236362', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000237671', 'ENSG00000275113'])
    Before_Marmoset_families.remove(['ENSG00000239862', 'ENSG00000250036'])
    Before_Marmoset_families.remove(['ENSG00000267335', 'ENSG00000267631'])
            
    
    Sg_GF = './Sg_GeneFamilies/'
    if not os.path.isdir(Sg_GF):
        os.mkdir(Sg_GF)

    Ensembl_GF = Sg_GF + 'noOutgroup/'
    if not os.path.isdir(Ensembl_GF):
        os.mkdir(Ensembl_GF)

    Before_Marmoset_folder = Ensembl_GF + 'Before_Marmoset/'
    if not os.path.isdir(Before_Marmoset_folder):
        os.mkdir(Before_Marmoset_folder)

    # output Before_Marmoset list
    with open(Before_Marmoset_folder + 'Before_Marmoset_List.txt', 'w+') as f:
        for family in Before_Marmoset_families:
            f.write('_'.join(family) + '\n')

    # Before Marmoset
    two_paralog_species_list = ['Gibbon', 'Gorilla', 'Macaque',
                                'Olive baboon', 'Orangutan', 'Marmoset']
    outgroup_species_list = []
    # same_ortholog = [family for family in Post_Tarsier_families if sum([ortholog_mapping[i]['Tarsier'] != '' and ortholog_mapping[i]['Bushbaby'] != '' for i in family])]
    # All post_Tarsier families have the same paralog present in outgroups

    outputFamilies(Before_Marmoset_families, Before_Marmoset_folder, ortholog_mapping, outgroup_species_list, two_paralog_species_list)    

