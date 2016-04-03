import urllib2, gzip, os, json
from Bio import SeqIO

#### http://pycogent.org/examples/query_ensembl.html
#from cogent.db.ensembl import HostAccount, Genome

if __name__ == '__main__':
##    species_to_taxid = {'Human':9606, 'Chimpanzee':9598,
##                        'Gorilla':9593, 'Orangutan':9601,
##                        'Macaque':9544}
##    species_to_protid = {species:[] for species in species_to_taxid}
##    metaphor_path = './metaPhors/'
##    for species in species_to_taxid:
##        gz_file = metaphor_path + str(species_to_taxid[species]) + '.txt.gz'
##        with gzip.open(gz_file, 'rb') as f:
##            f.readline()
##            for line in f:
##                content = line.split()
##                #species_to_protid[species].append([content[1], content[3]])

######################################################################################
#  Use Duplicated Genes Database to select Human 2-paralog gene families
######################################################################################

    DGD_tab_file = './DuplicatedGenesDatabase/dgd_Hsa_all_v71.tsv.gz'
    ## http://dgd.genouest.org/
    ### Tab info
    ##:chr          group_id        NB_Genes        start   end     strand  ENS_ID  Name
    ## Description     CCDS    Ensembl Human Gene      UCSC Stable ID          UniGene         UniParc         UniProt transcript name         UniProtKB Gene Name     UniProtKB/Swiss-Prot    UniProtKB/TrEMBL
    group_id_to_ENS = dict()
    ENS_to_TrEMBL = dict()
    ENS_to_name = dict()
    with gzip.open(DGD_tab_file, 'rb') as f:
        f.readline()
        for line in f:
            content = line.split()
            group_id = content[1]
            if int(content[2]) == 2: ## NB_Genes == 2
                ENS_to_TrEMBL[content[6]] = content[-1]
                ENS_to_name[content[6]] = content[7]
                if not group_id in group_id_to_ENS:
                    group_id_to_ENS[group_id] = [content[6]]
                else:
                    group_id_to_ENS[group_id].append(content[6])

    TrEMBL_to_ENS = {ENS_to_TrEMBL[ENS]:ENS for ENS in ENS_to_TrEMBL}
    name_to_ENS = {ENS_to_name[ENS]:ENS for ENS in ENS_to_name}

######################################################################################
#  Use OrthoMAM v.9 ortholog mapping to map to other species
######################################################################################
    OrthoMaMv9_align_folder = './OrthoMaMv9_align/'
    all_files = os.listdir(OrthoMaMv9_align_folder)
    all_ENS = [file_name.split('_')[0] for file_name in all_files]
    ENS_to_FastaFile = {file_name.split('_')[0]:file_name for file_name in all_files}

######################################################################################
#  Filter all 2-paralog gene families that both show up in all_ENS
######################################################################################
    select_pair = [group_id_to_ENS[group_id] for group_id in group_id_to_ENS if all([(ENS in all_ENS) for ENS in group_id_to_ENS[group_id]])]


######################################################################################
#  Now try to see if other species all have only 2 paralogs by metaPhors database
######################################################################################

    metaphor_id_file = 'ftp://phylomedb.org/metaphors/release-201601/id_conversion.txt.gz'
    ENS_to_metaPhor_file = './ENS_to_metaPhor.txt'
    ENS_to_metaPhor = dict()
    if not os.path.isfile(ENS_to_metaPhor_file):        
        with gzip.open(metaphor_id_file, 'rb') as f:
            f.readline()
            for line in f:
                content = line.split()
                if content[0] in TrEMBL_to_ENS:
                    print content[0]
                    ENS_to_metaPhor[TrEMBL_to_ENS[content[0]]] = content[2]
    else:
        with open(ENS_to_metaPhor_file, 'rb') as f:
            for line in f:
                content = line.split()
                ENS_to_metaPhor[content[0]] = content[1]

    metaPhor_to_ENS = {ENS_to_metaPhor[ENS]:ENS for ENS in ENS_to_metaPhor}
    
    #######
    # Now use HTML parser to search for species
    #######
    selected_metaPhor_id = []
    if not os.path.isfile('./selected_metaPhor_id.txt'):
        target_species_list = [u'Homo sapiens', u'Pan troglodytes', u'Gorilla gorilla', u'Pongo abelii', u'Macaca mulatta']
        for meta_id in metaPhor_to_ENS:
            print 'Query metaPhors id: ', meta_id
            website_query = "http://betaorthology.phylomedb.org/wsgi/query/getGenes?term=" + ENS_to_name[metaPhor_to_ENS[meta_id]]
            ### This method does not work for all genes, but some of them
            handle = urllib2.urlopen(website_query)
            species_to_metaid = json.loads(handle.readline())
            species_list = [i[u'species'] for i in species_to_metaid]
            #print species_list
            if all([target_species in species_list for target_species in target_species_list]):
                print "It is proved!"
                selected_metaPhor_id.append(meta_id)
    else:
        with open('./selected_metaPhor_id.txt', 'rb') as f:
            for line in f:
                selected_metaPhor_id.append(line.split()[0])
    
    #######
    # Now make intersection between selected meta_id and group_id
    #######

    selected_group_id_list = [group_id for group_id in group_id_to_ENS if all([(ENS_to_metaPhor[ENS] in selected_metaPhor_id) if ENS in ENS_to_metaPhor.keys() else False for ENS in group_id_to_ENS[group_id] ])]
    

######################################################################################
#  Now construct fasta files
######################################################################################
    
    # Use gene name for file names and Gene family names
    GeneFamily_Folder = './GeneFamilies/'
    ENS_fasta_species_list = ['Pongo', 'Macaca', 'Gorilla', 'Pan', 'Homo']
    with open('./GeneFamilyList.txt', 'w+') as f:    
        for group_id in selected_group_id_list:
            ENS_1 = group_id_to_ENS[group_id][0]
            ENS_2 = group_id_to_ENS[group_id][1]
            paralog1 = ENS_to_name[ENS_1]
            paralog2 = ENS_to_name[ENS_2]
            f.write(paralog1 + '_' + paralog2 + '\n')
##            if ENS_1 in ENS_to_FastaFile and ENS_2 in ENS_to_FastaFile:
##                f.write(paralog1 + '_' + paralog2 + '\n')
##                
##                with open(GeneFamily_Folder + paralog1 + '_' + paralog2 + '.fasta', 'w+') as g:
##                    with open('./OrthoMaMv9_align/' + ENS_to_FastaFile[ENS_1], 'rb') as fasta_file:
##                        species_to_fasta = dict()
##                        to_write = False
##                        for line in fasta_file:
##                            if line[0] == '>':
##                                species = line.split()[0][1:]
##                            else:
##                                sequence = line.split()[0]
##                                species_to_fasta[species] = sequence
##
##                    for species in ENS_fasta_species_list:
##                        g.write('>' + species + '_' + paralog1 + '\n')
##                        g.write(species_to_fasta[species] + '\n')
##
##                
##                    with open('./OrthoMaMv9_align/' + ENS_to_FastaFile[ENS_2], 'rb') as fasta_file:
##                        species_to_fasta = dict()
##                        to_write = False
##                        for line in fasta_file:
##                            if line[0] == '>':
##                                species = line.split()[0][1:]
##                            else:
##                                sequence = line.split()[0]
##                                species_to_fasta[species] = sequence
##
##                    for species in ENS_fasta_species_list:
##                        g.write('>' + species + '_' + paralog2 + '\n')
##                        g.write(species_to_fasta[species] + '\n')
##                
