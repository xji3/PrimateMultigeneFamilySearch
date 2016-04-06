import urllib2, gzip, os, json, pickle, sys

sys.path.insert(0, '../metaphors_api/')
import dbClient

from Bio import SeqIO

#### http://rest.ensembl.org/documentation/info/sequence_id
import requests, sys

#### http://pycogent.org/examples/query_ensembl.html
#from cogent.db.ensembl import HostAccount, Genome


def getENSsequence(ENS_id, seq_type = 'CDS'):
    server = "http://rest.ensembl.org"
    ext = "/sequence/id/" + ENS_id + '?type=' + seq_type
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    return r.text



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
    remove_group_id = [group_id for group_id in group_id_to_ENS if not all([(ENS in all_ENS) for ENS in group_id_to_ENS[group_id]])]
    for group_id in remove_group_id:
        del group_id_to_ENS[group_id]

######################################################################################
#  Output a list of less accurate gene families first
######################################################################################
    GeneFamily_Folder = './GeneFamilies/'
    ENS_fasta_species_list = ['Pongo', 'Macaca', 'Gorilla', 'Pan', 'Homo']
    with open('./GeneFamilyList_LessAccurate.txt', 'w+') as f:
        for group_id in group_id_to_ENS:
            ENS_1, ENS_2 = group_id_to_ENS[group_id]
            paralog1 = ENS_to_name[ENS_1]
            paralog2 = ENS_to_name[ENS_2]
            f.write(paralog1 + '_' + paralog2 + '\n')
            with open(GeneFamily_Folder + paralog1 + '_' + paralog2 + '.fasta', 'w+') as g:
                #f.write(paralog1 + '_' + paralog2 + '\n')
                with open('./OrthoMaMv9_align/' + ENS_to_FastaFile[ENS_1], 'rb') as fasta_file:
                    species_to_fasta = dict()
                    to_write = False
                    for line in fasta_file:
                        if line[0] == '>':
                            species = line.split()[0][1:]
                        else:
                            sequence = line.split()[0]
                            species_to_fasta[species] = sequence

                for species in ENS_fasta_species_list:
                    g.write('>' + species + '_' + paralog1 + '\n')
                    g.write(species_to_fasta[species] + '\n')

            
                with open('./OrthoMaMv9_align/' + ENS_to_FastaFile[ENS_2], 'rb') as fasta_file:
                    species_to_fasta = dict()
                    to_write = False
                    for line in fasta_file:
                        if line[0] == '>':
                            species = line.split()[0][1:]
                        else:
                            sequence = line.split()[0]
                            species_to_fasta[species] = sequence

                for species in ENS_fasta_species_list:
                    g.write('>' + species + '_' + paralog2 + '\n')
                    g.write(species_to_fasta[species] + '\n')


######################################################################################
#  Now try to see if other species all have only 2 paralogs by metaPhors database
######################################################################################
    
    m = dbClient.metaphors()
    
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
    
##    #######
##    # Now use HTML parser to search for species
##    #######
##    selected_metaPhor_id = []
##    if not os.path.isfile('./selected_metaPhor_id.txt'):
##        target_species_list = [u'Homo sapiens', u'Pan troglodytes', u'Gorilla gorilla', u'Pongo abelii', u'Macaca mulatta']
##        for meta_id in metaPhor_to_ENS:
##            print 'Query metaPhors id: ', meta_id
##            website_query = "http://betaorthology.phylomedb.org/wsgi/query/getGenes?term=" + ENS_to_name[metaPhor_to_ENS[meta_id]]
##            ### This method does not work for all genes, but some of them
##            handle = urllib2.urlopen(website_query)
##            species_to_metaid = json.loads(handle.readline())
##            species_list = [i[u'species'] for i in species_to_metaid]
##            #print species_list
##            if all([target_species in species_list for target_species in target_species_list]):
##                print "It is proved!"
##                selected_metaPhor_id.append(meta_id)
##    else:
##        with open('./selected_metaPhor_id.txt', 'rb') as f:
##            for line in f:
##                selected_metaPhor_id.append(line.split()[0])
    #######
    # Now use metaphors API to search for orthologs in species
    #######
    selected_metaPhor_id = []
    target_species_list  = [9544L, 9593L, 9598L, 9601L, 9606L]
    for meta_id in metaPhor_to_ENS:
        print 'Query metaPhors id: ', meta_id
        gene_name = ENS_to_name[metaPhor_to_ENS[meta_id]]
        protid = m.get_metaid(gene_name)
        if protid:
            orthologs = m.get_orthologs(protid)
            if all([species in orthologs for species in target_species_list]):
                selected_metaPhor_id.append(meta_id)


    
    #######
    # Now make intersection between selected meta_id and group_id
    #######

    selected_group_id_list = [group_id for group_id in group_id_to_ENS if all([(ENS_to_metaPhor[ENS] in selected_metaPhor_id) if ENS in ENS_to_metaPhor.keys() else False for ENS in group_id_to_ENS[group_id] ])]

######################################################################################
#  Now get meta_id of other species and then link back to ENS id
######################################################################################

    if not os.path.isfile('./Homo_meta_id_to_others.p'):
        target_species_list = [u'Homo sapiens', u'Pan troglodytes', u'Gorilla gorilla', u'Pongo abelii', u'Macaca mulatta']
        Homo_meta_id_to_others = dict()
        for meta_id in selected_metaPhor_id:
            print 'Query metaPhors id: ', meta_id
            Homo_meta_id_to_others[meta_id] = dict()
            website_query = "http://betaorthology.phylomedb.org/wsgi/query/getGenes?term=" + ENS_to_name[metaPhor_to_ENS[meta_id]]
            ### This method does not work for all genes, but some of them
            handle = urllib2.urlopen(website_query)
            species_to_metaid = json.loads(handle.readline())
            for item in species_to_metaid:
                species = item[u'species']
                if species in target_species_list:
                    Homo_meta_id_to_others[meta_id][str(species)] = str(item[u'metaid'])
        pickle.dump( Homo_meta_id_to_others, open('./Homo_meta_id_to_others.p', 'w+'))
    else:
        Homo_meta_id_to_others = pickle.load(open('./Homo_meta_id_to_others.p', 'rb'))
    
    metaphor_id_file = '/Users/Xiang/Downloads/id_conversion.txt.gz'
    ENS_to_metaPhor_file = './ENS_to_metaPhor.txt'
    new_meta_id = [Homo_meta_id_to_others[meta_id][species] for species in ['Homo sapiens', 'Pan troglodytes', 'Pongo abelii', 'Macaca mulatta', 'Gorilla gorilla']
                   for meta_id in Homo_meta_id_to_others.keys()]
    meta_to_TrEMBL = dict()
    meta_to_TrEMBL_file = './meta_to_TrEMBL.txt'
    if not os.path.isfile(meta_to_TrEMBL_file):        
        with gzip.open(metaphor_id_file, 'rb') as f:
            f.readline()
            for line in f:
                content = line.split()
                if content[2] in new_meta_id:
                    print content[2]
                    meta_to_TrEMBL[content[2]] = content[0]
    else:
        meta_to_TrEMBL = pickle.load(open(meta_to_TrEMBL_file, 'rb'))

    TrEMBL_to_Ensembl_file = './TrEMBL_to_Ensembl_mapping.txt'
    ### http://www.uniprot.org/mapping/
    TrEMBL_to_Ensembl = dict()
    with open(TrEMBL_to_Ensembl_file, 'rb') as f:
        f.readline()
        for line in f:
            content = line.split()
            TrEMBL_to_Ensembl[content[0]] = content[1]


######################################################################################
#  Now construct fasta files
######################################################################################
    
    # Use gene name for file names and Gene family names
    GeneFamily_Folder = './GeneFamilies/'
    ENS_fasta_species_list = ['Pongo', 'Macaca', 'Gorilla', 'Pan', 'Homo']
    with open('./GeneFamilyList.txt', 'w+') as f:    
        for group_id in selected_group_id_list:
            print group_id
            ENS_1 = group_id_to_ENS[group_id][0]
            ENS_2 = group_id_to_ENS[group_id][1]
            paralog1 = ENS_to_name[ENS_1]
            paralog2 = ENS_to_name[ENS_2]
            meta_id_1 = ENS_to_metaPhor[ENS_1]
            meta_id_2 = ENS_to_metaPhor[ENS_2]
            
            f.write(paralog1 + '_' + paralog2 + '\n')
                
            if ENS_1 in ENS_to_FastaFile and ENS_2 in ENS_to_FastaFile:
                with open(GeneFamily_Folder + paralog1 + '_' + paralog2 + '.fasta', 'w+') as g:
                    #f.write(paralog1 + '_' + paralog2 + '\n')
                    with open('./OrthoMaMv9_align/' + ENS_to_FastaFile[ENS_1], 'rb') as fasta_file:
                        species_to_fasta = dict()
                        to_write = False
                        for line in fasta_file:
                            if line[0] == '>':
                                species = line.split()[0][1:]
                            else:
                                sequence = line.split()[0]
                                species_to_fasta[species] = sequence

                    for species in ENS_fasta_species_list:
                        g.write('>' + species + '_' + paralog1 + '\n')
                        g.write(species_to_fasta[species] + '\n')

                
                    with open('./OrthoMaMv9_align/' + ENS_to_FastaFile[ENS_2], 'rb') as fasta_file:
                        species_to_fasta = dict()
                        to_write = False
                        for line in fasta_file:
                            if line[0] == '>':
                                species = line.split()[0][1:]
                            else:
                                sequence = line.split()[0]
                                species_to_fasta[species] = sequence

                    for species in ENS_fasta_species_list:
                        g.write('>' + species + '_' + paralog2 + '\n')
                        g.write(species_to_fasta[species] + '\n')

                    continue

##            else:
##                print ENS_1, ENS_2
##                with open('./manual/'  + paralog1 + '_' + paralog2 + '_info.txt', 'w+') as h:
##                    h.write('\t'.join(['Species', 'Gene', 'TrEMBL id', 'metaPhors id',  'Ensembl id', '\n']))
##                    for species in ['Homo sapiens', 'Pan troglodytes', 'Pongo abelii', 'Macaca mulatta', 'Gorilla gorilla']:  # Homo_meta_id_to_others[Homo_meta_id_to_others.keys()[0]].keys()
##                        for ENS_id in [ENS_1, ENS_2]:
##                            meta_id = ENS_to_metaPhor[ENS_id]
##                            paralog = ENS_to_name[ENS_id]
##                            species_meta_id = Homo_meta_id_to_others[meta_id][species]
##                            species_TrEMBL_id = meta_to_TrEMBL[species_meta_id]
##                            species_ENS_id = TrEMBL_to_Ensembl[species_TrEMBL_id]
##
##                            h.write('\t'.join([species, paralog, species_TrEMBL_id, species_meta_id, species_ENS_id, '\n']))



##                    with open(GeneFamily_Folder + paralog1 + '_' + paralog2 + '.fasta', 'w+') as g:
##                        for species in ['Homo sapiens', 'Pan troglodytes', 'Pongo abelii', 'Macaca mulatta', 'Gorilla gorilla']:
##                            g.write('>' + species.split()[0] + '_' + paralog1 + '\n')
##                            g.write('\n')
##                        for species in ['Homo sapiens', 'Pan troglodytes', 'Pongo abelii', 'Macaca mulatta', 'Gorilla gorilla']:
##                            g.write('>' + species.split()[0] + '_' + paralog2 + '\n')
##                            g.write('\n')

                    
##                if meta_id_1 in Homo_meta_id_to_others and meta_id_2 in Homo_meta_id_to_others:
##                    f.write(paralog1 + '_' + paralog2 + '\n')
##                    print group_id
##                    for species in ['Homo sapiens', 'Pan troglodytes', 'Pongo abelii', 'Macaca mulatta', 'Gorilla gorilla']:  # Homo_meta_id_to_others[Homo_meta_id_to_others.keys()[0]].keys()
##                        species_meta_id = Homo_meta_id_to_others[meta_id_1][species]
##                        species_TrEMBL_id = meta_to_TrEMBL[species_meta_id]
##                        species_ENS_id = TrEMBL_to_Ensembl[species_TrEMBL_id]
##                        
##                        g.write('>' + species.split()[0] + '_' + paralog1 + '\n')
##                        g.write(getENSsequence(species_ENS_id) + '\n')
##
##
##                    for species in ['Homo sapiens', 'Pan troglodytes', 'Pongo abelii', 'Macaca mulatta', 'Gorilla gorilla']:  # Homo_meta_id_to_others[Homo_meta_id_to_others.keys()[0]].keys()
##                        species_meta_id = Homo_meta_id_to_others[meta_id_2][species]
##                        species_TrEMBL_id = meta_to_TrEMBL[species_meta_id]
##                        species_ENS_id = TrEMBL_to_Ensembl[species_TrEMBL_id]
##                        
##                        g.write('>' + species.split()[0] + '_' + paralog1 + '\n')
##                        g.write(getENSsequence(species_ENS_id) + '\n')
##
                

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
