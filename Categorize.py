import os, subprocess

def get_appearance(ENS_file, species_list):
    assert(os.path.isfile(ENS_file))
    show_up_species = []
    with open(ENS_file) as f:
        for line in f:
            if line[0] == '>':
                show_up_species.append(line[1:-1])
    return {species:(species in show_up_species) for species in species_list}

if __name__ == '__main__':
    family_list = []
    with open('./GeneFamilyList_LessAccurate_ENS_id.txt', 'rb') as f:
        for line in f:
            family_list.append(line[:-1].split('_'))

######################################################################################
#  Use OrthoMAM v.9 ortholog mapping to map to other species
######################################################################################
    OrthoMaMv9_align_folder = './OrthoMaMv9_align/'
    all_files = os.listdir(OrthoMaMv9_align_folder)
    all_ENS = [file_name.split('_')[0] for file_name in all_files]
    ENS_to_FastaFile = {file_name.split('_')[0]:file_name for file_name in all_files}

######################################################################################
#  Now record gene appearance in each species
######################################################################################

    species_list = ['Homo', 'Pan', 'Gorilla', 'Pongo', 'Macaca', 'Callithrix', 'Tarsius', 'Otolemur',
                    'Tupaia', 'Mus']

    Total_count = dict()
    for pair in family_list:
        paralog1_ENS, paralog2_ENS = pair
        ENS_file_1 = OrthoMaMv9_align_folder + ENS_to_FastaFile[paralog1_ENS]
        ENS_file_2 = OrthoMaMv9_align_folder + ENS_to_FastaFile[paralog2_ENS]

        paralog1_appearance = get_appearance(ENS_file_1, species_list)
        paralog2_appearance = get_appearance(ENS_file_2, species_list)

        Total_count['_'.join(pair)] = {species:(paralog1_appearance[species] + paralog2_appearance[species]) for species in species_list}

    with open('./Total_count.txt', 'w+') as f:
        f.write('\t'.join(species_list) + '\t pair \n')
        for key in Total_count:
            f.write('\t'.join([str(Total_count[key][species]) for species in species_list]) + '\t' + key + '\n')


