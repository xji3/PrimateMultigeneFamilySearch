import os

if __name__ == '__main__':
    case_list = ['Post_NWM', 'Post_OWM', 'Post_Gibbon', 'Post_Oran']
    #case_list = ['Before_Marmoset']
    inner_path = './withOutgroup/'
    #inner_path = './noOutgroup/'

    post_dup_list = ['N1', 'N2', 'N3', 'N4']
    #post_dup_list = ['N1']
    for case_iter in range(len(case_list)):
        case = case_list[case_iter]
        post_dup = post_dup_list[case_iter]
        pairs = []
        all_pairs = inner_path + case + '/' + case + '_list.txt'
        all_pairs = inner_path + case + '/' + case + '_list.txt'
        with open(all_pairs, 'r') as f:
            for line in f.readlines():
                pairs.append(line.replace('\n','').split('_'))

        model = 'MG94'

        sh_line = 'sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
        IGC_bash_file = './' + case + '_' + model + '_Primate_IGC.sh'
        
        with open(IGC_bash_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                file_size = os.stat(inner_path + case + '/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta').st_size
                if file_size > 1000:
                    f.write(sh_line + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh \n')
                    with open('./ShFiles/' + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh', 'w+') as g:
                        g.write('#!/bin/bash' + '\n')
                        g.write('python Run.py --model ' + model + ' --paralog1 ' + paralog[0]
                                + ' --paralog2 ' + paralog[1] + ' --no-force --no-clock --post_dup ' + post_dup
                                + ' --aln_folder ' + case + '/ ' + '\n')


        sh_line = 'sbatch -o ./Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_'
        IGC_bash_file = './Force_' + case + '_' + model + '_Primate_IGC.sh'
        
        with open(IGC_bash_file, 'w+') as f:
            f.write('#!/bin/bash' + '\n')
            for paralog in pairs:
                file_size = os.stat(inner_path + case + '/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta').st_size
                if file_size > 1000:
                    f.write(sh_line + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh \n')
                    with open('./ShFiles/Force_' + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh', 'w+') as g:
                        g.write('#!/bin/bash' + '\n')
                        g.write('python Run.py --model ' + model + ' --paralog1 ' + paralog[0]
                                + ' --paralog2 ' + paralog[1] + ' --force --no-clock --post_dup ' + post_dup
                                + ' --aln_folder ' + case + '/ ' + '\n')

