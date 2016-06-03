import os

if __name__ == '__main__':
    pairs = []
    all_pairs = './Post_Mouse/Post_Mouse_list.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    model = 'MG94'

    sh_line = 'sbatch -o Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
    IGC_bash_file = './' + model + '_Primate_IGC.sh'
    
    with open(IGC_bash_file, 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for paralog in pairs:
            f.write(sh_line + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh \n')
            with open('./ShFiles/' + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run.py --model ' + model + ' --paralog1 ' + paralog[0]
                        + ' --paralog2 ' + paralog[1] + ' --no-force --no-clock' + '\n')


    sh_line = 'sbatch -o Primate-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Force_'
    IGC_bash_file = './Force_' + model + '_Primate_IGC.sh'
    
    with open(IGC_bash_file, 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for paralog in pairs:
            f.write(sh_line + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh \n')
            with open('./ShFiles/Force_' + '_'.join(paralog) + '_' + model + '_nonclock' + '.sh', 'w+') as g:
                g.write('#!/bin/bash' + '\n')
                g.write('python Run.py --model ' + model + ' --paralog1 ' + paralog[0]
                        + ' --paralog2 ' + paralog[1] + ' --force --no-clock' + '\n')

