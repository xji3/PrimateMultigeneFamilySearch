from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]
    Force = None
    alignment_file = './' + args.aln_folder + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = './Primate_Tree.newick'
    if args.force:
        if args.model == 'MG94':
            Force = {5:0.0}
        elif args.model == 'HKY':
            Force = {4:0.0}
    else:
        Force = None
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = args.model, Force = Force, clock = args.clock, post_dup = args.post_dup)
    test.get_mle(True, True, 0, 'BFGS')
    test.get_individual_summary(summary_path = './Summary/')
    test.get_SitewisePosteriorSummary(summary_path = './Summary/')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required = True, help = 'Substitution Model')
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
    parser.add_argument('--post_dup', required = True, help = 'First post duplication node name')
    parser.add_argument('--aln_folder', required = True, help = 'Alignment folder')
    
    
    main(parser.parse_args())


##    paralog = ['ENSG00000099194', 'ENSG00000145284']
##    Force = None
##    model = 'MG94'
##    alignment_file = './Post_Mouse/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##    newicktree = './Primate_Tree.newick'
##    if Force:
##        if model == 'MG94':
##            Force = {5:0.0}
##        elif model == 'HKY':
##            Force = {4:0.0}
##    else:
##        Force = None
##    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = False)
##    test.get_mle(True, True, 0, 'BFGS')
##    test.get_individual_summary(summary_path = './Summary/')
##    test.get_SitewisePosteriorSummary(summary_path = './Summary/')

