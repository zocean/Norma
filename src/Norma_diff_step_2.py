#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 28 Mar 2019 13:02:07

import os,sys,argparse
from progressbar import ProgressBar
from bx.bbi.bigwig_file import BigWigFile
'''import custom function/class'''
from utility import *

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--res',type=int,dest="res",help="resolution in kb")
    p.add_argument('--genome',type=str,dest="genome",help="chromosome size file")
    p.add_argument('--bw_list',type=str,dest="bw_list",nargs="+",help="bigwig files")
    p.add_argument('--label',type=str,dest="label",nargs="+",help="label associated with bigwig files")
    p.add_argument('--mode',type=str,dest="mode",nargs="+",help="annotation mode")
    p.add_argument('--output',type=str,dest="output",help="output file")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def make_genome_window(genome_size, res):
    '''
    cut genome into windows
    '''
    win_list = []
    for chrom in sorted(genome_size.keys()):
        # if chromosome size < resolution
        if genome_size[chrom] < res:
            win_list.append(Region(chrom, 0, genome_size[chrom]))
        # cut chromosome into windows
        else:
            bin_size = int(genome_size[chrom])/int(res) + 1
            start = 0
            stop = res
            for nn in range(bin_size):
                # take care the last window
                if stop > genome_size[chrom]:
                    stop = genome_size[chrom]
                win = Region(chrom, start, stop)
                win_list.append(win)
                start = stop
                stop += res
    return win_list

def main():
    global args
    args = parse_arg()
    args.res = args.res*1000
    # check parameters
    try:
        assert len(args.bw_list) == len(args.label)
    except AssertionError:
        print >>sys.stderr, "number of bigwig file and number of label must be matched"
        exit(1)
    print >>sys.stderr, "# check parameters done"
    # parse the bigwig file
    anno_list = []
    mode_list = args.mode
    label_list = args.label
    for nn in range(len(args.bw_list)):
        anno_list.append(BigWigFile(open(args.bw_list[nn])))
    print >>sys.stderr, "# load data done"
    # build the table
    genome_size = load_genome_size(args.genome)
    win_list = make_genome_window(genome_size, args.res)
    # make the annotation
    print >>sys.stderr, "# begin annotation"
    progress = ProgressBar()
    for nn in progress(range(len(win_list))):
        win = win_list[nn]
        for mm in range(len(anno_list)):
            win.get_anno(anno_list[mm], label_list[mm], mode_list[mm], genome_size)
    progress.finish()
    # report 
    print >>sys.stderr, "# begin report"
    win_list = sorted(win_list, key = lambda x: (x.chrom, x.start))
    fout = open(args.output, 'w')
    print >>fout, win_list[0].header()
    progress = ProgressBar()
    for nn in progress(range(len(win_list))):
        win = win_list[nn]
        is_na = True
        for label in win.label:
            if win.anno[label] != 'NA':
                is_na = False
                break
        if is_na:
            continue
        else:
            print >>fout, win.write()
    progress.finish()
    fout.close()
    
if __name__=="__main__":
    main()
