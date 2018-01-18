#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 18 Jan 2018 15:15:14

import os,sys,argparse
from progressbar import ProgressBar
import pysam
from bx.bbi.bigwig_file import BigWigFile
import tabix
import random
'''import custom function/class'''
from TFBS_Evo.my_utility import *
from utility import *

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bed',type=str,dest="bed",help="bed file")
    p.add_argument('--update',dest="update",action="store_true",help="if set, will update the file assumeing header exists")
    p.add_argument('--exclude',type=str,dest="exclude",nargs="+",help="exclude region in tabix format")
    p.add_argument('--anno',type=str,dest="anno",nargs="+",help="annotation file list")
    p.add_argument('--mode',type=str,dest="mode",nargs="+",help="mode list")
    p.add_argument('--label',type=str,dest="label",nargs="+",help="label list")
    p.add_argument('--genome',type=str,dest="genome",help="genome fasta file")
    p.add_argument('--genome_size',type=str,dest="genome_size",help="genome size file")
    p.add_argument('--output',type=str,dest="output",help="output file name")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def LoadExclude(file_list):
    exclude_list = []
    if file_list is None:
        return exclude_list
    for filename in file_list:
        exclude_list.append(tabix.open(filename))
    return exclude_list

def OptionsCheck():
    text = ""
    text += "# DHS peak file: %s\n" % (args.bed)
    if args.update:
        text += "# update existing annotation: True\n"
    else:
        text += "# update existing annotation: False\n"
    text += "# genome fasta file: %s\n" % (args.genome)
    text += "# genome chromosome size file: %s\n" % (args.genome_size)
    if args.exclude is not None:
        text += "# exclude region file: " + ','.join(filename for filename in args.exclude) + '\n'
    text += "# output file: %s\n" % (args.output)
    for nn in range(len(args.anno)):
        anno = args.anno[nn]
        label = args.label[nn]
        mode = args.mode[nn]
        text += "# Annotaion label: %s\n" % (label)
        text += "#\t\tAnnotation mode: %s\n" % (mode)
        text += "#\t\tAnnotation file: %s\n" % (anno)
    print >>sys.stderr, text

def Main():
    global args
    args=ParseArg()
    # check parameter
    if args.anno is None and args.label is None and args.mode is None:
        args.anno = []; args.label = []; args.mode = []
        pass
    else:
        assert len(args.anno) == len(args.label)
        assert len(args.anno) == len(args.mode)
    OptionsCheck()
    # load genome fasta/size
    genome = pysam.Fastafile(args.genome)
    genome_size = LoadGenomeSize(args.genome_size)
    # load exclude tabix list
    if args.exclude is not None:
        print >>sys.stderr, "load exclude regions"
        exclude_list = LoadExclude(args.exclude)
    else:
        exclude_list = []
    # load bed region
    if args.update:
        region_list = LoadRegionAnno(args.bed)
    else:
        region_list = LoadRegion(args.bed)
    # load annotation files into tabix, bigwig, fasta or bam object
    anno_list = []
    label_list = args.label
    mode_list = args.mode
    for nn in range(len(mode_list)):
        mode = mode_list[nn]
        if mode in ['coverage', 'coverage_per', 'count', 'anno', 'nearby', 'tfbs']:
            try:
                anno_list.append(tabix.open(args.anno[nn]))
            except:
                print >>sys.stderr, args.anno[nn]
                exit(1)
        elif mode in ['signal_mean', 'signal_sum']:
            anno_list.append(BigWigFile(open(args.anno[nn])))
        elif mode in ['dnase_count','dnase_ave']:
            anno_list.append(SamHandler(pysam.Samfile(args.anno[nn])))
        elif mode in ['gc', 'gc_count']:
            anno_list.append(pysam.Fastafile(args.anno[nn]))
        elif mode in ['size']:
            anno_list.append(args.anno[nn])
        else:
            print >>sys.stderr, "Unknown mode: %s" % (mode)
            exit(1)
    # create center region
    progress = ProgressBar()
    for nn in progress(range(len(region_list))):
        region = region_list[nn]
        region.region_init(exclude_list, genome_size)
    # annotate 
    print >>sys.stderr, "Annotation Region"
    progress = ProgressBar()
    for nn in progress(range(len(region_list))):
        region = region_list[nn]
        for nn in range(len(anno_list)):
            anno = anno_list[nn]
            label = label_list[nn]
            mode = mode_list[nn]
            region.get_anno(anno, label, mode, genome_size)
            try:
                region.get_anno(anno, label, mode, genome_size)
            except:
                print >>sys.stderr, "[anno] skip %s:%d-%d" % (region.chrom, region.start, region.stop)
                continue
    print >>sys.stderr, "Annotation bed done"
    # write to output file
    if args.output == args.bed: # re-annotate existing file
        fout = WriteToFile(args.output+'.tmp')
    else:
        fout = WriteToFile(args.output)
    print >>fout, region_list[0].header()
    region_list = sorted(region_list, key = lambda bed: (bed.chrom, bed.start))
    for region in region_list:
        try:
            print >>fout, region.write()
        except:
            print >>sys.stderr, "[report] skip %s:%d-%d" % (region.chrom, region.start, region.stop)
            continue
    fout.close()
    # program finished 
    if args.output == args.bed:
        os.system("mv %s %s" % (args.output + '.tmp', args.output)) 

if __name__=="__main__":
    Main()
