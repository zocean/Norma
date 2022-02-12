#!/usr/bin/python
# Programmer : Yang Zhang
# Contact: yzhan116@illinois.edu
# Last-modified: 23 Jun 2021 05:11:58 PM

import os,sys,argparse
from math import log
import numpy as np
import pysam
import re

def parse_argument():
    ''' This Function Parse the Argument '''
    p = argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog = 'Library dependency : pysam')
    p.add_argument('-v','--version',action = 'version', version = '%(prog)s 0.5')
    p.add_argument('-r','--res',dest = "res", type = int, default = 10, help = "wig file resultion(bp)")
    p.add_argument('-w','--win',dest = "win", type = int, default = 1000,  help = "Sliding window size(bp)")
    p.add_argument('-e','--exp',dest = "bamexp", type = str, help = "pulldown group bam file")
    p.add_argument('-c','--con',dest = "bamcon", type = str, help = "control group bam file")
    p.add_argument('-g','--genome',dest = "genome", type = str, help = "genome chromosome size file (used for wigToBigWig)")
    p.add_argument('-o','--output',dest = "output", type = str, help = "the prefix of output wig files without .wig, eg. ") 
    p.add_argument('--wig2bw',dest = "wig2bw", type = str, help = "program location (full path eg. /home/zocean/wigToBigWig of wigToBigWig program, 'sys' means wigToBigWig is the system path")
    if len(sys.argv)<2:
        sys.exit(p.print_help())
    return p.parse_args()

##########
# Utility
##########

def logging(text):
    print >>sys.stderr, "Logging: " + text

def warning(text):
    print >>sys.stderr, "Warning: " + text

def error(text):
    print >>sys.stderr, "Error: " + text

############
# Functions
############

def filterPick(list,filter):
    '''
    return the items match the filter pattern
    '''
    return [ ( l ) for l in list for m in (filter(l),) if m]

def delete_key(dic, key_list):
    for key in key_list:
        del dic[key]
    return dic

def BamToBin(bamfile,binsize):
    '''
    This function count read number ( based on the largest mappable region on the genome) into bins (same size as resolution, with minimal value 1bp)
    '''
    # creat a hash dictionary for chromosome name in the header of bamfile
    bin_hash = {}
    # read the indexed bam file using pysam
    samfile = pysam.Samfile(bamfile,"rb")
    # get number of mapped reads
    total_mapped = samfile.mapped
    # read chromosome name (chr1 etc.) and the length into a dictionary
    chr_table = {}
    for nn in range(len(samfile.references)):
        chr_table[samfile.references[nn]] = samfile.lengths[nn]
    # remove user defined chromosome
    chr_removed_list = filterPick(chr_table.keys(),chr_filter)
    chr_table = delete_key(chr_table,chr_removed_list)
    # initialization count table
    for chrom,length in chr_table.items():
        bin_num = length/binsize
        if length%binsize != 0:
            bin_num += 1
        # put 0 at in bins, innitiation bin count
        bin_hash[chrom] = [0.000000 for row in range(bin_num)]
    # look through every read in bamfile
    n = 0
    last_percent = 0
    for alignment in samfile:
        if alignment.tid < 0:
            continue
        pos = alignment.pos + alignment.alen
        chr_name = samfile.getrname(alignment.tid)
        if chr_table.get(chr_name, None) is None:
            continue
        try:
            bin_hash[chr_name][pos/binsize] += 1.0
        except KeyError:
            warning("read position %d is larger than chromosome size %d" % (pos,chr_table[chr_name]))
            continue
        # add processing information
        n += 1
        percent = float(n)/float(total_mapped)*100
        if int(percent) - last_percent == 2:
            print >>sys.stderr, "="*int(int(percent)/100.0*80.0) + '>' + '{:.2f}'.format(percent) + '%' + '\r',
            last_percent = int(percent)
    print >>sys.stderr, ""
    samfile.close()
    return bin_hash

def BinToWin(bin_hash,step):
    '''
    This function roll up read count from bins into windows. The second argument "step" tells the function how many bins would be merged together
    '''
    # read chromosome name from hash dictionary
    chr_list = bin_hash.keys()
    # create a hash dictionary for read count in window
    win_hash = {}
    if step == 1:
        for chrom in chr_list:
            win_hash[chrom] = list(bin_hash[chrom])
    else:
        for chrom in chr_list:
            # initialization window hash dictionary
            win_no = len(bin_hash[chrom]) - step + 1
            win_hash[chrom] = [0 for n in range(win_no)]
            count_list = bin_hash[chrom]
            buffer = sum(count_list[0:step])
            win_hash[chrom][0] = buffer
            for j in range(1,win_no):
                win_hash[chrom][j] = buffer + count_list[j-1+step] - count_list[j-1]
                buffer = win_hash[chrom][j] # update buffer
    return win_hash

def CalAve(table):
    value_list = table.values()
    array = [ aa for dd in value_list for aa in dd if aa is not None]
    return np.nanmean(array)

def WinNormalize(exp_win,con_win):
    '''
    normalize pulldown with input
    '''
    chrs_A = exp_win.keys()
    chrs_B = con_win.keys()
    # check chromosome between pulldown sample and control sample
    if sorted(chrs_A) != sorted(chrs_B):
        print "The chromosome list between two bam files are different, are you sure continue? (Yes/No)"
        go_on = raw_input("> ")
        if go_on == "No":
            sys.exit("interrupt by user")
        elif go_on == "Yes":
            pass
        else:
            sys.exit("Unknown answer")
    # get average number of reads of window
    con_ave_win = CalAve(con_win)
    chrom_list = list(set(chrs_A) & set(chrs_B))
    norm_win = {}
    for chrom in chrom_list:
        norm_win[chrom] = [0 for nn in range(len(exp_win[chrom]))]
        out_list = norm_win[chrom]
        exp_list = exp_win[chrom]
        con_list = con_win[chrom]
        for nn in range(len(out_list)):
            if con_list[nn] > 1e-6 and exp_list[nn] > 1e-6:
                out_list[nn] = (exp_list[nn])*(con_ave_win/(con_list[nn]))
            else:
                out_list[nn] = None # no reads in exp and control 
    return norm_win

def GetRatio(norm_win):
    '''
    calculate ratio
    '''
    table = {}
    ave_win = float(CalAve(norm_win))
    for chrom in norm_win.keys():
        table[chrom] = []
        for value in norm_win[chrom]:
            if value is not None:
                table[chrom].append(float(value)/ave_win)
            else:
                table[chrom].append(1.0) # chromosome gap region
    return table

def WriteWig(output_file,hash_table,genome_table,resolution,step):
    '''
    write to wig file for comparison
    '''
    win_half = int(resolution*step/2)
    chrs = hash_table.keys()
    chrs.sort()
    fo = open(output_file+".wig","w")
    #wigTobigwig do not support the track line
    #print >>fo, "track type=wiggle_0 name=\""+ output_file + "\" description=\"" + output_file + "\" visibility=full autoScale=on color=50,150,255"
    for ii in range(len(chrs)):
        chr_list = hash_table[chrs[ii]]
        is_end = False
        print >>fo, "variableStep chrom=" + chrs[ii] + " span=" + str(resolution)
        for jj in range(len(chr_list)):
            if is_end:
                break
            if jj*resolution + win_half + resolution > genome_table[chrs[ii]]: # last window
                new_span = genome_table[chrs[ii]] - jj*resolution - win_half
                print >>fo, "variableStep chrom=" + chrs[ii] + " span=" + str(new_span)
                is_end = True
            print >>fo, "%d\t%.6f" % (jj*resolution + win_half, log(chr_list[jj],2))
    fo.flush()
    fo.close()

def WriteBedGraph(output_file,hash_table,genome_table,resolution,step):
    '''
    write to bedgraph file for comparison
    '''
    win_half = int(resolution*step/2)
    chrs = hash_table.keys()
    chrs.sort()
    fo = open(output_file+".bedgraph","w")
    for ii in range(len(chrs)):
        chr_list = hash_table[chrs[ii]]
        for jj in range(len(chr_list)):
            if jj*resolution + win_half > genome_table[chrs[ii]]:
                warning("window position %d is larger than chromosome size %d for chromosome %s" % (jj*resolution + win_half, genome_table[chrs[ii]], chrs[ii]))
                continue
            if jj*resolution + win_half + resolution > genome_table[chrs[ii]]: # last window
                new_span = jj*resolution + win_half + resolution - genome_table[chrs[ii]]
                print >>fo, "%s\t%s\t%s\t%.6f" % (chrs[ii], jj*resolution+win_half, jj*resolution+win_half+new_span, chr_list[jj])
            else:
                print >>fo, "%s\t%s\t%s\t%.6f" % (chrs[ii], jj*resolution+win_half, jj*resolution+win_half+resolution, chr_list[jj])
    fo.flush()
    fo.close()

def LoadGenome(filename):
    genome = {}
    try:
        fin = open(filename,"r")
    except IOError:
        error("Can't open file: %s" % (filename))
        exit()
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        size = int(row[1])
        genome[chrom] = size
    return genome

def ReportOptions():
    text = "# TSA-seq_normalize.py version: 0.2\n"
    text += "# experimental bam file: %s\n" % (args.bamexp)
    text += "# control bam file: %s\n" % (args.bamcon)
    text += "# size of sliding window: %d\n" % (args.win)
    text += "# resolution of output wig files: %d\n" % (args.res)
    text += "# genome chromosome size file: %s\n" % (args.genome)
    text += "# output file name prefix: %s" % (args.output)
    print >>sys.stderr, text

def Main():
    global args
    global chr_filter 
    chr_filter = re.compile('(random|chrM|hap|Un)').search
    args = parse_argument()
    # report parameters
    ReportOptions()
    # check whether the arguments are in possible range
    if (args.res > (args.win/2)) and args.res != args.win:
        sys.exit("resolution can not be larger than half size of sliding window")
    elif (args.res < 1):
        sys.exit("resultion can not be less than 1")
    else:
        if args.win % args.res == 0:
            step = args.win/args.res
        else:
            error("window size (dividend) must be divided by resolution (divisor) with no remaindar")
            exit(1)
    # load genome size
    genome_table = LoadGenome(args.genome)
    # pipeline begin
    logging("*************** analysis begin ***************") 
    # Step 1 put control bam reads into bins and aggregate into windows
    logging("Step1: count read number in control bam file")
    logging("Step1.1 put reads into bins.")
    Con_bins = BamToBin(args.bamcon, args.res)
    logging("Step1.2 roll up from bins to windows.")
    Con_win = BinToWin(Con_bins, step)
    # Step 2 put pulldown bam reads in bins and aggregate into windows
    logging("Step2: count read number in pulldown bam file")
    logging("Step2.1 put reads into bins.")
    Exp_bins = BamToBin(args.bamexp, args.res)
    logging("Step2.2 roll up from bins to windows.")
    Exp_win = BinToWin(Exp_bins, step)
    # release memory
    del Exp_bins
    del Con_bins
    # Step 3 Calculate Reads per bin adjusted for input control DNA variations for pulldown sample
    logging("Step3 normalize pulldown read count using input control")
    #WriteBedGraph(args.output+"_Exp", Exp_win, genome_table, args.res, step)
    #WriteBedGraph(args.output+"_Con", Con_win, genome_table, args.res, step)
    norm_win = WinNormalize(Exp_win, Con_win)
    # release memory
    del Exp_win
    del Con_win
    # Step 4 get Enrichment ratio
    logging("Step4 get enrichment ratio in log2 scale")
    ratio_win = GetRatio(norm_win)
    # Step 5 write to wig files
    logging("Step5 write ")
    logging("## Step5 Write log2(ratio) normalized read count to wig file")
    WriteWig(args.output,ratio_win,genome_table,args.res,step)
    if args.wig2bw is not None:
        logging("Step6 convert wig to bigwig")
        if args.wig2bw == 'sys':
            os.system("wigToBigWig %s %s %s" % (args.output+".wig", args.genome, args.output+".bw"))
        else:
            os.system(args.wig2bw + " %s %s %s" % (args.output+".wig", args.genome, args.output+".bw")) 
    #print "## Step6.2 write experimental group read count to wig file"
    #print "## Step6.3 write control group read count to wig file"
    #finish analysis
    logging("*************** analysis done ***************")

if __name__=="__main__":
    Main()
