#!/usr/bin/python
# Programmer : zocean
# Date: 
# Last-modified: 27 Jun 2017 11:08:35 PM

import os,sys,argparse
import math 
import numpy as np
from scipy import exp2
from bx.bbi.bigwig_file import BigWigFile
from TSA_utility import *

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bw',type=str,dest="bigwig",help="bigwig file")
    p.add_argument('-g','--genome',type=str,dest="genome",help="human genome file (two column, first col is chromosome name, second chromosome is chromosome size")
    p.add_argument('-w','--window',type=int,dest="window",help="window size used to aggregate data")
    p.add_argument('--smooth',action="store_true",help="if set, will smooth data, smooth is strongly recommended")
    p.add_argument('-o','--output',type=str,dest="output",help="output folder, output file will be written into that folder")
    p.add_argument('-n','--name',type=str,dest="name",help="output file prefix, bed file and wig file will be generated")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def Smooth(x,window_len=21,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s = np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def Main():
    global args
    args=ParseArg()
    bw = BigWigFile(open(args.bigwig))
    CheckFolderExist(args.output)
    fout = WriteToFile(args.output + '/' + args.name + '.bed')
    wout = WriteToFile(args.output + '/' + args.name + '.wig')
    genome = LoadGenome(args.genome)
    if args.smooth:
        logging("Options: turn on smooth mode")
    for chrom in SortGenome(genome):
        chrom_size = genome[chrom]
        logging("Process: %s\t%d" % (chrom,chrom_size))
        array = bw.get_as_array(chrom,0,chrom_size)
        invalid = np.isnan(array)
        array[invalid] = 0
        agg_array = []
        start = 0
        stop = args.window
        for nn in range(int(math.ceil(len(array)/float(args.window)))):
            if stop >= len(array):
                stop = len(array)
                agg_array.append(np.mean(array[start:stop]))
                break
            agg_array.append(np.mean(array[start:stop]))
            start += args.window
            stop += args.window
        agg_array = np.array(agg_array)
        if args.smooth:
            smooth_array = Smooth(agg_array)
        else:
            smooth_array = agg_array
        print >>wout, "variableStep chrom=%s span=%d" % (chrom,args.window)
        for nn,value in enumerate(smooth_array):
            if nn == 0: 
                print >>fout, "%s\t0\t%d\t%.6f" % (chrom,(nn+1)*args.window,float(value))
                print >>wout, "%d\t%.6f" % (nn+1,value) 
            elif nn == len(smooth_array) - 1:
                print >>fout, "%s\t%d\t%d\t%.6f" % (chrom,nn*args.window,chrom_size,float(value))
                print >>wout, "variableStep chrom=%s span=%d" % (chrom,chrom_size-((nn)*args.window))
                print >>wout, "%d\t%.6f" % (nn*args.window+1,float(value))
            else:
                print >>fout, "%s\t%d\t%d\t%.6f" % (chrom,nn*args.window,(nn+1)*args.window,float(value))
                print >>wout, "%d\t%.6f" % (nn*args.window+1,float(value))
    fout.flush()
    wout.flush()
    wig2bw = "wigToBigWig -clip %s %s %s" % ( args.output + '/' + args.name + '.wig', args.genome, args.output + '/' + args.name + '.bw')
    os.system(wig2bw)
    logging("Finish: TSA_smooth DONE!!!")

if __name__=="__main__":
    Main()
