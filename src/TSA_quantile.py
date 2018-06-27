#!/usr/bin/python
# Programmer : zocean
# Date: 
# Last-modified: 01 Nov 2014 04:57:03 PM

import os,sys,argparse
import numpy as np
import tempfile
from TSA_utility import *

"""
This script will read wig file and huamn gap region file (N), then it will calculate the quantile of the genome wide value and report the segment bed file with different quantile range
"""

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-w','--wig',type=str,dest="wig",help="wig file")
    p.add_argument('-q','--quantile',type=int,dest="quantile",help="quantile segment number, 10 means divide 100%% into 10 segment")
    p.add_argument('--quantilelist',type=str,dest="quantilelist",help="quantile list, eg, 0,10,20,25,45,50,100, the first one must be 0 and the last one must be 100 and should be in an ascending order, the delimiter must be comma")
    p.add_argument('-g','--gap',type=str,dest="gap",help="gap region file")
    p.add_argument('-o','--output',type=str,dest="output",help="output file name")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def LoadGap(fin_name):
    table = {}
    for line in ReadFromFile(args.gap):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])+1 # convert to 1-base
        stop = int(row[2])
        if chrom not in table.keys():
            table[chrom] = []
        table[chrom].append((start,stop))
    return table

def FindQuantile(wig_file,quantile_num,gap,is_list=False):
    """find quantile segment"""
    chrom = None
    span = None
    start = None
    array = []
    num = 1
    for line in ReadFromFile(wig_file):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed" % (num))
        num += 1
        if row[0] == 'variableStep': # new track chrom begin
            chrom = row[1].split('=')[1]
            span = int(row[2].split('=')[1])
            start = None
        else:
            is_overlap = False
            try:
                start = int(row[0])
                value = float(row[1])
            except ValueError:
                warning("Warning: Read wig error. Unknown line format for line:%s" % (line.strip()))
            for mm in range(len(gap[chrom])):
                if start + span <= gap[chrom][mm][0]:
                    continue
                if start > gap[chrom][mm][1]:
                    continue
                if float(min((start + span),gap[chrom][mm][1]) - gap[chrom][mm][0])/float(span) < 0.75:
                    continue
                is_overlap = True
            if not is_overlap:
                array.append(value)
    if is_list is False:
        quantile_list = (range(0,100,100/quantile_num) + [100])[1:]
    else:
        quantile_list = sorted([int(item.strip()) for item in quantile_num.split(',')])
        if 0 not in quantile_list or 100 not in quantile_list:
            error("0 and 100 must exist in quantilelist")
            exit(1)
        quantile_list = quantile_list[1:]
    return CombineList(quantile_list, np.percentile(array,quantile_list))

def CombineList(list_a,list_b):
    result = []
    if len(list_a) != len(list_b):
        error("Can't combine two lists with different length")
        exit(1)
    for nn in range(len(list_a)):
        result.append((list_a[nn],list_b[nn]))
    return result

def WriteToSeg(out,wig_file,gap,quantile):
    num = 1
    for line in ReadFromFile(wig_file):
        if line.strip().startswith('#') or line.strip() == '' or line.strip() == '\n':
            continue
        row = line.strip().split()
        if row[0] == 'track':
            continue
        if num % 10000 == 0:
            logging("Process: %d line processed" % (num))
        num += 1
        if row[0] == 'variableStep':
            chrom = row[1].split('=')[1]
            span = int(row[2].split('=')[1])
            start = None
        else:
            is_overlap = False
            try:
                start = int(row[0])
                value = float(row[1])
            except ValueError:
                warning("Parsing wig file error: Unknown line format for line: %s" % (line.strip()))
            for mm in range(len(gap[chrom])):
                if start + span <= gap[chrom][mm][0]:
                    continue
                if start > gap[chrom][mm][1]:
                    continue
                if float(min((start + span),gap[chrom][mm][1]) - gap[chrom][mm][0])/float(span) < 0.75:
                    continue
                print >>out, chrom + '\t' + str(start-1) + '\t' + str(start+span-1) + '\t' + 'N' + '\t' + '{:.6f}'.format(value)
                is_overlap = True
                break
            if is_overlap:
                continue
            else:
                for nn in range(len(quantile)):
                    if value < quantile[nn][1]:
                        break
                print >>out, chrom + '\t' + str(start-1) + '\t' + str(start+span-1) + '\t' + str(quantile[nn][0]) + '\t' + '{:.6f}'.format(value)

def MergeSegment(filename,outname):
    """Merge adjacent bed region with same TSA group"""
    out = WriteToFile(outname)
    group = None
    chrom = None
    for line in ReadFromFile(filename):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        if row[0] != chrom:
            if chrom is not None:
                print >>out, "%s\t%s\t%s\t%s\t%.6f" % (chrom,start,stop,group,np.mean(value))
            chrom = row[0]
            start = row[1]
            stop = row[2]
            group = row[3]
            value = [float(row[4])]
        else:
            if row[3] != group:
                print >>out, "%s\t%s\t%s\t%s\t%.6f" % (chrom,start,stop,group,np.mean(value))
                chrom = row[0]
                start = row[1]
                stop = row[2]
                group = row[3]
                value = [float(row[4])]
            else:
                stop = row[2]
                value.append(float(row[4]))
    print >>out, "%s\t%s\t%s\t%s\t:%.6f" % (chrom,start,stop,group,np.mean(value))

def Main():
    global args
    args=ParseArg()
    out = WriteToFile(args.output)
    gap_table = LoadGap(args.gap)
    if args.quantile is not None:
        quantile_seg = FindQuantile(args.wig,args.quantile,gap_table)
    elif args.quantilelist is not None:
        quantile_seg = FindQuantile(args.wig,args.quantilelist,gap_table,True)
    else:
        error("quantile or quantilelist must be set")
        exit(1)
    logging("Results: Dividing upper bound is " + '\t'.join('{:.3f}'.format(item[1]) for item in quantile_seg))
    WriteToSeg(out,args.wig,gap_table,quantile_seg)
    out.flush()
    SortBed(args.output)
    mergename = "".join(args.output.split('.')[:-1]) + "_merge.bedgraph"
    MergeSegment(args.output,mergename)
    CreateTabix(args.output,args.output + '.bgz','bed')
    CreateTabix(mergename,mergename + '.bgz', 'bed')
    logging("Finish: TSA_quantile DONE!!!")

if __name__=="__main__":
    Main()
