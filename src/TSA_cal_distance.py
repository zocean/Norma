#!/usr/bin/python
# Programmer : zocean
# Date: 
# Last-modified: 13 Feb 2019 15:12:18

import os,sys,argparse
import math
import tempfile

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-b','--bw',type=str,dest="bw",help="big wig file")
    p.add_argument('-y0',type=float,dest="y0",help="y0")
    p.add_argument('-A',type=float,dest="A",help="A")
    p.add_argument('-R0',type=float,dest="R0",help="R0")
    p.add_argument('-g','--genome',type=str,dest="genome",help="genome size file")
    p.add_argument('-o','--output',type=str,dest="output",help="output bigwig file name")
    p.add_argument('--bw2wig',dest = "bw2wig", type = str, help = "path of bigWigToWig (full path eg. /home/zocean/bigWigToWig of bigWigToWig program, 'sys' means bigWigToWig is the system path")
    p.add_argument('--wig2bw',dest = "wig2bw", type = str, help = "path of wigToBigWig (full path eg. /home/zocean/wigToBigWig of wigToBigWig program, 'sys' means wigToBigWig is the system path")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

"""Basic functions to handle IO"""
def ReadFromFile(fin_name):
    try:
        fin = open(fin_name,"r")
    except IOError:
        error("Can't open file:%s to read" % (fin_name))
        exit(1)
    return fin

def WriteToFile(fout_name):
    if fout_name == 'stdout' or not fout_name:
        fout = sys.stdout
    else:
        try:
            fout = open(fout_name,"w")
        except IOError:
            error("Can't open file:%s to write. Use stdout instead." % (fout_name))
            fout = sys.stdout
    return fout

def Main():
    global args
    args=ParseArg()
    out = tempfile.NamedTemporaryFile(delete=False)
    wig_file = tempfile.NamedTemporaryFile(delete=False)
    if args.bw2wig == 'sys':
        os.system("bigWigToWig %s %s" % (args.bw, wig_file.name))
    else:
        os.system("%s %s %s" % (args.bw2wig, args.bw, wig_file.name))
    for line in ReadFromFile(wig_file.name):
        if line.strip().startswith('#') or line.strip() == '':
            print >>out, line.strip()
        row = line.strip().split()
        if row[0] in ['variableStep','track','fixedStep']:
            print >>out, line.strip()
        else:
            pos = row[0]
            if row[1] in ['0','-0']: # skip gap region
                print >>out, pos + '\t' + '0'
                continue
            value = float(row[1])
            if pow(2,value) - args.y0 < 1e-6:
                dis = max(0.0,math.log(1e-6/args.A)/args.R0)
            else:
                dis = max(0.0,math.log((pow(2,value) - args.y0)/args.A)/args.R0)
            print >>out, pos + '\t' + '{:.4f}'.format(dis)
    out.flush()
    if args.wig2bw == 'sys':
        os.system("wigToBigWig %s %s %s" % (out.name, args.genome, args.output))
    else:
        os.system("%s %s %s %s" % (args.wig2bw, out.name, args.genome, args.output))
    os.system("rm %s" % (wig_file.name))
    os.system("rm %s" % (out.name))
 
if __name__=="__main__":
    Main()
