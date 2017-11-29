#!/home/yzhan116/Software/Python-2.7.5/python-2.7.5
# Programmer : zocean
# Date: 
# Last-modified: 30 Nov 2015 04:08:47

import os,sys,argparse

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',type=str,dest="input",help="input filename (TSA-seq_anno output file)")
    p.add_argument('-o','--output',type=str,dest="output",help="output distance file used for R")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def ReadFromFile(filename):
    try:
        fin = open(filename,'r')
    except IOError:
        error("Can't open file: %s to read" % (filename))
        exit(1)
    return fin

def WriteToFile(filename):
    try:
        fout = open(filename,'w')
    except IOError:
        warning("Can't open file: %s to write, use stdout instead" % (filename))
        fout = sys.stdout
    return fout

def Main():
    global args
    args=ParseArg()
    out = WriteToFile(args.output)
    print >>out, "chrom\tR\tdiff\tabs"
    num = 0
    header = {}
    for line in ReadFromFile(args.input):
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        if num == 0: # load header
            for nn in range(len(row)):
                header[row[nn]] = nn
            num += 1
            continue
        default = float(row[header['default']])
        chrom = row[header['chrom']]
        for item in sorted(header.keys()):
            if item in ['chrom','start','stop','default']:
                continue
            else:
                value = float(row[header[item]]) - default
                print >>out, "%s\t%s\t%.6f\t%.6f" % (chrom,item.replace('R',''), value, abs(value))
 
if __name__=="__main__":
    Main()
