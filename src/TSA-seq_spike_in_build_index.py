#!/home/yzhan116/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu 
# Last-modified: 18 May 2016 21:20:50

import os,sys,argparse

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-f','--fasta',type=str,dest="fasta",help="build bowtie1(2) index")
    p.add_argument('-m','--mode',type=str,dest="mode",help="bowtie1 or bowtie2")
    p.add_argument('--bowtie',type=str,dest="bowtie",help="location of bowtie-build program version should be compatiable with option '--mode'")
    p.add_argument('-o','--output',type=str,dest="output",help="output index file base")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def Main():
    global args
    args=ParseArg()
    if not os.path.isfile(args.fasta):
        print >>sys.stderr, "Can't open fasta file: %s\n" % (args.fasta)
        exit(1)
    if args.mode == 'bowtie1':
        if args.bowtie is not None:
            os.system("%s %s %s" % (args.bowtie, args.fasta, args.output))
        else:
            os.system("bowtie-build %s %s" % (args.fasta, args.output))
    elif args.mode == 'bowtie2':
        if args.bowtie is not None:
            os.system("%s %s %s" % (args.bowtie, args.fasta, args.output))
        else:
            os.system("bowtie2-build %s %s" % (args.fasta, args.output))
    
if __name__=="__main__":
    Main()
