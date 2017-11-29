#!/home/yzhan116/Software/Python-2.7.5/python-2.7.5
# Programmer : zocean
# Date: 
# Last-modified: 19 May 2016 14:57:36

import os,sys,argparse
import pysam
import pybedtools 

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--pulldown',type=str,dest="pulldown",help="pulldown bam file")
    p.add_argument('--control',type=str,dest="control",help="no primary antibody control bam file")
    p.add_argument('-m','--mode',type=str,dest="mode",help="get ratio in which mode, 'spike' means reads were mapped to spike-in bowtie index, 'genome' means reads were mapped to species whole genome")
    p.add_argument('--pos',type=str,dest="pos",help="if mode is 'genome', user need to provide spike-in location in six-column version bed file")
    p.add_argument('--tot_pulldown',type=int,dest="tot_pulldown",help="total number of reads (including unmapped reads) in pulldown sample, use this option to bypass the procedure to count the reads number from bam file, espeically when unmapped reads were removed from bam file")
    p.add_argument('--tot_control',type=int,dest="tot_control",help="total number of reads (including unmapped reads) in control sample")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

def GetRatioSpike():
    pulldown = pysam.Samfile(args.pulldown)
    control = pysam.Samfile(args.control)
    if args.tot_pulldown is not None:
        tot_pulldown = args.tot_pulldown
        tot_control = args.control
    else:
        tot_pulldown = pulldown.mapped + pulldown.unmapped
        tot_control = control.mapped + control.unmapped
    print >>sys.stderr, "Total number of reads in pulldown sample: %d" % (tot_pulldown)
    print >>sys.stderr, "Total number of reads in control sample: %d" % (tot_control)
    spike_pulldown = pulldown.mapped
    spike_control = control.mapped
    print >>sys.stderr, "Total number of reads mapped in spike-in region: %d in pulldown sample" % (spike_pulldown)
    print >>sys.stderr, "Total number of reads mapped in spike-in region: %d in control sample" % (spike_control)
    ratio = float(spike_control)/float(spike_pulldown)*float(tot_pulldown)/float(tot_control)
    print >>sys.stderr, "Ratio is %.6f" % (ratio)
    pulldown.close()
    control.close()

def GetRatioGenome():
    pulldown = pysam.Samfile(args.pulldown)
    control = pysam.Samfile(args.control)
    if args.tot_pulldown is not None:
        tot_pulldown = args.tot_pulldown
        tot_control = args.control
    else:
        tot_pulldown = pulldown.mapped + pulldown.unmapped
        tot_control = control.mapped + control.unmapped
    print >>sys.stderr, "Total number of reads in pulldown sample: %d" % (tot_pulldown)
    print >>sys.stderr, "Total number of reads in control sample: %d" % (tot_control)
    # get spike-in read within bed file
    pulldown_bam = pybedtools.BedTool(args.pulldown)
    control_bam = pybedtools.BedTool(args.control)
    spike_pulldown = pulldown_bam.intersect(args.pos,f=0.5).count()
    spike_control = control_bam.intersect(args.pos,f=0.5).count()
    print >>sys.stderr, "Total number of reads mapped in spike-in in pulldown sample: %d" % (spike_pulldown)
    print >>sys.stderr, "Total number of reads mapped in spike-in in control sample: %d" % (spike_control)
    ratio = float(spike_control)/float(spike_pulldown)*float(tot_pulldown)/float(tot_control)
    print >>sys.stderr, "Ratio is %.6f" % (ratio)
    pulldown.close()
    control.close()
    pybedtools.cleanup()

def CheckOptions():
    if args.mode == 'genome':
        if args.pos is None:
            print >>sys.stderr, "In 'genome' mode, user must provide spike-in location as three to six column bed file"
            exit(1)
    elif args.mode == 'spike':
        pass
    else:
        print >>sys.stderr, "Unknown mode: %s" % (args.mode)
    if args.tot_pulldown is not None and args.tot_control is None:
        print >>sys.stderr, "If 'tot_pulldown' is given, 'tot_control' must also be given"
        exit(1)
    if args.tot_control is not None and args.tot_pulldown is None:
        print >>sys.stderr, "If 'tot_control' is given, 'tot_pulldown' must also be given"
        exit(1)

def Main():
    global args
    args=ParseArg()
    CheckOptions()
    if args.mode == 'spike':
        GetRatioSpike()
    elif args.mode == 'genome':
        GetRatioGenome()
    else:
        print >>sys.stderr, "Unknown mode: %s" % (args.mode)
    
if __name__=="__main__":
    Main()
