#!/usr/bin/python
# Programmer : zocean
# Last-modified: 01 Feb 2023 04:15:16 PM

import os,sys,argparse
from math import log
import pysam
import re

def parse_argument():
    ''' This Function Parse the Argument '''
    p = argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog = 'Library dependency : pysam')
    p.add_argument('-v','--version',action = 'version', version = '%(prog)s 0.5')
    p.add_argument('-r','--res',dest = "resolution",type = int, default = 10, help = "wig file resultion(bp)")
    p.add_argument('-w','--win',dest = "window",type = int, default = 1000,  help = "Sliding window size(bp)")
    p.add_argument('-e','--exp',dest = "BamExp",type = str, help = "experimental group bam file")
    p.add_argument('-c','--con',dest = "BamCon",type = str, help = "control group bam file")
    p.add_argument('-l','--length',dest = "readlength",type = int, default = 100, help = " sequencing read length")
    p.add_argument('-o','--wig_name',dest = "wig_name",type = str, help = "the name of wig file without .wig") 
    if len(sys.argv)<2:
        sys.exit(p.print_help())
    return p.parse_args()

def filterPick(list,filter):
    '''
    return the items match the filter pattern
    '''
    return [ ( l ) for l in list for m in (filter(l),) if m]

def delete_key(dic, key_list):
    for key in key_list:
        del dic[key]
    return dic

def BamToBin(bamfile,binsize,norm_factor):
    '''
    This function count read number ( based on the mid of a read) into bins (an analogical word to resolution, with minimal number of 1bp.
    '''
    # creat a hash dictionary for chromosome name in the header of bamfile
    bin_hash = {}
    # read the indexed bam file using pysam
    samfile = pysam.Samfile(bamfile,"rb")
    # read chromosome name (chr1 et.al) and the length into a dict
    chr_dict = {}
    for i in range(len(samfile.references)):
        chr_dict[samfile.references[i]] = samfile.lengths[i]
    # remove user define chromosome list pattern 
    chr_removed_list = filterPick(chr_dict.keys(),chr_filter)
    chr_dict = delete_key(chr_dict,chr_removed_list)
    chrs = chr_dict.keys()
    lengths = []
    for i in range(len(chr_dict)):
        lengths.append(chr_dict[chrs[i]])
    # innitiation total number of bins
    total_bins_num = 0
    # normalization factor goes here
    for i in range(len(chrs)):
        bins = lengths[i]/binsize
        if lengths[i]%binsize != 0:
            bins += 1
        total_bins_num += bins
        # put 0 at in bins, innitiation bin count
        bin_hash[chrs[i]] = [0.000000 for row in range(bins)]
    # look through every read in bamfile
    n = 0
    for alignment in samfile:
        if alignment.tid < 0:
            continue
        pos = alignment.pos + alignment.alen
        chr_name = samfile.getrname(alignment.tid)
        if chr_name in chrs:
            bin_hash[chr_name][pos/binsize] += 1.0/norm_factor
        # add processing information
        n += 1
        if n % 100000 == 0:
            print ".",
	    #print "> %d reads processed\r" % n,
    print ""
    samfile.close()
    return bin_hash

def BinToWin(bin_hash,step):
    '''
    This function roll up read count from bins into window. The second argument "step" tells the function how many bins would be merged together
    '''
    # read chromosome name from hash dictionary
    chrs = bin_hash.keys()
    # create a hash dictionary for read count in window
    win_hash = {}
    for i in range(len(chrs)):
        # innitiation window hash dictionary
        win_no = len(bin_hash[chrs[i]]) - step + 1
        #print "#bin is %d, #step is %d, #win is %d" %(len(bin_hash[chrs[i]]), step, win_no)
        win_hash[chrs[i]] = [0 for n in range(win_no)]
        chr_list = bin_hash[chrs[i]]
        buffer = sum(chr_list[0:step])
        win_hash[chrs[i]][0] = buffer
        for j in range(1,win_no):
            win_hash[chrs[i]][j] = buffer + chr_list[j-1+step] - chr_list[j-1]
            # update buffer
            buffer = win_hash[chrs[i]][j]
        print ".",
	#print "> %d chrosome done\n" % (i+1),
    print ""
    return win_hash

def WinCompare(first_hash,second_hash):
    '''
    calculate log2(ratio) between experimental count to control count table
    '''
    chrs_A = first_hash.keys()
    chrs_B = second_hash.keys()
    if len(chrs_A) != len(chrs_B):
        print "The chromsome number between two bam files are different, are you sure continue? (Yes/No)"
        go_on = raw_input("> ")
        if go_on == "No":
            sys.exit("interupt by user")
    chrs = list(set(chrs_A) & set(chrs_B))
    compare_hash = {}
    for i in range(len(chrs)):
        compare_hash[chrs[i]] = [0 for row in range(len(first_hash[chrs[i]]))]
        chr_list = compare_hash[chrs[i]]
        A_list = first_hash[chrs[i]]
        B_list = second_hash[chrs[i]]
        for j in range(len(chr_list)):
            chr_list[j] = ((A_list[j] + 0.01) / (B_list[j] + 0.01))
        print ".",
    print ""
    return compare_hash

def WriteWig(output_file,hash_table,resolution,step):
    '''
    write to wig file for comparison
    '''
    chrs = hash_table.keys()
    chrs.sort()
    fo = open(output_file+".wig","w")
    #wigTobigwig do not support the track line
    #print >>fo, "track type=wiggle_0 name=\""+ output_file + "\" description=\"" + output_file + "\" visibility=full autoScale=on color=50,150,255"
    for i in range(len(chrs)):
        chr_list = hash_table[chrs[i]]
        print >>fo, "variableStep chrom=" + chrs[i] + " span=" + str(resolution)
        for j in range(len(chr_list)):
            #if chr_list[j] != 0.0:
            print >>fo, "%s\t%.4f" % (str(j*resolution + int(resolution*step/2)),log(chr_list[j],2))
    fo.close()
    # write chromesome name and length to file ( can be used as parameter of wigTobigwig
    chrs = hash_table.keys()
    fo = open(output_file+"_chrom.sizes","w")
    for i in range(len(chrs)):
	sizes = len(hash_table[chrs[i]])*resolution + int(resolution*step)
        chrom = chrs[i]
        print >>fo, "%s\t%d" % (chrom, sizes)
    fo.close()
    # convert wig2bigwig
    os.system("wigToBigWig %s %s %s" % (output_file+'.wig', output_file+"_chrom.sizes", output_file+'.bw'))

def WriteWig_input(output_file,hash_table,resolution,step):
    '''
    write to wig file for input
    '''
    chrs = hash_table.keys()
    chrs.sort()
    fo = open(output_file+".wig","w")
    #wigTobigwig do not support the track line
    #print >>fo, "track type=wiggle_0 name=\""+ output_file + "\" description=\"" + output_file + "\" visibility=full autoScale=on color=50,150,255"
    for i in range(len(chrs)):
        chr_list = hash_table[chrs[i]]
        print >>fo, "variableStep chrom=" + chrs[i] + " span=" + str(resolution)
        for j in range(len(chr_list)):
            if chr_list[j] > 0.0:
                print >>fo, "%s\t%.4f" % (str(j*resolution + int(resolution*step/2)),log(chr_list[j],2))
	    else:
		print >>fo, "%s\t%.4f" % (str(j*resolution + int(resolution*step/2)),-10)
    fo.close()
    # write chromesome name and length to file ( can be used as parameter of wigTobigwig
    chrs = hash_table.keys()
    fo = open(output_file+"_chrom.sizes","w")
    for i in range(len(chrs)):
	sizes = len(hash_table[chrs[i]])*resolution + int(resolution*step)
        chrom = chrs[i]
        print >>fo, "%s\t%d" % (chrom, sizes)
    fo.close()

def BamMetrics(bamfile,readlength):
    '''
    This function will get the number of mappable read and the percentage of genome covered by at least one read (this number is not accurate but an appopriation in order to increase the computation speed
    '''
    samfile = pysam.Samfile(bamfile,"rb")
    #coverage_table = BamToBin(bamfile,readlength,1)
    #chrs = coverage_table.keys()
    chr_dict = {}
    for i in range(len(samfile.references)):
        chr_dict[samfile.references[i]] = samfile.lengths[i]
    # remove user define chromosome list pattern
    chr_removed_list = filterPick(chr_dict.keys(),chr_filter)
    chr_dict = delete_key(chr_dict,chr_removed_list)
    chrs = chr_dict.keys()
    lengths = []
    for i in range(len(chr_dict)):
        lengths.append(chr_dict[chrs[i]])   
    metrics = {}
    metrics["mappable_reads"] = 0
    metrics["#chromosome"] = 0
    metrics["#chromosome"] = 0
    metrics["genome_size"] = 0
    idx_stats = samfile.get_index_statistics()
    for idx in idx_stats:
        if idx.contig in chrs:
            metrics["mappable_reads"] += idx.mapped
            metrics["#chromosome"] += 1
            metrics["genome_size"] += chr_dict[idx.contig]
    #change normlization method
    #estimate_coverage = 0
    #for i in range(len(chrs)):
    #    chr_list = coverage_table[chrs[i]]
    #    for j in range(len(chr_list)):
    #        if chr_list[j] > 0:
    #            estimate_coverage += 1
    #metrics["mappable_bp"] = estimate_coverage * readlength 
    samfile.close()
    return metrics

def Main():
    global args
    global chr_filter 
    chr_filter = re.compile('(random|chrM|hap|Un)').search
    args = parse_argument()
    # check whether the arguments are in possible range
    if (args.resolution > (args.window/2)):
        sys.exit("resolution can not be larger than half size of sliding window")
    elif (args.resolution < 1):
        sys.exit("resultion can not be less than 1")
    else:
        step = args.window/args.resolution
    # analysis begin
    print "*************** analysis begin ***************"
    # Step 1 get experimental group bam alignment metics
    print "## Step1 get experimental group bam alignment metric."
    metrics_Exp = BamMetrics(args.BamExp, args.readlength)
    norm_factor_Exp = float(metrics_Exp["mappable_reads"]) * float(args.window) / float(metrics_Exp["genome_size"])
    print "## Experimental group: normalization factor: %f" % norm_factor_Exp
    print "## Experimental group: %d read is mappable" % metrics_Exp["mappable_reads"]
    print "## Experimental group: genome size is %d" % metrics_Exp["genome_size"]
    #print "## Experimental group: %d base pair is covered by at least 1 read" % metrics_Exp["mappable_bp"]

    # Step 2 get control group bam alignment metics
    print "\n## Step2 get control group bam alignment metric."
    metrics_Con = BamMetrics(args.BamCon, args.readlength)
    norm_factor_Con = float(metrics_Con["mappable_reads"]) * float(args.window) / float(metrics_Con["genome_size"])
    print "## Control group: normalization factor: %.6f" % norm_factor_Con
    print "## Control group: %d read is mappable" % metrics_Con["mappable_reads"]
    print "## Control group: genome size is %d" % metrics_Con["genome_size"]
    #print "## Control group: %d base pair is covered by at least 1 read" % metrics_Con["mappable_bp"]
    
    # Step 3 put experimental group reads in bins and count reads in windows
    print "\n## Step3 count experimental group read number using sliding window."
    print "## Step3.1 put reads into bins."
    Exp_bins = BamToBin(args.BamExp, args.resolution, norm_factor_Exp)
    print "## Step3.2 roll up from bins to windows."
    Exp_window = BinToWin(Exp_bins, step)
    print "                     "
    del Exp_bins

    # Step 4 put control group reads in bins and count reads in windows
    print "## Step4 count experimental group read number using sliding window."
    print "## Step4.1 put reads into bins."
    Con_bins = BamToBin(args.BamCon, args.resolution, norm_factor_Con)
    print "## Step4.2 roll up from bins to windows."
    Con_window = BinToWin(Con_bins, step)
    print "                     "
    del Con_bins

    # Step 5 compare two tables
    print "## Step5 Compare normalized read count between two bam files"
    compare_window = WinCompare(Exp_window, Con_window)
    print "                     "
    #del Exp_window
    #del Con_window

    # Step 6 write to wig files
    print "## Step6.1 Write log2(ratio) normalized read count to wig file"
    WriteWig(args.wig_name,compare_window,args.resolution,step)
    #print "## Step6.2 write experimental group read count to wig file"
    #WriteWig_input(args.BamExp.split('.')[0],Exp_window,args.resolution,step)
    #print "## Step6.3 write control group read count to wig file"
    #WriteWig_input(args.BamCon.split('.')[0],Con_window,args.resolution,step)
    
    #finish analysis
    print "*************** analysis done ***************"
if __name__=="__main__":
    Main()


