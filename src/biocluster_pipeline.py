#!/usr/bin/python
# Programmer : Yang Zhang 
# Contact: zocean636@gmail.com
# Last-modified: 24 Jan 2019 15:20:08

import os,sys,argparse

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--conf',type=str,dest="conf",help="configure file")
    p.add_argument('--dry_run',dest="dry_run",action="store_true",help="set this parameter if just want to test environment. No real script will be procssed")
    if len(sys.argv) < 2:
        print p.print_help()
        exit(1)
    return p.parse_args()

class Run(object):
    def __init__(self):
        # global parameter
        self.genome_size = None
        self.bowtie2_index = None
        self.wig2bigwig = None
        self.norma = None
        # tool specific parameter
        self.bowtie_opt = None
        self.norma_opt = None
        # experiment specific parameter
        self.exp_name = None
        self.fastq_pulldown = None
        self.fastq_input = None
        self.label_pulldown = None
        self.label_input = None
        # output file
        self.output_folder = None
        self.out_bam_pulldown = None
        self.out_bam_input = None
        self.out_bowtie_log_pulldown = None
        self.out_bowtie_log_input = None
        self.out_bam_pulldown_rmdup = None
        self.out_bam_input_rmdup = None
        self.out_norma_output = None
        self.out_norma_log = None
    def build(self, conf):
        # check required parameters
        for parameter in ['genome_size', 'bowtie2_index', 'wig2bigwig', 'norma', 'exp_name', 'fastq_pulldown', 'fastq_input', 'label_pulldown', 'label_input', 'output_folder']:
            if conf.get(parameter, None) is None:
                print >>sys.stderr, "%s parameter not found" % (parameter)
                exit(1)
        # run initiation
        self.genome_size = conf['genome_size']
        self.bowtie2_index = conf['bowtie2_index']
        self.wig2bigwig = conf['wig2bigwig']
        self.norma = conf['norma']
        if conf.get('bowtie_opt', None) is not None and conf['bowtie_opt'] != "":
            self.bowtie_opt = conf['bowtie_opt']
        if conf.get('norma_opt', None) is not None and conf['norma_opt'] != "":
            self.norma_opt = conf['norma_opt']
        self.exp_name = conf['exp_name']
        self.fastq_pulldown = conf['fastq_pulldown'].split(',')
        self.fastq_input = conf['fastq_input'].split(',')
        self.label_pulldown = conf['label_pulldown']
        self.label_input = conf['label_input']
        # output 
        self.output_folder = conf['output_folder']
        if not os.path.isdir(self.output_folder):
            os.makedirs(self.output_folder)
        self.out_bam_pulldown = os.path.join(self.output_folder, '%s.bam' % (self.label_pulldown))
        self.out_bam_input = os.path.join(self.output_folder, '%s.bam' % (self.label_input))
        self.out_log_bowtie_pulldown = os.path.join(self.output_folder, 'log_bowtie_%s.txt' % (self.label_pulldown))
        self.out_log_bowtie_input = os.path.join(self.output_folder, 'log_bowtie_%s.txt' % (self.label_input))
        self.out_bam_pulldown_rmdup  = os.path.join(self.output_folder, '%s.rmdup.bam' % (self.label_pulldown))
        self.out_bam_input_rmdup = os.path.join(self.output_folder, '%s.rmdup.bam' % (self.label_input))
        self.out_norma_output = os.path.join(self.output_folder, self.exp_name)
        self.out_log_norma = os.path.join(self.output_folder, 'log_norma_%s' % (self.exp_name))
    def pipeline(self, dry_run = False):
        # 
        print >>sys.stderr, "# Start Norma pipeline for experiment: %s" % (self.exp_name)
        print >>sys.stderr, "# Step 1: Align the pulldown fastq to the reference genome" 
        cmd = self.__run_bowtie(self.bowtie2_index, self.fastq_pulldown, self.bowtie_opt, self.out_bam_pulldown, self.out_log_bowtie_pulldown)
        if dry_run:
            print >>sys.stderr, cmd
        else:
            os.system(cmd)
        print >>sys.stderr, "# Step 1: Alignment done: check %s for running log" % (self.out_log_bowtie_pulldown) 
        print >>sys.stderr, ""
        #
        print >>sys.stderr, "# Step 2: PCR duplicates removal for pulldown" 
        cmd = self.__run_rmdup(self.out_bam_pulldown, self.out_bam_pulldown_rmdup)
        if dry_run:
            print >>sys.stderr, cmd
        else:
            os.system(cmd)
        print >>sys.stderr, ""
        #
        print >>sys.stderr, "# Step 3: Align the input fastq to the reference genome"
        cmd = self.__run_bowtie(self.bowtie2_index, self.fastq_input, self.bowtie_opt, self.out_bam_input, self.out_log_bowtie_input)
        if dry_run:
            print >>sys.stderr, cmd
        else:
            os.system(cmd)
        print >>sys.stderr, "# Step 3: Aligment done: check %s for running log" % (self.out_log_bowtie_input)
        print >>sys.stderr, ""
        #
        print >>sys.stderr, "# Step 4: PCR duplicates removal for input"
        cmd = self.__run_rmdup(self.out_bam_input, self.out_bam_input_rmdup)
        if dry_run:
            print >>sys.stderr, cmd
        else:
            os.system(cmd)
        print >>sys.stderr, ""
        #
        print >>sys.stderr, "# Step 5: Run Norma to get the TSA-seq signal"
        cmd = self.__run_norma(self.norma, self.out_bam_pulldown_rmdup, self.out_bam_input_rmdup, self.out_norma_output, self.out_log_norma, self.norma_opt, self.wig2bigwig, self.genome_size)
        if dry_run:
            print >>sys.stderr, cmd
        else:
            os.system(cmd)
        print >>sys.stderr, "# Step 5: Norma done: check %s for running log" % (self.out_log_norma)
    def __run_bowtie(self, genome_index, fastq_list, other_opt, output_file, log_file):
        if other_opt is not None:
            cmd =  "bowtie2 %s -x %s -U %s 2>%s | samtools view -S -bh - | samtools sort -o %s" % (other_opt, genome_index, ' '.join(fastq_list), log_file, output_file)
        else:
            cmd =  "bowtie2 -x %s -U %s 2>%s | samtools view -bS " % (genome_index, ' '.join(fastq_list), log_file, output_file) 
        cmd += '\n' + "samtools index %s" % (output_file)
        return cmd
    def __run_rmdup(self, input_bam, output_bam):
        cmd = "samtools rmdup -s %s %s" % (input_bam, output_bam)
        cmd += '\n' + "samtools index %s" % (output_bam)
        return cmd
    def __run_norma(self, norma_script, pulldown_bam, input_bam, output, log, other_opt, wig2bigiwg, genome_size):
        if other_opt is not None:
            cmd = "%s %s -g %s -e %s -c %s --wig2bw %s -o %s 2>&1 >%s" % (norma_script, other_opt, genome_size, pulldown_bam, input_bam, wig2bigiwg, output, log)
        else:
            cmd = "%s -g %s -e %s -c %s --wig2bw %s -o %s 2>&1 >%s" % (norma_script, other_opt, genome_size, pulldown_bam, input_bam, wig2bigiwg, output, log)
        return cmd

def parse_conf(filename):
    fin = open(filename, 'r')
    table = {}
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split('=')
        table[row[0].strip()] = row[1].strip()
    fin.close()
    return table

def main():
    global args
    args = parse_arg()
    # parse the configure table
    conf = parse_conf(args.conf)
    print >>sys.stderr, "# parse parameters done"
    # build Run
    TSA_seq_run = Run()
    TSA_seq_run.build(conf)
    print >>sys.stderr, "# build run done"
    # run
    print >>sys.stderr, "# run pipeline"
    TSA_seq_run.pipeline(args.dry_run)
    
if __name__=="__main__":
    main()
