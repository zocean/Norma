#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 19 Feb 2018 19:33:30

import os,sys,argparse
import math
import numpy as np
import tempfile
from operator import attrgetter
import pysam
import pyfasta
import pybedtools
from itertools import product
from TFBS_Evo.my_utility import *

##############
# Custom Class
##############

class BED(object):
    def __init__(self, chrom, start, stop):
        '''basic bed ojbect'''
        # public
        self.chrom = chrom
        self.start = start
        self.stop = stop
    def __repr__(self):
        return "%s\t%d\t%d" % (self.chrom, self.start, self.stop)

class Region(BED):
    def __init__(self, chrom, start, stop):
        '''bed or tfbs region'''
        # public
        BED.__init__(self, chrom, start, stop)
        self.mid = (self.start+self.stop)/2
        self.verbose = 0
        # private
        self.center = [BED(chrom, start, stop)] # list of BED object for initial region after we remove regions overlapping with the excluding region
        self.anno = {} # annotation table
        self.label = [] # label list
    def region_init(self, exclude_tabix_list, genome_size):
        '''
        get the corrected center region
        --exclude_tabix_list    exluding region in tabix format, removed from center regions
        '''
        self.center = subtract_bed(self.chrom, self.start, self.stop, exclude_tabix_list)
    def get_anno(self, anno, label, mode, genome_size):
        '''
        anno       python annotation object, eg. tabix object, bigwig object etc.
        label      annotation name
        mode       treatment for annotation object
        genome_size     chromosome size table
        '''
        self.__add_label(label.split(':')[0])
        if mode == 'coverage':
            self.anno[label] = self.__get_coverage(anno)
        elif mode == 'coverage_per':
            self.anno[label] = self.__get_coverage_per(anno)
        elif mode == 'signal_mean':
            self.anno[label] = self.__get_mean_signal(anno)
        elif mode == 'signal_sum':
            self.anno[label] = self.__get_sum_signal(anno)
        elif mode == 'count':
            self.anno[label] = self.__get_center_count(anno)
        elif mode == 'anno':
            self.anno[label] = self.__get_anno(anno)
        elif mode == 'dnase_count':
            self.anno[label] = self.__get_dnase_cov(anno)
        elif mode == 'dnase_ave':
            self.anno[label] = self.__get_dnase_ave(anno)
        elif mode == 'gc':
            self.anno[label] = self.__get_gc(anno)
        elif mode == 'gc_count':
            self.anno[label] = self.__get_gc_count(anno)
        elif mode == 'size':
            self.anno[label] = self.__get_center_size(anno)
        elif mode == 'nearby':
            self.anno[label] = self.__get_nearby_count(anno)
        elif mode == 'exclude':
            self.anno[label] = self.__get_exclude_label(anno)
        elif mode == 'tsi_peak':
            self.anno[label] = self.__get_tsi_peak(anno)
        elif mode == 'tsi_signal':
            self.anno[label] = self.__get_tsi_signal(anno)
        elif mode == 'tfbs':
            label_name,tfbs_list = label.split(':')
            self.anno[label_name] = self.__get_tfbs_anno(anno, tfbs_list.split(','))
        else:
            print >>sys.stderr, "Unknown mode: %s" % (mode)
            exit(1)
    def __add_label(self, label):
        if label in self.label:
            if self.verbose > 1:
                print >>sys.stderr, "Colname: %s already exists. Update it" % (label)
            else:
                pass # label already there, will not change it
        else:
            self.label.append(label)
    def __get_coverage(self, tabix_dbi):
        '''get the number of base pair of region intersecting with annotation region'''
        o_array = []
        for bed in self.center:
            try:
                result = tabix_dbi.query(bed.chrom, bed.start, bed.stop)
            except:
                continue
            if result is not None:
                for overlap in result:
                    o_start = int(overlap[1])
                    o_stop = int(overlap[2])
                    for nn in range(max(o_start, bed.start), min(o_stop, bed.stop)):
                        o_array.append(nn)
        o_sum = len(set(o_array))
        return str(o_sum)
    def __get_coverage_per(self, tabix_dbi):
        '''get the percentage of DHS intersecting with annotation region'''
        o_array = []
        total_size = 0
        for bed in self.center:
            total_size += bed.stop - bed.start
            try:
                result = tabix_dbi.query(bed.chrom, bed.start, bed.stop)
            except:
                continue
            if result is not None:
                for overlap in result:
                    o_start = int(overlap[1])
                    o_stop = int(overlap[2])
                    for nn in range(max(o_start, bed.start), min(o_stop, bed.stop)):
                        o_array.append(nn)
        o_sum = len(set(o_array))
        if o_sum == 0:
            return '{:.6f}'.format(0.0)
        if total_size > 0:
            return '{:.6f}'.format(float(o_sum)/float(total_size))
        else:
            return 'NA'
    def __get_mean_signal(self, bw):
        '''get the mean signal value within region'''
        array = np.array([])
        for bed in self.center:
            signal = bw.get_as_array(bed.chrom, bed.start, bed.stop)
            if signal is not None:
                array = np.append(array, signal)
        if len(array) < 1: # resolve the warning issue with numpy when array with zero element
            return 'NA'
        value = np.nanmean(array)
        if np.isnan(value):
            return 'NA'
        else:
            return '{:.6f}'.format(value)
    def __get_sum_signal(self, bw):
        '''get sum signal over region'''
        array = np.array([])
        for bed in self.center:
            signal = bw.get_as_array(bed.chrom, bed.start, bed.stop)
            if signal is not None:
                array = np.append(array, signal)
        value = np.nansum(array)
        if np.isnan(value):
            return 'NA'
        else:
            return '{:.6f}'.format(value)
    def __get_center_count(self, tabix_dbi):
        '''get the number of region overlapping region'''
        num = 0
        for bed in self.center:
            try:
                result = tabix_dbi.query(bed.chrom, bed.start, bed.stop)
                for overlap in result:
                    num += 1
            except:
                continue
        return str(num)
    def __get_anno(self, tabix_dbi):
        '''get the segment feature (4th column is group label) in tabix that most overlapping with region'''
        result = tabix_dbi.query(self.chrom, self.start, self.stop)
        if result is None:
            return 'NA'
        else:
            o_len = 0
            o_group = None
            for overlap in result:
                o_start = int(overlap[1])
                o_stop = int(overlap[2])
                name = overlap[3]
                size = min(o_stop, self.stop) - max(o_start, self.start)
                if size > o_len:
                    o_len = size
                    o_group = name
            if o_group is None:
                o_group = 'NA'
            return o_group
    def __get_bam_cov(self, bam):
        '''
        get the number of read count within region
        --bam   SamHandler object below
        '''
        cut_count = 0
        for bed in self.center:
            cut_for, cut_rev = bam.get_count(bed.chrom, bed.start, bed.stop)
            cut_count += sum(cut_for) + sum(cut_rev)
        return str(cut_count)
    def __get_bam_ave(self, bam):
        '''
        get the average read count within region
        --bam   SamHandler object below
        '''
        cut_count = 0
        total_size = 0
        for bed in self.center:
            total_size += bed.stop - bed.start
            cut_for, cut_rev = bam.get_count(bed.chrom, bed.start, bed.stop)
            cut_count += sum(cut_for) + sum(cut_rev)
        if total_size < 1:
            return 'NA'
        else:
            return '{:.6f}'.format(float(cut_count)/float(total_size))
    def __get_gc(self, genome):
        '''
        get the gc content within DHS
        --genome    pysam or pyfasta object
        '''
        table = {'A':0,'C':0,'G':0,'T':0}
        for bed in self.center:
            seq = get_seq_from_bed(bed.chrom, bed.start, bed.stop, genome)
            for nt in ['A','C','G','T']:
                table[nt] += seq.count(nt)
        total = sum(table.values())
        if total < 1:
            return 'NA'
        else:
            return '{:.6f}'.format(float(table['C']+table['G'])/float(total))
    def __get_gc_count(self, genome):
        '''
        get the gc count within DHS
        --genome    pysam or pyfasta object
        '''
        table = {'A':0,'C':0,'G':0,'T':0}
        for bed in self.center:
            seq = get_seq_from_bed(bed.chrom, bed.start, bed.stop, genome)
            for nt in ['A','C','G','T']:
                table[nt] += seq.count(nt)
        return str(table['C']+table['G'])
    def __get_center_size(self, anno):
        '''
        get region size
        --anno  will not be used at all
        '''
        size = 0
        for bed in self.center:
            size += (bed.stop-bed.start)
        return str(size)
    def __get_exclude_label(self, anno):
        '''
        get whether the region is included in filtered region
        --anno  a dict contains the position (0-base, with chrom) information, if value of a key is 0 this key(pos) will not be considered
        '''
        for bed in self.center:
            if bed.stop - bed.start == 1:
                if bed.start in anno[bed.chrom]:
                    return str(anno[bed.chrom][bed.start])
                else:
                    return '0'
            else:
                for pos in range(bed.start,bed.strop):
                    if pos in anno[bed.chrom]:
                        if anno[bed.chrom][pos] == 0:
                            return '0'
                    else:
                        return '0'
        return '1'
    def __get_tsi_peak(self, tabix_list):
        '''
        get tissue specificity of region
        --tabix_list    a list of tabix object
        '''
        array = []
        for peak_tabix in tabix_list:
            if int(self.__get_center_count(peak_tabix)) > 0 + 1e-6:
                array.append(1.0)
            else:
                array.append(0.0)
        return self.cal_entropy(array)
    def __get_tsi_signal(self, bw_list):
        '''
        get tissue specificity of region
        --bw_list       a list of bigwig object
        '''
        array = []
        for bw in bw_list:
            value_ave = self.__get_mean_signal(bw)
            if value_ave == 'NA':
                continue
            else:
                array.append(float(value_ave))
        return self.cal_entropy(array)
    def __get_tfbs_anno(self, tabix_dbi, tfbs_list):
        try:
            result = tabix_dbi.query(self.chrom, self.start, self.stop)
            for overlap in result:
                name = overlap[3]
                if name in tfbs_list:
                    return '1'
        except:
            return '0'
        return '0'
    def add_value(self, label, value):
        self.__add_label(label)
        self.anno[label] = value
    def overlap(self, tabix_dbi, marker):
        '''if region overlap with region with name (4th column) equals to marker return True else return False'''
        try:
            result = tabix_dbi.query(self.chrom, self.start, self.stop)
            for overlap in result:
                name = overlap[3]
                if name in marker:
                    return True
        except:
            return False
        return False
    def write(self, label_list = None):
        assert len(self.anno) == len(self.label)
        if label_list is None: # use label list order by default
            return "%s\t%d\t%d\t%d" % (self.chrom, self.start, self.stop, (self.start+self.stop)/2) + '\t' + '\t'.join(self.anno[label.split(':')[0]] for label in self.label)
        else: # support re-order
            return "%s\t%d\t%d\t%d" % (self.chrom, self.start, self.stop, (self.start+self.stop)/2) + '\t' + '\t'.join(self.anno[label.split(':')[0]] for label in label_list)
    def header(self,label_list = None):
        if label_list is None:
            return "#chrom\tstart\tstop\tmid\t" + '\t'.join(label for label in self.label)
        else:
            return "#chrom\tstart\tstop\tmid\t" + '\t'.join(label for label in label_list)
    def write_to_pickle(self, fout_file):
        return
    def load_from_pickle(self, fin_file):
        return

class SamHandler(object):
    '''inheritate from pysam object:'''
    def __init__(self,samfile):
        # public
        self.samfile = samfile
        # private
        self.genome_table = {}
        self.tid2name = {}
        # init function
        self.get_genome()
    def get_genome(self):
        assert len(self.samfile.references) == len(self.samfile.lengths)
        for nn in range(len(self.samfile.references)):
            chrom = self.samfile.references[nn]
            tid = self.samfile.gettid(chrom)
            size = self.samfile.lengths[nn]
            self.genome_table[chrom] = int(size)
            self.tid2name[tid] = chrom
    def get_count(self, chrom, start, stop):
        '''get the cleavage count (can also be viewed as count of reads) of a region'''
        if not isinstance(chrom,str) or not isinstance(start, int) or not isinstance(stop,int):
            error("chromosomen name (%s) must be string, start (%s) and stop (%s) position must be integer" % (chrom,start,stop))
            exit(1)
        cut = {'+':[0]*(stop-start), '-':[0]*(stop-start)}
        for alignread in self.samfile.fetch(str(chrom),start,stop):
            #tid_id = alignread.rname
            # get cut site
            is_rev = False
            if alignread.is_reverse:
                is_rev = True
                cut_site = int(alignread.aend)
            else:
                cut_site = int(alignread.pos) - 1
            # take care of start or end of chromosome
            cut_site = max(0, cut_site)
            cut_site = min(self.genome_table[chrom]-1, cut_site)
            # skip read outside given region
            if cut_site < start or cut_site > stop -1 :
                continue
            # update cut_table
            if is_rev:
                cut['-'][cut_site-start] += 1
            else:
                cut['+'][cut_site-start] += 1
        return cut['+'], cut['-']

class Gene(object):
    '''Genome annotation object'''
    def __init__(self,chrom,start,stop,strand,exon_list,cds_list,genename,geneid):
        # public
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand
        self.exon = exon_list
        self.cds = cds_list
        self.name = genename
        self.geneid = geneid
        # private
        self.promoter = {}
        self.downstream = {}
        self.tss = None
        self.tes = None
        self.utr5 = None
        self.utr3 = None
    def __repr__(self):
        return "%s\t%d\t%d\t%s\t%s" % (self.chrom,self.start,self.stop,self.strand,self.name) + '\t' + ','.join(str(item[0])+'-'+str(item[1]) for item in self.exon)
    def getpromoter(self,up,down,genome_size):
        key = "%d_%d" % (up,down)
        if self.strand == '+':
            self.promoter[key] = (max(self.start-up,0), min(self.start-down,genome_size[self.chrom]))
        elif self.strand == '-':
            self.promoter[key] = (max(self.stop+down,0), min(self.stop+up,genome_size[self.chrom]))
        else:
            error("Unknown strand: %s" % (self.strand))
    def getintron(self):
        self.intron = []
        for nn in range(len(self.exon)-1):
            intron_start = self.exon[nn][1]
            intron_stop = self.exon[nn+1][0]
            assert intron_start < intron_stop
            self.intron.append((intron_start,intron_stop))
    def gettss(self):
        if self.strand == '+':
            self.tss = self.start
        else:
            self.tss = self.stop - 1
    def gettes(self):
        if self.strand == '+':
            self.tes = self.stop - 1
        else:
            self.tes = self.start
    def getdownstream(self,size,genome_size):
        key = str(size)
        if self.strand == '+':
            self.downstream[key] = (self.stop, min(self.stop+size,genome_size[self.chrom]))
        elif self.strand == '-':
            self.downstream[key] = (max(0,self.start-size),self.start)
        else:
            error("Unknown strand: %s" % (self.strand))
            exit(1)
    def getutr(self):
        if len(self.cds) < 1:
            return
        assert self.exon[0][0] <= self.cds[0][0]
        assert self.exon[-1][1] >= self.cds[-1][1]
        if self.exon[0][0] < self.cds[0][0]:
            self.utr5 = (self.exon[0][0],self.cds[0][0])
        else:
            pass
        if self.exon[-1][1] > self.cds[-1][1]:
            self.utr3 = (self.cds[-1][1],self.exon[-1][1])
        else:
            pass
        if self.strand == '-': # if gene in reverse strand swap two utr
            tmp = self.utr3
            self.utr3 = self.utr5
            self.utr5 = tmp
        elif self.strand == '+':
            pass
        else:
            error("Unknown strand: %s" % (self.strand))
            exit(1)
    def build(self,promoter_set,downstream_set,genome_size): # get private variable value
        self.getintron()
        self.gettss()
        self.gettes()
        self.getutr()
        if promoter_set is not None:
            up,down = promoter_set
            self.getpromoter(up,down,genome_size)
        if downstream_set is not None:
            self.getdownstream(downstream_set,genome_size)
 
#################
# Custom Function
#################

def logging(text):
    print >>sys.stderr, "Logging: " + text

def warning(text):
    print >>sys.stderr, "Warning: " + text

def error(text,code=None):
    if code is not None:
        print >>sys.stderr, "Error Num:%d %s" % (code,text)
    else:
        print >>sys.stderr, "Error: %s" % (text)

def check_folder_exist(folder,is_cover=True):
    """Check the existence of folder, it not exist will create it"""
    if not os.path.exists(folder):
        os.mkdir(folder)
    else:
        if is_cover:
            warning("output folder %s already exists, may overwrite the content inside it." % (folder))
        else:
            error("output folder %s already exists, will exit." % (folder))
            exit(1)

def check_file_exist(filename,is_cover=True):
    """Check the existence of folder"""
    if not os.path.isfile(filename):
        pass
    else:
        if is_cover:
            warning("output file %s already exsits, will overwrite the content inside it." % (filename))
        else:
            error("output file %s already exists, will exit." % (filename))
            exit(1)

def read_from_file(fin_name):
    try:
        fin = open(fin_name,"r")
    except IOError:
        error("Can't open file:%s to read." % (fin_name))
        exit(1)
    return fin

def write_to_file(fout_name,is_cover=True):
    if fout_name == 'stdout' or not fout_name:
        fout = sys.stdout
        warning("use stdout to report result")
    else:
        if check_file_exist(fout_name):
            if is_cover:
                warning("file:%s already exist. Old file will be covered by new file." % (fout_name))
                try:
                    fout = open(fout_name,"w")
                except IOError:
                    warning("Can't open file:%s to write. Use stdout instead." % (fout_name))
                    fout = sys.stdout
            else:
                error("file:%s already exist" % (fout_name))
                exit(1)
        else:
            try:
                fout = open(fout_name,"w")
            except IOError:
                warning("Can't open file:%s to write. Use stdout instead." % (fout_name))
                fout = sys.stdout
    return fout
 
def load_region(filename):
    '''load bed file into a simple list'''
    region_list = []
    fin = read_from_file(filename)
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])
        stop = int(row[2])
        region_list.append(Region(chrom, start, stop))
    fin.close()
    return region_list

def load_region_anno(filename):
    '''load bed region annotation file, the first row is the header containing the annotation label for each anno data'''
    region_list = []
    fin = read_from_file(filename)
    for line in fin:
        if line.strip() == '':
            continue
        row = line.strip().split()
        if line.strip().startswith('#'): # header
            if len(row) > 3: # with additional column
                label_list = row[3:]
            continue
        chrom = row[0]
        start = int(row[1])
        stop = int(row[2])
        region = Region(chrom,start,stop)
        region.label = label_list
        assert len(row) - 3 == len(label_list)
        for nn in range(len(label_list)):
            region.anno[label_list[nn]] = row[3+nn]
        region_list.append(region)
    fin.close()
    return region_list

def load_gene(bed_filename):
    '''load 12 column gene annotation downloaded form UCSC'''
    fin = read_from_file(bed_filename)
    genelist = []
    for line in fin: 
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])
        stop = int(row[2])
        name = row[3]
        strand = row[5]
        cds_start = int(row[6])
        cds_stop = int(row[7])
        exon_len = [int(item) for item in row[10].rstrip(',').split(',')]
        exon_start = [int(item) for item in row[11].rstrip(',').split(',')]
        try:
            assert len(exon_len) == len(exon_start)
        except:
            warnning("exon count and exon start is not same for gene:%s"  %(name))
            continue
        exon_list = []
        cds_list = []
        for nn in range(len(exon_start)):
            exon_start_pos = start + exon_start[nn]
            exon_stop_pos = exon_start_pos + exon_len[nn]
            exon_list.append((exon_start_pos,exon_stop_pos))
            if cds_start <= exon_stop_pos and cds_stop >= exon_start_pos:
                cds_start_pos = max(cds_start,exon_start_pos)
                cds_stop_pos = min(cds_stop,exon_stop_pos)
                if cds_start_pos < cds_stop_pos:
                    cds_list.append((cds_start_pos,cds_stop_pos))
        genelist.append(Gene(chrom,start,stop,strand,exon_list,cds_list,name,'NA'))
    genelist_sorted = sorted(genelist, key = attrgetter('chrom','start'))
    fin.close()
    return genelist_sorted

def merge_segment(filename):
    '''merge the segment in a 4 column bed file, assuming the bed has been sorted'''
    fin = read_from_file(filename)
    tmpout = tempfile.NamedTemporaryFile(delete=False)
    tmp_chrom = None
    tmp_start = None
    tmp_stop = None
    tmp_group = None
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])
        stop = int(row[2])
        group = row[3]
        if tmp_start is None: # first line
            tmp_chrom = chrom
            tmp_start = start
            tmp_stop = stop
            tmp_group = group
        elif tmp_chrom != chrom: # new chromsome
            print >>tmpout, "%s\t%d\t%d\t%s" % (tmp_chrom,tmp_start,tmp_stop,tmp_group)
            tmp_chrom = chrom
            tmp_start = start
            tmp_stop = stop
            tmp_group = group
        else: # same chromsome
            try:
                assert tmp_stop == start # should be adjacent with each other
            except AssertionError:
                print >>sys.stderr, chrom, tmp_stop, start
                print >>sys.stderr, "Nearby regions are not adjacent. Is the bed file sorted."
                exit(1)
            if tmp_group == group: # same group name extend previous stored segment
                tmp_stop = stop 
            else: # different segment, update it
                print >>tmpout, "%s\t%d\t%d\t%s" % (tmp_chrom,tmp_start,tmp_stop,tmp_group)
                tmp_start = start
                tmp_stop = stop
                tmp_group = group
    fin.close()
    print >>tmpout, "%s\t%d\t%d\t%s" % (tmp_chrom,tmp_start,tmp_stop,tmp_group)
    tmpout.flush()
    tmpout.close()
    os.system("cp %s %s" % (tmpout.name,filename))
    os.unlink(tmpout.name)

def load_genome_size(filename):
    table = {}
    fin = read_from_file(filename)
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        table[row[0]] = int(row[1])
    return table

def get_seq_from_bed(chrom, start, stop, genome):
    '''
    get seq from a genomic region
    --chrom     chromosome name
    --start     start position (0-base)
    --stop      stop position (0-base)
    --genome    pysam.Fastafile or pyfasta.Fasta object
    '''
    if isinstance(genome, pysam.Fastafile):
        seq = str(genome.fetch(str(chrom), start, stop)).upper()
    elif isinstance(genome, pyfasta.Fasta):
        seq = str(genome[chrom][start:stop]).upper()
    else:
        print >>sys.stderr, "Unknown genome fasta class"
        exit(1)
    return seq

def subtract_bed(chrom, start, stop, exclude_tabix_list):
    '''return bed region list that not within exclude tabix file'''
    if len(exclude_tabix_list) < 1 or exclude_tabix_list is None or (len(exclude_tabix_list) ==1 and exclude_tabix_list[0] is None):
        return [BED(chrom,start,stop)]
    region_list = []
    remove = {}
    pos = [start, stop]
    # mark region should be removed
    for exclude_tabix in exclude_tabix_list:
        result = exclude_tabix.query(chrom, start, stop)
        for overlap in result:
            o_chrom = overlap[0]
            o_start = int(overlap[1])
            o_stop = int(overlap[2])
            exclude_start = max(o_start, start)
            exclude_stop = min(o_stop, stop)
            remove[(exclude_start, exclude_stop)] = True
            pos.append(exclude_start)
            pos.append(exclude_stop)
        continue
        try:
            result = exclude_tabix.query(chrom, start, stop)
            for overlap in result:
                o_chrom = overlap[0]
                o_start = int(overlap[1])
                o_stop = int(overlap[2])
                exclude_start = max(o_start, start)
                exclude_stop = min(o_stop, stop)
                remove[(exclude_start, exclude_stop)] = True
                pos.append(exclude_start)
                pos.append(exclude_stop)
        except:
            pass
    # boundary list
    pos = sorted(list(set(pos)))
    # report
    for nn in range(len(pos)-1):
        start = pos[nn]
        stop = pos[nn+1]
        if (start, stop) in remove:
            continue
        else:
            is_remove = False
            for key in remove.keys():
                remove_start = key[0]
                remove_stop = key[1]
                if start >= remove_start and stop <= remove_stop:
                    is_remove = True
                    break
            if is_remove:
                continue
            else:
                region_list.append(BED(chrom, start, stop))
    return region_list

def sort_bed(filename,is_cover):
    """sort bed file by chrom then by start"""
    tmpout = tempfile.NamedTemporaryFile(delete=False)
    genome = {}
    fin = read_from_file(filename)
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        genome[chrom] = 0
    fin.close()
    sorted_chrom = sort_genome(genome)
    table = {}
    fin = read_from_file(filename)
    for line in fin:
        if line.strip().startswith('#') or line.strip() == '':
            continue
        row = line.strip().split()
        chrom = row[0]
        start = int(row[1])
        stop = int(row[2])
        other = "\t".join(item for item in row[3:])
        try:
            table[chrom].append([chrom,start,stop,other])
        except KeyError:
            table[chrom] = [[chrom,start,stop,other]]
    sorted_table = {}
    for chrom in table.keys():
        sorted_table[chrom] = sorted(table[chrom],key=lambda bed:bed[1])
    for chrom in sorted_chrom:
        for bed in sorted_table[chrom]:
            print >>tmpout, bed[0] + '\t' + str(bed[1]) + '\t' + str(bed[2]) + '\t' + bed[3]
    tmpout.flush()
    fin.close()
    if is_cover:
        os.system("cp %s %s" % (tmpout.name,filename))
        os.unlink(tmpout.name)
        return filename
    else:
        return tmpout.name

def create_tabix(filename,zippedname,mode):
    os.system("bgzip -c %s >%s" % (filename,zippedname))
    os.system("tabix -p %s %s" % (mode,zippedname))

def test():
    #TODO
    return
    
if __name__=="__main__":
    test()
