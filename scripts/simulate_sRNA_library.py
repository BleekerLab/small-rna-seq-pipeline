#!/usr/bin/env python

from string import digits
import sys
import subprocess
import os
from random import randint
from random import shuffle
from math import log10
from time import time
import datetime
################################################################################
#FUNCTIONs

def percent_complete(total_reads, read_count):
    #Function that produces on-screen percentage marker for tracking progress
    #Disabled if --suppress_percent present in options

    return
    if opt_suppress == False:

    #Reads through the counts as a percent of total
    #Reports with a [=============]34% type status bar

        if total_reads % 100 == read_count:
            p = read_count / ( total_reads / 20 )
            q = int(float(read_count) / float(total_reads) * 100)
            sys.stdout.write('\r')
            # the exact output you're looking for:
            sys.stdout.write("[%-20s] %d%% " % ('='*p, q))
            sys.stdout.flush()
        elif read_count == 1:
            p = 0
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%% " % ('='*p, 0))
            sys.stdout.flush()
        else:
            pass
    else:
        pass

def validate_path(path):
    #this function is called to validate the existance of an os path

    if os.path.exists(path) == False:
        print "Fatal:  Invalid path specified.\n%s\n" % (path)
        sys.exit()

def validate_genome():
    #Function to validate the genome integrity and components

    #check if the file exists
    if os.path.exists(opt_g) == True:
        print opt_g, "--> Found"
    else:
        print opt_g, "--> Fatal: genome not found"
        sys.exit()

    #This is an example of a bowtie index for a given genome
    #Expected to be: genome_location.1.ebwt
    bowtie_index = '%s.1.ebwt' % (opt_g)

    #checks for the ebwt existance
    if os.path.exists(bowtie_index) == True:
        print ".ebwt index Found"
    else:
        print "FATAL: (bowtie-build) .ebwt index not found!"
        print "confirm the existence of .ebwt index files"
        print "(searched for: %s)" % (bowtie_index)
        print ""
        print "bowtie index must be constructed using bowtie-build in this fashion:"
        print "bowtie-build [./genome.fa] [genome.fa]"
        sys.exit()

    #checks for the existance of the .fai samtools index
    s_fai = opt_g + ".fai"
    if os.path.exists(s_fai) == True:
        print ".fai index found"
    else:
        print "Failed to open expected fai index file %s.\nAttempting to create using samtools faidx ..." % (s_fai)

    #if this fails, samtools will attempt to build it's own
        subprocess.call(["samtools","faidx", opt_g])

        if os.path.exists(s_fai) == True:
            print "successful!"
        else:
            print "FAILED"
            sys.exit()

def miRNA_bins():
    # builds a list of miRNA bins from an annotation file

    print "\nAssigning miRNA bins from: " + opt_mir
    bins = {}
    mir_loc_list = []

    #opening the annotation, reading by line
    with open(opt_mir, 'r') as f:
        for line in f:
            if line[0] != '#':
                if line[0].translate(None, digits) == '':
                    if 'miRNA_primary_transcript' in line:
                        line = line.split('\t')
                        #storing the start and stop locations for every loci in a list

                        start = "chr%s:%s" % (line[0], line[3])
                        stop  = "chr%s:%s" % (line[0], line[4])
                        mir_loc_list.append(start)
                        mir_loc_list.append(stop)



                else:
                    if 'miRNA_primary_transcript' in line:
                        line = line.split('\t')
                        #storing the start and stop locations for every loci in a list

                        start = "%s:%s" % (line[0], line[3])
                        stop  = "%s:%s" % (line[0], line[4])

                        if chrom_id_dict['genome_type'] == 'something_upper':

                            if start.lower() == start:

                                start = "Chr%s" % (start.strip('chr'))
                                stop = "Chr%s" % (stop.strip('chr'))

                        mir_loc_list.append(start)
                        mir_loc_list.append(stop)


    #with the list of all mir chromosomal locations
    for i in mir_loc_list:

        #converts to a genomic bin location -
        #genomic coordinates, sorting each site into two bins which overlap by half
        try:
            key = int(locus_index(i).gen_coord() / bin_size)
        except:
            print 'FATAL: Chromosome name from .gff3 not found in genome.  Are these formatted the same?'
            print i, '<- format in mir file'
            print chrom_id_dict
            sys.exit()

        #saving a dictionary of these bins, with the key representing the bin_size sized bin
        #bins[3] represents miRNAs
        bins[key] = 0,0,0,1


    print "Done!"
    return(bins)

def bowtie():
    #runs bowtie on the template sRNA library, storing as a .bam

    proc_start = time()

    #opening save-file
    f = file(opt_bam, 'w')

    # print bow_p, 'bowp'
    # print opt_g, 'optg'
    # print opt_f, 'optf'
    #process for bowtie alignment
    p1_call = "bowtie -f --all -m 1000 --best --strata -v 0 -p %s --sam %s %s" % (bow_p, opt_g, opt_f)
    print '\nCall for bowtie: \n%s' % (p1_call)
    p1 = subprocess.Popen([p1_call], shell = True, stdout=subprocess.PIPE, bufsize = 1)
    #process for saving as a .bam
    p2 = subprocess.Popen(["samtools", "view", "-b", '-'], stdin = p1.stdout, stdout = f, bufsize = 1)

    #waits till the end then kills the processes
    p2.wait()
    f.close()

    proc_end = time()

    print "Bowtie + samtools --> Time elapsed: %s minutes" % (round((proc_end-proc_start) / 60, 2))

def other_bins(bins):
    #assigning bins other than miRNA

    print "\nAssigning Phased, siRNA and garbage RNA"
    # print 'Totaling length'
    # proc_start = time()

    # this quantify the bin amounts for phased, si and garbage RNA
    # acts very similarly to miRNA_bins()

    #this variable represents a black-hole for unwanted output
    devnull = open(os.devnull, 'wb')

    #process counting the lines of the incoming bam, for percentage use
    # p = subprocess.Popen("samtools view %s | wc -l" % (opt_bam), stdout = subprocess.PIPE, stderr = devnull, shell = True)
    # total = int(p.communicate()[0].strip())
    # proc_end = time()
    # print 'total lines:', total
    #
    # print "Totaling --> Time elapsed: %s minutes" % (round((proc_end-proc_start) / 60, 2))


    print "\nReading and Placing..."

    assign_count = 1
    proc_start = time()

    #process to pipe .bam file
    p = subprocess.Popen(["samtools view %s" % (opt_bam)], stdout = subprocess.PIPE, stderr = devnull, shell = True, bufsize = -1)

    #loop reads through all of the lines, 1 at a time
    while True:

        try:
            line = p.stdout.readline()
            #print line

            #these variables and function calls calculate the on-screen percentage
            assign_count += 1
            if assign_count%1000000 == 0:
                if assign_count < 10000000:
                    print assign_count
                else:
                    if assign_count%10000000 == 0:
                        print assign_count

            # percent_complete(total, assign_count)

            #getting the chromosomal position of a read in the sample sam
            chr_pos = "%s:%s" % (line.split('\t')[2], line.split('\t')[3])
            seq = line.split('\t')[9]


            #key represents a bin number, calculated by rounding from the genomic coord. divided by 150 nt bin
            key = int(locus_index(chr_pos).gen_coord() / bin_size)

            #if this key doesn't exist yet, it is instantiated
            try:
                bins[key]
            except:
                bins[key] = 0,0,0,0


            #reads length 23 or 24 nt are assigned as candidate siRNAs
            if len(seq) == 23 or len(seq) == 24:

                bins[key] = bins[key][0] + 1 ,bins[key][1],bins[key][2], bins[key][3]

            #reads length 21 are assigned as tasiRNA, as miRNA containing bins have already been quantified
            if len(seq) == 21:
                bins[key] = bins[key][0], bins[key][1] + 1,bins[key][2], bins[key][3]

            #reads between 10 and 21 nt or more than 24 nt are deemed as "garbage"
            if len(seq) < 21 and len(seq) > 10 or len(seq) > 24:
                bins[key] = bins[key][0], bins[key][1], bins[key][2] + 1, bins[key][3]

        #this loop is ended by a key failing to exist, a sign that there are no more reads to count
        except IndexError:
            p.kill()
            break

######################

    proc_end = time()

    print "Bin formation --> Time elapsed: %s minutes" % (round((proc_end-proc_start) / 60, 2))

    print ""

    return(bins)

def bin_quantify(align_bins):
    #function to quantify the number of bins for each sRNA type, using differing bin-depths

    print 'Quantifying bins'
    si_bin_info = {}
    pha_bin_info = {}
    garb_bin_info = {}
    mi_bin_info = {}
    list_of_dicts = [si_bin_info, pha_bin_info, garb_bin_info, mi_bin_info]
    list_of_keys = ['alignment_count', 'bin_count', 'gr_eq_5', 'gr_eq_4', 'gr_eq_3', 'gr_eq_2', 'gr_eq_1']

    #instantiates keys for each type of sRNA
    for dicts in list_of_dicts:
        for keys in list_of_keys:
            dicts[keys] = 0

    #quantifies the number of reads in a bin, using cutoffs of 5-1 for depth
    #stores in the *_bin_info dictionaries
    #calls bin_info to solve this
    for k in align_bins:
        bin = align_bins[k]

        bin_info(bin, 0, si_bin_info)
        bin_info(bin, 1, pha_bin_info)
        bin_info(bin, 2, garb_bin_info)
        bin_info(bin, 3, mi_bin_info)

    #identifies minimum depth required for the simulation to complete
    #this is only important with si and tasiRNA
    #calls depth_calculator to solve this
    list_of_keys.append('min_bin_depth')
    si_bin_info['min_bin_depth'] = depth_calculator(si_bin_info, len(het_ns))
    pha_bin_info['min_bin_depth'] = depth_calculator(pha_bin_info, len(tasi_ns))
    garb_bin_info['min_bin_depth'] = 'Na'
    mi_bin_info['min_bin_depth'] = 'Na'

    #builds a dictionary of all of the gathered data
    #values are stored in a tuple 0-si, 1-pha, 2-garb, 3-mir
    out = {}
    for keys in list_of_keys:
        out[keys] = (list_of_dicts[0][keys], list_of_dicts[1][keys], list_of_dicts[2][keys], list_of_dicts[3][keys])

    return(out)

def bin_info(bin, index, bin_info):
    #function to calculate the number of reads per bin, given 5-1 depth cut-offs
    #index indicates which type of sRNA is of interest

    #totals non-mir depths
    total_in_bin = bin[0] + bin[1] + bin[2]
    value = bin[index]

    #checks to see if this sRNA is the 80% majority
    if value > total_in_bin * 0.8:
        TorF = True
    else:
        TorF = False

    #if it is, then tallies amounts for different depth cut-offs
    #if miRNAs are present in the bin, this bin is discluded from the others
    if value != 0 and TorF == True:
        bin_info['bin_count'] += 1
        bin_info['alignment_count'] += value
    if value >= 5 and TorF == True and bin[3] == 0:
        bin_info['gr_eq_5'] += 1
    if value >= 4 and TorF == True and bin[3] == 0:
        bin_info['gr_eq_4'] += 1
    if value >= 3 and TorF == True and bin[3] == 0:
        bin_info['gr_eq_3'] += 1
    if value >= 2 and TorF == True and bin[3] == 0:
        bin_info['gr_eq_2'] += 1
    if value >= 1 and TorF == True and bin[3] == 0:
        bin_info['gr_eq_1'] += 1

def depth_calculator(info, req_loci):
    #checks to see what depth is required to satisfy the simulated library size
    #1.2 x the exact amount needed is the cutoff
    #if the bin quantity is not sufficient with depth 1, than a pseudo-annotation CANNOT be constructed

    if info['gr_eq_5'] > int(req_loci * 1.2):
        return(5)
    elif info['gr_eq_4'] > int(req_loci * 1.2):
        return(4)
    elif info['gr_eq_3'] > int(req_loci * 1.2):
        return(3)
    elif info['gr_eq_2'] > int(req_loci * 1.2):
        return(2)
    elif info['gr_eq_1'] > int(req_loci * 1.2):
        return(1)
    else:
        print '\nError:  Sample library depth is insufficient to generate simulated library'
        sys.exit()

def target_list():
    #function to assign bins to each sRNA type

    totals = align_bins['totals']
    align_bins.pop('totals', None)

    total_bins = align_bins['num_of_bins']
    align_bins.pop('num_of_bins', None)

    print ""

    si_list = []
    pha_list = []
    garb_list = []
    mi_list = []

    #sorts through bins, assigning each bin to only 1 type of sRNA
    #bins containing an miRNA read are automatically chosen for mir
    #all other sRNA types need an 80% majority

    for k in align_bins:
        #total of non miRNA depths
        total_in_bin = align_bins[k][0] + align_bins[k][1] + align_bins[k][2]
        #siRNA
        if align_bins[k][0] >= min_bin_depth[0] and align_bins[k][0] > total_in_bin * 0.8 and align_bins[k][3] == 0:
            si_list.append(int(k*bin_size))
        #tasiRNA
        elif align_bins[k][1] >= min_bin_depth[1] and align_bins[k][1] > total_in_bin * 0.8 and align_bins[k][3] == 0:
            pha_list.append(int(k*bin_size))
        #garbage
        elif align_bins[k][2] >= 5 and align_bins[k][3] == 0:
            garb_list.append(int(k*bin_size))
        #miRNA
        elif align_bins[k][3] != 0:
            mi_list.append(int(k*bin_size))

    #quantifying to find out which sRNA type has the most bins
    si =  len(si_list)
    pha = len(pha_list)
    garb = len(garb_list)
    mi = len(mi_list)
    most = max(len(si_list), len(pha_list), len(garb_list), len(mi_list))

    #writes 0's to the shorter lists to make them all the same length
    for r in range(0,most):
        try:
            si_list[r]
        except:
            si_list.append(0)
        try:
            pha_list[r]
        except:
            pha_list.append(0)
        try:
            garb_list[r]
        except:
            garb_list.append(0)
        try:
            mi_list[r]
        except:
            mi_list.append(0)


    #combines the lists into a single list of a 4-var tuple (si, tasi, garb, mi)
    z = zip(si_list, pha_list, garb_list, mi_list)

    #overwriting any files with the location [opt_p]
    with open(opt_p, "w") as f:
        f.write("")

    #writes the pseudo-annotation, with an informative header
    with open(opt_p, "w+") as f:
        for k in z:
            try:
                pseudo
                f.write(pseudo)
            except:
                try:
                    loc = opt_bam
                except:
                    loc = opt_f

                print "Writing Pseudo annotation:", opt_p
                info_header = '''#PA location: %s
#Source lib: %s
#Source MIR: %s
#Minimum Bin depth: (1-5 alignments, with lower values indicating lower confidence)
#\t%s siRNA
#\t%s tasiRNA
#Number of candidate bins:
#\t%s siRNA
#\t%s tasiRNA
#\t%s non-small RNA
#\t%s miRNA
#siRNA\tphaRNA\tnon-small\tmiRNA\n''' % (opt_g, loc, opt_mir, min_bin_depth[0], min_bin_depth[1], si, pha, garb, mi)

                f.write(info_header)

            pseudo = "%s\t%s\t%s\t%s\n" % (k[0],k[1],k[2],k[3])




def chrom_id(opt_g):
    # reads from the genome's .fai file, and spits out a dictionary where:
    # "chromosome name" is the key
    # Chr_start, Chr_end is the value, in respect to their whole genome locations
    with open("%s.fai" % (opt_g), "r") as g:
        chrom_id_dict = {}
        for line in g:
            chrom_key = str(line.split("\t")[0])#.lower()
            chrom_id_dict[chrom_key] = int(line.split("\t")[2]), int(line.split("\t")[2]) + int(line.split("\t")[1])

    GT = 'string'
    for key in chrom_id_dict:
        if key.translate(None, digits) == '':
            GT = 'numeric'

    if GT == 'string':
        for key in chrom_id_dict:
            if key.lower() == key:
                GT = 'all_lower'
            else:
                GT = 'something_upper'


    chrom_id_dict['genome_type'] = GT

    return(chrom_id_dict)

class locus_index:
    #this class is assuming that opt_g.fai and chrom_id_dict are defined objects

    def __init__(self, locus_input):
        #this is if the locus_input came from a sam file
        if "\t" in locus_input:
            self.chr  = locus_input.split("\t")[2]#.lower()
            self.chr_pos = locus_input.split("\t")[3]

        #this is if the locus_input came from a chromosomal loc
        # format: chr5:4557424
        if ":" in locus_input and "-" not in locus_input:
            self.chr  = locus_input.split(":")[0]#.lower()
            self.chr_pos = locus_input.split(":")[1]

        #this is if the locus_input came from a chromosomal interval
        # format: chr5:4557424-4557430
        if ":" in locus_input and "-" in locus_input:
            self.chr  = locus_input.split(":")[0]#.lower()
            self.chr_pos = locus_input.split(":")[1].split("-")[0]

        #this is if the locus_input came from a genomic loc, as a tuple
        if type(locus_input) == int:
            self.gen_pos = int(locus_input)

        if type(locus_input) == tuple:
            self.gen_pos = int(locus_input[0])
            self.gen_end   = int(locus_input[1])

    #reads coordinates in chromosomal form, and returns them as a single variable genomic coordinate
    def gen_coord(self):
        if self.chr == "*":
            coord = 0
        else:
            coord = int(chrom_id_dict[self.chr][0]) + int(self.chr_pos)
        return (int(coord))

    #reads coordinates in genomic form, and returns as a chromosomal choordinate
    def chr_coord(self):
        for w in chrom_id_dict:
            if self.gen_pos > chrom_id_dict[w][0] and self.gen_pos < chrom_id_dict[w][1]:
                try:
                    coord = "%s:%s-%s" % (w, self.gen_pos - chrom_id_dict[w][0], self.gen_end - chrom_id_dict[w][0])
                except:
                    coord = "%s:%s-%s" % (w, self.gen_pos - chrom_id_dict[w][0], self.gen_pos + interval - chrom_id_dict[w][0])
                return (coord)



def get_bins(RNA_type):
    #reads the pseudo-annotation, returning a list of genomic positions for a single type of sRNA
    #these positions indicate the first nt of a bin

    genome_positions = []

    if "miRNA" in RNA_type:
        index = 3
    if "siRNA" in RNA_type:
        index = 0
    if "tasiRNA" in RNA_type:
        index = 1
    if "garbage" in RNA_type:
        index = 2
    with open(opt_p, "r") as f:
        for line in f:
            if line[0] != '#':
                bin_loc = line.split('\t')[index]
                if int(bin_loc) != 0:
                    genome_positions.append(bin_loc)

    return(genome_positions)


def get_a_locus(genome_positions, interval):
    #randomly chooses bin from possible bins, returning the (start,end) of the locus
    #as well as the specific number for the bin

    rand = randint(2,len(genome_positions)-1)
    chosen_bin = genome_positions[rand]
    locus_start = int(chosen_bin) + randint(1,bin_size-1)
    locus_end = locus_start + interval
    locus = (locus_start, locus_end)
    return((locus, chosen_bin))

def occupied_checker(occupied_list, locus, locus_type):
    #checks for loci overlap, returns 1 if occupied, the loci if not

    locus_start = locus[0]
    locus_end   = locus[1]

    #looking for this in a list of occupied coordinates
    for k in occupied_list:

        #if the locus coordinates fall within an occupied space,
        if (k[0] < locus_start and k[1] > locus_start) or (k[0] < locus_end and k[1] > locus_end - 50):
            return(1)

        #otherwise returns the loci and type
        else:
            return((locus_start, locus_end, locus_type))

def get_a_seq(gen_loc):
    #returns the seq for a total loci region

    chr_location = locus_index(gen_loc).chr_coord()
    # try:
    #     #if 'chr' in chr_location:
    #     #    chr_location = 'Chr'.join(chr_location.split('chr'))
	#     pass
    #
    #
    # except TypeError:
    #     print ''
    #     print chr_location, 'chr_location'
    #     print gen_loc, 'gen_loc'
    #     return(None)
    #     # sys.exit()
    #     # pass

    #process calling the samtools faidx tool, which will give the precise sequence for this range


    if chrom_id_dict['genome_type'] == 'numeric':
        chr_location = chr_location.strip('Chr')


    # print chrom_id_dict['genome_type']
    # print chr_location

    #print chr_location

    try:
        p1 = subprocess.Popen(["samtools","faidx", opt_g, chr_location], stdout = subprocess.PIPE)
        seq = p1.communicate()[0]
        seq = "".join(seq.split("\n")[1:])
    except TypeError:
        return("")

    return(seq)

def filter(seq):
    #filtering out easily identifiable highly-reptetive regions

    seq = seq.upper()


    # Unknown reads
    # ex. ACTGCTAGTCGANNNTATANNCGACG
    if "N" in seq:
        return(True)


    # DiNucleotide Repeats
    # ex.  ATATATATATATATATATATATATAT
    evens = list(seq)[::2]
    odds  = list(seq)[1::2]
    all = evens + ["null"] + odds
    num = 0
    longest = 0
    old = ''
    for i in all:
        if i == old:
            num += 1
            if num > longest:
                longest = num
        if i != old:
            num = 0
        if num == "null":
            old = ''
            num = 0
        old = i

    #if the dinucleotide repeat is 16 bp in length
    if longest > 8:
        return(True)



    #returning false if no bad sequences found
    return(False)

def final_read_check(read):
    #true means it's good, false means bad!

    if read < 20:
        return(False)
    if read > 24:
        return(False)

    return(True)

def get_nts(reads2get, entries2get):
    #returns a list of loci, with the value being their read depth

    logmax = log10(reads2get)
    loginterval = logmax / entries2get
    reads_done = 0
    this_log = 0
    out = []

    #building this list as 'out', with the depths defined by a logarithmic function
    while reads_done < reads2get:
        reads = 1 + (10**this_log)
        if reads_done + reads > reads2get:
            reads = reads2get - reads_done
        out.append(int(reads))
        this_log   += loginterval
        reads_done += reads

    #output is shuffled to randomize order
    shuffle(out)
    return(out)

def errorify(perfect):
    #generates errors based on opt_e value
    #simulates sequencing errors/mutations

    pick = float(randint(1,1e5) / 1e5)

    try:
        #if an error is picked, a single base will be mutated to one of the 3 others
        if pick < opt_e:

            mutated = randint(0,len(perfect)-1)


            old_base = perfect[mutated]
            new_base = get_new_base(old_base)

            after = len(perfect) - mutated - 1
            MD = '%(mutated)s%(new_base)s%(after)s' % locals()



            s = list(perfect)
            s[mutated] = new_base.lower()
            err_read = "".join(s)
            err_read = (MD,err_read)

        else:
            err_read = (0,perfect)
        return(err_read)
    except ValueError:
        return("error")

def get_new_base(old_base):
    #selects a new base for a given error simulation
    #can't choose itself

    new_base = old_base
    while new_base in old_base:
        #print "while"
        pick_base = randint(1,4)
        if pick_base == 1:
            new_base = "A"
        if pick_base == 2:
            new_base = "T"
        if pick_base == 3:
            new_base = "G"
        if pick_base == 4:
            new_base = "C"
    return(new_base)

def pick_a_strand():
    #MIR
    #chooses + or - strand, randomly. returns + or -

    pick = randint(1,2)
    if pick == 1:
        return('+')
    if pick == 2:
        return('-')

def pick_an_arm(locus):
    #MIR
    #chooses the orientation of the miRNA locus, returning the coordinates for each feature

    loc_start = int(locus[0])
    left_start  = loc_start + 17
    left_stop   = loc_start + 17 + 20
    right_start = loc_start + 17 + 20 + 48
    right_stop  = loc_start + 17 + 20 + 48 + 20

    pick = randint(1,2)
    if pick == 1:
        return( (right_start, right_stop, left_start, left_stop) )
    if pick == 2:
        return( (left_start, left_stop, right_start, right_stop) )

def get_mir_mod(start, stop, strand, modtype):
    #MIR
    #applies modification type to start and stop loci, in a strand dependent manner

    if "+" in strand:
        start = start + modtype[0]
        stop  = stop  + modtype[1]
    if "-" in strand:
        start = start - modtype[1]
        stop  = stop  - modtype[0]
    return( (start,stop) )

def get_mir_modseq(mod, locus, seq, strand):
    #MIR
    #grabs the exact sequence for a given modification

    loc_start = locus[0]
    loc_stop  = locus[1]
    if "+" in strand:
        mod_start = mod[0]
        mod_stop  = mod[1] + 1
    if "-" in strand:
        mod_start = mod[0]
        mod_stop  = mod[1] + 1

    #offset = mod_start = loc_start
    #length = mod_stop - mod_start + 1
    modseq_f = seq[mod_start-loc_start:mod_stop-loc_start]

    if "-" in strand:
        modseq = rev_comp(modseq_f)

    else:
        modseq = modseq_f

    return (modseq)

def get_mir_perfects(locus, seq, strand, mir_star):
    #MIR
    #this function is the parent of several, generating a dictionary
    #of all possible miRNA variants for a given locus
    #key: ["modification_type:mir_or_star"]
    #val: (start, stop, seq)
    #returns this dictionary, as a complete list of possibilities

    mir_start = mir_star[0]
    mir_stop  = mir_star[1]
    hash = {}


    mir_or_star = "miR"
    modtype = (0,0)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (0,-1)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (0,-2)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (1,1)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (1,0)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-1,-1)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-1,-2)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-2,-2)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-2,-3)
    mod    = get_mir_mod(mir_start, mir_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)


    star_start = mir_star[2]
    star_stop  = mir_star[3]

    mir_or_star = "star"

    modtype = (0,0)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (0,-1)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (0,-2)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (1,1)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (1,0)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-1,-1)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-1,-2)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-2,-2)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    modtype = (-2,-3)
    mod    = get_mir_mod(star_start, star_stop, strand, modtype)
    modseq = get_mir_modseq(mod, locus, seq, strand)
    hash["%s%s:%s" % (mir_or_star, modtype[0], modtype[1])] = (mod[0], mod[1] ,modseq)

    return(hash)

def rev_comp(seq):
    #returns the exact reverse complement of an input sequence

    out = []
    inp = list(seq)
    for k in reversed(inp):
        #upper case
        if k == "A":
            out.append("T")
        if k == "T":
            out.append("A")
        if k == "G":
            out.append("C")
        if k == "C":
            out.append("G")

        #lower case
        if k == "a":
            out.append("t")
        if k == "t":
            out.append("a")
        if k == "g":
            out.append("c")
        if k == "c":
            out.append("g")
    out = "".join(out)
    return (out)

def sample_mir_keys():
    #MIR
    #This randomly chooses what modification types have taken place on a read
    #and whether a read is mir or star.  This is based on a probability table.
    #miR with no modifications is an 80% chance, for example

    pick = randint(0,1000)
    if pick <  600:
        key = "miR0:0"
    if pick >= 600 and pick < 800:
        key = "star0:0"
    if pick >= 800 and pick < 840:
        key = "miR0:-1"
    if pick >= 840 and pick < 850:
        key = "miR0:-2"
    if pick >= 850 and pick < 890:
        key = "miR1:1"
    if pick >= 890 and pick < 900:
        key = "miR1:0"
    if pick >= 900 and pick < 920:
        key = "miR-1:-1"
    if pick >= 920 and pick < 925:
        key = "miR-1:-2"
    if pick >= 925 and pick < 945:
        key = "miR-2:-2"
    if pick >= 945 and pick < 950:
        key = "miR-2:-3"
    if pick >= 950 and pick < 960:
        key = "star0:-1"
    if pick >= 960 and pick < 965:
        key = "star0:-2"
    if pick >= 965 and pick < 975:
        key = "star1:1"
    if pick >= 975 and pick < 980:
        key = "star1:0"
    if pick >= 980:
        key = "star-1:-1"
    return(key)

def get_tasi_perfects(locus, seq):
    #TASI
    #builds a dictionary of possible tasi perfects
    #since these are phased, the strand, size and phase number are stored in the key
    #val = (start, stop, seq)

    start = locus[0]
    stop  = locus[1]
    hash  = {}

    for phase in range(0,6):

        #top strand first

        offset = 4 + (21 * phase)

        #21 mer
        tasi_seq = seq[offset : offset + 21]
        tasi_start = start + offset
        tasi_stop  = start + offset + 21 - 1
        hash["top_%s_21" % (phase)] = (tasi_start, tasi_stop, tasi_seq)

        #20 mer
        tasi_seq = seq[offset : offset + 20]
        tasi_start = start + offset
        tasi_stop  = start + offset + 20 - 1
        hash["top_%s_20" % (phase)] = (tasi_start, tasi_stop, tasi_seq)

        #22 mer
        tasi_seq = seq[offset : offset + 22]
        tasi_start = start + offset
        tasi_stop  = start + offset + 22 - 1
        hash["top_%s_22" % (phase)] = (tasi_start, tasi_stop, tasi_seq)


        #now the bottom strand

        offset = 2 + (21*phase)

        #21 mer
        tasi_seq = rev_comp(seq[offset : offset + 21])
        tasi_start = start + offset
        tasi_stop  = start + offset + 21 - 1
        hash["bot_%s_21" % (phase)] = (tasi_start, tasi_stop, tasi_seq)

        #20 mer
        tasi_seq = rev_comp(seq[offset : offset + 20])
        tasi_start = start + offset
        tasi_stop  = start + offset + 20 - 1
        hash["bot_%s_20" % (phase)] = (tasi_start, tasi_stop, tasi_seq)

        #22 mer
        tasi_seq = rev_comp(seq[offset : offset + 22])
        tasi_start = start + offset
        tasi_stop  = start + offset + 22 - 1
        hash["bot_%s_22" % (phase)] = (tasi_start, tasi_stop, tasi_seq)

    return(hash)

def sample_tasi_keys():
    #TASI
    #This randomly chooses the phase, strand and read length.  Length is based on a probability table.
    #80% chance of 21 nt length, for example
    #pure random for strandedness and phase

    strand_pick = randint(0,1000)

    if strand_pick < 500:
        strand = "top"
    else:
        strand = "bot"

    phase_pick = randint(0,5)
    size_pick  = randint(0,1000)

    if size_pick < 800:
        size = 21
    if size_pick >= 800 and size_pick < 900:
        size = 20
    else:
        size = 22

    key = "%s_%s_%s" % (strand, phase_pick, size)

    return(key)

def get_het_perfects(locus, seq):
    #     start = locus[0]
    #     stop  = locus[1]
    #     hash  = {}



    perfects = {}



    for length in range(35,61):

        for r in range(0,len(seq)):
            transcript = seq[r : r + length]



            if len(transcript) == length:


                for len_sRNA in range(21,25):

                    start = locus[0] + r
                    stop  = start + len_sRNA - 1
                    perfects['%sbase_%slen_L_%smer_top' % (r, length, len_sRNA)] = (start, stop, transcript[0:len_sRNA])
                    perfects['%sbase_%slen_L_%smer_bot' % (r, length, len_sRNA)] = (start, stop, rev_comp(transcript[0:len_sRNA]))

                    stop  = locus[0] + r + length - 1
                    start = stop - len_sRNA + 1
                    perfects['%sbase_%slen_R_%smer_top' % (r, length, len_sRNA)] = (start, stop, transcript[-len_sRNA:])
                    perfects['%sbase_%slen_R_%smer_bot' % (r, length, len_sRNA)] = (start, stop, rev_comp(transcript[-len_sRNA:]))

    return (perfects)

def get_het_seq(key, seq, locus):
    # key = (pos_pick, tran_len, end, size, strand)

    transcript = seq[key[0] : key[0] + key[1]]

    if key[2] == "L":
        start = locus[0] + key[0]
        stop  = start + key[3] - 1

        if key[4] == "top":
            return(start, stop, transcript[0:key[3]])
        else:
            return(start, stop, rev_comp(transcript[0:key[3]]))

    else:
        stop  = locus[0] + key[0] + key[1] - 1
        start = stop - key[3] + 1

        if key[4] == "top":
            return(start, stop, transcript[-key[3]:])
        else:
            return(start, stop, rev_comp(transcript[-key[3]:]))






def sample_het_keys(loci_size):
    #HET
    #This randomly chooses the strand, position and read length.
    #strandedness and position in loci are pure random
    #length is decreasing likelihood from 24 to 21 nt

    # picking top or bottom strand
    strand_pick = randint(0,1)
    if strand_pick == 1:
        strand = "top"
    else:
        strand = "bot"

    # picking right or left end of primary transcript
    end_pick = randint(0,1)
    if end_pick == 1:
        end = "R"
    else:
        end = "L"

    # picking size of sRNA read: 21 - 24 bp, weighted towards 24
    size_pick = randint(0, 1000)
    if size_pick < 900:
        size = 24
    elif size_pick >= 900 and size_pick < 950:
        size = 23
    elif size_pick >= 950 and size_pick < 980:
        size = 22
    else:
        size = 21

    # picking length of primary transcript: 35-60 bp, weighted towards smaller by log10 scale
    # 10^1.544 = 35 and 10^1.778 = 60
    while True:
        tran_len = float(randint(1544, 1778))
        tran_len = int(10 ** (tran_len / 1000))
        if tran_len > 34 and tran_len < 61:
            break



    # picking location in loci for primary transcript formation
    while True:
        pos_pick = randint(0, loci_size)

        if pos_pick + tran_len < loci_size:
            break


    key = (pos_pick, tran_len, end, size, strand)
    return(key)

def get_het_interval():
    grab = randint(230, 300)
    grab = float(grab) / 100
    return(int(10 ** grab))


def RNA_assignment(RNA_type, loci_depth):

    #This function performs most of the calls for building the simulated library
    #RNA_type = 'miRNA', 'tasiRNA', or 'siRNA'
    #interval = loci size, based on type
    #loci_depth is a list of loci, with the value being the number of
    #reads to simulate for that site (generated by get_nts)

    #gathers up all possible positions for this type
    genome_positions = get_bins(RNA_type)



    locus_n = 0
    read_count = 0
    #total is used in generating percentage bar
    total_reads = sum(loci_depth)

    print "\nSampling %s" % (RNA_type)
    print total_reads, "reads to assign"
    print len(loci_depth), "loci to assign"
    print len(genome_positions), "available bins"

    #iterates through the list of depths
    for depth in loci_depth:
        while True:

            if RNA_type == 'siRNA':
                interval = get_het_interval()
            elif RNA_type == 'miRNA':
                interval = 125
            elif RNA_type == 'tasiRNA':
                interval = 140

            #grab a locus from possible positions
            gotten = get_a_locus(genome_positions, interval)
            locus = gotten[0]
            chosen_bin = gotten[1]

            #look for overlaps with occupied bins
            check = occupied_checker(occupied_list, locus, RNA_type)

            #grabs the sequence and checks for filters
            seq = get_a_seq(locus)
            if filter(seq) == False:

                if seq != None:
                    if len(seq) > 3:

                        locus_n += 1

                        #did it return an overlap w/ an occupied bin?
                        if check == 1:
                        #yes, try again
                            pass
                        else:
                        #adds this sequence to the occupieds and removes from possible choices
                            occupied_list.append(check)
                            genome_positions.remove(chosen_bin)

                            ##################################
                            # siRNA Simulation

                            if RNA_type == 'siRNA':

                                #for the number of reads in a loci
                                for i in range(1, depth+1):

                                    read_count += 1
                                    # percent_complete(total_reads, read_count)

                                    #chooses a key for this read
                                    key = sample_het_keys(interval)

                                    #grabs the sequence, start and stop for this read.
                                    het_seq = get_het_seq(key, seq, locus)

                                    if "top" in key:
                                        het_strand = "+"
                                    if "bot" in key:
                                        het_strand = "-"

                                    #loci-info output
                                    OUTS = "HET_%s\t%s\t.\t%s\t%s\n" % (locus_n, locus_index(locus).chr_coord(), depth, len(seq))



                                    #fasta-file output (sequences errorified)
                                    err_seq = errorify(het_seq[2])
                                    if err_seq == 'error':
                                        pass
                                    else:
                                        OUTF = ">HET_%s_%s_%s_%s_%s_loci%sbp\n%s\n" % (locus_n, i, locus_index(het_seq).chr_coord(), het_strand, err_seq[0], len(seq), err_seq[1])

                                        if len(err_seq[1]) > 15:
                                            #writes fasta
                                            with open(fasta_output,"a+") as f:
                                                f.write(OUTF)

                                #writes info
                                with open(info_output,"a+") as f:
                                    f.write(OUTS)

                                break


                            ##################################
                            # tasiRNA Simulation

                            if RNA_type == 'tasiRNA':

                                #building list of perfects
                                tasi_perfects = get_tasi_perfects(locus, seq)

                                #for the number of reads in a loci
                                for i in range(1, depth+1):

                                    read_count += 1
                                    percent_complete(total_reads, read_count)

                                    #chooses a key for this read
                                    key = sample_tasi_keys()
                                    if "top" in key:
                                        tasi_strand = "+"
                                    if "bot" in key:
                                        tasi_strand = "-"

                                    #loci-info output
                                    OUTS = "TAS_%s\t%s\t.\t%s\n" % (locus_n, locus_index(locus).chr_coord(), depth)

                                    #fasta-file output (sequences errorified)
                                    length = len(tasi_perfects[key][2])
                                    if final_read_check(length) is True:
                                        err_seq = errorify(tasi_perfects[key][2])
                                        OUTF = ">TAS_%s_%s_%s_%s_%s\n%s\n" % (locus_n, i, locus_index(tasi_perfects[key]).chr_coord(), tasi_strand, err_seq[0], err_seq[1])
                                        #writes fasta
                                        with open(fasta_output,"a+") as f:
                                            f.write(OUTF)
                                            pass

                                #writes info
                                with open(info_output,"a+") as f:
                                    f.write(OUTS)
                                break

                            ##################################
                            # miRNA Simulation

                            if RNA_type == 'miRNA':

                                #chooses the characteristics of this loci (strand/orientation)
                                mir_strand = pick_a_strand()
                                mir_star = pick_an_arm(locus)

                                #builds a list of perfect reads
                                mir_perfects = get_mir_perfects(locus, seq, mir_strand, mir_star)


                                #loci-info output
                                OUTS = "MIRNA_%s\t%s\t%s\t%s\n" % (locus_n, locus_index(locus).chr_coord(), mir_strand, depth)

                                #writes loci-info
                                with open(info_output,"a+") as f:
                                    f.write(OUTS)

                                #for the number of reads in a loci
                                for i in range(1, depth+1):

                                    read_count += 1
                                    percent_complete(total_reads, read_count)

                                    #chooses a key for this read
                                    key = sample_mir_keys()

                                    #fasta-file output (sequence errorified)
                                    err_seq = errorify(mir_perfects[key][2])
                                    if err_seq == 'error':
                                        pass
                                    else:
                                        OUTF = ">MIRNA_%s_%s_%s_%s_%s\n%s\n" % (locus_n, i, locus_index(mir_perfects[key]).chr_coord(), mir_strand, err_seq[0], err_seq[1])

                                        #writes fasta
                                        with open(fasta_output,"a+") as f:
                                            f.write(OUTF)

                                break
    return(fasta_output)


def find_mmap_reads(fa_loc):

    devnull = open(os.devnull, 'wb')

    proc_start = time()



    fa_loc_stripped = fa_loc.strip('.fa')

    max_output = '.%(fa_loc_stripped)s_mmap.fa' % locals()
    al_output  = '.%(fa_loc_stripped)s_unique.fa' % locals()
    un_output  = '.%(fa_loc_stripped)s_unaligned.fa' % locals()

    #process for bowtie alignment
    p1_call = "bowtie -f --sam -m 1 --max %s --al %s --un %s -p %s -v 1 --best --strata %s %s" % (max_output, al_output, un_output, bow_p, opt_g, fa_loc)


    print '\nQuantifying multimappers with the following call: \n%s' % (p1_call)
    p_bowtie = subprocess.Popen([p1_call], shell = True, stdout=devnull, bufsize = 1)

    p_bowtie.wait()

    print''
    print 'writing:', fa_loc
    with open(fa_loc, 'w') as output_f:
        if os.path.isfile(max_output) == True:
            with open(max_output, 'r') as f:
                for line in f:
                    if line[0] == '>':
                        output_f.write('%s_mmap\n' % (line.strip()))
                    else:
                        output_f.write(line)
        if os.path.isfile(al_output) == True:
            with open(al_output, 'r') as f:
                for line in f:
                    if line[0] == '>':
                        output_f.write('%s_uniq\n' % (line.strip()))
                    else:
                        output_f.write(line)
        if os.path.isfile(un_output) == True:
            with open(un_output, 'r') as f:
                for line in f:
                    if line[0] == '>':
                        output_f.write('%s_nomap\n' % (line.strip()))
                    else:
                        output_f.write(line)

    proc_end = time()

    print "Multimap Identification --> Time elapsed: %s minutes" % (round((proc_end-proc_start) / 60, 2))

################################################################################
ver_num = 2.0
program_name = sys.argv[0]

usage = """\n%s: %s

Simulate plant small RNA-seq data

Dependencies:
python 2.6.6
samtools 1.1
bowtie 1.0

Usage:
sim_sRNA_library.py [options] -f sample_lib.fasta -m species_specific_miRNA_annotation.gff3 -g indexed_reference_genome
  --- OR ---
sim_sRNA_library.py [options] -b sample_lib_alignment_(all).bam -m species_specific_miRNA_annotation.gff3 -g indexed_reference_genome
  --- OR ---
sim_sRNA_library.py [options] -p psuedo_annotation.txt -g indexed_reference_genome

Required:
-f : fasta file of template sRNA library, already adapter trimmed.  This will be the basis for producing a psuedo_annotation.
-s : sam file of aligned sRNA library.  For proper use of the program, this should be produced using: bowtie --all --best --strata -v 0
-m : annotation of known mature miRNAs for the species of interest, gff3 format.  Available from miRBase.org
-p : location of psuedo_annotation file (.txt), generated by this program

Options:
-o : output identifier for psuedo_annotation and simulated library generation. Default: random.
-v : print version number and quit
-h : print help message and quit
-r : desired total number of reads, in millions. Default: 5
-e : per-read probability of a single nt sequencing error (substitution). Default: 1E-4
--suppress_percent : when called, percent progress bars will not be displayed.
--bowtie_p : number of processing cores for bowtie. Increases speed with higher RAM usage. Default: 1

""" % (program_name, ver_num)


print '\n--Initializing--'

if "-h" in sys.argv or len(sys.argv) == 1:
    print usage
    sys.exit()

if "-v" in sys.argv:
    print "%s: %s\n" % (program_name, ver_num)


if "-r" in sys.argv:
    opt_r = sys.argv[int(sys.argv.index('-r'))+1]
    if float(opt_r) != int(opt_r):
        print "Fatal: Option -r must be an integer.\n%s\n" % (usage)
        sys.exit()
    elif int(opt_r) < 1 or int(opt_r) > 5:
        print "Fatal: Option -r must be an between 1 and 5.\n"
        sys.exit()
    else:
        opt_r = int(opt_r) * 1000000
else:
    opt_r = 5000000


if "-e" in sys.argv:
    opt_e = sys.argv[int(sys.argv.index('-e'))+1]
    if float(opt_e) >= 0 and float(opt_e) < 1:
        opt_e = float(opt_e)
    else:
        print "Fatal: Invalid value for option -e.\n%s\n" % (str(usage))
else:
    opt_e = 0.0001

if "-o" in sys.argv:
    opt_o = sys.argv[int(sys.argv.index('-o'))+1]

else:
    value = datetime.datetime.now()
    mins = ''.join(str(value).split(' ')[-1].split('.')[0].split(':'))
    ms = str(value).split(' ')[-1].split('.')[1][0:2]
    opt_o = ''.join((mins,ms))
    print "Output ID: not specified, assigning -->", opt_o

if '-p' not in sys.argv:
    if '-m' not in sys.argv:
        print "Fatal: Annotation of miRNA sequences required for generation of psuedo-annotation.\n%s\n" % (usage)
        sys.exit()

    if '-b' not in sys.argv:
        if '-f' not in sys.argv:
            print "Fatal: Alignment (-b, .bam) or reads (-f, .fasta) required for generation of psuedo-annotation.\n%s\n" % (usage)
            sys.exit()

if '-p' in sys.argv:
    opt_p = sys.argv[int(sys.argv.index('-p'))+1]
    validate_path(opt_p)

    run_type = 'p'
else:
    opt_p = 'pseudo_%s.txt' % (opt_o)
    if "-b" in sys.argv:
        opt_bam = sys.argv[int(sys.argv.index('-b'))+1]
        validate_path(opt_bam)
        run_type = 'b'

    elif "-f" in sys.argv:
        opt_f = sys.argv[int(sys.argv.index('-f'))+1]
        validate_path(opt_f)
        run_type = 'f'
        opt_bam = 'sim_%s.bam' % (opt_o)
    if "-m" in sys.argv:
        opt_mir = sys.argv[int(sys.argv.index('-m'))+1]
        validate_path(opt_mir)




if "-g" in sys.argv:
    opt_g = sys.argv[int(sys.argv.index('-g'))+1]
    validate_path(opt_g)
else:
    print "Fatal: Must specify a reference genome location.  This genome must be the same for psuedo annotation and simulated library construction.\n%s\n" % (usage)
    sys.exit()

if "--suppress_percent" in sys.argv:
    opt_suppress = True
else:
    opt_suppress = False

if "--bowtie_p" in sys.argv:
    bow_p = sys.argv[int(sys.argv.index('--bowtie_p'))+1]
else:
    bow_p = 1
if "-bowtie_p" in sys.argv:
    print "unrecognized command: '-bowtie_p'  (did you mean --bowtie_p?)"
    sys.exit()

#checking for required software
try:
    subprocess.check_call(["which","samtools"], stdout=subprocess.PIPE)
    print "samtools installation: Found"
except:
    print "FATAL: samtools check failed. samtools must be installed to use this script.\n%s\n" % (usage)
    sys.exit()
try:
    subprocess.check_call(["which","bowtie"], stdout=subprocess.PIPE)
    print "bowtie installation: Found"
except:
    print "FATAL: bowtie check failed. samtools must be installed to use this script.\n%s\n" % (usage)
    sys.exit()

print '\n--Settings--'
print 'simulator.py version:', ver_num
print 'Output ID:', opt_o
print 'Err. rate:', opt_e
print 'Read count:', opt_r
print 'Suppress perc.:', opt_suppress
print 'Bowtie_p:', bow_p


if run_type == 'p':
    print 'Run-type:   p (using pseudo_annotation)'
    print 'Psuedo_Annotation: ', opt_p
if run_type == 'f':
    print 'Run-type:   f (aligning fasta-formatted template and building psuedo_annotation)'
    print 'Template lib:      ', opt_f
    print 'MIR annotation:    ', opt_mir
if run_type == 'b':
    print 'Run-type:   b (using alignment file to build psuedo_annotation)'
    print 'Template alignment:', opt_bam
    print 'MIR annotation:    ', opt_mir

################################################################################
# Function Calls



start = time()

#bin size has been non-optionally set to 150 nt
bin_size = 150


print '\n--Validating Genome--'
#validating genome and related indexes
validate_genome()
#creating a dictionary to translate chromosomal coords to genomic coords
chrom_id_dict = chrom_id(opt_g)

#calculating the number of reads for each sRNA type
mir_reads =  int(opt_r * 0.3)
tasi_reads = int(opt_r * 0.05)
het_reads =  opt_r - tasi_reads - mir_reads

#generating lists of loci, carrying the value of how many reads should
#come from that loci.  These lists are shuffled so the depths are not
#sequential, but in random order.
mir_ns = get_nts(mir_reads, 100)
tasi_ns = get_nts(tasi_reads, 20)
het_ns = get_nts(het_reads, 10000)

# if f(asta) is chosen, the program will build a pseudo-ann based off of a new bowtie alignment
if 'f' in run_type:
    print "\n--Building Pseudo-Annotation--"
    print "Running bowtie for template read library"
    print "Calling: bowtie --all --best --strata -v 0 [genome] [fasta]"
    print "(this will take some time, due to --all option)"
    print "Increasing --bowtie_p from default: 1 can improve speed, for systems"
    print "with more RAM"
    bowtie()

    align_bins = miRNA_bins()
    align_bins = other_bins(align_bins)
    bin_info = bin_quantify(align_bins)
    align_bins['totals'] = bin_info['bin_count']
    align_bins['num_of_bins'] = bin_info['alignment_count']
    min_bin_depth = bin_info['min_bin_depth']
    target_list()

# if b(am) is chosen, the program will build a pseudo-ann based off an old alignment
# for correct results, this alignment must be performed as described in fasta-runtype
if 'b' in run_type:
    print "Bam-file of template library alignments: %s --> Found" % (opt_bam)
    align_bins = miRNA_bins()
    align_bins = other_bins(align_bins)
    bin_info = bin_quantify(align_bins)
    align_bins['totals'] = bin_info['bin_count']
    align_bins['num_of_bins'] = bin_info['alignment_count']
    min_bin_depth = bin_info['min_bin_depth']
    target_list()

# if p(seudo-annotation) is chosen, the program skips this building step
if 'p' in run_type:
    #this prints off header information from annotation
    print "\nPseudo-Annotation file:  %s --> Found" % (opt_p)
    print "Header info:"
    with open(opt_p, 'r') as file:
        for line in file:
            if line[0] == '#':
                print line[:-1]
            else:
                break


#Initiates new output files, writing over old ones
output_ID = opt_o

info_output = "./sim_" + output_ID + ".txt"
with open(info_output, "w") as f:
    f.write("")
fasta_output = "./sim_" + output_ID + ".fa"
with open(fasta_output, "w") as f:
    f.write("")

#initializing a list of occupied loci
occupied_list = [(0,1)]

#calling the actual simulator function
output_fasta_location = RNA_assignment('miRNA', mir_ns)
output_fasta_location = RNA_assignment('tasiRNA', tasi_ns)
output_fasta_location = RNA_assignment('siRNA', het_ns)

find_mmap_reads(output_fasta_location)



print ""

end = time()

print "Simulation complete!"
print "total simulator --> Time elapsed: %s minutes" % (round((end-start) / 60, 2))




sys.exit()

#program complete


README = '''SYNOPSIS
    sim_sRNA_library.py - v0.3

    simulation small-RNA seq libraries based on a template library

    Copyright (C) 2015  Nathan R. Johnson; Michael J. Axtell

AUTHORS
    Nathan R. Johnson, Penn State University, jax523@gmail.com
    Michael J. Axtell, Penn State Universtiy, mja18@psu.edu

LICENSE
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


DEPENDENCIES (all in PATH)
    python 2.6.6
    samtools 1.1
    bowtie 1.0

USAGE
    sim_sRNA_library.py [options] -f sample_lib.fasta -m species_specific_miRNA_annotation.gff3 -g indexed_reference_genome
      --- OR ---
    sim_sRNA_library.py [options] -b sample_lib_alignment_(all).bam -m sim_sRNA_library_specific_miRNA_annotation.gff3 -g indexed_reference_genome
      --- OR ---
    simulator.py [options] -p psuedo_annotation.txt -g indexed_reference_genome

REQUIRED
    -f : fasta file of template sRNA library, already adapter trimmed.  This will be the basis for producing a psuedo_annotation
    -s : sam file of aligned sRNA library.  For proper use of the program, this should be produced using: bowtie --all --best --strata -v 0
    -m : annotation of known mature miRNAs for the species of interest, gff3 format.  Available from miRBase.org
    -p : location of psuedo_annotation file (.txt), generated by this program

OPTIONS
    -o : output identifier for psuedo_annotation and simulated library generation. Default: random.
    -v : print version number and quit
    -h : print help message and quit
    -r : desired total number of reads, in millions. Default: 5
    -e : per-read probability of a single nt sequencing error (substitution). Default: 1E-4
    --suppress_percent : when called, percent progress bars will not be displayed.

METHODS     

    Approach

        This method uses real sRNA-seq data as a template to generate simulated
        sequencing data of user-selectable number of reads. This has been
        developed to test the accuracy of small RNA alignment methods, as the
        known origin of a read allows a researcher to know the rate of
        misaligned reads.

    Genome

        There is no strict requirement for genomic masking in this program, but
        it is recommended to use non- or soft-masked genomes, as sRNAs may come
        from highly repetititve portions of the genome. It is imparative that
        the reference genome and mir annotation (.gff3) have the same convention
        for chromosome names.

        Additionally, the genome must be indexed for bowtie, using the bowtie-
        build module. This program expects bowtie indecies to be present in the
        reference genome folder, looking for the presence of the following file
        as proof:

        reference:   Osativa.fa ebwt proof:  Osativa.fa.1.ebwt

    Template annotation - miRNA

        Identification of candidate locations for miRNAs comes from annotation
        data in the form of a .gff3 file. These data are available for many
        species at mirbase.org.

    Template library - si, tasi, non-sRNA

        sRNA-seq data in fasta format is used as the template input. This
        library should already be adapter trimmed. The library is then aligned
        using bowtie --all --best --strata -v 0, reporting all alignments of
        every read and saving it as a .bam, binary alignment file. This file may
        be used for subsequent runs of the program to save alignment time.

        This template is used in the identification of candidate locations for
        heterochromatic siRNA, tasiRNA and garbage RNA (fragmented or non-sRNA).

    Pseudo-annotation construction

        Using the alignment and annotation data from the previous steps, the
        program constructs a list of candidate 'bins', 150 nt in length. All
        reads falling within a bin will be stored as a depth value for their
        given type. These types are defined simply by size, with siRNA
        encompassing 23 and 24 nt reads, tasiRNA encompassing 21 nt reads, and
        garbage RNA encompassing reads between 15 and 20 nt, as well as over 24
        nt in length.

        Once depths have been completely assigned, each bin will be chosen for
        only one type of small RNA.

        To be chosen for a given type, that type must have an 80%% majority in
        that bin, giving expression level importance.  Any bins containing a
        miRNA read will be automatically chosen as miRNA, therefore discluded
        from the other types.  A minimum depth is also required for a bin to be
        chosen for a type. Min-depth is 1-5 reads deep, and will be selected
        automatically by the program, to choose the most strict requirement that
        still allows enough loci for the later simulation steps.  If a minimum
        depth of 1 does not allow enough loci for later steps, the program deems
        this library to be of too low quality to be an adequate template.

        Once all bins have been chosen for a type, their genomic locations are
        output to a .txt file known as a pseudo-annotation, which will be used
        for subsequent steps.  This file may be reused for repeated simulations,
        as it is non-probabalistic for a given library.

    Read numbers and types

        30%% of the simulated reads will come from roughly 100 MIRNA loci, 5%%
        of the simulated reads will come from roughly 20 tasiRNA/phased
        secondary siRNA loci, and the remaining 65%% from roughly 10,000
        heterochromatic siRNA loci.

        The abundance of reads from each locus is distributed on a log-linear
        scale (e.g., plotting the log10 of read number as a function of
        abundance rank yields a straight line).

    MIRNA simulation

        Valid MIRNA loci are selected from annotated MIRNA loci, based on the
        'pseudo-annotation'. The locus has a hypothetical size of 125nts.

        The strand of the MIRNA precursor is randomly selected, as is the arm
        from which the mature miRNA and star come from. There is no actual
        hairpin sequence necessarily present at MIRNA loci .. it is only the
        pattern of reads that is being simulated.

        The mature miRNA and mature miRNA* are defined as 'master' positions.
        The left-most 'master' position in a locus is a 21-mer starting at
        position 17 of the locus. The right-most 'master' position is a 21 mer
        at position 85. Assignment of the arm (i.e., whether the left-most or
        right-most 'master' position is the miRNA or miRNA*) is random at each
        locus.

        Once a locus has been found and 'master' positions defined, each read is
        simulated according to the following probabilities. In the following
        list, "miR" means mature miRNA, "star" means miRNA*. The numbers after
        each indicate the offset at 5' and 3' ends relative to the master
        positions. So, "miR0:1" means the mature miRNA sequence, starting at the
        master 5' end, and ending 1 nt after the master 3' end.

        60%%  miR0:0, 20%% star0:0,, 4%% miR0:-1, 1%% miR0:-2, 4%% miR1:1, 1%%
        miR1:0, 2%% miR-1:-1, 0.5%% miR-1:-2, 2%% miR-2:-2, 0.5%% miR-2:-3, 1%%
        star0:-1, 0.5%% star0:-2, 1%% star1:1, 0.5%% star1:0, 2%% star-1:-1.

        Sampling of simulated MIRNA-derived reads continues until the required
        number of reads for a particular locus is recovered.

    tasiRNA/phased siRNA simulation

        TAS loci are selected from bins chosen in the 'pseudo-annotation'. Each
        locus has a nominal size of 140nts.

        Each locus is simulated to be diced in 6 21 nt phases. At each phasing
        position, 21mers are the dominant size, with 20mers and 22mers being
        less frequent .. the 20 and 22nt variants vary in their 3' positions
        relative to the 'master' 21nt RNAs.

        Once a locus has been identified, and all possible 20, 21, and 22mers
        charted, each read is simulated according to the following
        probabilities:

        Strand of origin is 50%% top, 50%% bottom.

        Phase position is equal chance for all (e.g. 1/6 chance for any
        particular phase location).

        80%% of the time, the 21mer is returned, 10%% of the time the 20mer, and
        10%% the 22 mer.

    Heterochromatic siRNA simulation

        Heterochromatic siRNA loci are selected from bins chosen in the 'pseudo-
        annotation'. Their nominal locus size is 200-1000 nts, logarithmically
        weighted to smaller loci.

        Precursors are simulated as 35-60nt regions, logarithmically weighted to
        smaller precursors. Eligible products from a loci are 21-24nt in length,
        coming from either strand or end of the precursor. eligible from these
        loci. Each simulated read is identified with the following
        probabilities:

        50%% chance for top of bottom strand origin.

        50%% chance for left or right end of precursor.=.

        90%% chance of a 24 mer, 5%% chance of a 23 mer, 3%% chance of a 22 mer,
        and 2%% chance of a 21 mer.

    Simulation of sequencing errors

        Regardless of small RNA type, once a given read has been simulated,
        there is a chance of introucing a single-nucleotide substitution
        relative to the reference genome at a randomly selected position in the
        read. The chance of introducing this error is given by option -e
        (Default is 1 in 10,000).

    No overlapping loci

        None of the simulated loci are allowed to have any overlap with each
        other.

    Filtering of sequences

        Simple filtering is applied to remove sequences that are ambiguous or
        easily identifiable as highly repetitive.  Ambiguous sequences are said
        to contain any "N" bases (non-ATGC). Dinucleotide repeats are filtered
        and defined as 8 uninterrupted pairs of the same 2 bases, ex.
        'ATATATATATATATAT'.  No other sequences are filtered, to avoid bias in
        selecting sequences.

OUTPUT   

    Files

        Four files are created in the working directory.  Option -o defines the
        [out] character.

        sim_[out].bam : A binary alignment file of the input template library.

        pseudo_[out].txt : A pseudo-annotation text file, identifying bins for
        simulated reads.

        sim_[out].txt : A tab-delimited text file giving the coordinates and
        names for each simulated locus.

        sim_[out].fa : A FASTA file of all of the simulated reads.

    Naming conventions

        The FASTA header of each simulated read encode the basic information
        about the read. Several different fields are separated by "_"
        characters. For instance the following read header
        ">MIRNA_2_2_chr7:123063894-123063914_+_0" means ..

        MIRNA: This came from a simulated MIRNA locus

        2: The 2nd simulated locus of this type

        2: Read number 2 from this locus

        chr7:123063894-123063914: The true origin of this read.

        +: The genomic strand of this read.

        0: The number of sequencing errors simulated into the read.'''
