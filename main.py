from maxentpy.maxent import load_matrix5, load_matrix3
from maxentpy import maxent_fast
from maxentpy.maxent_fast import load_matrix
from itertools import product
import os
import sys
import argparse
import numpy as np
import pybedtools
import gzip

global mat5
global mat3

mat5 = load_matrix5()
mat3 = load_matrix3()

def parser_args(args):
    parser = argparse.ArgumentParser(prog="novasplice")
    parser.add_argument('-v', '--vcf', help="Full path to the vcf file being used", type=str, required=True)
    parser.add_argument('-r', '--reference', help="Full path to the reference genome being used", type=str, required=True)
    parser.add_argument('-g', '--gtf', help="Full path to the reference gtf being used", type=str, required=True)
    parser.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is working directory under /novasplice_output", type=str, required=False, default="./novasplice_output")
    return parser.parse_args(args)

def randomlygenerate_threess():
    choices = ["A", "C", "G", "T"]
    return ''.join(np.random.choice(choices, size=23))

def generate_all_kmers(k):
    return product(["A", "C", "G", "T"], repeat=k)

def compute_five_score(dna):
    if len(dna) != 9:
        print("ERROR INCORRECT LENGTH: %s" % dna)
    return maxent_fast.score5(dna, matrix=mat5)

def compute_three_score(dna):
    if len(dna) != 23:
        print("ERROR INCORRECT LENGTH: %s" % dna)
    return maxent_fast.score3(dna, matrix=mat3)

def generate_splicingbedfile_fromgtf(gtf, output):
    with gzip.open(gtf, 'rt') as g:
        with open(os.path.join(output, "splice-site.bed"), 'w') as bed:
            for lines in g:
                line = lines.rstrip().split()
                if line[2] != "exon": continue
                exstart = int(line[3])
                exend = int(line[4])
                chrom = line[0]
                bed.write("%s\t%s\t%s\n" % (chrom, exstart-3, exstart+6))
                bed.write("%s\t%s\t%s\n" % (chrom, exstart-21, exstart+2))
                bed.write("%s\t%s\t%s\n" % (chrom, exend-3, exend+6))
                bed.write("%s\t%s\t%s\n" % (chrom, exend-21, exend+2))

def generate_fastafile_frombed(output, ref):
    bedfile = pybedtools.BedTool(os.path.join(output, "splice-site.bed"))
    fasta = pybedtools.BedTool(ref)
    bedfile = bedfile.sequence(fi=fasta)
    return bedfile.seqfn

def main():
    args = parser_args(sys.argv[1:])

    #create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    basepath = args.output

    # identify canonical splicing regions
    generate_splicingbedfile_fromgtf(args.gtf, basepath)
    fasta = generate_fastafile_frombed(basepath, args.reference)
    for j in open(fasta):
        print(j)


    """
    with open(os.path.join(basepath, "five_score_summary"), 'a') as five_file:
        for string in generate_all_kmers(9):
            dna = ''.join(string)
            score = compute_five_score(dna)
            five_file.write("%s\t%s\n" % (dna, score))
    with open(os.path.join(basepath, "three_score_summary"), 'a') as three_file:
        for i in range(300000):
            dna = randomlygenerate_threess()
            score = compute_three_score(dna)
            three_file.write("%s\t%s\n" % (dna, score))
    """

if __name__ == "__main__":
    main()
