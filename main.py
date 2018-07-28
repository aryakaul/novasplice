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
from novasplice.hisat2_extract_exons import extract_exons

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

def extract_exon_boundaries(gtf, output):
    with gzip.open(gtf, 'rt') as g:
        extract_exons(g, output)

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

def generate_variantbedfile_fromvcf(vcf, output):
    with open(vcf, 'r') as vcf:
        with open(os.path.join(output, "variant-site.bed"), 'w') as new:
            for lines in vcf:
                line = lines.rstrip()
                if line.startswith("#"): continue
                line = line.split()
                chrom = line[0]
                pos = int(line[1])
                variant_name = line[2]
                for x in range(9):
                    st = pos-9+x
                    en = pos+x
                    name = variant_name + "-5ss-" + str(x)
                    string = "%s\t%s\t%s\t%s\n" % (chrom, st, en, name)
                    new.write(string)
                for y in range(23):
                    st = pos-23+y
                    en = pos+y
                    name = variant_name + "-3ss-" + str(y)
                    string = "%s\t%s\t%s\t%s\n" % (chrom, st, en, name)
                    new.write(string)

def generate_splicingbed_withexonbound(output):
    with open(os.path.join(output, "exon-boundaries"), 'r') as exons:
        with open(os.path.join(output, "splice-site.bed"), 'w') as bed:
            for lines in exons:
                line = lines.rstrip().split()
                chrom = line[0]
                exst = int(line[1])
                exen = int(line[2])
                direct = line[3]
                bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exst-20, exst+3, direct))
                bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exen-2, exen+7, direct))

def generate_fastafile_frombed(output, ref):
    bedfile = pybedtools.BedTool(os.path.join(output, "splice-site.bed"))
    fasta = pybedtools.BedTool(ref)
    bedfile = bedfile.sequence(fi=fasta, s=True)
    return bedfile.seqfn

def main():
    args = parser_args(sys.argv[1:])

    #create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    basepath = args.output

    # extract and score reference fasta from file
    extract_exon_boundaries(args.gtf, args.output)
    generate_splicingbed_withexonbound(args.output)
    fasta = generate_fastafile_frombed(args.output, args.reference)
    generate_variantbedfile_fromvcf(args.vcf, args.output)
    sys.exit(2)
     

    x = {}
    for j in open(fasta):
        j = j.rstrip()
        if len(j) == 0 or j.startswith(">"): print(j)
        elif len(j) == 9:
            fivescore = compute_five_score(j)
            print("%s\t%s" % (j, fivescore))
            if j[3:5] not in x: x[j[3:5]] = [0,0]
            x[j[3:5]][0] += 1
        elif len(j) == 23:
            threescore = compute_three_score(j)
            print("%s\t%s" % (j, threescore))
            if j[18:20] not in x: x[j[18:20]] = [0,0]
            x[j[18:20]][1] += 1
        else:
            print("ERROR")
            sys.exit(2)
    sys.exit(2)


if __name__ == "__main__":
    main()
