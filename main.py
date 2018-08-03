import glob
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
    parser.add_argument('-p', '--percent', help="Lower bound percent to call novel splice site", type=float, required=False, default=0.05)
    parser.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is working directory under /novasplice_output", type=str, required=False, default="./novasplice_output")
    parser.add_argument('-l', '--libraryname', help="Name of the final file novasplice outputs with predictions", type=str, required=False, default="novasplice_predictions")
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

def generate_variantbedfile_fromvcf(vcf, output, donor):
    if donor:
        name = os.path.join(output, "variant-site-donorsites.bed")
    else:
        name = os.path.join(output, "variant-site-acceptorsites.bed")
    with open(vcf, 'r') as vcf:
        with open(name, 'w') as new:
            for lines in vcf:
                line = lines.rstrip()
                if line.startswith("#"): continue
                line = line.split()
                chrom = line[0]
                pos = int(line[1])
                variant_name = line[0]+'/'+line[1]+'/'+line[2]
                ref = line[3]
                alt = line[4]
                if donor:
                    for x in range(9):
                        st = pos-9+x
                        en = pos+x
                        name = variant_name + "-5ss-" + str(x) + "-" + ref + "->" + alt
                        string = "%s\t%s\t%s\t%s\n" % (chrom, st, en, name)
                        new.write(string)
                else:
                    for y in range(23):
                        st = pos-23+y
                        en = pos+y
                        name = variant_name + "-3ss-" + str(y) + "-" + ref + "->" + alt
                        string = "%s\t%s\t%s\t%s\n" % (chrom, st, en, name)
                        new.write(string)

def generate_splicingbed_withexonbound(output):
    with open(os.path.join(output, "exon-boundaries.bed"), 'r') as exons:
        with open(os.path.join(output, "splice-site.bed"), 'w') as bed:
            for lines in exons:
                line = lines.rstrip().split()
                chrom = line[0]
                exst = int(line[1])
                exen = int(line[2])
                direct = line[3]
                bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exst-20, exst+3, direct))
                bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exen-2, exen+7, direct))
    a = pybedtools.BedTool(os.path.join(output, "splice-site.bed"))
    b = a.sort()
    return b

def generate_fastafile_frombed(ref, bed):
    bedfile = pybedtools.BedTool(bed)
    fasta = pybedtools.BedTool(ref)
    bedfile = bedfile.sequence(fi=fasta, s=True, name=True)
    bedfile = bedfile.sequence(fi=fasta, name=True)
    return bedfile.seqfn

def generate_variantfastafile(fasta, output):
    with open(fasta, 'r') as f:
        with open(output, 'w') as g:
            for lines in f:
                line = lines.rstrip().split()
                if len(line) == 0: continue
                if line[0].startswith(">"):
                    ref = line[0].split('>')
                    loc = int((ref[1].split('-'))[2])
                    reference_allele = (ref[1].split('-'))[3]
                    alternate_allele = ref[2].split('(')[0]
                    g.write(''.join(line)+"\n")
                else:
                    old_dna = line[0]
                    new_dna = list(old_dna)
                    index = len(old_dna)-1-loc
                    if new_dna[index] != reference_allele:
                        print("ERROR")
                        print(new_dna)
                        print(index)
                        print(loc)
                        sys.exit(2)
                    new_dna[index] = alternate_allele
                    new_dna = "".join(new_dna)
                    g.write(new_dna+"\n")

def extract_canonical_score(inputvcf, output, ref, donor):
    nametoscore = {}
    with open(inputvcf, 'r') as iv:
        with open(os.path.join(output, "intermedbed.bed"), 'w') as n:
            for lines in iv:
                line = lines.rstrip().split()
                name = line[0]+"/"+line[1]+"/"+line[2]
                #nametoscore[name] = -50
                n.write(line[0]+"\t"+line[11]+"\t"+line[12]+"\t"+name+"\n")
    fasta = pybedtools.BedTool(ref)
    bedfile = pybedtools.BedTool(os.path.join(output, 'intermedbed.bed'))
    bedfile = bedfile.sequence(fi=fasta, s=True, name=True)
    with open(bedfile.seqfn, 'r') as f:
        for lines in f:
            if lines.startswith(">"):
                line = lines.split(">")[1]
                name = line.split("(")[0]
            else:
                if len(lines) == 0: continue
                seq = lines.rstrip()
                if donor: score = compute_five_score(seq)
                else: score = compute_three_score(seq)
                nametoscore[name] = score
                name = ""
    return nametoscore

def compare_scores(variantsitesfasta, canonicalscoredict, percent, output, donor, lib):
    sorted_dict = {}
    with open(variantsitesfasta, 'r') as f:
        for lines in f:
            line = lines.rstrip()
            if len(line) == 0: continue
            if line.startswith(">"):
                name = line.split(">")[1].split("-")[0]
                sortedname = line.split(">")[1]
                canonicalscore = canonicalscoredict[name]
            else:
                novelss = line
                if donor: novscore = compute_five_score(novelss)
                else: novscore = compute_three_score(novelss)
                lb = canonicalscore-(percent*canonicalscore)
                if novscore>canonicalscore or novscore>lb:
                    print("Potential splice site found!")
                    print("Novel SS: %s\t Score: %s\t Canonical Score: %s" % (novelss, novscore, canonicalscore))
                    sorted_dict[sortedname] = (novelss, novscore, canonicalscore)
    with open(os.path.join(output, lib), 'a') as output:
        for key,value in sorted(sorted_dict.items(), key=lambda kv:kv[1][1], reverse=True):
            output.write("%s\t%s\t%s\t%s\n" % (key, value[0], value[1], value[2]))


def main():
    args = parser_args(sys.argv[1:])

    #create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    basepath = args.output

    # generate a bed file of each canonical splice site
    extract_exon_boundaries(args.gtf, os.path.join(args.output, "exon-boundaries.bed"))
    canonicalsplicingbed = generate_splicingbed_withexonbound(args.output)

    # for every variant in vcf, find the closest upstream and downstream canonical splice site associated with the variant
    vcf = pybedtools.BedTool(args.vcf)
    vcf.sort()
    closest_upstream = vcf.closest(canonicalsplicingbed, D="ref", id=True).moveto(os.path.join(args.output, "close-up.bed"))
    closest_downstream = vcf.closest(canonicalsplicingbed, D="ref", iu=True).moveto(os.path.join(args.output, "close-down.bed"))

    #for every variant, compute the set of 9 possible donor sites with that variant
    generate_variantbedfile_fromvcf(args.vcf, args.output, True)
    fastaref = generate_fastafile_frombed(args.reference, os.path.join(args.output, "variant-site-donorsites.bed"))
    generate_variantfastafile(fastaref, os.path.join(args.output, "variant-site-donorsites.fa"))
    nametoscore = extract_canonical_score(os.path.join(args.output, 'close-up.bed'), args.output, args.reference, True)
    compare_scores(os.path.join(args.output, "variant-site-donorsites.fa"), nametoscore, args.percent, args.output, True, args.libraryname)


    #for every variant, compute the set of 23 possible acceptor sites with that variant
    generate_variantbedfile_fromvcf(args.vcf, args.output, False)
    fastaref = generate_fastafile_frombed(args.reference, os.path.join(args.output, "variant-site-acceptorsites.bed"))
    generate_variantfastafile(fastaref, os.path.join(args.output, "variant-site-acceptorsites.fa"))
    nametoscore = extract_canonical_score(os.path.join(args.output, 'close-down.bed'), args.output, args.reference, False)
    compare_scores(os.path.join(args.output, "variant-site-acceptorsites.fa"), nametoscore, args.percent, args.output, False, args.libraryname)

    toRem = []
    toRem += (glob.glob(os.path.join(args.output, '*.bed')))
    toRem += (glob.glob(os.path.join(args.output, '*.fa')))
    for i in (toRem):
        os.remove(i)

if __name__ == "__main__":
    main()
