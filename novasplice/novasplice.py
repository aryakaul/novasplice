#!/usr/bin/env python
import time
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

global mat5
global mat3

mat5 = load_matrix5()
mat3 = load_matrix3()

def parser_args(args):
    parser = argparse.ArgumentParser(prog="novasplice")
    parser.add_argument('-v', '--vcf', help="Full path to the sorted vcf file being used", type=str, required=False)
    parser.add_argument('-vz', '--zippedvcf', help="Full path to the sorted zipped vcf file being used", type=str, required=False)
    parser.add_argument('-r', '--reference', help="Full path to the reference genome being used", type=str, required=False)
    parser.add_argument('-rz', '--zippedreference', help="Full path to the zipped reference genome being used", type=str, required=False)
    parser.add_argument('-b', '--bed', help="Full path to the reference exon boundary bed file being used", type=str, required=True)
    parser.add_argument('-c', '--chrlens', help="Full path to the chromosome length file being used", type=str, required=True)
    parser.add_argument('-p', '--percent', help="Lower bound percent to call novel splice site", type=float, required=False, default=0.05)
    parser.add_argument('-o', '--output', help="Path to the output folder to dump simdigree's output to. Default is working directory under /novasplice_output", type=str, required=False, default="./novasplice_output")
    parser.add_argument('-i', '--intermediate', help="Path to output folder that will hold intermediate files generated, not specific to the provided vcf. Especially useful when running NovaSplice on a large number of VCFs that all come from the same reference and make use of the same --bed option.", required=False, default="./novasplice-temp")
    parser.add_argument('-t', '--temp', help="Full path to an alternative directory to use for temp files. Default is /tmp", type=str, required=False)
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

def generate_variantbedfile_fromclosest(closest, bed, donor):
    with open(closest, 'r') as closebed:
        with open(bed, 'w') as new:
            printed = False
            for lines in closebed:
                line = lines.rstrip()
                if line.startswith("#"): continue
                line = line.split()
                chrom = line[0]
                pos = int(line[1])
                variant_name = line[0]+'/'+line[1]+'/'+line[2]+","+line[-7]+"_"+line[-6]+"_"+line[-5]
                ref = line[3]
                alt = line[4]
                if len(alt) > 1 or len(ref) > 1:
                    if not printed:
                        print("VCF file contains Non-SNPs. NovaSplice will continue, but it'll be slower because of these")
                        printed = True
                    continue
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

def generate_variantbedfile_fromvcf(vcf, bed, donor, zipped):
    if zipped:
        vcf = gzip.open(vcf, 'rt')
    else:
        vcf = open(vcf, 'r')
    with open(bed, 'w') as new:
        printed = False
        for lines in vcf:
            line = lines.rstrip()
            if line.startswith("#"): continue
            line = line.split()
            chrom = line[0]
            pos = int(line[1])
            variant_name = line[0]+'/'+line[1]+'/'+line[2]
            ref = line[3]
            alt = line[4]
            if len(alt) > 1 or len(ref) > 1:
                if not printed:
                    print("VCF file contains Non-SNPs. NovaSplice will continue, but it'll be slower because of these")
                    printed = True
                continue
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
    vcf.close()

def filter_snps_only(feature):
    if len(feature[3])==1 and len(feature[4])==1:
        return True
    return False

def generate_splicingbed_withexonbound(ss, chrlens, bed):
    with open(bed, 'r') as exons:
        with open(ss, 'w') as bed:
            with open(chrlens, 'r') as lens:
                file_contents = lens.read().split()[2:]
                lensdict = {}
                for ctr,content in enumerate(file_contents):
                    if ctr % 2 == 0:
                        chrom = content
                        lensdict[chrom] = 0
                    else:
                        lensdict[chrom] += int(content)
                for lines in exons:
                    line = lines.rstrip().split()
                    chrom = line[0]
                    exst = int(line[1])
                    exen = int(line[2])
                    if exst-20 < 0 or exen-3 < 0: continue
                    if exst+3 > lensdict[chrom] or exen+20 > lensdict[chrom]: continue
                    direct = line[5]
                    if direct == "+":
                        bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exst-20, exst+3, direct))
                        bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exen-3, exen+6, direct))
                    elif direct == "-":
                        bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exst-6, exst+3, direct))
                        bed.write("%s\t%s\t%s\t.\t.\t%s\n" % (chrom, exen-3, exen+20, direct))
    a = pybedtools.BedTool(ss)
    a = a.sort().moveto(ss)

def generate_fastafile_frombed(ref, bed):
    bedfile = pybedtools.BedTool(bed)
    fasta = pybedtools.BedTool(ref)
    bedfile = bedfile.sequence(fi=fasta, s=True, name=True)
    return bedfile.seqfn

def writevcfheader(path, vcf, zipped):
    if zipped:
        vcf_file = gzip.open(vcf,'rt')
    else:
        vcf_file = open(vcf, 'r')

    with open(path, 'w') as out:
        for lines in vcf_file:
            if lines.startswith('#'):
                out.write(lines)
            else:
                return

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
                    if new_dna[index].lower() != reference_allele.lower():
                        print("ERROR")
                        print(new_dna)
                        print(index)
                        print(reference_allele)
                        print(alternate_allele)
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
                name = line[0]+"/"+line[1]+"/"+line[2]+","
                diff = line[-7:]
                if diff[0] == ".": continue
                name += diff[0]+"_"+diff[1]+"_"+diff[2]
                n.write(diff[0]+"\t"+diff[1]+"\t"+diff[2]+"\t"+name+"\n")
    fasta = pybedtools.BedTool(ref)
    bedfile = pybedtools.BedTool(os.path.join(output, 'intermedbed.bed'))
    bedfile = bedfile.sequence(fi=fasta, s=True, name=True)
    with open(bedfile.seqfn, 'r') as f:
        for lines in f:
            if lines.startswith(">"):
                line = lines.split(">")[1]
                canonical_name = line.split("(")[0].split(",")[1]
            else:
                if len(lines) == 0: continue
                seq = lines.rstrip()
                if "N" in seq or "n" in seq: continue
                if donor:
                    #TODO why do I need this check?
                    if len(seq) != 9: continue
                    score = compute_five_score(seq)
                else:
                    #TODO why do I need this check?
                    if len(seq) != 23: continue
                    score = compute_three_score(seq)
                nametoscore[canonical_name] = score
                name = ""
    return nametoscore

def generate_canonicalss_score(fastaname, novelsplice, donor):
    alternate_allele = fastaname.split(">")[1].split("(")[0]
    reference_allele = fastaname.split("-")[3]
    loc = len(novelsplice) - int(fastaname.split("-")[2]) - 1
    if alternate_allele != novelsplice[loc]: print("ERROR")
    canonsplice = list(novelsplice)
    canonsplice[loc] = reference_allele
    canonsplice = "".join(canonsplice)
    if donor: score = compute_five_score(canonsplice)
    else: score = compute_three_score(canonsplice)
    return score


def compare_scores(variantsitesfasta, canonicalscoredict, percent, output, donor, lib):
    sorted_dict = {}
    with open(variantsitesfasta, 'r') as f:
        for lines in f:
            line = lines.rstrip()
            if len(line) == 0: continue
            if line.startswith(">"):
                name = line.split(">")[1].split("-")[0]
                canonicalss_loc = line.split(">")[1].split("-")[0].split(",")[1]
                sortedname = line[1:]
                try:
                    canonicalscore = canonicalscoredict[canonicalss_loc]
                except:
                    if donor: word = "upstream"
                    else: word = "downstream"
                    print("%s does not have a closest %s canonical splice site" % (name, word))
                    canonicalscore = -1
            else:
                novelss = line
                if "N" in novelss or "n" in novelss: continue
                if canonicalscore < 0: continue
                canonscore = generate_canonicalss_score(sortedname, novelss, donor)
                if canonscore > canonicalscore: continue
                if donor: novscore = compute_five_score(novelss)
                else: novscore = compute_three_score(novelss)
                lb = canonicalscore-(percent*canonicalscore)
                if novscore>canonicalscore or novscore>lb:
                    print("Potential splice site found!")
                    print("Novel SS: %s\t Novel Score: %s\t Score before variant: %s\tCanonical Score: %s" % (novelss, novscore, canonscore, canonicalscore))
                    sorted_dict[sortedname] = (novelss, novscore, canonscore, canonicalscore, canonicalss_loc)
    with open(os.path.join(output, lib), 'a') as output:
        for key,value in sorted(sorted_dict.items(), key=lambda kv:kv[1][1], reverse=True):
            output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (key, value[0], value[1], value[2], value[3], value[4]))

def merge_two_files(bigvcf, vcfheader, newoutput):
    with gzip.open(newoutput, 'wt') as newfile:
        with open(vcfheader,'r') as header:
            headercontent = header.read()
            newfile.write(headercontent)
        os.remove(vcfheader)
        for lines in bigvcf:
            newfile.write(str(lines).rstrip()+'\n')

def main():
    args = parser_args(sys.argv[1:])

    if args.temp:
        pybedtools.helpers.set_tempdir(args.temp)

    if not args.vcf and not args.zippedvcf:
        print("ERROR. VCF required, please use -v or -vz")
        if not args.zippedvcf.endswith(".gz"):
            print("ERROR. --vz used with non-gzipped file (must end with '.gz')")
        if args.vcf.endswith(".gz"):
            print("ERROR. --vz option should be used, gzipped file detected")
        sys.exit(2)
    if not args.reference and not args.zippedreference:
        print("ERROR. Reference fasta required, please use -r or -rz")
        sys.exit(2)

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if not os.path.exists(args.intermediate):
        os.makedirs(args.intermediate)

    # Go through the exon-boundaries bed file, and generate a new file
    #  called $INTERMED/splice-site.bed. This bed file contains the
    #  organism's canonical splice site coordinates.
    canonical_splicesites_bedfile = os.path.join(args.intermediate, "splice-site.bed")
    if os.path.exists(canonical_splicesites_bedfile):
        print("Bed file containing splice sites is already generated. Moving on!")
    else:
        print("Generating a splicing bed file using exon boundaries...")
        start = time.time()
        generate_splicingbed_withexonbound(canonical_splicesites_bedfile, args.chrlens, args.bed)
        end = time.time()
        print("Finished generating. Time took %s" % (end-start))

    # Find the subset of the VCF which contains non-coding variants
    #    use this VCF for NovaSplice calculations.
    print("Subsetting VCF to noncoding variants and SNPs")
    start = time.time()
    pybedtools.cleanup()
    vcfheader_file=os.path.join(args.output, "vcf-header")
    if args.vcf:
        writevcfheader(vcfheader_file, args.vcf, False)
        vcf = pybedtools.BedTool(args.vcf)
    else:
        writevcfheader(vcfheader_file, args.zippedvcf, True)
        vcf = pybedtools.BedTool(args.zippedvcf)
    exon_bounds = pybedtools.BedTool(args.bed).sort()
    subset_vcf_location = os.path.join(args.output, "coding-excludedvariants.vcf.gz")
    subsetvcf = vcf.intersect(exon_bounds, v=True, sorted=True).filter(filter_snps_only)
    merge_two_files(subsetvcf, vcfheader_file, subset_vcf_location)
    end = time.time()
    print("Finished subsetting. Time took %s" % (end-start))

    # We now find the closest upstream/downstream canonical splice sites
    #    non-coding variant in the subsetted VCF
    print("Finding closest canonical splice sites to each non-coding variant")
    start = time.time()
    pybedtools.cleanup()
    subset_vcf = pybedtools.BedTool(subset_vcf_location)
    canon_bed = pybedtools.BedTool(canonical_splicesites_bedfile)
    subset_vcf.closest(canon_bed, D="b", id=True, io=True, output=os.path.join(args.output,"close-up.bed"))
    subset_vcf.closest(canon_bed, D="b", iu=True, io=True, output=os.path.join(args.output, "close-down.bed"))
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    # DONOR SPECIFIC
    # For every variant, compute the set of 9 possible donor sites with that variant. The
    #    output is a bed file that has 9 entries for every variant.
    variantsite_location = os.path.join(args.output,"variant-site-donorsites.bed")
    print("Generating a variant bed file from vcf...")
    start = time.time()
    #generate_variantbedfile_fromvcf(subset_vcf_location, variantsite_location, True, True)
    generate_variantbedfile_fromclosest(os.path.join(args.output,"close-down.bed"), variantsite_location, True)
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    # Extract the sequences of every variant and store them in a file. Note that this
    #    file is not saved anywhere except $TMP and is solely based on the reference
    print("Generating a fasta file from variant bed file...")
    start = time.time()
    if args.reference:
        fastaref = generate_fastafile_frombed(args.reference, variantsite_location)
    else:
        fastaref = generate_fastafile_frombed(args.zippedreference, variantsite_location)
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    # Mutate the fasta file with the variants found in the vcf file
    print("Mutating fasta file per SNPs in VCF...")
    start = time.time()
    generate_variantfastafile(fastaref, os.path.join(args.output, "variant-site-donorsites.fa"))
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    # Create a dictionary that relates a given variant to the closest upstream canonical 
    #   splice site's score
    print("Scoring canonical and novel splice-sites...")
    start = time.time()
    if args.reference:
        nametoscore = extract_canonical_score(os.path.join(args.output, 'close-down.bed'), args.output, args.reference, True)
    else:
        nametoscore = extract_canonical_score(os.path.join(args.output, 'close-down.bed'), args.output, args.zippedreference, True)
    with open(os.path.join(args.output, args.libraryname), 'w') as out:
        out.write("Novel SS\tNovel Score\tScore before variant\tClosest Canonical Score\tLocation of closest canonical ss\n")
    compare_scores(os.path.join(args.output, "variant-site-donorsites.fa"), nametoscore, args.percent, args.output, True, args.libraryname)
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))
    print("Finished analysis for donor sites. Now doing same for acceptor sites.")

    #ACCEPTOR SITES


    variantsite_location = os.path.join(args.output,"variant-site-acceptorsites.bed")
    print("Generating a variant bed file from vcf...")
    start = time.time()
    #generate_variantbedfile_fromvcf(subset_vcf_location, variantsite_location, False, True)
    generate_variantbedfile_fromclosest(os.path.join(args.output,"close-up.bed"), variantsite_location, False)
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    print("Generating a fasta file from variant bed file...")
    start = time.time()
    if args.reference:
        fastaref = generate_fastafile_frombed(args.reference, variantsite_location)
    else:
        fastaref = generate_fastafile_frombed(args.zippedreference, variantsite_location)
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    print("Mutating fasta file per SNPs in VCF...")
    start = time.time()
    generate_variantfastafile(fastaref, os.path.join(args.output, "variant-site-acceptorsites.fa"))
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    print("Scoring canonical and novel splice-sites...")
    start = time.time()
    if args.reference:
        nametoscore = extract_canonical_score(os.path.join(args.output, 'close-up.bed'), args.output, args.reference, False)
    else:
        nametoscore = extract_canonical_score(os.path.join(args.output, 'close-up.bed'), args.output, args.zippedreference, False)
    compare_scores(os.path.join(args.output, "variant-site-acceptorsites.fa"), nametoscore, args.percent, args.output, False, args.libraryname)
    end = time.time()
    print("Finished generating. Time took %s" % (end-start))

    print("Finished! Removing intermediate files now...")

    toRem = []
    toRem += (glob.glob(os.path.join(args.output, '*.bed')))
    toRem += (glob.glob(os.path.join(args.output, '*.vcf.gz')))
    toRem += (glob.glob(os.path.join(args.output, '*.fa')))
    for i in (toRem):
        os.remove(i)

if __name__ == "__main__":
    main()
