import argparse
import sys
import random
import os

def parse_args(args):
    parser = argparse.ArgumentParser(description="Check help flag")
    parser.add_argument("-l", "--length", type=int, required=True, help= "Length of random sequences to generate")
    parser.add_argument("-n", "--number", type=int, required=True, help="Number of random sequences to generate")
    parser.add_argument("-f", "--fastafile", required=False, action='store_true', help="Print fasta file")
    parser.add_argument("-r", "--rna", required=False, action='store_true', help="Use RNA instead of DNA")
    return parser.parse_args()

def generateSeq(length, noOfSeqs=1, fastafile=False, rna=False):
    generated = []
    for seqs in range(noOfSeqs):
        seq = ""
        for nt in range(length):
            seq += rollDie(rna)
        generated.append(seq)
    if fastafile:
        seqctr = 0
        with open("./random_len-"+str(length)+"_no-"+str(noOfSeqs)+".fasta", 'w') as filein:
            for seq in generated:
                filein.write(">p%s\n" % (seqctr))
                filein.write(seq)
            seqctr += 1
    return generated

def rollDie(rna):
    x = random.random()
    dice = [0.25, 0.5, 0.75, 1]
    sides = ["A", "C", "G", "T"]
    if rna:
        sides = ["A", "C", "G", "U"]
    if 0<x<=dice[0]:
        return sides[0]
    for i in range(1, len(dice)):
        if dice[i-1]<x<=dice[i]:
            return sides[i]

def outputGenerated(generated, rna=False):
    ntctr = {"A":0,"G":0,"C":0,"T":0}
    if rna:
        ntctr = {"A":0, "G":0, "C":0, "U":0}
    totctr = 0
    for i in generated:
        for nt in i:
            ntctr[nt]+=1
            totctr+=1
        print(i)
    print("Observed nucleotide frequencies in set:")
    print("Adenine: %s out of %s (%s%%)" % (ntctr["A"], totctr, (100*float(ntctr["A"]/totctr))))
    print("Guanine: %s out of %s (%s%%)" % (ntctr["G"], totctr, (100*float(ntctr["G"]/totctr))))
    print("Cytosine: %s out of %s (%s%%)" % (ntctr["C"], totctr, (100*float(ntctr["C"]/totctr))))
    if rna:
        print("Uracil: %s out of %s (%s%%)" % (ntctr["U"], totctr, (100*float(ntctr["U"]/totctr))))
    else:
        print("Thymine: %s out of %s (%s%%)" % (ntctr["T"], totctr, (100*float(ntctr["T"]/totctr))))

def main():
    args = parse_args(sys.argv[1:])
    gen = generateSeq(args.length, args.number, args.fastafile, args.rna)
    outputGenerated(gen, args.rna)

if __name__ == "__main__":
    main()
