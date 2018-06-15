import argparse
import sys


def parse_args(args):
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument("-i", "--fastaFile", help="Direct path to the fasta file being analyzed.", required=True)
    parser.add_argument("-o", "--outputFile", help="Direct path to the desired output file.", required=False, default="./NovaSplice-predictions.txt")
    return parser.parse_args()

def readInFasta(fastapath):
    with open(fastapath, 'r') as fastain:
        fastadict = {}
        for lines in fastain:
            line = lines.rstrip()
            if line.startswith(">"):
                currID = line
                fastadict[currID]=""
            else:
                fastadict[currID]+=line
    return fastadict

def generateMutagenesis(dnaseq):
    mutatedNTs = {}
    nts = ['A','C','G','U']
    for ntindex in range(len(dnaseq)):
        for nt in nts:
            deepcopy = dnaseq[:]
            if nt != dnaseq[ntindex]:
                deepcopy[ntindex] = nt
                mutatedNTs[("").join(deepcopy)]=[ntindex, str(dnaseq[ntindex])+"->"+str(nt)]
    return mutatedNTs

def donorScore(dnaseq, index):
    if (index-1)<0:
        postDNA = dnaseq[index:]
    else:
        postDNA = dnaseq[index-1:]
    if len(postDNA) < 5:
        return 0
    if postDNA[0:2] == "GU":
        postpostDNA = postDNA[2:]
        for ntidx in range(len(postpostDNA)):
            if postpostDNA[ntidx] == "A":
                remainingLength = len(postpostDNA)-ntidx
                

def main():
    args = parse_args(sys.argv[1:])
    seqDict = readInFasta(args.fastaFile)
    for seqids in seqDict:
        mutatedNTs = generateMutagenesis(list(seqDict[seqids]))
        for mutations in mutatedNTs:
            donorScore(mutations, mutatedNTs[mutations][0])


if __name__=="__main__":
    main()
