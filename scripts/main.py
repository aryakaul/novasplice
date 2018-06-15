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

def predictSplice(dnaseq):
    


def main():
    args = parse_args(sys.argv[1:])
    seqDict = readInFasta(args.fastaFile)
    for ids in seqDict:
        storage = predictSplice(seqDict[ids])

if __name__=="__main__":
    main()
