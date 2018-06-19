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
    #extract 2 bp around mutation
    if (index-1)<0:
        postDNA = dnaseq[index:]
    else:
        postDNA = dnaseq[index-1:]
    
    #discard edge case of strands being too short
    if len(postDNA) < 5:
        return 0

    #check if mutation leads to donor splice site
    if postDNA[0:2] == "GU":
        print("Potential donor site found at index %s" % index)
        postpostDNA = postDNA[2:]

        #look for branch point
        for ntidx in range(len(postpostDNA)):
            if postpostDNA[ntidx] == "A":
                print("Potential branch point found at index %s" % ntidx)
                branchpointpost = postpostDNA[ntidx+1:]
                
                #look for polypyr tail
                score, end = computePolyPyrProb(branchpointpost)
                
                if score == 1:
                    print("Potential poly-pyrimidine tail found at index %s" % end)
                    polypyrpost = branchpointpost[end:]
                    
                    #look for acceptor splice site
                    if polypyrpost[0:2] == "AG":
                        print("Potential acceptor site found")
                        print("Splice site found!")
                        return 1
    return 0

def computePolyPyrProb(polypyr):
    if len(polypyr)==0:
        return 0,0
    uracilctr=0
    pyrctr=0
    for i in range(len(polypyr)):
        if polypyr[i] == "C" or polypyr[i]=="U":
            pyrctr+=1
            if polypyr[i]=="U":
                uracilctr+=1
        else:
            if uracilctr != 0 and pyrctr > 13:
                return 1,i
            else:
                return 0,i

def branchScore(dnaseq, index):

    #check if it can be a branch point
    if dnaseq[index]=="A":
        print("Potential branch point found")
        score, end = computePolyPyrProb(dnaseq[index+1:])
        
        #check polypyr tail exists
        if score == 1:
            polypyrpost = dnaseq[index+1:][end:]
            
            #check if acceptor site exists
            if polypyrpost[0:2] == "AG":
                print("Potential acceptor site found")
                prebranch = dnaseq[:index]
                
                #search for donor site
                for i in range(len(prebranch)-1):
                    if prebranch[i:i+2]=="GU":
                        print("Potential donor site found")
                        print("Splice site found!")
                        return 1 
    return 0

def acceptorScore(dnaseq, index):
    
    #find sequence before potential acceptor splice site
    extractPrevDNA = dnaseq[:index+1]
    if extractPrevDNA[len(extractPrevDNA)-2:len(extractPrevDNA)] == "AG":
        print("Potential acceptor splice site found")
        polyPyrTailPlus = extractPrevDNA[:len(extractPrevDNA)-2][::-1]
        score, end = computePolyPyrProb(polyPyrTailPlus)
        
        #check for polypyr tail
        if score == 1:
            print("Potential poly-pyrimidine tail found at index %s" % end)
            polypyrpre = polyPyrTailPlus[end:]
            
            #check for branch point
            if polypyrpre[0] == "A":
                print("Potential branch point found")
                branchpre = polypyrpre[1:]
                
                #check for acceptor site
                for i in range(len(branchpre)-1):
                    if branchpre[i:i+2] == "UG":
                        print("Potential donor site found")
                        print("Splice site found!")
                        return 1
    return 0

def main():
    print(branchScore("GUAUUUUUUUUUUUUUUUUUAG", 2))
    sys.exit(2) 
    args = parse_args(sys.argv[1:])
    seqDict = readInFasta(args.fastaFile)
    for seqids in seqDict:
        mutatedNTs = generateMutagenesis(list(seqDict[seqids]))
        for mutations in mutatedNTs:
            print("%s potential mutation being studied" % (mutations))
            print("Potential donor sites being analyzed")
            score = donorScore(mutations, mutatedNTs[mutations][0])
            mutatedNTs[mutations].append(score) 
            print("Potential branch points being analyzed")
            score = donorScore(mutations, mutatedNTs[mutations][0])
            print("-------------------------\n")

if __name__=="__main__":
    main()
