import argparse
import sys
import os

def parse_args(args):
    parser = argparse.ArgumentParser(description="Check the help flag")
    parser.add_argument("-i", "--fastaFile", help="Direct path to the fasta file being analyzed.", required=True)
    parser.add_argument("-o", "--outputPath", help="Direct path to the desired output folder.", required=False, default="./NovaSplice/")
    return parser.parse_args()

def readInFasta(fastapath):
    with open(fastapath, 'r') as fastain:
        fastadict = {}
        for lines in fastain:
            line = lines.rstrip()
            if line.startswith(">"):
                currID = line.split(">")[1]
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
    numberOfSpliceSites = 0
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
                        numberOfSpliceSites+=1
    return numberOfSpliceSites

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
    return 0,i

def branchScore(dnaseq, index):
    numberOfSpliceSites = 0

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
                        numberOfSpliceSites+=1
    return numberOfSpliceSites

def acceptorScore(dnaseq, index):
    numberOfSpliceSites = 0

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
                        numberOfSpliceSites+=1
    return numberOfSpliceSites

def outputSpliceSumm(mutatedNTs, outputpath):
    indexDict = {}
    ntidxlookup = {'A':0, 'C':1, 'G':2, 'U':3}
    for mutations in mutatedNTs:
        totalSplice = sum(mutatedNTs[mutations][2:5])
        index = mutatedNTs[mutations][0]
        ntandref = mutatedNTs[mutations][1]
        newnt = ntandref.split(">")[1]
        ref = ntidxlookup[ntandref.split("-")[0]]
        ntidx = ntidxlookup[newnt]
        if index not in indexDict:
            indexDict[index] = [0,0,0,0]
        indexDict[index][ntidx] += totalSplice
        indexDict[index][ref] = 'R'
    with open(outputpath, 'w') as f:
        f.write("LOCATION\tA\tC\tG\tT\n")
        for i in indexDict:
            f.write("%s\t%s\t%s\t%s\t%s\n" % (i, indexDict[i][0], indexDict[i][1], indexDict[i][2], indexDict[i][3]))

def main():
    args = parse_args(sys.argv[1:])
    seqDict = readInFasta(args.fastaFile)
    if not os.path.isdir(args.outputPath):
        os.makedirs(args.outputPath)
    for seqids in seqDict:
        mutatedNTs = generateMutagenesis(list(seqDict[seqids]))
        for mutations in mutatedNTs:
            print("%s potential mutation being studied" % (mutations))
            print("Potential donor sites being analyzed")
            donorsplice = donorScore(mutations, mutatedNTs[mutations][0])
            mutatedNTs[mutations].append(donorsplice)

            print("Potential branch points being analyzed")
            branchsplice = branchScore(mutations, mutatedNTs[mutations][0])
            mutatedNTs[mutations].append(branchsplice)

            print("Potential acceptor sites being analyzed")
            acceptsplice = acceptorScore(mutations, mutatedNTs[mutations][0])
            mutatedNTs[mutations].append(acceptsplice)
            print("-------------------------\n")
        summaryoutputpath=os.path.join(args.outputPath,seqids+".summary")
        outputSpliceSumm(mutatedNTs, summaryoutputpath)
if __name__=="__main__":
    main()
