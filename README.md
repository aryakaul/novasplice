# NovaSplice
##### Developed by Arya Kaul (last updated on 2018-06-25)
NovaSplice is a tool to determine locations of potential novel splicing given mutagenesis at those locations.


## Requirements
* Python3

## Installation

To use NovaSplice, first ensure that the requirements are correctly installed. Then run:
```
git clone https://github.com/aryakaul/novasplice.git
cd novasplice/scripts
python3 main.py --help
```

If all requirements are correctly installed you should see the following message:
```
usage: main.py [-h] -i FASTAFILE [-o OUTPUTPATH]

Check the help flag

optional arguments:
    -h, --help            show this help message and exit
    -i FASTAFILE, --fastaFile FASTAFILE 
        Direct path to the fasta file being analyzed.
    -o OUTPUTPATH, --outputPath OUTPUTPATH
        Direct path to the desired output folder.
```

## Usage

### NovaSplice
Use NovaSplice by running:

```
python3 main.py -i ./FASTAFILE -o ./OUTPUTFOLDER/
```

The -i flag is required, and if no argument is specified for -o, then it will create a folder in the working directory with name `NovaSplice` and dump output files there. The bulk of the output is directed to stdout, and if that interests you, you should pipe the output of the command to a new file.

If I pass NovaSplice a FASTA file with the following content:
```
>ex0
GGAUUUUUUUUUUUUUUUAGU
```

An example of the first 13 lines of stdout is reproduced below:
```
Index of mutation: 0 Ref->Alt:G->A
        Potential donor sites being analyzed
        Potential branch points being analyzed
                Potential branch point found
        Potential acceptor sites being analyzed
Index of mutation: 0 Ref->Alt:G->C
        Potential donor sites being analyzed
        Potential branch points being analyzed
        Potential acceptor sites being analyzed
Index of mutation: 0 Ref->Alt:G->U
        Potential donor sites being analyzed
        Potential branch points being analyzed
        Potential acceptor sites being analyzed
```
As can be seen, NovaSplice searches every possible mutation event for the three potential sites explained below. The file `ex0.summary` in the created output folder has the following format:

|LOCATION|A|C|G|T|
|--------|-|-|-|-|
|0       |2|2|3|R|

An 'R' would denote the fact that this is the reference at this location. The numbers would correspond to the # of all possible splice sites generated at that position and with that nucleotide substitution. Note that multiple splice sites are possible at a given location and nucleotide since the distance between the donor splice site and the branch point is arbitrary. 

## Description
The following is how NovaSplice's algorithm currently works:
1. Given some arbitrary RNA sequence of length *n*, generate every possible single nucleotide mutagenesis event. This leads to a set of possible mutated sequences of size ``3^{n}``, where each sequence has length *n*. Denote this set *R*
2. For *s* in *R*:
    * Check if mutation in *s* -> donor splice site
    * Check if mutation in *s* -> branch point
    * Check if mutation in *s* -> acceptor splice site
3. Output those sequences and mutation events that lead to a splice site.

A splice site is defined as a stretch of sequences having each of the following requirements:
* A 'GU' sequence
* A possible branch point some arbitrary position upstream of the 'GU'
* A polypyrimidine tail following the potential branch point consisting of >13 pyrimidines and having more than 50% Uracil content
* A 'AG' immediately after the polypyrimidine tail

As of the writing of this. This procedure outputs a binary score (1/0) with no stochastic modeling. This is the next step in NovaSplice's development.
