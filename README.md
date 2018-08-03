# NovaSplice
##### Developed by Arya Kaul (last updated on 2018-08-02)
NovaSplice is a tool to determine locations of potential novel splicing given mutagenesis at those locations.


## Requirements
* Python3
* maxentpy
* pybedtools

## Usage

```
novasplice -v $VCF_FILE -r $REFERENCE_FASTA_FILE -g $REFERENCE_GTF_FILE 
```

Beyond those required arguments, Novasplice also offers 

```
-l, --libraryname : the name of the output novasplice predictions file. Default is 'novasplice-predictions'
-p, --percent : float representing the lowerbound percent difference Novasplice calls a novel splice region at. i.e. if maxent scores the canonical site as 100, setting a percent of 0.1 means anything above 90 (100-0.1*100) gets called as a novel splice site. Default is 0.05
-o, --output : output directory to pipe results to
```

## Description
The following is how Novasplice works:
1. For every variant, *V*, in the provided VCF:
    1. Identify and score the closest canonical splice site to *V*, denote it *C*
    2. Score(*C*) = *S*
    3. Generate the set of all possible splice sites, *PS*, containing *V*
    4. For every splice site, *P'*, in *PS*:
        1. Score(*P'*) = *S'*
        2. If *S'* > *S* or within some user-defined bound less than *S*:
            * Report possible splice site
