[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/novasplice/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/novasplice/badges/downloads.svg)](https://anaconda.org/bioconda/novasplice)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/novasplice/badges/license.svg)](https://anaconda.org/bioconda/novasplice)

## Installation
If you already have `bioconda` installed, simply run:
```
conda install -c bioconda novasplice
```

Afterwards, run:
```
novasplice --help
```

If you don't have `bioconda`, you may run the following commands to get it:
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Find more information [at the bioconda website](https://bioconda.github.io).

## Usage

```
novasplice -v/-vz $VCF_FILE -r/-rz $REFERENCE_FASTA_FILE -b $EXON_BED_FILE -c $CHR_LENS
```

Beyond those required arguments, Novasplice also offers 

```
-l, --libraryname : the name of the output novasplice predictions file. Default is 'novasplice-predictions'
-p, --percent : float representing the lowerbound percent difference Novasplice calls a novel splice region at. i.e. if maxent scores the canonical site as 100, setting a percent of 0.1 means anything above 90 (100-0.1*100) gets called as a novel splice site. Default is 0.05
-o, --output : output directory to pipe results to. If it doesn't exist, create it.
-i, --intermediate : intermediate directory to dump intermediate files. If the folder exists already, use the files within it. Useful for multiple runs with similar VCFs
-t, --temp : temp directory to use. 
```

### Software Requirements
* Python3
* maxentpy
* pybedtools

*Note*: All requirements will automatically be installed via conda

## Software Workflow
The following is how Novasplice works:
1. For every variant, *V*, in the provided VCF:
    1. Identify and score the closest canonical splice site to *V*, denote it *C*
    2. Score(*C*) = *S*
    3. Generate the set of all possible splice sites, *PS*, containing *V*
    4. For every splice site, *P'*, in *PS*:
        1. Score(*P'*) = *S'*
        2. If *S'* > *S* or within some user-defined bound less than *S*:
            * Report possible splice site

## In-Depth Workflow

1. We use the bed file provided in the `-b` option to create a new bed file. This bed file is called `splice-sites.bed`. What is `splice-sites.bed`? `splice-sites.bed` contains the regions corresponding to the 5' and 3' splice sites. Maxentscan, the splice site prediction tool used by NovaSplice (at the time of this writing), requires a 9 bp sequence to score a potential 5' splice site (3 bases in exon + 6 bases in intron) and a 23 bp sequence to score a potential 3' splice site (20 bases in intron + 3 bases in exon). The procedure to generate this file is merely to loop through every exon within the `exon-boundaries.bed` file and write two new bed entries corresponding to that exon's 3/5' splice site. Note that this is dependent upon the directionality of the strand. If a file with this name already exists under the specified intermediate directory, then this step is skipped.

3. Take the given sorted VCF, and run 
```
bedtools intersect -a $VCF -b `exon-boundaries.bed` -v --header --sorted > subset.vcf
```

What does this command do? We generate a new VCF, called `coding-excludedvariants.vcf.gz` that contains the subset of variants that do not fall within an exon. Why is this important? It'll speed up the computation, and NovaSplice is designed to be used for the identification of novel splice sites in non-coding regions of DNA. As I've said before, if this file exists; we SKIP.

4. The next steps are a series of `bedtools closest` commands. This functionality might get more refined as NovaSplice is developed; however, currently it makes the following important assumption: **If a novel splice site occurs within a non-coding region of DNA, that novel splice site will only compete with the closest upstream and downstream canonical splice sites.** While I can make a logical defense of this assumption within my head, this does not appear (to me) to be immediately obvious. Irregardless, here is what I do:
    1. **Case 1** - I have a variant on the + strand, and it induces a novel 5' splice site. In this case, I run `bedtools closest --ignore-upstream`. The result is the location of the closest canonical 5' splice site to my variant of interest. This gets appended to the file `close-down.bed`.
    2. **Case 2** - I have a variant on the + strand, and it induces a novel 3' splice site. In this case, I run `bedtools closest --ignore-downstream`. The result is the location of the closest canonical 3' splice site to my variant of interest. This gets appended to the file `close-up.bed`
    3. **Case 3** - I have a variant on the - strand, and it induces a novel 5' splice site. In this case, I run `bedtools closest --ignore-upstream`. The result is the location of the closest canonical 5' splice site to my variant of interest. This gets appended to the file `close-down.bed`.
    4. **Case 4** - I have a variant on the - strand, and it induces a novel 3' splice site. In this case, I run `bedtools closest --ignore-downstream`. The result is the location of the closest canonical 3' splice site to my variant of interest. This gets appended to the file `close-up.bed`.

So what is the result of all this? I get two files: `close-up.bed` and `close-down.bed`. If everything went well, `close-up.bed` should only contain 5' splice sites, and `close-down.bed` should only contain 3' splice sites. For some reason which I don't fully understand, *occasionally* I see 5' splice sites in my 3' splice site file, and other times I see 3' splice sites in my 5' splice sites.

5. We then generate a bed file called `variant-site-donorsites.bed` and `variant-site-acceptorsites.bed` these files consist of 9/23 entries for every variant in the VCF file and denote all possible splice sites that could have been generated with the variant. 

6. We then derive the reference sequence from these positions and mutate the sequence according to the variants found, these are stored in `variant-site-donorsites.fa` and `variant-site-acceptorsites.fa`

7. We compute 3 scores using maxentscan. The first is the score of the closest canonical splice site to each variant. The next is the score of the reference splice site region before the variant in the VCF was introduced, and the last is the score of the reference splice site region after the variant was introduced. If the canonical score is < 0, we skip. If the score of the reference region before the variant was > canonical score we skip. If the score of the reference region after the variant was > canonical score and > score before the variant, we report it.

8. Final file contains the sorted list of all novel splice sites sorted by the final score of the reference region after the variant.

### Acknowledgments
This tool was developed by [Arya Kaul](aryakaul.github.io) in consultation with the Sunyaev Lab at the Harvard Medical School. 
