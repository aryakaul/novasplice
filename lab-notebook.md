# NovaSplice
#### Arya Kaul

#### 2018-06-15
Google drive folder has been created with relevant papers. [Link is here](https://drive.google.com/drive/folders/14_M6-YHQTlVhZLqZwtEI6FyJ9Phcl88U?usp=sharing).

In aggregate, they suggest that the best way to identify intronic regions that could present themselves as novel splicing regions is to identify those regions that contain the closest relationship to donor/acceptor splice sites and then artificially mutate them. For every mutation event, we could measure the likelihood of it being a splice site using existing methods and assign a probability of this mutation causing a novel splicing event. Unfortunately, this method would be computationally burdensome and its complexity would increase if we wanted to consider >1 mutagenesis events.

That being said, is there a more elegant way to do this? I don't (currently) see one.

Let's try to implement the naive solution first. This code is found in scripts/main.py.

Goal?
Given an arbitrary FASTA file, let's assign possible mutation events to each nucleotide and its associated potential for causing a splicing event. 

To this end, I have written a function to parse and read into memory an arbitrary FASTA file. 
I then wrote a function that given a dna sequence, generates every possible *single* nucleotide mutagenesis event and outputs a list of them. For every generated nucleotide mutagenesis event, we associate the location of mutagenesis and the REF->ALT nucleotide change.

Now, I need to build a function that can take a nucleotide sequence, and an associated index. And assign a score (likelihood function) to that nucleotide index to the probability of that index being involved in a splicing event.

**UPDATE** - I think it would be more simple to split the initial function into 3 separate scores. One score for the mutation leading to a donor splice site, one score for the mutation leading to a branch point adenine, and one score for the mutation leading to a acceptor splice site.

#### 2018-06-18
I've successfully completed the `donorScore` function, and cleaned up some logic in the beginning. I've also added more thorough commenting throughout.

