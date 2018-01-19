# vdjRec version 0.1
TCR amino acid sequence generation and recombination probability estimation in R. Some R skills are required to use advanced features.  

## About
This software allows to calculate recombination probability of TCR sequence, simulate TCR repertoire, and identify candidate condition-associated T-cell receptor sequences using only patient cohort. For details of the approach, see [preprint](https://www.biorxiv.org/content/early/2017/09/27/195057).  

## Software requirements
Any OS where R is available (Linux, OS X, Windows), however parallel computing is currently not available on Windows.  

## Installation

1. Install R distribution of choice (i.e. from [R Core team](https://cloud.r-project.org/), or [Microsoft R Open](https://mran.microsoft.com/open/) )
2. Install BioStrings package from bioconductor: open R console and execute following commands: 
```R
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```
3. Download this repository and unzip demo.zip.


## Quick start
#### Sequence generation

Currently only TCRB generation is supported (more loci would be available soon!). 

User need to specify V, J and desired number of recombination events. 

*Note*: _translate_ and _inframe_only_ options discard out-of-frame sequences from the output.

```R
source("generation.R")
#Generate both inframe and out-of-frame nucleotide sequences:
#all 100 recombination events are in the output 
gen_beta(100,V="TRBV9",J="TRBJ2-1",translate=F,inframe_only=F) 

#Lets generate some in frame TCR beta amino acid sequences:
#Only in frame sequences are in the output (appox. third of recombination events), output is translated
gen_beta(100,V="TRBV9",J="TRBJ2-1",translate=T,inframe_only=T)

#Lets generate some in frame TCR beta nucleotide sequences:
#Only in frame sequences are in the output (appox. third of recombination events), output is NOT translated
gen_beta(100,V="TRBV9",J="TRBJ2-1",translate=F,inframe_only=T)
```
Immunoseq format uses different V and J segment naming than IMGT, use supplied conversion tables if V and J names are in immunoseq format: 
```R
gen_beta(100,V=Vv["TCRBV09-01"],J=Jv["TCRBJ02-01"],translate=F,inframe_only=F)
```

#### Generative probability estimation
How to estimate generative probability for CDR3 amino acid sequence? 
Generate a lot of TCR aminoacid sequences and count, how many of them correspond to our sequences of interest!

To save memory, sequence generation is performed in _iter_ batches of _nrec_ size, and may be performed in parallel on several _cores_, so total number of simulated sequences is _nrec_ \*_iter_\*_cores_
```R
source("generation.R")
#load sample file with TCRB CDR3 amino acid sequences
CDR3s<-read.csv2("demo/TCRBV05-01_TCRBJ01-01.csv")
#estimate counts
CDR3s_p<-estimate_pgen_aa(CDR3s,iter=3,nrec=5e5,V="TRBV5-1",J="TRBJ1-1",colname="CDR3.amino.acid.sequence")#specify column name with CDR3 aa seqs.
#we could easily convert counts to probs
CDR3s_p$Pgen=(CDR3s_p$sim_num)/(5e5*3)

#We could do more with more cores (currently OS X and Linux). 
#Note, that now total number of generated sequences is nrec*iter*cores = 5e5*3*2.
CDR3s_p<-estimate_pgen_aa(CDR3s,iter=3,nrec=5e5,cores=2,V="TRBV5-1",J="TRBJ1-1",colname="CDR3.amino.acid.sequence")
```

To check for contamination, we could also save generated nucleotide variants for few sequences of interest (_target_ parameter): 
```R
source("generation.R")
CDR3s<-read.csv2("demo/TCRBV07-06_TCRBJ01-04.csv")
#function now returns list with data and variants of target
CDR3s_p<-estimate_pgen_aa_div(CDR3s,iter=3,nrec=5e5,V="TRBV7-6",J="TRBJ1-4",targets="CASSLAPGATNEKLFF")
CDR3s_p$variants

#as expected, all have same amino acid sequence as target
translate_cdr(CDR3s_p$variants)
#table with counts is now in data element of the list: 
sum(CDR3s_p$data$sim_num) #number of simulated amino acid sequences corresponding on data
```

#### Identification of condition-associated clonotypes

This is most advanced part.
Note, that all analysis is needed to be done for each VJ combination separately.

We would need two things: 
1. Table with CDR3 amino acid sequence and columns indicating presence (or read count) in our donors. Typically, this is a result of merge of sample dataframes by CDR3 amino acid sequence. See _demo\/TCRBV05\-01\_TCRBJ01\-01.csv_ for example.
2. Table with sizes of repertoire of given donors VJ combination (= number of unique CDR3 beta clonotypes with this VJ combination), and column indicating if you want to use this sample for current analysis (normally this indicates, that sample corresponds to patient cohort). See _demo\/samples\_TCRBV05\-01\_TCRBJ01\-01.csv_ for example.

Table 1. example: 

CDR3.amino.acid.sequence | Donor1.Read.count | Donor2.Read.count | Donor3.Read.count
------------------------ | ----------------- | ----------------- | -----------------
CASSLAPGATNEKLFF         | 0                 | 103               | 108
CASSLPGTGEKLFF           | 1                 | 0                 | 2

Table 2. example: 

sample | count | analysis
-----  | ----- | --------
Donor1 | 3040  | TRUE
Donor2 | 5304  | TRUE
Donor3 | 1356  | FALSE

Below is a real world example of analysis for TRBV7-6 TRBJ1-4 combination from Emerson et al, Nature genetics, 2017.

```R
source("generation.R")
source("analysis.R")

#load sharing data for TRBV7-6 TRBJ1-4 combination for data from Emerson et al, Nature genetics, 2017.
CDR3s<-read.csv2("demo/TCRBV07-06_TCRBJ01-04.csv")

#lets generate 2e9 sequences to estimate generation probability. 
#This takes some time, approx 2 hours on 8-core intel i7 processor. 
#For demo purposes skip it, and load precomputed table (see below)
CDR3s_p<-estimate_pgen_aa(CDR3s,iter=500,nrec=5e5,cores=8, V="TRBV7-6",J="TRBJ1-4")

#Analogous command for single core (use on windows), this simulation takes approx 16 hours of time, so do not run it.
CDR3s_p<-estimate_pgen_aa(CDR3s,iter=500*8,nrec=5e5,cores=1, V="TRBV7-6",J="TRBJ1-4")

#To save time, we could load precomputed table:
CDR3s_p<-read.csv2("demo/res_TCRBV07-06_TCRBJ01-04.csv")
#read sample sheet
sizes<-read.csv2("demo/samples_TCRBV07-06_TCRBJ01-04.csv")

CDR3s_p<-do_analysis(CDR3s_p,total=2e9,sizevec=sizes$count,indvec=sizes$analysis)
```

This outputs the same table with the additional columns, such as p-value and effect size:

```R
#filter clones with 0 in silico rearrangements:
CDR3s_p<-CDR3s_p[CDR3s_p$sim_num>0,]

#do multiple testing correction and output significant results only with few columns:
CDR3s_p[p.adjust(CDR3s_p$pval_post,method="holm")<0.05,c("sim_num","CDR3.amino.acid.sequence","donors","P_post","pval_post","effect_size")]

```

Here is short final output description. See methods section in the [manuscript](https://www.biorxiv.org/content/early/2017/09/27/195057) for detailed description. 

Field          | Description
-----          | -----------
_sim\_num_     | raw number of simulated TCR amino acid sequences having this CDR3AA
_ML_           | ML estimate of probability of observing sequence from data (P\_data)
_donors_       | number of donors in selected cohort having this CDR3AA
_P\_post_      | theoretical probability to observe sequence, estimated from recombination model
_pval\_post_   | p\-value (not corrected for multiple testing)
_effect\_size_ | log10(effect size) is log10(_ML_)\-log10(_P\_post_)

#### Experiment design
For details see _Designing the experiment_ section in the manuscript. Idea is to make a simulation, to check, if clone with given _P\_data_ and _P\_post_, which is _q_ times lower than _P\_data_ would be found in cohort of size _n_, sequencing depths vector _nvec_, significance threshold _thres_ by our method. One could do _niter_ simulations, and function would return number of simulations when clone is below significant threshold, and also number of donors with clone. Function return list of number of significant tests for each _P\_data_ value, _power_ , and number of donors with sequence in each simulation, _sizes_.  
Let's do a quick example:
```R
source("analysis.R")

#do niter=100 simulations for clone with effect size q=5,
#and cohort size=30 for given values of pdata:
tst_q5<-do_power_analysis(thres = 0.0001,niter = 100,q=5,n=30,nvec=rep(1e3,n),pdata=10^-seq(7,2,length.out = 18))
#lets plot the results!
#Number of significant test depending of clone Pdata
plot(10^-seq(7,2,length.out = 18),tst_q5$power,log="x",type="l",xlab="Pdata'",ylab="# significant results (out of 100)",ylim=c(0,100))

#Average (over 100 simulations for each Pdata) number of donors with sequence:
plot(10^-seq(7,2,length.out = 18),rowSums(tst10c$sizes)/100,log="x",type="l",xlab="Pdata'",ylab="# of donors with sequence")

```
