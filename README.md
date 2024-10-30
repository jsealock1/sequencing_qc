# Example Code for Conducting Sequencing QC using Hail

## Installing Hail
Hail documentation: [hail.is](https://hail.is/) <br/>
Hail requires Python3 and Java 11 JRE <br/>
Install Hail: pip install hail <br/>
More detailed instructions: https://hail.is/docs/0.2/getting_started.html <br/>

## Input and Output
Input: Hail can handle a variety of genomic datasets, including VCFs, gVCFs, PLINK bed/bim/fam, and bgen files. More information can be found here: https://hail.is/docs/0.2/methods/impex.html <br/>
<br/>
Output: The final script here outputs a Hail matrix table. Additionally, several intermediate files are generated: <br/>
- variant qc'ed matrix table (step 1) <br/>
- LD pruned variant hail table (step 2) <br/>
- pairwise relatedness hail table (step 2) <br/>
- hail table of unrelated samples (step 2) <br/>
- hail table of PC 1-10 (step 3) <br/>
- hail table of samples passing outlier filtering (step 4) <br/>
- hail table of samples passing sex check (step 5) <br/>

## QC Steps
This repo contains 6 scripts to run basic quality control filtering using Hail: <br/>
01_variant_qc.py : applies filters for genotypes, removes low-complexity regions, and filters to exome regions (if necessary) <br/>
02_relatedness_filtering.py : calculates pairwise kinship and filters to a maximally independent set <br/>
03_ancestry.py : conducts PCA and assigns ancestry <br/>
04_sample_qc.py : uses the gnomad package to calculate stratified sample qc metrics and conduct filtering within ancestry and cohort <br/>
05_sex_check.py : calculates F-statistic on X-chromosme to assign sex and filter out samples with discordant sex and gender <br/>
06_output_qced_data.py : filters to passing samples from each step to output qc'ed matrix table <br/>

## Important Information
All of the quality control steps outlined here can be conducted directly in python using the Hail package. 
However, the code here is provided as a starting point for QC. 
We suggest analysts visualize variant/sample QC metrics and adjust QC parameters accordingly. 


