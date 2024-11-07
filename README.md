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
- high quality site matrix table for sample qc (step 1) <br/>
- pairwise relatedness hail table (step 2) <br/>
- hail table of unrelated samples (step 2) <br/>
- hail table of PC 1-10 (step 3) <br/>
- hail table of samples passing outlier filtering (step 4) <br/>
- hail table of samples passing sex check (step 5) <br/>
- matrix table for running GATK's VQSR (step 6) <br/>

## How to use this repo
This repository is meant to act as a template for basic quality control steps using Hail. All of the quality control steps outlined here can be conducted directly in python using the Hail package. <br/>
The scripts are separated by different QC steps, with each step building off the previous. <br/>
The order of the steps are denoted by the numerical prefix of each files. <br/>
The beginning of each file contains input/output file paths. For each step, specify your own files to run the Hail code on your dataset. <br/>
**The code here is provided as a starting point for QC. Every dataset will require unique QC, therefore, we suggest analysts visualize variant QC metrics, sample QC metrics, PCA plots, and relatedness distributions to adjust QC parameters accordingly.**

## QC Steps
This repo contains 6 scripts to run basic quality control filtering using Hail: <br/>
01_high_quality_sites.py : applies filters for genotypes, removes low-complexity regions, filters to exome regions (if necessary), and filters to common and highly called sites to conduct sample QC (steps 2-5) <br/>
02_relatedness_filtering.py : calculates pairwise kinship and filters to a maximally independent set <br/>
03_ancestry.py : conducts PCA and assigns ancestry using the gnomad python package <br/>
04_sample_qc.py : uses the gnomad package to calculate stratified sample qc metrics and conduct filtering within ancestry and cohort <br/>
05_sex_check.py : calculates F-statistic on X-chromosme to assign sex and filter out samples with discordant sex and gender <br/>
06_output_mt_for_vqsr.py : filters to passing samples and outputs a matrix table of high-quality samples in order to run GATK's VQSR (more info: https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR) <br/>
07_output_qced_data.py : filter to passing samples, removes low-complexity regions, filters to exome regions, and removes sites failing VQSR to produce a QC'ed matrix table <br/>



## Help
Help forum for troubleshooting and general Hail questions: https://discuss.hail.is/ <br/>
Questions about this repo: jsealock@broadinstitute.org 

