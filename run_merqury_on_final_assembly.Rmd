---
title: "Run Merqury on Final BFR Assembly"
author: "Jason Toy"
date: "1/26/2023"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Activate Merqury (v1.3)
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/merqury
```

## Test merqury
```{bash eval = FALSE}
Rscript /hb/groups/bernardi_lab/programs/merqury/share/merqury/plot/plot_spectra_cn.R --help
```

## Create working directory
```{bash eval = FALSE}
mkdir /hb/groups/bernardi_lab/jason/BFR_merqury

cd /hb/groups/bernardi_lab/jason/BFR_merqury/
```

## Run Meryl to get kmer counts
```{bash eval = FALSE}
meryl k=21 count /hb/home/jatoy/bfren_genome/illumina_reads/*trim_paired.fq.gz output BFR_k21.meryl threads=24

```

## Run Merqury
```{bash eval = FALSE}
REF=/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta

/hb/groups/bernardi_lab/programs/merqury/bin/merqury.sh /hb/groups/bernardi_lab/jason/BFR_merqury/BFR_k21.meryl $REF BFR_k21_merq
```

## Plot copy number spectra plots
```{bash eval = FALSE}
Rscript /hb/groups/bernardi_lab/programs/merqury/share/merqury/plot/plot_spectra_cn.R -f BFR_k21_merq.BFR_ref_nuc_mt_final_v2.spectra-cn.hist -o BFR_k21_merq.BFR_ref_nuc_mt_final_v2.spectra-cn
```
