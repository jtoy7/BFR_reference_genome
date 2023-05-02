---
title: "MinIon genome assembly and analysis - guppy v5.0.15"
author: "Jason Toy"
date: "9/29/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Set-Up Directory

Move or upload all .fast5 output files to your own home sub-directory. My directories are ```/hb/home/jatoy/bfren_genome/minion/runs/cell_*_fast5```. (separate directories for cell_12 and cell_17)

Create two more directories in ```~/bfren_genome/minion/runs/``` called ```cell_12_fastq/``` and ```cell_17_fastq/```. We will direct the output files from guppy to these directories.
\  

## Base Calling

Run guppy (v5.0.15) to convert .fast5 files to .fastq files (base calling). This can be done automatically during the sequencing run, but the run will take longer, and if the computer crashes during that time, you may lose data. Recording the data in the .fast5 file format saves time and is safer. The guppy software then does the base calling and makes a new .fastq file for each .fast5 file. This way you can also run newer basecalling models on the same sequencing data as more accurate models are released.
\  

Create the following slurm script file and save it as "guppy_v5_0_15.slurm"
```{bash eval = FALSE}
#!/bin/bash
#SBATCH --job-name=guppy_2021-10-01
#SBATCH --output=guppy_2021-10-01.out
#SBATCH --error=guppy_02021-10-01.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jatoy@ucsc.edu
#SBATCH --partition=96x24gpu4
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=24

# Above are all the specifics SLURM needs to run your specific job.

# --ntasks-per-node will be used in doParallel.R to specify the number
# of cores to use on the machine. Using 24

# now load programs needed

module load guppy


#run the command:


guppy_basecaller -c dna_r9.4.1_450bps_sup.cfg --num_callers 24 --device cuda:0 --input_path /hb/home/jatoy/bfren_genome/minion/runs/cell_12_fast5 --save_path /hb/home/jatoy/bfren_genome/minion/runs/cell_12_fastq

# --device option says, "Use the first GPU"
# dna_r9.4.1_450bps_sup.cfg is the "super accurate" bonito model for r9 flow cells

```
\  

Submit the job:
```{bash eval = FALSE}
sbatch guppy_v5_0_15.slurm
```

When the job has completed, you should have a .fastq corresponding to each .fast5 file.
\  
  
  
Now repeat with run 17 data:
Edit or create a new slurm script file and save it as "guppy_v5_0_15_run17.slurm"
```{bash eval = FALSE}
#!/bin/bash
#SBATCH --job-name=guppy_2021-10-01_run17
#SBATCH --output=guppy_2021-10-01_run17.out
#SBATCH --error=guppy_02021-10-01_run17.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jatoy@ucsc.edu
#SBATCH --partition=96x24gpu4
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=24

# Above are all the specifics SLURM needs to run your specific job.

# --ntasks-per-node will be used in doParallel.R to specify the number
# of cores to use on the machine. Using 24

# now load programs needed

module load guppy


#run the command:


guppy_basecaller -c dna_r9.4.1_450bps_sup.cfg --num_callers 24 --device cuda:0 --input_path /hb/groups/bernardi_lab/MinIon/BFR_17 --save_path /hb/home/jatoy/bfren_genome/minion/runs/cell_17_fastq

# --device option says, "Use the first GPU"
# dna_r9.4.1_450bps_sup.cfg is the "super accurate" bonito model for r9 flow cells

```

\  

Now concatenate all the .fastq files (from both sequencing runs) into a single file called ```bfren_ref_bigfile.fastq```:
```{bash eval = FALSE}
cd ~/minion/runs

cat cell_12_fastq/pass/*.fastq cell_17_fastq/pass/*.fastq > bfren_ref_bigfile.fastq
```
\  

## Evaluation of Sequencing Run

Use NanoStat (v1.5.0) to summarize the sequencing run and calculate stats. I did this in an interactive session:

```{bash eval = FALSE}

#request a node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

#ssh into your allocated node
ssh $SLURM_NODELIST


#run commands
module load miniconda3.9            #if you don't run this first, NanoStat won't load properly

conda activate /hb/groups/bernardi_lab/programs/nanostat_v1.5.0

cd ~/bren_genome/minion/runs/

NanoStat --fastq bfren_ref_bigfile.fastq --threads 24 --name bfren_ref_bigfile.fastq_statfastq_report
```
\  
  
bfren_ref_bigfile.fastq_statfastq_report:
```
General summary:        
Mean read length:                 2,979.6
Mean read quality:                   14.4
Median read length:               2,072.0
Median read quality:                 14.5
Number of reads:              5,608,949.0
Read length N50:                  4,850.0
STDEV read length:                 2,955.3
Total bases:             16,712,372,110.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:	5608949 (100.0%) 16712.4Mb
>Q7:	5608949 (100.0%) 16712.4Mb
>Q10:	5608667 (100.0%) 16712.0Mb
>Q12:	4764274 (84.9%) 14351.0Mb
>Q15:	2330888 (41.6%) 7254.2Mb
Top 5 highest mean basecall quality scores and their read lengths
1:	37.7 (2036)
2:	35.0 (662)
3:	32.1 (274)
4:	28.6 (675)
5:	28.2 (54)
Top 5 longest reads and their mean basecall quality score
1:	87499 (14.5)
2:	83825 (11.0)
3:	80049 (15.0)
4:	77894 (14.8)
5:	77796 (13.9)
```
\  

Compare with results from just run 12, using guppy v3:
bfren_r_bigfile.fastq_statfastq_report:
```
General summary:        
Mean read length:                2,197.2
Mean read quality:                  10.8
Median read length:                998.0
Median read quality:                11.3
Number of reads:             2,080,084.0
Read length N50:                 4,654.0
Total bases:             4,570,356,427.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:	2011943 (96.7%) 4464.6Mb
>Q7:	1903219 (91.5%) 4267.0Mb
>Q10:	1504444 (72.3%) 3438.1Mb
>Q12:	700977 (33.7%) 1631.8Mb
>Q15:	8653 (0.4%) 7.1Mb
Top 5 highest mean basecall quality scores and their read lengths
1:	20.9 (431)
2:	19.8 (1878)
3:	19.6 (288)
4:	19.3 (985)
5:	19.2 (270)
Top 5 longest reads and their mean basecall quality score
1:	86361 (12.0)
2:	84389 (8.2)
3:	82997 (9.6)
4:	79414 (12.0)
5:	77315 (11.8)
```
\  

Check for fastq errors
First make sure your concatenated fastq files have correctly merged using Fastq Utils (v0.25.1).
```{bash eval = FALSE}

#request a node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

#ssh into your allocated node
ssh $SLURM_NODELIST


#load miniconda
module load miniconda3.9

#activate fastq_utils
conda activate /hb/groups/bernardi_lab/programs/fastq_utils_v0.25.1

#change working directory
cd ~/bfren_genome/minion/runs

#run command
fastq_info bfren_ref_bigfile.fastq.gz

```
\  

fastq_info output:
```
fastq_utils 0.25.1
DEFAULT_HASHSIZE=39000001                                                                 
Scanning and indexing all reads from bfren_ref_bigfile.fastq.gz                           
Read name provided with no suffix                                                         
5600000Scanning complete.                                                                                                                                                           
Reads processed: 5608949                                                                  
Memory used in indexing: ~1025 MB                                                         
------------------------------------                                                      
Number of reads: 5608949                                                                  
Quality encoding range: 34 123                                                            
Quality encoding: sanger                                                                  
Read length: 27 87499 2072                                                                
OK   
```
\  

Run FastQC (v0.11.7) on bigfile to get more detailed QC info on the reads:
```{bash eval = FALSE}

module load fastqc

module load java    #may not be necessary anymore, but if you get a java error, try loading this module


fastqc  bfren_ref_bigfile.fastq.gz -t 24 -o ~/bfren_genome/minion/runs/fastqc/
```
\  

If all looks good, move on to next step to remove NanoPore adapters
\  


## Scaffold Assembly

First, compress the "bigfile" fasta using gzip (the next step requires a zipped file as input).

```{bash eval = FALSE}
gzip -v bfren_ref_bigfile.fastq
```
\  

Then, run PoreChop (v2.4) to remove the NanoPore adaptor sequences from your reads.

```{bash eval = FALSE}
#access node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST

#load miniconda module
module load miniconda3.9

#activate PoreChop
conda activate /hb/groups/bernardi_lab/programs/porechop_v2.4

#change to working directory
cd ~/bfren_genome/minion/runs/

#run PoreChop
porechop -i bfren_ref_bigfile.fastq.gz -o bfren_ref_bigfile_pc2.fastq.gz --threads 24
```
\  
  
  
Porechop output:
```
Trimming adapters from read ends
	   SQK-NSK007_Y_Top: AATGTACTTCGTTCAGTTACGTATTGCT
	SQK-NSK007_Y_Bottom: GCAATACGTAACTGAACGAAGT
       1D2_part_2_start: CTTCGTTCAGTTACGTATTGCTGGCGTCTGCTT
	     1D2_part_2_end: CACCCAAGCAGACGCCAGCAATACGTAACT

5,608,949 / 5,608,949 (100.0%)

5,300,176 / 5,608,949 reads had adapters trimmed from their start (185,081,259 bp removed)
3,482,449 / 5,608,949 reads had adapters trimmed from their end (55,258,957 bp removed)


Splitting reads containing middle adapters
5,608,949 / 5,608,949 (100.0%)

31,973 / 5,608,949 reads were split based on middle adapters


Saving trimmed reads to file
pigz not found - using gzip to compress

Saved result to /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc.fastq.gz 
```
\  

Run fastq_info again on porechopped reads to check for errors in new fastq file:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/fastq_utils_v0.25.1

cd ~/bfren_genome/minion/runs

fastq_info bfren_ref_bigfile_pc.fastq
```
\  


fastq_info output:
```
fastq_utils 0.25.1
DEFAULT_HASHSIZE=39000001
Scanning and indexing all reads from bfren_ref_bigfile_pc.fastq
Read name provided with no suffix
5600000Scanning complete.

Reads processed: 5614527
Memory used in indexing: ~1026 MB
------------------------------------
Number of reads: 5614527
Quality encoding range: 34 123
Quality encoding: sanger
Read length: 13 87449 2029
OK    
```
\  

  
Rerun NanoStat on porechopped reads to summarize the sequencing run and calculate stats.
\  

First unzip bfren_ref_bigfile_pc.fastq.gz
\  

```{bash eval = FALSE}
gunzip bfren_ref_bigfile_pc.fastq.gz
```
\  

Then run NanoStat again:
```{bash eval = FALSE}

#request a node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

#ssh into your allocated node
ssh $SLURM_NODELIST


#run commands
module load miniconda3.9            #if you don't run this first, NanoStat won't load properly

conda activate /hb/groups/bernardi_lab/programs/nanostat_v1.5.0

cd ~/bren_genome/minion/runs/

NanoStat --fastq bfren_ref_bigfile_pc.fastq --threads 24 --name bfren_ref_bigfile_pc.fastq_statfastq_report
```
\  


NanoStat Results for bfren_ref_bigfile_pc.fastq:
```
General summary:         
Mean read length:                  2,932.8
Mean read quality:                    14.6
Median read length:                2,029.0
Median read quality:                  14.7
Number of reads:               5,612,978.0
Read length N50:                   4,830.0
STDEV read length:                 2,947.6
Total bases:              16,461,999,922.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:	5612977 (100.0%) 16462.0Mb
>Q7:	5612965 (100.0%) 16462.0Mb
>Q10:	5610436 (100.0%) 16459.2Mb
>Q12:	4821562 (85.9%) 14193.0Mb
>Q15:	2529345 (45.1%) 7398.0Mb
Top 5 highest mean basecall quality scores and their read lengths
1:	90.0 (169)
2:	90.0 (499)
3:	90.0 (169)
4:	89.5 (233)
5:	51.2 (261)
Top 5 longest reads and their mean basecall quality score
1:	87449 (14.5)
2:	83790 (11.0)
3:	80001 (15.0)
4:	77865 (14.8)
5:	77769 (13.9)
```
\  

Run NanoPlot on the reads to look at read size/quality distributions
```{bash eval = FALSE}
#request a node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST

module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/nanoplot

NanoPlot -t 24 --fastq  bfren_ref_bigfile_pc.fastq --plots hex dot --outdir nanoplot
```
\  

Run FastQC again after runnging PoreChop and compare the results:
``` {bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST

module load fastqc

fastqc bfren_ref_bigfile_pc.fastq --outdir fastqc
```
\  


Now create two filtered sets of reads from the "bigfile": one for longer reads of lower quality and one for shorter reads of higher quality.

```{bash eval = FALSE}
#access node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST


#load miniconda and activate NanoFilt
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/nanofilt


#make sure you're in the correct working directory
cd ~/minion/runs/


#run NanoFilt
gunzip -c bfren_ref_bigfile_pc.fastq.gz | NanoFilt -l 1000 | gzip > bfren_r_bigfile_pc_l000.fastq.gz

gunzip -c bfren_ref_bigfile_pc.fastq.gz | NanoFilt -l 500 | gzip > bfren_r_bigfile_pc_500.fastq.gz


#or if the bigfile_pc.fastq file is already unzipped
NanoFilt bfren_ref_bigfile_pc.fastq -l 1000 | gzip > bfren_ref_bigfile_pc_l000.fastq.gz
NanoFilt bfren_ref_bigfile_pc.fastq -l 500 | gzip > bfren_ref_bigfile_pc_500.fastq.gz
```
\  

Align long reads into scaffolds using two different approaches: once using only the reads > 1000bp, and once using all reads > 500 bp

```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/wtdbg2


wtdbg2 -x ont -g 600m -t 24 -L 1000 -i bfren_ref_bigfile_pc_l000.fastq.gz -fo bfr_ont_1000 #reads > 1000bp

wtdbg2 -x ont -g 600m -t 24 -L 500 -i bfren_ref_bigfile_pc_500.fastq.gz -fo bfr_ont_500 #reads >500bp
```
\  

Now make a consensus sequence of the scaffolds.

```{bash eval = FALSE}
wtpoa-cns -t 24 -i bfr_ont_1000.ctg.lay.gz -fo bfr_ont_1000.ctg.fa

wtpoa-cns -t 24 -i bfr_ont_500.ctg.lay.gz -fo bfr_ont_500.ctg.fa
```
\  

Get stats on your new assemblies and compare them.

```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/assembly_stats

assembly-stats bfr_ont_1000.ctg.fa > bfr_ont_1000.ctg.fa_assemblystats 

assembly-stats bfr_ont_500.ctg.fa > bfr_ont_500.ctg.fa_assemblystats 
```
\  

The assembly that used all reads >500bp was more contiguous (1,008 scaffolds) with a greater average scaffold length (588,787.14 bp) than the assembly that used only reads >1000bp (1,033 scaffolds, average length 575,385.33), so from here on, we will use the more inclusive subset (reads >500 bp) and it's corresponding assembly.

Run NanoStat for the final filtered dataset:
```{bash eval = FALSE}

#request a node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

#ssh into your allocated node
ssh $SLURM_NODELIST


#run commands
module load miniconda3.9            #if you don't run this first, NanoStat won't load properly

conda activate /hb/groups/bernardi_lab/programs/nanostat_v1.5.0

cd ~/bren_genome/minion/runs/

NanoStat --fastq bfren_ref_bigfile_pc_500.fastq --threads 24 --name bfren_ref_bigfile_pc_500.fastq_statfastq_report
```
\  


NanoStat Results for bfren_ref_bigfile_pc_500.fastq:
```
General summary:         
Mean read length:                  3,268.7
Mean read quality:                    14.6
Median read length:                2,418.0
Median read quality:                  14.7
Number of reads:               4,966,516.0
Read length N50:                   4,886.0
STDEV read length:                 2,973.0
Total bases:              16,234,043,863.0
Number, percentage and megabases of reads above quality cutoffs
>Q5:	4966515 (100.0%) 16234.0Mb
>Q7:	4966503 (100.0%) 16234.0Mb
>Q10:	4965319 (100.0%) 16231.7Mb
>Q12:	4276732 (86.1%) 14000.6Mb
>Q15:	2250354 (45.3%) 7299.6Mb
Top 5 highest mean basecall quality scores and their read lengths
1:	37.7 (2036)
2:	35.0 (662)
3:	28.6 (675)
4:	27.7 (506)
5:	27.4 (651)
Top 5 longest reads and their mean basecall quality score
1:	87449 (14.5)
2:	83790 (11.0)
3:	80001 (15.0)
4:	77865 (14.8)
5:	77769 (13.9)
```
\  


# Run Busco (v5.2.2) to assess assembly completeness
\  

```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i ~/bfren_genome/minion/runs/redbean/minlen_500/bfr_ont_500.ctg.fa -o busco -l actinopterygii_odb10 -m genome --cpu 24
```
\  

BUSCO results
```
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

	***** Results: *****

	C:95.3%[S:94.4%,D:0.9%],F:1.9%,M:2.8%,n:3640	   
	3471	Complete BUSCOs (C)			   
	3437	Complete and single-copy BUSCOs (S)	   
	34	Complete and duplicated BUSCOs (D)	   
	69	Fragmented BUSCOs (F)			   
	100	Missing BUSCOs (M)			   
	3640	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 5.34c21f2
```

## Polish Scaffold Assembly with minimap2 and Racon

To polish, we will first align the same MinIon reads we just used for scaffolding (all reads >500 bp) to the scaffold assembly we just created using minimap2 (v2.17-r941).

```{bash eval = FALSE}
#access node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST


#load miniconda module and activate minimap2
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/minimap2


#make sure you're in the correct working directory
cd ~/minion/runs/


#run minimap2
minimap2 -ax map-ont -t 24 redbean/minlen_500/bfr_ont_500.ctg.fa bfren_ref_bigfile_pc_500.fastq.gz > bfr_ont_500_mmap.sam

# -a instructs minimap2 to output in .sam format
# -x map-ont sets the tuning parameters to presets optimized for ONT reads
```
\  

Now that the reads are mapped to the draft assembly, we will build a new consensus assembly using Racon (v1.4.13).

```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/racon

#run racon
racon -u -m 5 -x -4 -g -8 -w 500 -t 24 bfren_ref_bigfile_pc_500.fastq.gz bfr_ont_500_mmap.sam redbean/minlen_500/bfr_ont_500.ctg.fa > bfr_ont_500_racon.fa

# usage = racon [options ...] <mapped sequences> <alignment .sam file> <target sequences>
# -u = same as --include-unpolished; output unpolished target (reference) sequences
# -m = same as --match <int>; default is 3; score for matching bases
# -x = same as --mismatch <int>; default is -5; score for mismatching bases
# -g = same as --gap <int>, default: -4, gap penalty (must be negative)
# -w = same as --window-length <int>; default is 500; size of window on which POA is performed
# -t = number of threads

```
\  

Run assembly-stats on the now polished genome.

```{bash eval = FALSE}
module load miniconda3

conda activate /hb/groups/bernardi_lab/programs/assembly_stats

#run assembly-stats
assembly-stats bfr_ont_500_racon.fa > bfr_ont_500_racon.fa_assemblystats
```
\  

Note: enter these stats and the original assembly stats into a spreadsheet and add a new entry after each round of polishing to keep track of how each round influences the quality of the genome.
\  

# BUSCO the once-polished assembly to see if completeness improved
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i ~/bfren_genome/minion/runs/bfr_ont_500_racon.fa -o busco_racon -l actinopterygii_odb10 -m genome --cpu 24
```
\  

BUSCO results
```
# BUSCO version is: 5.2.2 
# The lineage dataset is: actinopterygii_odb10 (Creation date: 2021-02-19, number of genomes: 26, number of BUSCOs: 3640)
# Summarized benchmarking in BUSCO notation for file /hb/home/jatoy/bfren_genome/minion/runs/bfr_ont_500_racon.fa
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

	***** Results: *****

	C:97.4%[S:96.7%,D:0.7%],F:1.0%,M:1.6%,n:3640	   
	3547	Complete BUSCOs (C)			   
	3520	Complete and single-copy BUSCOs (S)	   
	27	Complete and duplicated BUSCOs (D)	   
	37	Fragmented BUSCOs (F)			   
	56	Missing BUSCOs (M)			   
	3640	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 5.34c21f2
```
\  

## Second NanoPore Polish
\  

We will now polish the once-polished assembly again with the same nanopore data as last time.  
First map the nanopore reads to the once-polished assembly (racon output) using minimap2:
```{bash eval = FALSE}
#access node
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST


#load miniconda module and activate minimap2
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/minimap2


#make sure you're in the correct working directory
cd ~/minion/runs/


#run minimap2
minimap2 -ax map-ont -t 24 bfr_ont_500_racon.fa bfren_ref_bigfile_pc_500.fastq.gz > bfr_ont_500_mmap_2.sam

# -a instructs minimap2 to output in .sam format
# -x map-ont sets the tuning parameters to presets optimized for ONT reads
```
\  

Now make the consensus sequences (Racon)
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/racon

#run racon
racon -u -m 5 -x -4 -g -8 -w 500 -t 24 bfren_ref_bigfile_pc_500.fastq.gz bfr_ont_500_mmap_2.sam bfr_ont_500_racon.fa > bfr_ont_500_racon_2.fa

# usage = racon [options ...] <mapped sequences> <alignment .sam file> <target sequences>
# -u = same as --include-unpolished; output unpolished target (reference) sequences
# -m = same as --match <int>; default is 3; score for matching bases
# -x = same as --mismatch <int>; default is -5; score for mismatching bases
# -g = same as --gap <int>, default: -4, gap penalty (must be negative)
# -w = same as --window-length <int>; default is 500; size of window on which POA is performed
# -t = number of threads

```
\  

```
[racon::Polisher::initialize] loaded target sequences 2.900655 s                          [racon::Polisher::initialize] loaded sequences 320.725576 s                               [racon::Polisher::initialize] loaded overlaps 87.950252 s                                 [racon::Polisher::initialize] aligning overlaps [====================] 91.493828 s        [racon::Polisher::initialize] transformed data into windows 52.285708 s                   [racon::Window::generate_consensus] warning: contig 398 might be chimeric in window 2!    [racon::Polisher::polish] generating consensus [====================] 909.857116 s        [racon::Polisher::] total = 1480.237076 s  
```

Run assembly-stats again:
```{bash eval = FALSE}
module load miniconda3

conda activate /hb/groups/bernardi_lab/programs/assembly_stats

#run assembly-stats
assembly-stats bfr_ont_500_racon_2.fa > bfr_ont_500_racon_2.fa_assemblystats
```
\  

Run BUSCO again:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i ~/bfren_genome/minion/runs/bfr_ont_500_racon_2.fa -o busco_racon_2 -l actinopterygii_odb10 -m genome --cpu 24
```
\  

BUSCO results:
```
# BUSCO version is: 5.2.2 
# The lineage dataset is: actinopterygii_odb10 (Creation date: 2021-02-19, number of genomes: 26, number of BUSCOs: 3640)
# Summarized benchmarking in BUSCO notation for file /hb/home/jatoy/bfren_genome/minion/runs/bfr_ont_500_racon_2.fa
# BUSCO was run in mode: genome
# Gene predictor used: metaeuk

	***** Results: *****

	C:97.7%[S:97.0%,D:0.7%],F:0.9%,M:1.4%,n:3640	   
	3555	Complete BUSCOs (C)			   
	3529	Complete and single-copy BUSCOs (S)	   
	26	Complete and duplicated BUSCOs (D)	   
	32	Fragmented BUSCOs (F)			   
	53	Missing BUSCOs (M)			   
	3640	Total BUSCO groups searched		   

Dependencies and versions:
	hmmsearch: 3.1
	metaeuk: 5.34c21f2
```
\  


## Polish with Short-Read (Illumina) Data
\  

First combine all read 1 files into one file and all read 2 files into another file:
```{bash eval = FALSE}

cd /hb/home/jatoy/bfren_genome/illumina_reads

cat BFR_ref_CSFP210002260-1a_H3MTYDSX2_L1_1.fq BFR_ref_CSFP210002260-1a_H55VJDSX2_L2_1.fq > BFR_ref_sr_1.fq

cat BFR_ref_CSFP210002260-1a_H3MTYDSX2_L1_2.fq BFR_ref_CSFP210002260-1a_H55VJDSX2_L2_2.fq > BFR_ref_sr_2.fq
```


# FastQC
Check read quality and assess level of adapter contamination for the combined files using FastQC (v0.11.7):
```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST

module load fastqc

cd /hb/home/jatoy/bfren_genome/illumina_reads

fastqc  BFR_ref_sr_1.fq --outdir fastqc
fastqc  BFR_ref_sr_2.fq --outdir fastqc
```
\  
These actually looked really clean, so I think Novogene trimmed them already.

# Trim reads using Trimmomatic (v0.39)
``` {bash eval = FALSE}
module load trimmomatic

java -jar /hb/software/apps/trimmomatic/gnu-0.39/trimmomatic-0.39.jar PE -threads 24 -phred33 -summary summary_stats BFR_ref_sr_1.fq BFR_ref_sr_2.fq  BFR_ref_sr_1_trim_paired.fq.gz BFR_ref_sr_1_trim_unpaired.fq.gz BFR_ref_sr_2_trim_paired.fq.gz BFR_ref_sr_2_trim_unpaired.fq.gz ILLUMINACLIP:/hb/software/apps/trimmomatic/gnu-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:2 TRAILING:2 MINLEN:25
```
\  

Options:  
`ILLUMINACLIP:TruSeq3-PE.fa:2:30:10` Remove adapters  
`LEADING:2` Remove leading low quality or N bases (below quality 3)  
`TRAILING:2` Remove trailing low quality or N bases (below quality 3)  
`MINLEN:25` Drop reads below the 36 bases long  
`SLIDINGWINDOW:4:15` Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15
\  

Trimming summary:
```
Input Read Pairs: 95669159
Both Surviving Reads: 95215363
Both Surviving Read Percent: 99.53
Forward Only Surviving Reads: 452066
Forward Only Surviving Read Percent: 0.47
Reverse Only Surviving Reads: 0
Reverse Only Surviving Read Percent: 0.00
Dropped Reads: 1730
Dropped Read Percent: 0.00
```
\  

# Run FastQC again and compare quality and adapter presence before and after trimming
Check read quality and assess level of adapter contamination for the combined files using FastQC (v0.11.7):
```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive

ssh $SLURM_NODELIST

module load fastqc

cd /hb/home/jatoy/bfren_genome/illumina_reads

fastqc  BFR_ref_sr_1_trim_paired.fq.gz --outdir fastqc
fastqc  BFR_ref_sr_2_trim_paired.fq.gz --outdir fastqc
```
\  

# BWA (Burrows-Wheeler Alignment)
First, create an index for the current reference (`bfr_ont_500_racon_2.fa`) using the bwa index command (v0.7.17-r1188):
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/bwa

bwa index ~/bfren_genome/minion/runs/bfr_ont_500_racon_2.fa
```
\  

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa.  
\  
Now map Illumina short reads to assembly:
```{bash eval = FALSE}
bwa mem -t 24 bfr_ont_500_racon_2.fa BFR_ref_sr_1_trim_paired.fq.gz BFR_ref_sr_2_trim_paired.fq.gz > BFR_ref_sr_bwa_aligned.sam
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -Sb -@ 24 -O BAM -o BFR_ref_sr_bwa_aligned.bam BFR_ref_sr_bwa_aligned.sam

# S = input format is auto-detected
# b = output format BAM
# @ = number of threads
# O = specify output format
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_sr_bwa_aligned_sorted.bam -O BAM -@ 24 BFR_ref_sr_bwa_aligned.bam
```
\  

Pilon requires the BAM input file to be indexed, so create the index file (.bai):
```{bash eval = FALSE}
samtools index BFR_ref_sr_bwa_aligned_sorted.bam
```


Now run Pilon (v1.23) to polish second Racon assembly with the short read BAM file:
```{bash eval = FALSE}
salloc --partition=256x44 --nodes=1 --time=96:00:00 --exclusive

module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/pilon

java -Xmx200G -jar /hb/groups/bernardi_lab/programs/pilon/share/pilon-1.23-2/pilon-1.23.jar --genome bfr_ont_500_racon_2.fa --bam BFR_ref_sr_bwa_aligned_sorted.bam  --output BFR_ref_pilon --outdir /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_1 --diploid --changes --vcf --tracks --threads 10 > /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_1/BFR_ref_pilon_out_256.log
```
\ 
I ran Pilon with the --diploid parameter. Remy did not do this, but according to Bruce Walker (the creator of Pilon), "Really, the only thing "--diploid" does is affect whether to report hererozygous SNPs and small indels. Since many of the initial Pilon applications were bacterial, by default it treats mixed evidence as ambiguous rather than as a heterozygous SNP or indel. Pilon's support for diploid calls isn't very sophisticated; in particular, it isn't able to generate multiple haplotypes through local reassembly (it's only calling heterzygosity via the alignment pileup information)."
\  

Now run assembly-stats on the new Pilon assembly:
```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/assembly_stats

assembly-stats BFR_ref_pilon.fasta > BFR_ref_pilon.fasta_assemblystats
```
\

Run BUSCO on Pilon assembly:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_1/BFR_ref_pilon.fasta -o busco_pilon_1 -l actinopterygii_odb10 -m genome --cpu 24
```
\  

## Second Pilon Polish
\  
# BWA (Burrows-Wheeler Alignment)
First, create an index for the current reference (`BFR_ref_pilon.fasta`) using the bwa index command (v0.7.17-r1188):
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/bwa

bwa index ~/bfren_genome/illumina_reads/pilon/polish_1/BFR_ref_pilon.fasta
```
\  

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa.  
\  

Now map Illumina short reads to assembly:
```{bash eval = FALSE}
bwa mem -t 24 ~/bfren_genome/illumina_reads/pilon/polish_1/BFR_ref_pilon.fasta ~/bfren_genome/illumina_reads/BFR_ref_sr_1_trim_paired.fq.gz ~/bfren_genome/illumina_reads/BFR_ref_sr_2_trim_paired.fq.gz > BFR_ref_sr_bwa_aligned_2.sam
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -Sb -@ 24 -O BAM -o BFR_ref_sr_bwa_aligned_2.bam BFR_ref_sr_bwa_aligned_2.sam

# S = input format is auto-detected
# b = output format BAM
# @ = number of threads
# O = specify output format
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_sr_bwa_aligned_sorted_2.bam -O BAM -@ 24 BFR_ref_sr_bwa_aligned_2.bam
```
\  

Pilon requires the BAM input file to be indexed, so create the index file (.bai):
```{bash eval = FALSE}
samtools index BFR_ref_sr_bwa_aligned_sorted_2.bam
```


Now run Pilon (v1.23) to polish the first Pilon assembly with the new short read BAM file:
```{bash eval = FALSE}
salloc --partition=256x44 --nodes=1 --time=96:00:00 --exclusive
ssh $SLURM_NODELIST


module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/pilon

java -Xmx200G -jar /hb/groups/bernardi_lab/programs/pilon/share/pilon-1.23-2/pilon-1.23.jar --genome ~/bfren_genome/illumina_reads/pilon/polish_1/BFR_ref_pilon.fasta --bam /hb/home/jatoy/bfren_genome/illumina_reads/bwa/polish_2/BFR_ref_sr_bwa_aligned_sorted_2.bam --output BFR_ref_pilon_2 --outdir /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_2 --diploid --changes --vcf --tracks --threads 10 > /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_2/BFR_ref_pilon_2_out.log
```
\  


Now run assembly-stats on the new Pilon assembly:
```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/assembly_stats

assembly-stats BFR_ref_pilon_2.fasta > BFR_ref_pilon_2.fasta_assemblystats
```
\

Run BUSCO on Pilon assembly:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_2/BFR_ref_pilon_2.fasta -o busco_pilon_2 -l actinopterygii_odb10 -m genome --cpu 24
```
\  

## BLAST polished assembly (v2.12.0):
```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=96:00:00 --exclusive
ssh $SLURM_NODELIST


module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/blast


#set path for blast databases
export PATH=$PATH:/hb/groups/bernardi_lab/programs/blast
export BLASTDB=/hb/groups/bernardi_lab/programs/blobtools2/nt

#run blast
blastn -db nt -query /hb/home/jatoy/bfren_genome/illumina_reads/pilon/polish_2/BFR_ref_pilon_2.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-20 -out /hb/home/jatoy/bfren_genome/BFR_ref_pilon_2_blast.out -num_threads 12

# qseqid = Query Seq-id
# staxids = unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
# std = default format specifiers = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
```
\ 
  

#### Note about BLAST parameters  
	
*-max_target_seqs* = Number of aligned sequences to keep. Use with report formats that do not have separate definition line and alignment sections such as tabular (all outfmt > 4). Not compatible with num_descriptions or num_alignments. Ties are broken by order of sequences in the database.

**IMPORTANT!** - When using -max_target_seqs, "BLAST returns the first N hits that exceed the specified E-value threshold, which may or may not be the highest scoring N hits. The invocation using the parameter ‘-max_target_seqs 1’ simply returns the first good hit found in the database, not the best hit as one would assume. Worse yet, the output produced depends on the order in which the sequences occur in the database. For the same query, different results will be returned by BLAST when using different versions of the database even if all versions contain the same best hit for this database sequence. Even ordering the database in a different way would cause BLAST to return a different ‘top hit’ when setting the max_target_seqs parameter to 1." (Shah et al. 2019 Bioinformatics)

*-max_hsps* = Maximum number of HSPs (alignments) to keep for any single query-subject pair. The HSPs shown will be the best as judged by expect value. This number should be an integer that is one or greater. If this option is not set, BLAST shows all HSPs meeting the expect value criteria. Setting it to one will show only the best HSP for every query-subject pair.

*-outfmt* = Output format. The format specified above is the default format expected by blobtools2 as input.

More details on BLAST parameters here:
https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/
\  
  
\  

## Map Illumina short reads once more to your final assembly for use in Blobtools2
\  

#### BWA (Burrows-Wheeler Alignment)
First, create an index for the final reference (`BFR_ref_pilon_2.fasta`) using the bwa index command (v0.7.17-r1188):
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/bwa

bwa index ~/bfren_genome/illumina_reads/pilon/polish_2/BFR_ref_pilon_2.fasta
```
\  

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa.  
\  

Now map Illumina short reads to assembly:
```{bash eval = FALSE}
bwa mem -t 24 ~/bfren_genome/illumina_reads/pilon/polish_2/BFR_ref_pilon_2.fasta ~/bfren_genome/illumina_reads/BFR_ref_sr_1_trim_paired.fq.gz ~/bfren_genome/illumina_reads/BFR_ref_sr_2_trim_paired.fq.gz > BFR_ref_sr_bwa_aligned_3.sam
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -Sb -@ 24 -O BAM -o BFR_ref_sr_bwa_aligned_3.bam BFR_ref_sr_bwa_aligned_3.sam

# S = input format is auto-detected
# b = output format BAM
# @ = number of threads
# O = specify output format
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_sr_bwa_aligned_sorted_3.bam -O BAM -@ 24 BFR_ref_sr_bwa_aligned_3.bam
```
\  
  
\  
  
## Remove contamination from assembly using Blobtools2 (v3.0.0)
\  
To run Blobtools2 you will need:
  
1. Your final assembly fasta file ```BFR_ref_pilon_2.fasta```
  
2. Your .bam coverage file for the final assembly (e.g., using bwa to map illumina sequences to your final assembly) ```BFR_ref_sr_bwa_aligned_sorted_3.bam```
  
3. Your BLAST output for the final assembly ```BFR_ref_pilon_2_blast.out```
  
4. (optional) The BUSCO output full table for your final assembly ```full_table.tsv```  
\  
  
Once you have all the files, running blobtools is relatively easy and fast.  
\  

Start an interactive session:
``` {bash eval = FALSE}

salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive --cpus-per-task=8
ssh $SLURM_NODELIST
```
\  

Now create your blobtools database and populate it with the coverage data, taxonomic information, and BUSCO scores:
``` {bash eval = FALSE}
#------- load programs --------------
module load miniconda3.9
conda activate /hb/groups/bernardi_lab/programs/blobtools2

#------- Run Commands ---------------
#Assign source directory variable to environment
SRC=/hb/home/jatoy/bfren_genome/

#Create a BlobDir
blobtools add --create --fasta $SRC/illumina_reads/pilon/polish_2/BFR_ref_pilon_2.fasta $SRC/blobtools2/BFR_ref_pilon_2_blob

#Add BLAST hits
blobtools add --threads 8 --hits $SRC/BFR_ref_pilon_2_blast.out --taxrule bestsumorder --taxdump /hb/groups/bernardi_lab/programs/blobtools2/taxdump $SRC/blobtools2/BFR_ref_pilon_2_blob

#Add mapping coverage
blobtools add --threads 8 --cov $SRC/illumina_reads/bwa/final_mapping_for_blobtools/BFR_ref_sr_bwa_aligned_sorted_3.bam $SRC/blobtools2/BFR_ref_pilon_2_blob

#Add BUSCO
blobtools add --threads 8 --busco $SRC/illumina_reads/pilon/busco_pilon_2/run_actinopterygii_odb10/full_table.tsv $SRC/blobtools2/BFR_ref_pilon_2_blob
```
\  

`--taxrule`: BlobTools2 assigns a putative taxonomic origin to each scaffold at 8 ranks from superkingdom to species based on aggregating hits in the --hits files. The --taxrule flag determines the rule to use when assigning BLAST hits to taxa. Two options are available: "bestsum" and "bestsumorder". Of the taxa represented in the BLAST hits, "bestsum" assigns the taxon for which the sum of bitscores is greatest across all files while the default "bestsumorder" taxrule uses the maximum bitscore in the first file and only uses hits from subsequent files for taxa without a hit in previous files.
\  
  
\  

### Download the BlobDir to desktop to view dataset (v3.1.0)
\
First download Miniconda installer (from Anaconda website) and install via WSL command line:
```{bash eval = FALSE}
bash Miniconda3-latest-Linux-x86_64.sh
```
\  

Now install blobtools2:
```{bash eval = FALSE}
#create directory for blobtookit
mkdir blobtoolkit
cd blobtoolkit

#create conda env with dependencies
conda create -n blobtools2 -y python=3.6 docopt pyyaml ujson pysam tqdm nodejs seqtk
conda activate blobtools2

#download taxdump database
mkdir -p taxdump
cd taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd ..

#download blobtools2
git clone https://github.com/blobtoolkit/blobtools2
git clone https://github.com/blobtoolkit/specification
git clone https://github.com/blobtoolkit/insdc-pipeline

#download and install blobtoolkit viewer
git clone https://github.com/blobtoolkit/viewer
cd viewer
npm install
cd ..

#add new directories to path
export PATH=~/blobtoolkit/blobtools2:~/blobtoolkit/specification:~/blobtoolkit/insdc-pipeline/scripts:$PATH

#now run pip to install everything
pip install blobtoolkit
```
\  


Download BFR_ref_pilon_2_blob folder to desktop and run blobtools view.

``` {bash eval = FALSE}
scp -r jatoy@hb.ucsc.edu:~/bfren_genome/blobtools2/BFR_ref_pilon_2_blob/ ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools
```
\  

Create interactive html page with the blobtools results:

```{bash eval = FALSE}
#activate conda env
conda activate blobtools2

#change to desired directory
cd ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools

#run view command
blobtools view --remote BFR_ref_pilon_2_blob
```
\  
  
Running ```blobtools view``` will output a local url that you can copy and paste into a web browser to view the figures and summary:  
\  

**Blob plot**
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_pilon_2/BFR_ref_pilon_2_blob.blob.circle.png)

\  

**Cumulative composition plot**
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_pilon_2/BFR_ref_pilon_2_blob.cumulative.png)

\  

**Snail plot**
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_pilon_2/BFR_ref_pilon_2_blob.snail.png)
\  

Save the data table from the "table" tab in the viewer as a .csv. Open in excel and sort by "bestsumorder phylum". This will allow you to see which contigs were identified as contaminant taxa.
\  

\  

### Remove contaminants

On the cluster, make a file with the list of contig numbers to KEEP in the assembly by copying and pasting into a new text file (`BFR_ref_pilon_2_ctg_clean.txt`).
```{bash eval = FALSE}
nano BFR_ref_pilon_2_ctg_clean.txt

head BFR_ref_pilon_2_ctg_clean.txt
```
```
ctg1_pilon_pilon
ctg2_pilon_pilon
ctg3_pilon_pilon
ctg4_pilon_pilon
ctg5_pilon_pilon
ctg6_pilon_pilon
ctg7_pilon_pilon
ctg8_pilon_pilon
ctg9_pilon_pilon
ctg10_pilon_pilon 
```

\  

Using samtools v(1.14), we will extract from the assembly only the contigs identified as Chordata or no hit. Contigs identified as Proteobacteria will be excluded.

```{bash eval = FALSE}
module load samtools

cd ~/bfren_genome/blobtools2

xargs samtools faidx ../illumina_reads/pilon/polish_2/BFR_ref_pilon_2.fasta < BFR_ref_pilon_2_ctg_clean.txt > BFR_ref_clean.fasta
```
\  

Check final assembly stats
```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/assembly_stats

assembly-stats BFR_ref_clean.fasta
```

```
stats for BFR_ref_clean.fasta
sum = 595959787, n = 1003, ave = 594177.26, largest = 12326090
N50 = 2589815, n = 55
N60 = 1896895, n = 82
N70 = 1185759, n = 122
N80 = 758445, n = 184
N90 = 406612, n = 292
N100 = 1298, n = 1003
N_count = 0
Gaps = 0
```
\
In total, blobtools removed **5 contigs** containing **27,511 bp**.
\  

### Run BUSCO (v5.2.2) on final assembly

```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=96:00:00 --exclusive
ssh $SLURM_NODELIST

module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i /hb/home/jatoy/bfren_genome/blobtools2/BFR_ref_clean.fasta -o busco_blob_clean -l actinopterygii_odb10 -m genome --cpu 24
```

```
C:98.1%[S:97.4%,D:0.7%],F:0.5%,M:1.4%,n:3640 
3572 Complete BUSCOs (C)
3547 Complete and single-copy BUSCOs (S)
25 Complete and duplicated BUSCOs (D) 
20 Fragmented BUSCOs (F) 
48 Missing BUSCOs (M) 
3640 Total BUSCO groups searched 
```

BUSCO scores were not changed by the removal of contaminant sequences.
\  

\  

### Reorder and rename contigs

I would now like to order my contigs from largest to smallest and rename contigs to a sequential order.

Right now my fasta file looks like this:
```{bash eval = FALSE}
grep ">" BFR_ref_clean.fasta | head
```

```
>ctg1_pilon_pilon
>ctg2_pilon_pilon
>ctg3_pilon_pilon
>ctg4_pilon_pilon
>ctg5_pilon_pilon
>ctg6_pilon_pilon
>ctg7_pilon_pilon
>ctg8_pilon_pilon
>ctg9_pilon_pilon
>ctg10_pilon_pilon
```

```{bash eval = FALSE}
grep ">" BFR_ref_clean.fasta | tail
```

```
>ctg911_pilon_pilon
>ctg940_pilon_pilon
>ctg928_pilon_pilon
>ctg945_pilon_pilon
>ctg935_pilon_pilon
>ctg957_pilon_pilon
>ctg943_pilon_pilon
>ctg833_pilon_pilon
>ctg963_pilon_pilon
>ctg993_pilon_pilon
```
\  

We will use `seqkit` (v2.1.0) to order contigs by length and rename

```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/seqkit

#order contigs by length
seqkit sort --by-length --reverse BFR_ref_clean.fasta > BFR_ref_clean_ordered.fasta

#rename contigs
seqkit replace --pattern '.+' --replacement 'Contig_{nr}' BFR_ref_clean_ordered.fasta > BFR_ref_clean_ordered_renamed.fasta

#where {nr} is the record number, starting from 1
#"." in regex means "any character except line break"
#"+" in regex means "occurring one or more times"
#so "replace --pattern '.+'" basically means "replace all characters in the contig name
```
\  

After sorting by length:
```{bash eval = FALSE}
grep ">" BFR_ref_clean_ordered.fasta | head
```

```
>ctg24_pilon_pilon
>ctg1_pilon_pilon
>ctg4_pilon_pilon
>ctg2_pilon_pilon
>ctg6_pilon_pilon 
>ctg28_pilon_pilon 
>ctg10_pilon_pilon 
>ctg36_pilon_pilon 
>ctg83_pilon_pilon 
>ctg21_pilon_pilon 
```
\  

After renaming:
```{bash eval = FALSE}
grep ">" BFR_ref_clean_ordered_renamed.fasta | head
```

```
>Contig_1
>Contig_2
>Contig_3
>Contig_4
>Contig_5
>Contig_6
>Contig_7
>Contig_8
>Contig_9
>Contig_10
```

```{bash eval = FALSE}
grep ">" BFR_ref_clean_ordered_renamed.fasta | tail
```

```
>Contig_994
>Contig_995
>Contig_996
>Contig_997
>Contig_998
>Contig_999
>Contig_1000
>Contig_1001
>Contig_1002
>Contig_1003
```
\  

The final genome assembly `BFR_ref_ordered_renamed.fasta` will be renamed to `BFR_ref_final.fasta`

```{bash eval = FALSE}
cp BFR_ref_clean_ordered_renamed.fasta BFR_ref_final.fasta

mv ../final_assembly
```
\  

Get final assembly stats (should be same as before, since all we did was reorder contigs and change names):
```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/assembly_stats/

assembly-stats BFR_ref_final.fasta
```

```
stats for BFR_ref_final.fasta
sum = 595959787, n = 1003, ave = 594177.26, largest = 12326090
N50 = 2589815, n = 55
N60 = 1896895, n = 82
N70 = 1185759, n = 122
N80 = 758445, n = 184
N90 = 406612, n = 292
N100 = 1298, n = 1003
N_count = 0
Gaps = 0
```
\  

Extract the lengths of each contig in the fasta file and sort from largest to smallest. Our file has just been sorted from largest to smallest so it should already show a decreasing contig length.
```{bash eval = FALSE}
### Using samtools (v1.14) ###

module load samtools

#index fasta
samtools faidx BFR_ref_final.fasta

#select only first two fields
cut -f1-2 BFR_ref_final.fasta.fai > BFR_ref_final_lengths
```

```{bash eval = FALSE}
### Using seqkit (v2.1.0) and including gc content ###

conda activate /hb/groups/bernardi_lab/programs/seqkit/ 

seqkit fx2tab --name --length --gc BFR_ref_final.fasta > BFR_ref_final_lengths_gc
```
\

Download both files
```{bash eval = FALSE}
scp jatoy@hb.ucsc.edu:~/bfren_genome/final_assembly/BFR_ref_final_lengths* ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/
```
\  

Plot contig lengths (from largest to smallest)
```{R eval = TRUE, results = "hide", warning = FALSE, message = FALSE}
#load data and packages

setwd("C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome")

library(tidyverse)

bfr <- read_tsv("BFR_ref_final_lengths_gc.tsv", col_names = c("contig", "length", "gc")) %>% 
  mutate(id = row_number())
```

```{R eval = TRUE}
ggplot(data = bfr) +
  geom_col(aes(x = id, y = length/1000000)) +
  theme_bw() +
  xlab("Contig") +
  ylab("Length (Mb)")
  
```
\  

Plot gc content by contig
```{R eval = TRUE}
ggplot(data = bfr) +
  geom_col(aes(x = id, y = gc)) +
  theme_bw() +
  xlab("Contig") +
  ylab("GC")
```
\  

Plot gc content by contig length
```{R eval = TRUE}
ggplot(data = bfr) +
  geom_point(aes(x = length/1000000, y = gc)) +
  theme_bw() +
  xlab("Contig length (Mb)") +
  ylab("GC")
```
\  


---


## Ragoo and Alignment to Chromosome-Level Genomes of Closely Related Species

To align/scaffold your new genome to existing, chromosome level genomes, use RagTag and the fasta file for the reference genome. Before running, make sure the reference file is in a directory to which you have access.
```{bash eval = FALSE}
module load miniconda3

conda activate /hb/groups/bernardi_lab/programs/ragtag

module deactivate python-3.6.2 #I had to unload the default loaded python on hb to get this to work

#run ragoo
ragtag.py scaffold -t 20 /hb/groups/bernardi_lab/jason/EJA_ccgp_genome/Embiotoca_jacksoni/fEmbJac1.NCBI.p_ctg.fasta.gz BFR_ref_final.fasta -o ./ragtag_BFR_EJA

```
\  

By default, this will create a new directory called ```ragtag_output```, inside of which will be a new assembly file called ```ragtag.scaffold.fasta```. Use the -o flag to change the output directory name. Manually change the new scaffolded assembly name to ```BFR_EJA_ragtag.fasta```.
\  

Check to make sure it worked properly by printing the names of new scaffolds that should've been created in the last step:
```{bash eval = FALSE}
grep “>” BFR_EJA_ragtag.fasta
```
\  

Get summary stats for scaffolding:
```{bash eval = FALSE}
cat ragtag.scaffold.stats
```

```
placed_sequences    placed_bp   unplaced_sequences    unplaced_bp   gap_bp    gap_sequences
773                 594231042   230                   1728745       63300     633 
```
/ 

Side Note: If you ever want to pull out a single chromosome, follow the example code shown below. You could also do this multiple times to pull out only the large chromosomes.
```{bash eval = FALSE}
module load samtools

samtools faidx BFR_EJA_ragtag.fasta mtDNA_RagTag -o mtDNA_BFR_EJA_ragtag.fasta
```
\  

Next, download the three genomes (your final clean assembly, the new ragtag-scaffolded assembly, and the reference species) to your local machine. They will need to be compressed. So gzip them if they are not.
```{bash eval = FALSE}
gzip -v BFR_ref_final.fasta
gzip -v BFR_EJA_ragtag.fasta

#BFR_EJA_RAGOO
scp jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly/BFR_ref_final.fasta.gz ~/winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome

scp jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly/ragtag_BFR_EJA/BFR_EJA_ragtag.fasta.gz ~/winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome

#EJA
scp jatoy@hb.ucsc.edu:/hb/groups/bernardi_lab/jason/EJA_ccgp_genome/Embiotoca_jacksoni/fEmbJac1.NCBI.p_ctg.fasta.gz ~/winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome

```
\  

Go to the D-GENIES website (http://dgenies.toulouse.inra.fr/) to visually map the assemblies against each other. We will do this twice to map both the independent BFR genome and the EJA-RagTag-mapped BFR genome against the EJA genome. D-GENIES uses minimap2 for mapping. "Target" will be the EJA reference and "Query" will be the BFR genome.  

This will take a few minutes to upload and run depending on your internet speed.  
  
Once the analysis is complete. It the site will email you to let you know you results are ready. Or you can just wait for it to finish. Once the plot is generated, you can click the "Sort contigs" button on the bottom right to rearrange your contigs on the y-axis so that it is easier to visually compare to the reference.  

Download the plot image, association table, and .paf file from the drop-down menu at the top and the summary figure/tsv from the "Summary" button on the bottom right.  
\

#### **BFR_ref_final (no scaffolding) -> EJA**

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/D-GENIES/no_scaffolding/map_BFR_ref_final_to_fEmbJac1.NCBI.p_ctg.png)

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/D-GENIES/legend_small.png)  
\  
Summary  
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/D-GENIES/no_scaffolding/summary_BFR_ref_final_to_fEmbJac1.NCBI.p_ctg.png)  
\
  
\

#### **BFR_EJA_ragtag -> EJA**

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/D-GENIES/ragtag_scaffolding_to_eja/map_BFR_EJA_ragtag_to_fEmbJac1.NCBI.p_ctg.png)  
\

Summary  
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/D-GENIES/ragtag_scaffolding_to_eja/summary_BFR_EJA_ragtag_to_fEmbJac1.NCBI.p_ctg.png) 

Note that the two mappings look very similar which is a good sign for the quality of our BFR genome. Both the EJA-scaffolded and unscaffolded BFR assemblies are useful. The scaffolded assembly provides a more contiguous assembly, but it assumes that the genome structure (e.g. synteny, number of chromosomes) of BFR is the same as in EJA. Since these two species are very closely related, this may indeed be the case here, but will not always be true. The advantage of the unscaffolded assembly is that is was assembled completely independently of a reference, so there are no assumptions made about the spatial relationships of contigs.  
\

\

---
## Generate final blob and snail plots for final cleaned assembly
\

### BLAST final assembly (v2.12.0):
```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=96:00:00 --exclusive
ssh $SLURM_NODELIST


module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/blast


#set path for blast databases
export PATH=$PATH:/hb/groups/bernardi_lab/programs/blast
export BLASTDB=/hb/groups/bernardi_lab/programs/blobtools2/nt

#unzip fasta
gzip -d BFR_ref_final.fasta.gz

#run blast
blastn -db nt -query /hb/home/jatoy/bfren_genome/final_assembly/BFR_ref_final.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-20 -out /hb/home/jatoy/bfren_genome/final_assembly/BFR_ref_final_blast.out -num_threads 18

# qseqid = Query Seq-id
# staxids = unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
# std = default format specifiers = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
```
\ 

### Map Illumina short reads once more to your final assembly for use in Blobtools2
\  

#### BWA (Burrows-Wheeler Alignment)
First, create an index for the final assembly (`BFR_ref_final.fasta`) using the bwa index command (v0.7.17-r1188):
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/bwa

bwa index ~/bfren_genome/final_assembly/BFR_ref_final.fasta
```
\  

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa.  
\  

Now map Illumina short reads to assembly:
```{bash eval = FALSE}
bwa mem -t 24 ~/bfren_genome/final_assembly/BFR_ref_final.fasta ~/bfren_genome/illumina_reads/BFR_ref_sr_1_trim_paired.fq.gz ~/bfren_genome/illumina_reads/BFR_ref_sr_2_trim_paired.fq.gz > ~/bfren_genome/final_assembly/BFR_ref_sr_bwa_aligned_final.sam
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -Sb -@ 24 -O BAM -o BFR_ref_sr_bwa_aligned_final.bam BFR_ref_sr_bwa_aligned_final.sam

# S = input format is auto-detected
# b = output format BAM
# @ = number of threads
# O = specify output format
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_sr_bwa_aligned_sorted_final.bam -O BAM -@ 24 BFR_ref_sr_bwa_aligned_final.bam
```
\  

Calculate average depth of coverage of short reads from .bam file:
```{bash eval = FALSE}
samtools depth BFR_ref_sr_bwa_aligned_sorted_final.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

```
Average = 47.0181
```
\  

Get summary of coverage and depth per chromosome:
```{bash eval = FALSE}
samtools coverage BFR_ref_sr_bwa_aligned_sorted_final.bam -o BFR_ref_sr_bwa_aligned_sorted_final_coverage_summary
```

Run BUSCO on final assembly (should be the same as the last BUSCO but we will run it again anyway:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i /hb/home/jatoy/bfren_genome/final_assembly/BFR_ref_final.fasta -o busco_final -l actinopterygii_odb10 -m genome --cpu 24
```
\  
  
\  
  
## Run Blobtools2 (v3.0.0) to get final plots
\  

Start an interactive session:
``` {bash eval = FALSE}

salloc --partition=128x24 --nodes=1 --time=48:00:00 --exclusive --cpus-per-task=8
ssh $SLURM_NODELIST
```
\  

Now create your blobtools database and populate it with the coverage data, taxonomic information, and BUSCO scores:
``` {bash eval = FALSE}
#------- load programs --------------
module load miniconda3.9
conda activate /hb/groups/bernardi_lab/programs/blobtools2

#------- Run Commands ---------------
#Assign source directory variable to environment
SRC=/hb/home/jatoy/bfren_genome/

#Create a BlobDir
blobtools add --create --fasta $SRC/final_assembly/BFR_ref_final.fasta $SRC/blobtools2/BFR_ref_final_blob

#Add BLAST hits
blobtools add --threads 8 --hits $SRC/final_assembly/BFR_ref_final_blast.out --taxrule bestsumorder --taxdump /hb/groups/bernardi_lab/programs/blobtools2/taxdump $SRC/blobtools2/BFR_ref_final_blob

#Add mapping coverage
blobtools add --threads 8 --cov $SRC/final_assembly/BFR_ref_sr_bwa_aligned_sorted_final.bam $SRC/blobtools2/BFR_ref_final_blob

#Add BUSCO
blobtools add --threads 8 --busco $SRC/final_assembly/busco_final/run_actinopterygii_odb10/full_table.tsv $SRC/blobtools2/BFR_ref_final_blob
```
\  

Download BFR_ref_pilon_2_blob folder to desktop and run blobtools view.

``` {bash eval = FALSE}
scp -r jatoy@hb.ucsc.edu:~/bfren_genome/blobtools2/BFR_ref_final_blob/ ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools
```
\  

Create interactive html page with the blobtools results:

```{bash eval = FALSE}
#activate conda env
conda activate blobtools2

#change to desired directory
cd ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools

#run view command
blobtools view --remote BFR_ref_final_blob
```

Save all figures and data tables.  
\

\

Final assembly blob plot (post-cleaning):
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/final_assembly/BFR_ref_final_blob.blob.circle.png)  
\  

\

Final assembly snailplot (post-cleaning):
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/final_assembly/BFR_ref_final_blob.snail.png)
\  

\  

