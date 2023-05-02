---
title: "Mitochondrial Genome Sequence Extraction and Nuclear Genome Correction"
author: "Jason Toy"
date: "7/5/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Run BLAST against EJA mitogenome and MitoFish database (downloaded from website) to get an idea of how many contigs may have mitochondrial contamination
\  

```{bash eval = FALSE}
makeblastdb -in /hb/groups/bernardi_lab/jason/EJA_ccgp_genome/Embiotoca_jacksoni/fEmbJac1.mitochondria.fasta.gz -dbtype nucl


blastn -query /hb/home/jatoy/bfren_genome/final_assembly/BFR_ref_final.fasta -db /hb/groups/bernardi_lab/jason/EJA_ccgp_genome/Embiotoca_jacksoni/fEmbJac1.mitochondria.fasta -task megablast -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 98.6 -num_threads 20 -soft_masking true -outfmt 7 > BFR_ref_to_EJA_mito.blastout
```
\  

```{bash eval = FALSE}
makeblastdb -in /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/mitofish_db_complete_partial_mitogenomes.fa -dbtype nucl

blastn -query /hb/home/jatoy/bfren_genome/final_assembly/BFR_ref_final.fasta -db /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/mitofish_db_complete_partial_mitogenomes.fa -task megablast -word_size 28 -best_inal_EJA_mthit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -evalue 0.0001 -perc_identity 98.6 -num_threads 20 -soft_masking true -outfmt 7 > BFR_ref_to_mitofish_db.blastout
```
\  

\  

## Geneious Assembly
\  

1. Map trimmed illumina forward and reverse reads and trimmed ONT reads (>500 bp) to EJA mitogenome using minimap2:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/minimap2


minimap2 -ax sr -t 24 -L fEmbJac1.mitochondria.fasta /hb/home/jatoy/bfren_genome/illumina_reads/BFR_ref_sr_1_trim_paired.fq.gz > BFR_sr_1_to_EJA_mt_mmap.sam

minimap2 -ax sr -t 24 fEmbJac1.mitochondria.fasta /hb/home/jatoy/bfren_genome/illumina_reads/BFR_ref_sr_2_trim_paired.fq.gz > BFR_sr_2_to_EJA_mt_mmap.sam

minimap2 -ax map-ont -t 24 fEmbJac1.mitochondria.fasta /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc_500.fastq.gz > BFR_ont_to_EJA_mt_mmap.sam
```
Where the first file is the template genome, the second is the fasta file for the reads you want to map.
\  

2. Remove non-mapped reads from sam file and convert to bam:
```{bash eval = FALSE}
module load samtools

samtools view -F 4 -b BFR_sr_1_to_EJA_mt_mmap.sam > BFR_sr_1_to_EJA_mt_mmapped.bam
samtools view -F 4 -b BFR_sr_2_to_EJA_mt_mmap.sam > BFR_sr_2_to_EJA_mt_mmapped.bam
samtools view -F 4 -b BFR_ont_to_EJA_mt_mmap.sam > BFR_ont_to_EJA_mt_mmapped.bam
```

3. Convert the bam file to fasta format:
```{bash eval = FALSE}
samtools fasta BFR_sr_1_to_EJA_mt_mmapped.bam > BFR_sr_1_to_EJA_mt_mmapped.fasta
samtools fasta BFR_sr_2_to_EJA_mt_mmapped.bam > BFR_sr_2_to_EJA_mt_mmapped.fasta
samtools fasta BFR_ont_to_EJA_mt_mmapped.bam > BFR_ont_to_EJA_mt_mmapped.fasta
```
\  

4. Download mapped read fasta files and EJA mitogenome file and import into Geneious
5. Using Geneious, map ONT reads to the EJA reference using "Map to Reference" function (default parameters). Extract consensus sequence trimmed to the beginning and end of the reference sequence.
6. Extract only paired reads from the mapped forward and reverse illumina reads (Sequence -> Set Paired Reads).
7. Map paired illumina reads to the ONT consensus sequence (Map to Reference). Extract consensus sequence trimmed to the beginning and end of the reference sequence.
8. Download EJA mitogenome reference with annotations (NC_029362) from "NCBI" folder in Geneious
9. Map BFR draft mitogenome to NC_029362. Set NC_029362 as the reference sequence and then transfer the annotations to the BFR sequence ("Live Annotate & Predict" Tab -> "Transfer Annotations").
10. Examine the new annotated sequence and look for frame shifts and stop codons in coding sequences (indicated by *). I identified 3 stop codons in coding sequences that resulted from spurious frame shifts introduced by ONT sequencing errors at mononucleotide repeats. These were manually corrected based off ONT read evidence (I mapped the ONT reads to these sequences):
        - Added 6th "C" to 5xC repeat in ND2 gene (site 4807)
        - Added 7th "G" to 6xG repeatWe in COXI gene (site 6157)
        - Added 6th "G" to 5xG repeat in ND4 gene (site 11,123)
        - There was also one "N" in the mitochondrial sequence. It was a similar case where Nanopore reads had sometimes called a 6th C at the end of a repeat and sometimes did not, leading to a consensus-called an N. There were some Illumina reads covering this site when mapped to the mitogenome, and each of them called a 6th "C", so the "N" was manually replaced with a "C"
11. Export corrected mitogenome sequence as fasta.
12. Upload fasta sequence to MitoAnnotator website. Run annotation pipeline. Download annotation results files.
\  


## Mitochondrial sequence removal from Nuclear genome assembly
\  

1. In Geneious, BLAST mitogenome (querry) against nuclear genome assembly (target database, BFR_ref_final.fasta)

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_mitogenome/BFR_mt_genome_to_BFR_ref_final_megablast.png)

2. I filtered the list of contigs that had BLAST hits by considering for removal only those contigs where BLAST hits to mitochondrial sequence made up >20% of the contig length (criteria/cutoffs may change based on your judgment). All other hits were assumed to be potential NUMTs. 
3. Remove mitochondrial sequences from Contigs. I removed mitochondrial sequence from a total of four small contigs (Contig 771, Contig 842, Contig 917, Contig 965; largest=11,036 bp). In each case, the validity of the remaining contig fragments was assessed by mapping long read data to each fragment. Fragments or portions of fragments with clear support from long read data were kept as new contigs. In total, this process led to the complete removal of two contigs (Contigs 771 and 842) from the nuclear assembly, the splitting of one contig into two (Contig 965), and the trimming of another (Contig 917).
\  
To map ONT reads to the contig fragments, I needed to create smaller subsets of reads that mapped to the original contigs. I did this using minimap2 in the same way that I did for the mitogenome assembly.
\  

### Contig 771
``` {bash eval = FALSE}
module load samtools

# Extract Contig_771 from BFR assembly
samtools faidx BFR_ref_final.fasta Contig_771 -o Contig_771.fasta 


module load miniconda3.9
conda activate /hb/groups/bernardi_lab/programs/minimap2

# Map ONT reads to Contig_771
minimap2 -ax map-ont -t 24 Contig_771.fasta /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc_500.fastq > BFR_ont_to_ctg_771_mmap.sam

# Remove non-mapped reads from sam file and convert to bam
samtools view -F 4 -b BFR_ont_to_ctg_771_mmap.sam > BFR_ont_to_ctg_771_mmapped.bam

# Convert the bam file to fasta format
samtools fasta BFR_ont_to_ctg_771_mmapped.bam > BFR_ont_to_ctg_771_mmapped.fasta
```
\  

### Contig 917
``` {bash eval = FALSE}
module load samtools

# Extract Contig_917 from BFR assembly
samtools faidx BFR_ref_final.fasta Contig_917 -o Contig_917.fasta 


module load miniconda3.9
conda activate /hb/groups/bernardi_lab/programs/minimap2

# Map ONT reads to Contig_917
minimap2 -ax map-ont -t 24 Contig_917.fasta /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc_500.fastq > BFR_ont_to_ctg_917_mmap.sam

# Remove non-mapped reads from sam file and convert to bam
samtools view -F 4 -b BFR_ont_to_ctg_917_mmap.sam > BFR_ont_to_ctg_917_mmapped.bam

# Convert the bam file to fasta format
samtools fasta BFR_ont_to_ctg_917_mmapped.bam > BFR_ont_to_ctg_917_mmapped.fasta
```
\  

### Contig 965
``` {bash eval = FALSE}
module load samtools

# Extract Contig_965 from BFR assembly
samtools faidx BFR_ref_final.fasta Contig_965 -o Contig_965.fasta 


module load miniconda3.9
conda activate /hb/groups/bernardi_lab/programs/minimap2

# Map ONT reads to Contig_965
minimap2 -ax map-ont -t 24 Contig_965.fasta /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc_500.fastq > BFR_ont_to_ctg_965_mmap.sam

# Remove non-mapped reads from sam file and convert to bam
samtools view -F 4 -b BFR_ont_to_ctg_965_mmap.sam > BFR_ont_to_ctg_965_mmapped.bam

# Convert the bam file to fasta format
samtools fasta BFR_ont_to_ctg_965_mmapped.bam > BFR_ont_to_ctg_965_mmapped.fasta
```
\  

I did not need to do this for Contig 842 because the entire contig was identified as mitochondrial sequence, so the whole contig was removed.
\  

In Contig 771, basepairs 1:9968 were identified as mitochondrial sequence, and the remaining fragment from basepair 9969:11036 had no ONT read support, so the whole contig was removed.
\  

In Contig 917, basepairs 1:1471 and 2368:5110 were identified as mitochondrial sequence. There was some level of ONT read support for the remaining fragment in the middle (bp 1472:2367), so this sequence was kept as a contig.

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_mitogenome/Congtig_917_fragment_ont_mapped_to_contig_917_extraction.png)
\  

In Contig 965, basepairs 1986:2779 were identified as mitochondrial sequence. This left two contig fragments, one from 1:1985 and one from 2780:3517. There was some level of ONT read support for the second fragment so this was kept in it's entirety as a new contig. The first fragment had ONT read support for only basepairs 1:597 (there was a clear break point in sequence disagreement after site 597 in the alignment of the ONT reads to the contig fragment), so only this portion of the fragment was kept as a new contig.

![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_mitogenome/Congtig_965_break_point_1_ont_mapped_to_contig_965_extraction_1.png)
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_mitogenome/Congtig_965_break_point_3_ont_mapped_to_contig_965_extraction_2.png)
\  

3. In Geneious, make a copy of the BFR_ref_final.fasta assembly. Remove and insert contigs from the assembly using the "Extract Sequences from List" and "Group Sequences into a List" functions (Sequences menu) to create the final, corrected genome assembly, both with the mitogenome included as a contig (named mtDNA) and without.


## Generate final blob and snail plots for mitochondria-extracted full assembly
\

### Upload MT-extracted assembly (without the MT genome) fasta and the mitogenome fasta to the hb cluster:
```{bash eval = FALSE}
scp ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_mitogenome/BFR_final_mt_corrected_no_mt.fasta jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly

scp ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_mitogenome/BFR_mitogenome_final.fasta jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly
```
\  

Using `nano`, manually rename the header for the mitogenome to "mtDNA":
```{bash eval = FALSE}
nano BFR_mitogenome_final.fasta
```

## Reorder contigs in nuclear assembly by length

We will use `seqkit` (v2.1.0) to order contigs by length and rename

```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/seqkit

#order contigs by length
seqkit sort --by-length --reverse --line-width 0 BFR_final_mt_corrected_no_mt.fasta > BFR_final_mt_corrected_no_mt_ordered.fasta

#rename contigs with new order numbers
seqkit replace --line-width 0 --pattern '.+' --replacement 'Contig_{nr}' BFR_final_mt_corrected_no_mt_ordered.fasta > BFR_final_mt_corrected_no_mt_ordered_renamed.fasta

#where {nr} is the record number, starting from 1
#"." in regex means "any character except line break"
#"+" in regex means "occurring one or more times"
#so "replace --pattern '.+'" basically means "replace all characters in the contig name
#--line-width 0 keeps the output as a single-line fasta
```
\  

After sorting by length:
```{bash eval = FALSE}
grep ">" BFR_final_mt_corrected_no_mt_ordered.fasta | tail
```

```
>Contig_997
>Contig_998
>Contig_999
>Contig_1000
>Contig_1001
>Contig_1002
>Contig_1003
>Contig_917_fragment
>Contig_965_fragment_2
>Contig_965_fragment_1 
```
\  

After renaming:
```{bash eval = FALSE}
grep ">" BFR_final_mt_corrected_no_mt_ordered_renamed.fasta | tail
```

```
>Contig_993
>Contig_994
>Contig_995
>Contig_996
>Contig_997
>Contig_998
>Contig_999
>Contig_1000
>Contig_1001
>Contig_1002
```
\  

Now add the mitogenome as a sequence to the end of the nuclear genome fasta:
```{bash eval = FALSE}
cat BFR_final_mt_corrected_no_mt_ordered_renamed.fasta BFR_mitogenome_final.fasta > BFR_ref_nuc_mt_final.fasta 
```
\  

The final genome assembly created is named `BFR_ref_nuc_mt_final.fasta`. Move it into it's own directory called `./final_fasta`
\  

```{bash eval = FALSE}
mkdir final_fasta
mv BFR_ref_nuc_mt_final.fasta ./final_fasta
```
\  

Get final assembly stats (should be similar to the last one, since we have the same number of contigs (1003) and only removed ~8,000 bp of sequence):
```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/assembly_stats/

cd ./final_fasta

assembly-stats BFR_ref_nuc_mt_final.fasta
```

```
stats for BFR_ref_nuc_mt_final.fasta
sum = 595951806, n = 1003, ave = 594169.30, largest = 12326090
N50 = 2589815, n = 55
N60 = 1896895, n = 82
N70 = 1185759, n = 122
N80 = 758445, n = 184
N90 = 406612, n = 292
N100 = 597, n = 1003
N_count = 0
Gaps = 0   
```
\  


## Final Genome Assembly Evaluation
\  

### BLAST MT-extracted final assembly (v2.12.0):
```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=96:00:00 --exclusive --mem=0
ssh $SLURM_NODELIST


module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/blast


#set path for blast databases
export PATH=$PATH:/hb/groups/bernardi_lab/programs/blast
export BLASTDB=/hb/groups/bernardi_lab/programs/blobtools2/nt

#run blast
blastn -db nt -query /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_ref_nuc_mt_final.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-20 -out /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/blobtools/BFR_ref_nuc_mt_final_blast.out -num_threads 18

# qseqid = Query Seq-id
# staxids = unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
# std = default format specifiers = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
```
\ 

### Map Illumina short reads to the MT-extracted assembly for use in Blobtools2
\  

#### BWA (Burrows-Wheeler Alignment)
First, create an index for the final assembly (`BFR_final_mt_corrected.fasta`) using the bwa index command (v0.7.17-r1188):
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/bwa

bwa index /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_ref_nuc_mt_final.fasta
```
\  

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa.  
\  

Now map Illumina short reads to assembly:
```{bash eval = FALSE}
bwa mem -t 24 /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_ref_nuc_mt_final.fasta ~/bfren_genome/illumina_reads/BFR_ref_sr_1_trim_paired.fq.gz ~/bfren_genome/illumina_reads/BFR_ref_sr_2_trim_paired.fq.gz > /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/blobtools/BFR_ref_nuc_mt_final_sr_bwa_aligned.sam
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -Sb -@ 24 -O BAM -o BFR_ref_nuc_mt_final_sr_bwa_aligned.bam BFR_ref_nuc_mt_final_sr_bwa_aligned.sam

# S = input format is auto-detected
# b = output format BAM
# @ = number of threads
# O = specify output format
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_nuc_mt_final_sr_bwa_aligned_sorted.bam -O BAM -@ 24 BFR_ref_nuc_mt_final_sr_bwa_aligned.bam
```
\  

Calculate average depth of coverage of short reads from .bam file:
```{bash eval = FALSE}
samtools depth BFR_ref_nuc_mt_final_sr_bwa_aligned_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

```
Average = 47.0657
```
\  

Get summary of coverage and depth per chromosome:
```{bash eval = FALSE}
samtools coverage BFR_ref_nuc_mt_final_sr_bwa_aligned_sorted.bam -o BFR_ref_nuc_mt_final_sr_bwa_aligned_sorted_coverage_summary
```
\  

### Run BUSCO on final assembly

It should be the same as the last BUSCO run, but we will run it again anyway:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_ref_nuc_mt_final.fasta -o busco_BFR_ref_nuc_mt_final -l actinopterygii_odb10 -m genome --cpu 24
```
\  
```
        --------------------------------------------------
        |Results from dataset actinopterygii_odb10        |
        --------------------------------------------------
        |C:98.1%[S:97.4%,D:0.7%],F:0.5%,M:1.4%,n:3640     |
        |3572   Complete BUSCOs (C)                       |
        |3547   Complete and single-copy BUSCOs (S)       |
        |25     Complete and duplicated BUSCOs (D)        |
        |20     Fragmented BUSCOs (F)                     |
        |48     Missing BUSCOs (M)                        |
        |3640   Total BUSCO groups searched               |
        --------------------------------------------------
```
  
\  
  
### Run Blobtools2 (v3.0.0) to get final plots
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
SRC=/hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/

#Create a BlobDir
blobtools add --create --fasta $SRC/BFR_ref_nuc_mt_final.fasta $SRC/blobtools/BFR_ref_nuc_mt_final_blob

#Add BLAST hits
blobtools add --threads 18 --hits $SRC/blobtools/BFR_ref_nuc_mt_final_blast.out --taxrule bestsumorder --taxdump /hb/groups/bernardi_lab/programs/blobtools2/taxdump $SRC/blobtools/BFR_ref_nuc_mt_final_blob

#Add mapping coverage
blobtools add --threads 18 --cov $SRC/blobtools/BFR_ref_nuc_mt_final_sr_bwa_aligned_sorted.bam $SRC/blobtools/BFR_ref_nuc_mt_final_blob

#Add BUSCO
blobtools add --threads 18 --busco $SRC/blobtools/busco_BFR_ref_nuc_mt_final/run_actinopterygii_odb10/full_table.tsv $SRC/blobtools/BFR_ref_nuc_mt_final_blob
```
\  

Download BFR_final_mt_corrected_blob folder to desktop and run blobtools view.

``` {bash eval = FALSE}
scp -r jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/blobtools/BFR_ref_nuc_mt_final_blob ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final
```
\  

Create interactive html page with the blobtools results:

```{bash eval = FALSE}
#activate conda env
conda activate blobtools2

#change to desired directory
cd ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final

#run view command
blobtools view --remote BFR_ref_nuc_mt_final_blob
```
\  

```
View dataset at http://localhost:8001/view/BFR_ref_nuc_mt_final_blob/dataset/BFR_ref_nuc_mt_final_blob/blob
```
\  

Save all figures and data tables.  
\

\

Final MT-corrected assembly blob plot (post-cleaning, post MT-correction):
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final/BFR_ref_nuc_mt_final_blob.blob.circle.png)  

\  

Final MT-corrected assembly snail plot (post-cleaning, post MT-correction):
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final/BFR_ref_nuc_mt_final_blob.snail.png)
\  

## Calculate Nanopore long read coverage

Map ONT long reads to assembly:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/minimap2

minimap2 -ax map-ont -t 24 -o BFR_ref_nuc_mt_final_ont_mmap.sam --sam-hit-only /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_ref_nuc_mt_final.fasta /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc_500.fastq.gz
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -F 4 -b BFR_ref_nuc_mt_final_ont_mmap.sam > BFR_ref_nuc_mt_final_ont_mmap.bam

# b = output format BAM
# @ = number of threads
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_nuc_mt_final_ont_mmap_sorted.bam -O BAM -@ 24 BFR_ref_nuc_mt_final_ont_mmap.bam
```
\  

Calculate average depth of coverage of short reads from .bam file:
```{bash eval = FALSE}
samtools depth BFR_ref_nuc_mt_final_ont_mmap_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
\  

```
Average =  26.5039
```
\  

Get summary of coverage and depth per chromosome:
```{bash eval = FALSE}
samtools coverage BFR_ref_nuc_mt_final_ont_mmap_sorted.bam -o BFR_ref_nuc_mt_final_ont_mmap_sorted_coverage_summary
```
\

\  

## Ragoo and Alignment to E. jacksoni Reference Genome

To align/scaffold your new genome to existing reference genomes, use RagTag and the fasta file for the reference genome. Before running, make sure the reference file is in a directory to which you have access.
```{bash eval = FALSE}
module load miniconda3

conda activate /hb/groups/bernardi_lab/programs/ragtag

module deactivate python-3.6.2 #I had to unload the default loaded python on hb to get this to work

#run ragoo
ragtag.py scaffold -t 20 /hb/groups/bernardi_lab/jason/EJA_ccgp_genome/Embiotoca_jacksoni/fEmbJac1.NCBI.p_ctg.fasta /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_ref_nuc_mt_final.fasta -o /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/final_fasta/BFR_EJA_ragtag

```
\  

By default, this will create a new directory called ```ragtag_output```, inside of which will be a new assembly file called ```ragtag.scaffold.fasta```. Use the -o flag to change the output directory name.
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
789                 594362913   214                   1588892       64800     648
```
\ 

Now map BFR_ref_nuc_mt_final to the EJA reference assembly using D-GENIES and look at the dot plot.
\  

## Split Contig 2 into two new Contigs
When comparing BFR_ref_nuc_mt_final to the EJA reference genome using D-GENIES, there is what appears to be a translocation. About 39% of Contig 2 (4.7 Mb) mapped to Scaffold 2 of the E. jacksoni assembly, while 61% (7.4 Mb) mapped to Scaffold 22. To investigate whether this was indeed a translocation or just a misassembly, I mapped ONT long-reads to Contig 2 to see if there was good read support for joining these two regions in the same contig in the first place. Contig 2 stops mapping to EJA Scaffold 2 at site 4,746,531 and starts mapping to Contig 22 at site 4,746,539, so I focused on this area. There were only 19 reads covering this point in the contig, and there was a clear breakpoint where reads spanning the junction at site 4,746,539 only mapped well to the contig on one side of the junction or the other. Because there seemed to be little long-read evidence for this junction, I split Contig 2 in between site 4,746,538 and 4,746,539 to create two new contigs (temporarily named Contig 2.1 and Contig 2.2). I then once a gain reran the summary analyses for the new assembly as below:
\  


## Final Genome Assembly Evaluation
\

### Upload MT-extracted assembly (without the MT genome) fasta and the mitogenome fasta to the hb cluster:
```{bash eval = FALSE}
scp /home/jason/winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/BFR_ref_nuc_final_v2.fasta jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly_v2
```
\  

## Reorder contigs in nuclear assembly by length

We will use `seqkit` (v2.1.0) to order contigs by length and rename

```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/seqkit

#order contigs by length
seqkit sort --by-length --reverse --line-width 0 BFR_ref_nuc_final_v2.fasta > BFR_ref_nuc_final_v2_ordered.fasta

#rename contigs with new order numbers
seqkit replace --line-width 0 --pattern '.+' --replacement 'Contig_{nr}' BFR_ref_nuc_final_v2_ordered.fasta > BFR_ref_nuc_final_v2_ordered_renamed.fasta

#where {nr} is the record number, starting from 1
#"." in regex means "any character except line break"
#"+" in regex means "occurring one or more times"
#so "replace --pattern '.+'" basically means "replace all characters in the contig name
#--line-width 0 keeps the output as a single-line fasta
```
\  

After sorting by length:
```{bash eval = FALSE}
grep ">" BFR_ref_nuc_final_v2_ordered.fasta | head -30
```

```
>Contig_1
>Contig_3
>Contig_4
>Contig_5
>Contig_6
>Contig_7
>Contig_8
>Contig_9
>Contig_10
>Contig_11
>Contig_12
>Contig_2.2
>Contig_13
>Contig_14
>Contig_15
>Contig_16
>Contig_17
>Contig_18
>Contig_19
>Contig_20
>Contig_21
>Contig_22
>Contig_23
>Contig_24
>Contig_25
>Contig_26
>Contig_2.1
>Contig_27
>Contig_28
>Contig_29
```
\  

After renaming:
```{bash eval = FALSE}
grep ">" BFR_ref_nuc_final_v2_ordered_renamed.fasta | head -30
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
>Contig_11
>Contig_12
>Contig_13
>Contig_14
>Contig_15
>Contig_16
>Contig_17
>Contig_18
>Contig_19
>Contig_20
>Contig_21
>Contig_22
>Contig_23
>Contig_24
>Contig_25
>Contig_26
>Contig_27
>Contig_28
>Contig_29
>Contig_30
```
\  

Now add the mitogenome as a sequence to the end of the nuclear genome fasta:
```{bash eval = FALSE}
cat BFR_ref_nuc_final_v2_ordered_renamed.fasta /hb/home/jatoy/bfren_genome/final_assembly/mt_genome/geneious_assembly/BFR_mitogenome_final.fasta > BFR_ref_nuc_mt_final_v2.fasta
```
\  

The final genome assembly created is named `BFR_ref_nuc_mt_final_v2.fasta`. Move it into it's own directory called `/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta`
\  

```{bash eval = FALSE}
mkdir final_fasta
mv BFR_ref_nuc_mt_final_v2.fasta ./final_fasta
```
\  

Get final assembly stats (should have 1004 contigs now instead of 1003, but should be the same total length):
```{bash eval = FALSE}
conda activate /hb/groups/bernardi_lab/programs/assembly_stats/

cd ./final_fasta

assembly-stats BFR_ref_nuc_mt_final_v2.fasta
```

```
stats for BFR_ref_nuc_mt_final_v2.fasta
sum = 595951806, n = 1004, ave = 593577.50, largest = 12326090
N50 = 2589815, n = 56
N60 = 1896895, n = 83
N70 = 1185759, n = 123
N80 = 758445, n = 185
N90 = 406612, n = 293
N100 = 597, n = 1004
N_count = 0
Gaps = 0
```
\  

## Generate final blob and snail plots for BFR_ref_nuc_mt_final_v2 full assembly
\  

### BLAST BFR_ref_nuc_mt_final_v2 assembly (v2.12.0):
```{bash eval = FALSE}
salloc --partition=128x24 --nodes=1 --time=96:00:00 --exclusive --mem=0
ssh $SLURM_NODELIST


module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/blast


#set path for blast databases
export PATH=$PATH:/hb/groups/bernardi_lab/programs/blast
export BLASTDB=/hb/groups/bernardi_lab/programs/blobtools2/nt

#run blast
blastn -db nt -query /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-20 -out /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/blobtools/BFR_ref_nuc_mt_final_v2_blast.out -num_threads 23

# qseqid = Query Seq-id
# staxids = unique Subject Taxonomy ID(s), separated by a ';'(in numerical order)
# std = default format specifiers = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
```
\ 

### Map Illumina short reads to the BFR_ref_nuc_mt_final_v2 assembly for use in Blobtools2
\  

#### BWA (Burrows-Wheeler Alignment)
First, create an index for the final assembly (`BFR_ref_nuc_mt_final_v2`) using the bwa index command (v0.7.17-r1188):
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/bwa

bwa index /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta
```
\  

bwa index will create new files (.amb, .ann, .bwt, .pac, .sa) that will need to be within the directory to map sequences using bwa.  
\  

Now map Illumina short reads to assembly:
```{bash eval = FALSE}
bwa mem -t 24 /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta ~/bfren_genome/illumina_reads/BFR_ref_sr_1_trim_paired.fq.gz ~/bfren_genome/illumina_reads/BFR_ref_sr_2_trim_paired.fq.gz > /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/blobtools/BFR_ref_nuc_mt_final_v2_sr_bwa_aligned.sam
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -Sb -@ 24 -O BAM -o BFR_ref_nuc_mt_final_v2_sr_bwa_aligned.bam BFR_ref_nuc_mt_final_v2_sr_bwa_aligned.sam

# S = input format is auto-detected
# b = output format BAM
# @ = number of threads
# O = specify output format
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_nuc_mt_final_v2_sr_bwa_aligned_sorted.bam -O BAM -@ 24 BFR_ref_nuc_mt_final_v2_sr_bwa_aligned.bam
```
\  

Calculate average depth of coverage of short reads from .bam file:
```{bash eval = FALSE}
samtools depth BFR_ref_nuc_mt_final_v2_sr_bwa_aligned_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

```
Average = 47.0657
```
\  

Get summary of coverage and depth per chromosome:
```{bash eval = FALSE}
samtools coverage BFR_ref_nuc_mt_final_v2_sr_bwa_aligned_sorted.bam -o BFR_ref_nuc_mt_final_v2_sr_bwa_aligned_sorted_coverage_summary
```
\  

### Run BUSCO on final assembly

It should be the same as the last BUSCO run, but we will run it again anyway:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/busco_v5.2.2

busco -i /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta -o busco_BFR_ref_nuc_mt_final_v2 -l actinopterygii_odb10 -m genome --cpu 24
```
\  
```
        --------------------------------------------------
        |Results from dataset actinopterygii_odb10        |
        --------------------------------------------------
        |C:98.1%[S:97.4%,D:0.7%],F:0.5%,M:1.4%,n:3640     |
        |3572   Complete BUSCOs (C)                       |
        |3547   Complete and single-copy BUSCOs (S)       |
        |25     Complete and duplicated BUSCOs (D)        |
        |20     Fragmented BUSCOs (F)                     |
        |48     Missing BUSCOs (M)                        |
        |3640   Total BUSCO groups searched               |
        --------------------------------------------------
```
  
\  
  
### Run Blobtools2 (v3.0.0) to get final plots
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
#Assign source directory variable in environment
SRC=/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/

#Create a BlobDir
blobtools add --create --fasta $SRC/BFR_ref_nuc_mt_final_v2.fasta $SRC/blobtools/BFR_ref_nuc_mt_final_v2_blob

#Add BLAST hits
blobtools add --threads 18 --hits $SRC/blobtools/BFR_ref_nuc_mt_final_v2_blast.out --taxrule bestsumorder --taxdump /hb/groups/bernardi_lab/programs/blobtools2/taxdump $SRC/blobtools/BFR_ref_nuc_mt_final_v2_blob

#Add mapping coverage
blobtools add --threads 18 --cov $SRC/blobtools/BFR_ref_nuc_mt_final_v2_sr_bwa_aligned_sorted.bam $SRC/blobtools/BFR_ref_nuc_mt_final_v2_blob

#Add BUSCO
blobtools add --threads 18 --busco $SRC/blobtools/busco_BFR_ref_nuc_mt_final_v2/run_actinopterygii_odb10/full_table.tsv $SRC/blobtools/BFR_ref_nuc_mt_final_v2_blob
```
\  

Download BFR_ref_nuc_mt_final_v2_blob folder to desktop and run blobtools view.

``` {bash eval = FALSE}
scp -r jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/blobtools/BFR_ref_nuc_mt_final_v2_blob ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final_v2
```
\  

Create interactive html page with the blobtools results:

```{bash eval = FALSE}
#activate conda env
conda activate blobtools2

#change to desired directory
cd ./winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final_v2

#run view command
blobtools view --remote BFR_ref_nuc_mt_final_v2_blob
```
\  

```
http://localhost:8001/view/BFR_ref_nuc_mt_final_v2_blob/dataset/BFR_ref_nuc_mt_final_v2_blob/blob
```
\  

Save all figures and data tables.  
\

\

Final v2 assembly blob plot (post-cleaning, post MT-correction):
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final_v2/BFR_ref_nuc_mt_final_v2_blob.blob.circle.png)  

\  

Final v2 assembly snail plot (post-cleaning, post MT-correction):
![](C:/Users/jason/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/blobtools/BFR_ref_nuc_mt_final_v2/BFR_ref_nuc_mt_final_v2_blob.snail.png)
\  

## Calculate Nanopore long read coverage

Map ONT long reads to assembly:
```{bash eval = FALSE}
module load miniconda3.9

conda activate /hb/groups/bernardi_lab/programs/minimap2

minimap2 -ax map-ont -t 24 -o BFR_ref_nuc_mt_final_v2_ont_mmap.sam --sam-hit-only /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta /hb/home/jatoy/bfren_genome/minion/runs/bfren_ref_bigfile_pc_500.fastq.gz
```
\  

Convert the SAM file to BAM format using samtools (v1.10):
```{bash eval = FALSE}
module load samtools

samtools view -F 4 -b BFR_ref_nuc_mt_final_v2_ont_mmap.sam > BFR_ref_nuc_mt_final_v2_ont_mmap.bam

# b = output format BAM
# @ = number of threads
# o = output file name
```
\  

Sort BAM file:
```{bash eval = FALSE}
samtools sort -o BFR_ref_nuc_mt_final_v2_ont_mmap_sorted.bam -O BAM -@ 24 BFR_ref_nuc_mt_final_v2_ont_mmap.bam
```
\  

Calculate average depth of coverage of long reads from .bam file:
```{bash eval = FALSE}
samtools depth BFR_ref_nuc_mt_final_v2_ont_mmap_sorted.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
\  

```
Average =  26.504
```
\  

Get summary of coverage and depth per chromosome:
```{bash eval = FALSE}
samtools coverage BFR_ref_nuc_mt_final_v2_ont_mmap_sorted.bam -o BFR_ref_nuc_mt_final_v2_ont_mmap_sorted_coverage_summary
```
\

\  

## Ragoo and Alignment to E. jacksoni Reference Genome

To align/scaffold your new genome to existing reference genomes, use RagTag and the fasta file for the reference genome. Before running, make sure the reference file is in a directory to which you have access.
```{bash eval = FALSE}
module load miniconda3

conda activate /hb/groups/bernardi_lab/programs/ragtag

module unload python-3.6.2 #I had to unload the default loaded python on hb to get this to work

#run ragoo
ragtag.py scaffold -t 20 /hb/groups/bernardi_lab/jason/EJA_ccgp_genome/Embiotoca_jacksoni/fEmbJac1.NCBI.p_ctg.fasta /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta -o /hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_EJA_ragtag

```
\  

By default, this will create a new directory called ```ragtag_output```, inside of which will be a new assembly file called ```ragtag.scaffold.fasta```. Use the -o flag to change the output directory name.
\  

Check to make sure it worked properly by printing the names of new scaffolds that should've been created in the last step:
```{bash eval = FALSE}
grep ">" ragtag.scaffold.fasta
```
\  

Get summary stats for scaffolding:
```{bash eval = FALSE}
cat ragtag.scaffold.stats
```

```
placed_sequences        placed_bp       unplaced_sequences      unplaced_bp     gap_bp  gap_sequences
790                     594362913       214                     1588892         64900   649
```
\  

## Download `BFR_ref_nuc_mt_final_v2.fasta` and map it to the EJA genome using D-GENIES
\  
```{bash eval = FALSE}
scp jatoy@hb.ucsc.edu:/hb/home/jatoy/bfren_genome/final_assembly_v2/final_fasta/BFR_ref_nuc_mt_final_v2.fasta ~/winhome/OneDrive/Documents/ucsc/projects/landscape_genomics/kelp_perch/reference_genome/
```
