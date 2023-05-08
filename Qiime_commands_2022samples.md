Qiime2 commands for 2022 samples
================
Gemma Clucas
2023-05-08

## 1. Import the data

It is saved on my solid state hardrive in the ANML folder.

    cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/MountainChickadees/2022
    conda activate qiime2-2021.4

    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/ANML/LaurenWhitenack_MountainChickadees/Plate34/reads  \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_plate1.qza
      
    for K in {1..1}; do
      qiime demux summarize \
        --i-data demux_Plate$K.qza \
        --o-visualization demux_Plate$K.qzv
    done

## 2. Trim primers using cutadapt

The primer sequences are:

ANML_F GGTCAACAAATCATAAAGATATTGG  
ANML_R GGWACTAATCAATTTCCAAATCC

### Trim 3’ ends first

At the 3’ end of the read, the primer will have been read through after
reading the ANML region. I need to be looking for the reverse complement
of the reverse primer in read 1 (—p-adapter-f) and the reverse
complement of the forward primer in R2 (—p-adapter-r).

R primer reverse complement: GGATTTGGAAATTGATTAGTWCC  
F primer reverse complement: CCAATATCTTTATGATTTGTTGACC

    for K in {1..1}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences demux_Plate$K.qza \
        --p-adapter-f GGATTTGGAAATTGATTAGTWCC \
        --p-adapter-r CCAATATCTTTATGATTTGTTGACC \
        --o-trimmed-sequences trimd_Plate$K.qza \
        --verbose > cutadapt_out_Plate$K.txt
    done

To see how much data passed the filter:

    grep "Total written (filtered):" cutadapt_out_Plate1.txt 

About 82% seems to have passed the filters for each sample, which is the
same as in the BTBW samples we sequenced.

### Trim 5’ ends of reads

All R1 should begin with the forward primer: GGTCAACAAATCATAAAGATATTGG
(25 bases).  
All R2 should begin with the reverse primer: GGWACTAATCAATTTCCAAATCC (23
bases).

Trim these with the following commands:

    for K in {1..1}; do
      qiime cutadapt trim-paired \
        --i-demultiplexed-sequences trimd_Plate$K.qza \
        --p-front-f GGTCAACAAATCATAAAGATATTGG \
        --p-front-r GGWACTAATCAATTTCCAAATCC  \
        --o-trimmed-sequences trimd2_Plate$K.qza \
        --verbose > cutadapt_out2_Plate$K.txt
    done

To see how much data passed the filter for each sample:

    grep "Total written (filtered):" cutadapt_out2_Plate1.txt 

88% this time.

## 3. Denoise with Dada2 START HERE

Using the same settings that seemed to work well for the BTBW samples.

    for K in {1..1}; do
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs trimd2_Plate$K.qza \
        --p-trunc-len-f 130 \
        --p-trunc-len-r 130 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-min-overlap 50 \
        --p-n-threads 8 \
        --o-representative-sequences rep-seqs_Plate$K \
        --o-table table_Plate$K \
        --o-denoising-stats denoise_Plate$K
    done


    for K in {1..1}; do  
      qiime metadata tabulate\
        --m-input-file denoise_Plate$K.qza\
        --o-visualization denoise_Plate$K.qzv
    done