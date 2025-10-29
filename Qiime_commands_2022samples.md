Qiime2 commands for 2022 and 2023 MOCH samples
================
Gemma Clucas
2023-05-08

## 1. Import the data

Just one plate for the 2022 samples.

    cd /Users/gemmaclucas/GitHub/Fecal_metabarcoding/Mountain_Chickadee_diets/2022
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

### 2023 samples

The data for both plates of samples from 2023 was combined into one
sequencing run at UNH, so it’s all in one folder.

    /Users/gc547/Dropbox/GitHub_copied/Fecal_metabarcoding/Mountain_Chickadee_diets/2023
    conda activate qiime2-amplicon-2024.10

    qiime tools import\
      --type 'SampleData[PairedEndSequencesWithQuality]'\
      --input-path /Volumes/Data_SS1/ANML/LaurenWhitenack_MountainChickadees/Plate1and2_2023/reads  \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt\
      --output-path demux_plate1.qza
      
    for K in {1..1}; do
      qiime demux summarize \
        --i-data demux_Plate$K.qza \
        --o-visualization demux_Plate$K.qzv
    done

Blanks look very clean. Lots of variability in read numbers between
samples though.

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

88% this time - normal.

### 2023 Samples

Exact same code as above. 82% kept after the first trim, 88% after the
second, as before.

## 3. Denoise with Dada2

Using the same settings that seemed to work well for the BTBW samples.
Truncating the reads to 130bp avoided reads being lost to the
low-quality filter (because read quality dropped after about 160bp as we
got to the end of the reads) and given that the amplicons are usually
about 185bp in length, we would actually expect an overlap of 75bp on
average, so this is plenty to ensure proper merging.

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

To view the rep-seqs and table:

    qiime feature-table tabulate-seqs \
      --i-data rep-seqs_Plate1.qza \
      --o-visualization rep-seqs_Plate1
      
    qiime feature-table summarize \
        --i-table table_Plate1.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_Plate1

### 2023 Samples

Exact same code as above. Still some reads in one of the extraction
blanks, but the number is very low compared to the samples.

## 4. Assign taxonomy

Using the same COI database and classifier that Devon O’Rourke trained.
Details
(here)\[<https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129>\].

I downloaded the database sequences (bold_anml_seqs.qza), the taxonomy
strings that go with each sequence (bold_anml_taxa.qza) and the trained
classifier (bold_anml_classifier.qza) into this folder. Note this only
runs with the 2019 version of QIIME, which is what he must have been
running when he made the classifier.

    conda activate qiime2-2019.10

    qiime feature-classifier classify-sklearn \
      --i-classifier bold_anml_classifier.qza \
      --i-reads rep-seqs_Plate1.qza \
      --o-classification taxonomy.qza
      
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv

### 2023 Samples

I couldn’t use the original classifier because I can’t run the 2019
version of QIIME on this new M3 laptop (and it’s not available on the
cluster either). So I trained a classifier using the same database of
sequences that Devon had trained the previous classifier on
(bold_anml_seqs.qza and bold_anml_taxa.qza), so that they are
comparable, and this works now.

    qiime feature-classifier classify-sklearn \
      --i-classifier classifier.qza \
      --i-reads rep-seqs_Plate1.qza \
      --o-classification taxonomy.qza
      
    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv

## 5. Make some barplots

Just to see what’s in the samples.

    qiime taxa barplot \
      --i-table table_Plate1.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_before_filtering.qzv

### 2023 Samples

Same as above. Massive diversity of insects. Also some rodent and bear
DNA!

## 6. Remove non-arthropod reads and drop mock community and blanks

    qiime taxa filter-table \
      --i-table table_Plate1.qza \
      --i-taxonomy taxonomy.qza \
      --p-include Arthropoda \
      --o-filtered-table table_Arthropoda.qza
      
    qiime feature-table summarize \
        --i-table table_Arthropoda.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_Arthropoda

After removing non-arthropod reads, one of the blanks (34-3) had some
arthropod reads, but it’s only 14 reads so it is very minimal.

I am also going to remove reads that were only classified to Arthropoda
and not beyond, since they are not giving us any useful information. So
that means I can use the `--p-include` flag to include anything that
also has a class designation, since they all have `Arthropoda;c_` in the
taxonomy string.

    qiime taxa filter-table \
      --i-table table_Arthropoda.qza \
      --i-taxonomy taxonomy.qza \
      --p-include "Arthropoda;c_" \
      --o-filtered-table table_ArthropodaC.qza
      
    qiime feature-table summarize \
        --i-table table_ArthropodaC.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_ArthropodaC

Drop the blanks and mock community and then filter the rep-seqs for just
the species that are in the table_ArthropodaC, so that I can send that
and sequence counts to Lauren, so that she has the unrarefied data.

    qiime feature-table filter-samples \
      --i-table table_ArthropodaC.qza \
      --m-metadata-file metadata.txt \
      --p-where "Species='MOCH'" \
      --o-filtered-table table_ArthropodaC_MOCH.qza

    qiime feature-table filter-seqs \
      --i-data rep-seqs_Plate1.qza \
      --i-table table_ArthropodaC_MOCH.qza \
      --o-filtered-data rep-seqs_ArthropodaC_MOCH.qza

    qiime metadata tabulate \
      --m-input-file rep-seqs_ArthropodaC_MOCH.qza \
      --m-input-file taxonomy.qza \
      --o-visualization sequence_taxonomy_2022ArthropodaMOCH.qzv

Also make a barplot to download the final dataset as a csv.

    qiime taxa barplot \
      --i-table table_ArthropodaC_MOCH.qza \
      --m-metadata-file metadata.txt \
      --i-taxonomy taxonomy.qza \
      --o-visualization barplot_ArthropodaC_MOCH

### 2023 Samples

Same as above but going straight to the `--p-include "Arthropoda;c_"`
step as the previous step (filtering for just arthropods) was redundant
with the addition of this one. BLANK-MOCH-1-1 has 4 reads, while
BLANK-MOCH-1-2 has 134 reads, so this is nothing compared to the average
depth of the samples.

    qiime taxa filter-table \
      --i-table table_Plate1.qza \
      --i-taxonomy taxonomy.qza \
      --p-include "Arthropoda;c_" \
      --o-filtered-table table_ArthropodaC.qza
      
    qiime feature-table summarize \
        --i-table table_ArthropodaC.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_ArthropodaC
