Qiime2 commands for 2022 samples
================
Gemma Clucas
2023-05-08

## 1. Import the data

It is saved on my solid state hardrive in the ANML folder.

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

To view the rep-seqs and table:

    qiime feature-table tabulate-seqs \
      --i-data rep-seqs_Plate1.qza \
      --o-visualization rep-seqs_Plate1
      
    qiime feature-table summarize \
        --i-table table_Plate1.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_Plate1

## 4. Assign taxonomy

I am going to use the same COI database and classifier that Devon
O’Rourke trained. Details
(here)\[<https://github.com/GemmaClucas/Hubbard-Brook-Warbler-Diets#5-coi-database>\].

I have copied the database sequences (bold_anml_seqs.qza), the taxonomy
strings that go with each sequence (bold_anml_taxa.qza) and the trained
classifier (bold_anml_classifier.qza) into this folder.

    conda activate qiime2-2019.10

    qiime feature-classifier classify-sklearn \
      --i-classifier bold_anml_classifier.qza \
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

## 6. Remove non-arthropod reads

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
arthropod reads, but it’s only 14 reads so it will be dropped when we
rarefy.

## 7. Calculate alpha-rarefaction curves

First we have to collapse the taxonomy to the species level.

    qiime taxa collapse \
      --i-table table_Arthropoda.qza \
      --i-taxonomy taxonomy.qza \
      --p-level 7 \
      --o-collapsed-table table_Arthropoda_collapsed.qza

    qiime diversity alpha-rarefaction \
      --i-table table_Arthropoda_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 500 \
      --p-max-depth 30000 \
      --o-visualization alpha-rarefaction-500-30000
      
    qiime diversity alpha-rarefaction \
      --i-table table_Arthropoda_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 500 \
      --p-max-depth 10000 \
      --o-visualization alpha-rarefaction-500-10000

At this point I looked into alternatives to rarefying and came across
the [SRS
plugin](https://forum.qiime2.org/t/q2-srs-qiime2-plugin-for-library-size-normalization-by-scaling-with-ranked-subsampling-srs/17661),
which is an alternative (and apparently better?) than rarefying.

    qiime srs SRScurve \
      --i-table table_Arthropoda.qza \
      --p-step 100 \
      --p-max-sample-size 10000 \
      --p-rarefy-comparison \
      --p-rarefy-comparison-legend \
      --p-rarefy-repeats 100 \
      --p-srs-color 'blue' \
      --p-rarefy-color '#333333' \
      --o-visualization SRScurve-plot_10000.qzv \
      --verbose

Looking at the rarefaction curves, SRS plots, and the [SRS shiny
app](https://vitorheidrich.shinyapps.io/srsshinyapp/), I think I am
going to rarefy to a depth of 2722 reads, since that will retain 90% or
more of the diversity in all of the samples and minimise the number of
samples we have to drop. The next sample has 1856 reads, which would
mean losing a lot of taxa from each sample.

Note that for Kim’s 3 samples that I included on this plate, a higher
depth will be needed as those appear to be more diverse.

## 8. Normalise table

I am going to use SRS rather than rarefying for normalising the samples.

    qiime srs SRS \
      --i-table table_Arthropoda.qza \
      --p-c-min 2722 \
      --o-normalized-table table_normalised2722.qza \
      --verbose
      
    qiime feature-table summarize \
        --i-table table_normalised2722.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_normalised2722

## 9. Remake barplots

    qiime taxa barplot \
      --i-table table_normalised2722.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_normalised2722.qzv

I downloaded the CSV from here, to send to Lauren after editing.
