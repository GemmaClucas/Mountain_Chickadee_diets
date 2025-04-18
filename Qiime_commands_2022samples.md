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

### 2023 samples

The data for both plates was combined at UNH, so it’s all in one folder
on my hard drive SSD_Data1.

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

Blanks look good. Lots of variability in read numbers between samples
though.

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

### 2023 Samples

Exact same code as above. 82% kept after the first trim, 88% after the
second, as before.

## 3. Denoise with Dada2

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

### 2023 Samples

Exact same code as above. Still some reads in one of the extraction
blanks, but the number is low compared to the samples.

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

### 2023 Samples

I can’t use the old classifier because I can’t run the old version of
qiime on this laptop I don’t think (and it’s not available on the
cluster). But I trained a classifier using this version for Tasos and
Kim’s projects, training it with the database that Devon had trained the
previous classifier on, so that’s good. Kim noted that the BOLD database
is a lot larger now, so if I wanted to train a classifier on a newer
version of the database then I’d have to re-run both the 2022 and 2023
samples with that.

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

### 2023 Samples

I am also going to remove reads that were only classified to Arthropoda
and not beyond, since they are not giving us any useful information. So
that means I can use the `--p-include` flag to include anything that
also has a class designation, since they all have `Arthropoda;c_` in the
taxonomy string.

BLANK-MOCH-1-1 has 4 reads, while BLANK-MOCH-1-2 has 134 reads, so this
is nothing compared to the average depth of the samples.

    qiime taxa filter-table \
      --i-table table_Plate1.qza \
      --i-taxonomy taxonomy.qza \
      --p-include "Arthropoda;c_" \
      --o-filtered-table table_ArthropodaC.qza
      
    qiime feature-table summarize \
        --i-table table_ArthropodaC.qza \
        --m-sample-metadata-file metadata.txt \
        --o-visualization table_ArthropodaC

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
      --o-visualization alpha-rarefactionC-500-30000
      
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

### 2023 Samples

The alpha rarefaction curves look like there’s more diversity in these
samples compared to 2022, with a plateau in diversity at around 7000
reads. Looking at the SRS curves (just by uploading the collapsed
arthropodC table to the SRS shiny app) then if I rarefy to the same
depth as the 2022 samples (2722 reads) then I will retain 52% of the
samples and 84% of the global species richness. If I rarefy to 7000
reads, then I only retain 43 samples out of 112, but I retain 86% of the
global species richness. So it is a trade-off between samples and depth.
But do I need to rarefy at all? Only to calculate diversity stats.

    qiime taxa collapse \
      --i-table table_ArthropodaC.qza \
      --i-taxonomy taxonomy.qza \
      --p-level 7 \
      --o-collapsed-table table_ArthropodaC_collapsed.qza

    qiime diversity alpha-rarefaction \
      --i-table table_ArthropodaC_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 500 \
      --p-max-depth 30000 \
      --o-visualization alpha-rarefactionC-500-30000
      
    qiime diversity alpha-rarefaction \
      --i-table table_ArthropodaC_collapsed.qza \
      --m-metadata-file metadata.txt \
      --p-min-depth 500 \
      --p-max-depth 10000 \
      --o-visualization alpha-rarefactionC-500-10000

Reasons for increased diversity in 2023 samples: more insects available?
Or technical due to doing these extractions by hand?

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

Have not run for 2023 samples after removing reads that could not be
classified to class level. I think I will send Lauren the un-rarefied
dataset and ask her what she wants to do with it.

## 9. Remake barplots

    qiime taxa barplot \
      --i-table table_normalised2722.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_normalised2722.qzv

I downloaded the CSV from here, to send to Lauren after editing.

### 2023 Samples

Since I have decided not to normalise, as above, this is just the
barplot before rarefying.

    qiime taxa barplot \
      --i-table table_ArthropodaC.qza \
      --i-taxonomy taxonomy.qza \
      --m-metadata-file metadata.txt \
      --o-visualization barplot_ArthropodaC.qzv

## Can I output the table with the sequences?

    qiime metadata tabulate \
      --m-input-file rep-seqs_Plate1.qza \
      --m-input-file taxonomy.qza \
      --o-visualization tabulated-feature-metadata.qzv

Yes this works!

Do the same for the 2022 samples:

    cd ../2022/

    qiime metadata tabulate \
      --m-input-file rep-seqs_Plate1.qza \
      --m-input-file taxonomy.qza \
      --o-visualization tabulated-feature-metadata.qzv

## 10. What if I run the MiFish taxonomy classification on the sequences?

Checking for contamination…

    ../../AP_SG_Penguins_2022/MiFish/mktaxa_singlethreaded.py \
      ../../AP_SG_Penguins_2022/MiFish/ncbi-refseqs-withHuman.qza \
      ../../AP_SG_Penguins_2022/MiFish/ncbi-taxonomy-withHuman.qza \
      rep-seqs_Plate1.qza

The mktaxa script ran for a few minutes but then threw this error
`TypeError: 'NoneType' object is not subscriptable` which apparently
happens when a function returns `None` but you try to access an index or
key from its return value.

I’m going to try running this on the Mandarte data instead to see
whether it works. I’ll just put the files in this folder for now though
since I don’t have a repo for Mandarte.

    cd ../Mandarte_diet_2024/

    ../../AP_SG_Penguins_2022/MiFish/mktaxa_singlethreaded.py \
      ../../AP_SG_Penguins_2022/MiFish/ncbi-refseqs-withHuman.qza \
      ../../AP_SG_Penguins_2022/MiFish/ncbi-taxonomy-withHuman.qza \
      rep-seqs_Plate1.qza

This failed with the same error so does that mean just nothing was
classified? Maybe I should try with a version of qiime on the server and
the original mktaxa.py script. The original mktaxa script just hung with
the multithreaded version and the single threaded version quit with the
same error that I got running it on my laptop.

I went through the MOCH Unassigned reads and a bunch of the reads that
were only assigned to Animalia, and none of them were fish.

I should do the same for the Mandarte dataset though.
