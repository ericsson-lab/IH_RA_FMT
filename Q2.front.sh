#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=211210_hypoxia       # name for the job
#SBATCH --cpus-per-task=5              # number of cores
#SBATCH --mem=150G                       # total memory
#SBATCH --nodes 2
#SBATCH --time 2-00:00:00                 # time limit in the form days-hours:minutes
#SBATCH --mail-user=zlmg2b@umsystem.edu    # email address for notifications
#SBATCH --mail-type=FAIL,END,BEGIN           # email types
#SBATCH --partition Lewis            # max of 1 node and 4 hours; use `Lewis` for larger jobs
#--------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

module load miniconda3
source activate qiime2-2021.8

mkdir sequences
mkdir Dada2
mkdir Dada2/visualization
mkdir phylogeny
mkdir taxonomy
mkdir transfer
mkdir transfer/sequences
mkdir transfer/Dada2
mkdir transfer/taxonomy

cd demux_seqs

for file in *R1*.gz; do [ -f "$file" ] || continue; mv -vf "$file" "${file//_*R1_001.fastq.gz/_R1.fastq.gz}"; done
for file in *R2*.gz; do [ -f "$file" ] || continue; mv -vf "$file" "${file//_*R2_001.fastq.gz/_R2.fastq.gz}"; done

cd ..

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ../workflow/feature-classifier/silva-138-99-seqs-515-806.qza \
  --i-reference-taxonomy ../workflow/feature-classifier/silva-138-99-tax-515-806.qza \
  --o-classifier ../workflow/feature-classifier/silva-138-99-classifier-515-806.qza

qiime tools import \
  --type "SampleData[PairedEndSequencesWithQuality]" \
  --input-format PairedEndFastqManifestPhred33 \
  --input-path ./manifest.csv \
  --output-path ./sequences/demux_seqs.qza

qiime demux summarize \
  --i-data ./sequences/demux_seqs.qza \
  --o-visualization ./sequences/demux_seqs.qzv

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ./sequences/demux_seqs.qza \
  --p-cores $SLURM_CPUS_ON_NODE \
  --p-adapter-f 'ATTAGAWACCCBDGTAGTCC' \
  --p-front-f 'GTGCCAGCMGCCGCGGTAA' \
  --p-adapter-r 'TTACCGCGGCKGCTGGCAC' \
  --p-front-r 'GGACTACHVGGGTWTCTAAT' \
  --p-discard-untrimmed \
  --p-no-indels \
  --o-trimmed-sequences ./sequences/trimmed_demux_seqs.qza

qiime demux summarize \
  --i-data ./sequences/trimmed_demux_seqs.qza \
  --o-visualization ./sequences/trimmed_demux_seqs.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./sequences/trimmed_demux_seqs.qza \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 150 \
  --o-table ./Dada2/dada2_table.qza \
  --o-representative-sequences ./Dada2/dada2_rep_seqs.qza \
  --o-denoising-stats ./Dada2/dada2_stats.qza \
  --p-n-threads $SLURM_CPUS_ON_NODE

qiime metadata tabulate \
  --m-input-file ./Dada2/dada2_stats.qza  \
  --o-visualization ./Dada2/visualization/dada2_stats.qzv

qiime feature-table filter-seqs \
  --i-data Dada2/dada2_rep_seqs.qza \
  --m-metadata-file Dada2/dada2_rep_seqs.qza \
  --p-where 'length(sequence) >= 249 AND length(sequence) <= 257' \
  --o-filtered-data ./Dada2/dada2_rep_seqs_filtered.qza

qiime feature-table filter-features \
  --i-table Dada2/dada2_table.qza \
  --m-metadata-file Dada2/dada2_rep_seqs_filtered.qza \
  --o-filtered-table Dada2/dada2_table_filtered.qza

qiime feature-table summarize \
  --i-table ./Dada2/dada2_table_filtered.qza \
  --m-sample-metadata-file ./metadata.txt \
  --o-visualization ./Dada2/visualization/dada2_table_filtered.qzv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./Dada2/dada2_rep_seqs_filtered.qza \
  --o-alignment ./phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment ./phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree ./phylogeny/unrooted-tree.qza \
  --o-rooted-tree ./phylogeny/rooted-tree.qza

qiime feature-classifier classify-sklearn \
  --i-reads ./Dada2/dada2_rep_seqs_filtered.qza \
  --i-classifier ./silva-138-99-classifier-515-806.qza \
  --o-classification ./taxonomy/taxonomy.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./Dada2/dada2_rep_seqs_filtered.qza \
  --o-alignment ./phylogeny/aligned-rep-seqs.qza \
  --o-masked-alignment ./phylogeny/masked-aligned-rep-seqs.qza \
  --o-tree ./phylogeny/unrooted-tree.qza \
  --o-rooted-tree ./phylogeny/rooted-tree.qza

qiime feature-classifier classify-sklearn \
  --i-reads ./Dada2/dada2_rep_seqs_filtered.qza \
  --i-classifier ./silva-138-99-classifier-515-806.qza \
  --o-classification ./taxonomy/taxonomy.qza

qiime tools export \
  --input-path ./taxonomy/taxonomy.qza \
  --output-path ./taxonomy/

qiime metadata tabulate \
  --m-input-file ./Dada2/dada2_rep_seqs_filtered.qza \
  --m-input-file ./taxonomy.qza \
  --o-visualization ./Dada2/visualization/dada2_rep_seqs_filtered_taxa.qzv

biom add-metadata \
  -i Dada2/feature-table.biom \
  -o Dada2/feature-table_taxa.biom \
  --observation-metadata-fp taxonomy/taxonomy.tsv \
  --observation-header="Feature ID,Taxon" \
  --sc-separated taxonomy

biom convert \
  -i Dada2/feature-table_taxa.biom \
  -o Dada2/dada2_feature-table_taxa.tsv \
  --to-tsv \
  --output-metadata-id=Taxon \
  --tsv-metadata-formatter=naive \
  --header-key=Taxon

qiime metadata tabulate \
  --m-input-file ./taxonomy/taxonomy.qza \
  --o-visualization ./taxonomy/taxonomy_table.qzv

qiime taxa barplot \
  --i-table ./core-metrics-results/rarefied_table.qza \
  --i-taxonomy ./taxonomy/taxonomy.qza \
  --m-metadata-file ./metadata.txt \
  --o-visualization ./taxonomy/taxa_barplot.qzv

cp ./sequences/demux_seqs.qzv ./transfer/sequences/demux_seqs.qzv
cp ./sequences/trimmed_demux_seqs.qzv ./transfer/sequences/trimmed_demux_seqs.qzv
cp ./Dada2/visualization/dada2_stats.qzv ./transfer/./Dada2/dada2_stats.qzv
cp ./Dada2/visualization/dada2_table_filtered.qzv ./transfer/Dada2/dada2_table_filtered.qzv
cp ./Dada2/visualization/dada2_rep_seqs_filtered_taxa.qzv \
   ./transfer/Dada2/visualization/dada2_rep_seqs_filtered_taxa.qzv
cp ./Dada2/dada2_feature-table_taxa.tsv ./transfer/Dada2/dada2_feature-table_taxa.tsv
cp ./taxonomy/taxonomy_table.qzv ./transfer/taxonomy/taxonomy_table.qzv
cp ./taxonomy/taxa_barplot.qzv ./transfer/taxonomy/taxa_barplot.qzv

echo "### Ending at: $(date) ###"
