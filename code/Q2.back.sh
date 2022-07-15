#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH --job-name=q2_back        # name for the job
#SBATCH --cpus-per-task=5              # number of cores
#SBATCH --mem=100G                       # total memory
#SBATCH --nodes 2
#SBATCH --time 1-00:00:00                 # time limit in the form days-hours:minutes
#SBATCH --mail-user=zlmg2b@umsystem.edu    # email address for notifications
#SBATCH --mail-type=FAIL,END,BEGIN           # email types
#SBATCH --partition Lewis            # max of 1 node and 4 hours; use `Lewis` for larger jobs
#--------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

module load miniconda3
source activate qiime2-2021.8

qiime diversity core-metrics-phylogenetic \
  --i-table ./Dada2/dada2_table_filtered.qza \
  --i-phylogeny ./phylogeny/rooted-tree.qza \
  --p-sampling-depth 44960
  --m-metadata-file ./metadata.txt \
  --p-n-jobs-or-threads $SLURM_CPUS_ON_NODE \
  --output-dir ./	

qiime tools export \
  --input-path ./core-metrics-results/rarefied_table.qza \
  --output-path ./core-metrics-results/

qiime tools export \
  --input-path ./Dada2/dada2_table_filtered.qza \
  --output-path ./Dada2/

biom add-metadata \
  -i core-metrics-results/feature-table.biom \
  -o core-metrics-results/feature-table_taxa.biom \
  --observation-metadata-fp taxonomy/taxonomy.tsv \
  --observation-header="Feature ID,Taxon" \
  --sc-separated taxonomy

biom convert \
  -i core-metrics-results/feature-table_taxa.biom \
  -o core-metrics-results/feature-table_taxa.tsv \
  --to-tsv \
  --output-metadata-id=Taxon \
  --tsv-metadata-formatter=naive \
  --header-key=Taxon

cp -r ./core-metrics-results ./transfer/

### ---------------- ANCOM --------------
# pseudocount 1
mkdir ancom

qiime composition add-pseudocount \
  --i-table ./core-metrics-results/rarefied_table.qza \
  --p-psuedocount 1 \
  --o-composition-table ./core-metrics-results/rarefied_table_pseudocount1.qza

qiime composition ancom \
  --i-table ./core-metrics-results/rarefied_table_pseudocount1.qza \
  --m-metadata-file ./metadata.txt \
  --m-metadata-column Group \
  --o-visualization ./ancom/ancom_group.qzv

qiime tools export \
  --input-path ./ancom/ancom_group.qzv
  --output-path ./ancom/

qiime metadata tabulate \
  --m-input-file ./ancom/ancom.tsv \
  --m-input-file ./kw_aldex2.tsv \
  --m-input-file ./taxonomy.tsv
  --m-input-file ./taxonomy/taxonomy.qza \
  --o-visualization ./ancom/metadata_taxnomy_aldex2_phylum_ancom.qzv

qiime tools export \
  --input-path ./ancom/metadata_taxonomy_aldex2_phylum_ancom.qzv \
  --output-path ./ancom/

qiime empress community-plot \
  --i-tree ./phylogeny/rooted-tree.qza \
  --i-feature-table ./core-metrics-results/rarefied_table/qza \
  --m-sample-metadata ./metadata.txt \
  --m-feature-metadata ./ancom/metadata.tsv \
  --o-visualization ./empress_phylum_aldex2_ancom_taxa.qzv



echo "### Ending at: $(date) ###"
