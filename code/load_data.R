library(tidyverse)
library(cowplot)
library(qiime2R)
library(vegan)
library(ggpubr)
library(rstatix)


metadata <- read_delim("data/metadata.txt")

table <- read_tsv("data/core-metrics-results_44959/rarefied-feature-table_taxa.tsv",
                  skip = 1) %>% 
  rename(., featureid =  `#OTU ID`)

feature.counts <- read_csv("data/Dada2/dada2-sample-frequency-detail.csv")
colnames(feature.counts) <- c("sampleid", "feature.counts")

taxonomy <- read_qza("data/taxonomy/taxonomy.qza")$data %>% 
  parse_taxonomy()




