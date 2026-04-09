library(tidyverse)
library(phyloseq)
library(mia)
library(scuttle)
library(vegan)
library(ggpubr)
library(ggthemes)
#remotes::install_version("matrixStats", version="1.1.0")

data_df <- read.csv("../input/barcode_sample_info.csv", quote="")

rownames(data_df) <- data_df$samplecode
data_df <- data_df %>%
  mutate(treatment = str_replace(treatment, "Control_No-phosphorus",
                                 "Control")) %>%
  mutate(treatment = str_replace(treatment, "Control_Phosphorus",
                                 "Phosphorus")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_No-phosphorus",
                                 "Glyphosate")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_Phosphorus",
                                 "Glyphosate+Phosphorus"))

Sample_metadata <- data_df %>%
  dplyr::select(c("plant", "part", "treatment", "plot", "glyphosate", "phosphorus", "timepoint")) %>%
  dplyr::filter(!(str_detect(rownames(data_df), "_NA_")))
Sample_metadata$code <- paste0(Sample_metadata$plant,Sample_metadata$part,Sample_metadata$treatment,Sample_metadata$timepoint)

##### combining the asv abundance, taxonomy and sample datas #####

ASV_abundance <- as.matrix(read.delim("../input/ASV_table.tsv", row.names=1, quote=""))

ASV_tax <- as.matrix(read.delim("../input/ASV_tax_species.tsv", row.names=1, quote=""))

OTU <- phyloseq::otu_table(ASV_abundance, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(ASV_tax)
samples = phyloseq::sample_data(Sample_metadata)

alldata <- phyloseq(OTU, TAX, samples)

alldata_tse <- mia::convertFromPhyloseq(alldata)

# Remove cyanobacteria, NAs and eukaryota
alldata_tse <- alldata_tse[rowData(alldata_tse)$Domain != "Eukaryota" & rowData(alldata_tse)$Phylum != "Cyanobacteria", ]

# Adding QC metrics
alldata_tse <- addPerCellQCMetrics(alldata_tse, use.altexps = TRUE)
alldata_tse_rel <- alldata_tse[, colData(alldata_tse)$sum >= 50]

# Relative abundances
alldata_tse_rel <- transformAssay(alldata_tse, method = "relabundance")

##### agglomeration #####

# Agglomerating to each taxonomic level
altExp(alldata_tse_rel,"Phylum") <- agglomerateByRank(alldata_tse_rel, "Phylum")
altExp(alldata_tse_rel,"Order") <- agglomerateByRank(alldata_tse_rel, "Order")
altExp(alldata_tse_rel,"Family") <- agglomerateByRank(alldata_tse_rel, "Family")
altExp(alldata_tse_rel,"Genus") <- agglomerateByRank(alldata_tse_rel, "Genus")

##### Phylum Plots #####

top_taxa_all <- getTop(altExp(alldata_tse_rel, "Phylum"), top = 13, assay_name = "relabundance") 

rowData(altExp(alldata_tse_rel, "Phylum"))$Phylum <- if_else(rowData(altExp(alldata_tse_rel, "Phylum"))$Phylum %in% top_taxa_all, rowData(altExp(alldata_tse_rel, "Phylum"))$Phylum, "Other")

df_all <- assay(altExp(alldata_tse_rel, "Phylum"), "counts") %>%
  merge(rowData(altExp(alldata_tse_rel, "Phylum")), by = 0) %>%
  pivot_longer(cols = 2:(which(colnames(.) == "Domain")-1), names_to = "Samples", values_to = "counts") %>%
  dplyr::left_join(rownames_to_column(data.frame(colData(alldata_tse_rel)), var = "Samples")) %>%
  group_by(plant, timepoint, part, treatment, Phylum) %>%
  summarise(counts = sum(counts)) %>%
  mutate(treatment = str_remove_all(treatment, "ontrol|lyphosate|hosphorus|\\+"))

pot <- df_all %>% filter(plant == "Potato") %>%
  ggplot(aes(x = treatment, y = counts, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(. ~ part+timepoint) +
  scale_fill_tableau(palette = "Tableau 20") +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme_pubclean()

faba <- df_all %>% filter(plant == "Faba-bean") %>%
  ggplot(aes(x = treatment, y = counts, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(. ~ part+timepoint) +
  scale_fill_tableau(palette = "Tableau 20") +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme_pubclean()

oat <- df_all %>% filter(plant == "Oat") %>%
  ggplot(aes(x = treatment, y = counts, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(. ~ part+timepoint) +
  scale_fill_tableau(palette = "Tableau 20") +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme_pubclean()



ggarrange(pot, faba, oat, labels = c("a","b", "c"), common.legend = T,
          ncol = 1)

ggsave(path = "../output/", filename = "barplots_phylum.pdf", width = 10, height = 12, device='pdf', dpi=300)


##### Order Plots #####

top_taxa_all <- getTop(altExp(alldata_tse_rel, "Order"), top = 14, assay_name = "relabundance")

rowData(altExp(alldata_tse_rel, "Order"))$Order <- if_else(rowData(altExp(alldata_tse_rel, "Order"))$Order %in% top_taxa_all, rowData(altExp(alldata_tse_rel, "Order"))$Order, "Other")

df_all <- assay(altExp(alldata_tse_rel, "Order"), "counts") %>%
  merge(rowData(altExp(alldata_tse_rel, "Order")), by = 0) %>%
  pivot_longer(cols = 2:(which(colnames(.) == "Domain")-1), names_to = "Samples", values_to = "counts") %>%
  dplyr::left_join(rownames_to_column(data.frame(colData(alldata_tse_rel)), var = "Samples")) %>%
  group_by(plant, timepoint, part, treatment, Order) %>%
  summarise(counts = sum(counts)) %>%
  mutate(treatment = str_remove_all(treatment, "ontrol|lyphosate|hosphorus|\\+"))

pot <- df_all %>% filter(plant == "Potato") %>%
  ggplot(aes(x = treatment, y = counts, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(. ~ part+timepoint) +
  scale_fill_tableau(palette = "Tableau 20") +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme_pubclean()

faba <- df_all %>% filter(plant == "Faba-bean") %>%
  ggplot(aes(x = treatment, y = counts, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(. ~ part+timepoint) +
  scale_fill_tableau(palette = "Tableau 20") +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme_pubclean()

oat <- df_all %>% filter(plant == "Oat") %>%
  ggplot(aes(x = treatment, y = counts, fill = Order)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(. ~ part+timepoint) +
  scale_fill_tableau(palette = "Tableau 20") +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme_pubclean()

ggarrange(pot, faba, oat, labels = c("a","b", "c"), common.legend = T,
          ncol = 1)

ggsave(path = "../output", filename = "barplots_order.pdf", width = 10, height = 12, device='pdf', dpi=300)
