library(tidyverse)
library(phyloseq)
library(mia)
library(scuttle)
library(ggpubr)
library(car)
library(withr)
#remotes::install_version("matrixStats", version="1.1.0")

data_df <- read.csv("../input/barcode_sample_info.csv", quote="")

rownames(data_df) <- data_df$samplecode
Sample_metadata <- data_df %>%
  dplyr::select(c("plant", "part", "treatment", "plot", "glyphosate", "phosphorus", "timepoint")) %>%
  dplyr::filter(!(str_detect(rownames(data_df), "_NA_")))
Sample_metadata$code <- paste0(Sample_metadata$plant,Sample_metadata$part,Sample_metadata$treatment,Sample_metadata$timepoint)

##### combining the asv abundance, taxonomy and sample data #####

ASV_abundance <- as.matrix(read.delim("../input/ASV_table.tsv", row.names=1, quote=""))

ASV_tax <- as.matrix(read.delim("../input/ASV_tax_species.tsv", row.names=1, quote=""))

OTU <- phyloseq::otu_table(ASV_abundance, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(ASV_tax)
samples = phyloseq::sample_data(Sample_metadata)

alldata <- phyloseq::phyloseq(OTU, TAX, samples)

alldata_tse <- mia::convertFromPhyloseq(alldata)

# Remove cyanobacteria, NAs and eukaryota
alldata_tse <- alldata_tse[rowData(alldata_tse)$Domain != "Eukaryota" & rowData(alldata_tse)$Phylum != "Cyanobacteria" & !is.na(rowData(alldata_tse)$Phylum), ]

# Adding QC metrics
alldata_tse <- scuttle::addPerCellQCMetrics(alldata_tse, use.altexps = TRUE)

alldata_tse_rel <- alldata_tse

# groups

groups <- as.data.frame(colData(alldata_tse_rel)) %>%
  remove_rownames() %>%
  select(plant, part, timepoint) %>%
  distinct() %>%
  arrange(plant, part, timepoint)

##### using addAlpha with rarefaction within the method #####

# this takes a while
data_filt <- list()
withr::with_seed(112,
  for(i in 1:nrow(groups)){
  pl <- groups[i,1]
  pa <- groups[i,2]
  tp <- groups[i,3]

  name <- paste(pl, pa, tp, sep = "_")

  print(name)

  data_filt[[name]] <- alldata_tse_rel[,colData(alldata_tse_rel)$plant == pl & colData(alldata_tse_rel)$part == pa & colData(alldata_tse_rel)$timepoint == tp]

  a <- as.numeric(sort(colSums(assay(data_filt[[name]], "counts")))[1])
  b <- as.numeric(sort(colSums(assay(data_filt[[name]], "counts")))[2])
  c <- as.numeric(sort(colSums(assay(data_filt[[name]], "counts")))[3])
  d <- median(sort(colSums(assay(data_filt[[name]], "counts"))))

  if(a/d >= 0.15){
    data_filt[[name]] <- addAlpha(data_filt[[name]],
             assay.type = "counts",
             index = "shannon",
             sample = a, na.rm = TRUE,
             niter=100)
  } else if(b/d >= 0.15) {
    data_filt[[name]] <- addAlpha(data_filt[[name]],
             assay.type = "counts",
             index = "shannon",
             sample = b, na.rm = TRUE,
             niter=100)
  } else if(c/d >= 0.15){
    data_filt[[name]] <- addAlpha(data_filt[[name]],
                                  assay.type = "counts",
                                  index = "shannon",
                                  sample = c, na.rm = TRUE,
                                  niter=100)
  } else{
    data_filt[[name]] <- c
  }

}
)


pot_filt <- list()
withr::with_seed(112,
                 for(i in 1){
                   pl <- "Potato"
                   pa <- "Leaf"


                   name <- paste(pl, pa, sep = "_")

                   print(name)

                   pot_filt[[name]] <- alldata_tse_rel[,colData(alldata_tse_rel)$plant == pl & colData(alldata_tse_rel)$part == pa]

                   a <- as.numeric(sort(colSums(assay(pot_filt[[name]], "counts")))[1])
                   b <- as.numeric(sort(colSums(assay(pot_filt[[name]], "counts")))[2])
                   c <- as.numeric(sort(colSums(assay(pot_filt[[name]], "counts")))[3])
                   d <- median(sort(colSums(assay(pot_filt[[name]], "counts"))))

                   if(a/d >= 0.15){
                     pot_filt[[name]] <- addAlpha(pot_filt[[name]],
                                                   assay.type = "counts",
                                                   index = "shannon",
                                                   sample = a, na.rm = TRUE,
                                                   niter=100)
                   } else if(b/d >= 0.15) {
                     pot_filt[[name]] <- addAlpha(pot_filt[[name]],
                                                   assay.type = "counts",
                                                   index = "shannon",
                                                   sample = b, na.rm = TRUE,
                                                   niter=100)
                   } else if(c/d >= 0.15){
                     pot_filt[[name]] <- addAlpha(pot_filt[[name]],
                                                   assay.type = "counts",
                                                   index = "shannon",
                                                   sample = c, na.rm = TRUE,
                                                   niter=100)
                   } else{
                     pot_filt[[name]] <- c
                   }

                 }
)

##### Wilcoxon #####

indices <- list()
alphares <- list()
alphares2 <- list()
levenes <- list()
for(name in names(data_filt)){
  
  tmp_data <- as.data.frame(colData(data_filt[[name]])) %>%
    mutate(treatment = str_replace(treatment, "Control_No-phosphorus",
                                   "C")) %>%
    mutate(treatment = str_replace(treatment, "Control_Phosphorus",
                                   "P")) %>%
    mutate(treatment = str_replace(treatment, "Glyphosate_No-phosphorus",
                                   "G")) %>%
    mutate(treatment = str_replace(treatment, "Glyphosate_Phosphorus",
                                   "GP"))
  
  indices[[name]] <- tmp_data
  alphares[[name]] <- rstatix::wilcox_test(shannon ~ treatment, data = tmp_data)
  alphares2[[name]] <- rstatix::wilcox_effsize(shannon ~ treatment, data = tmp_data)
  levenes[[name]] <- leveneTest(shannon ~ as.factor(treatment), data = tmp_data)
  
}


##### Comparisons #####

alphares <- bind_rows(alphares, .id = "samplecode")
alphares2 <- bind_rows(alphares2, .id = "samplecode")
levenes <- bind_rows(levenes, .id = "samplecode") %>%
  remove_rownames() %>%
  filter(!is.na(`Pr(>F)`))

alphares_full <- dplyr::full_join(alphares, alphares2,
                                  by = join_by(samplecode, .y., group1, group2, n1, n2)) %>%
  select(-c("magnitude", "p.adj.signif", ".y.")) %>%
  mutate(samplecode = str_replace_all(samplecode, "-|_", " ")) %>%
  arrange(samplecode)


##### Plots #####

keepers <- alphares_full %>%
  filter(p.adj <= 0.05 & group1 == "C") %>%
  select(samplecode) %>%
  distinct() %>% pull(samplecode)

assumption_set <- levenes %>%
  filter(`Pr(>F)` >= 0.05) %>%
  select(samplecode) %>%
  mutate(samplecode = str_replace_all(samplecode, "-|_", " ")) %>%
  distinct() %>% pull(samplecode)

keepers <- keepers[keepers %in% assumption_set]

# contains "Oat Leaf Early"    "Potato Leaf Late"  "Potato Tuber Late"

pot_leaf_late <- indices$Potato_Leaf_Late
pot_leaf_late$treatment <- factor(pot_leaf_late$treatment, levels = c("C", "G", "GP", "P"))

pot_leaf_late_test <- alphares_full %>%
  filter(samplecode == "Potato Leaf Late") %>%
  mutate(y.position = c(4, 4.5, 4.75, 5, 5.25, 5.5))

A <- ggboxplot(pot_leaf_late, x = "treatment", y = "shannon") + 
  stat_pvalue_manual(pot_leaf_late_test, label = "p.adj", hide.ns = T) +
  ylab("Shannon index") +
  xlab("Treatment")


pot_tuber_late <- indices$Potato_Tuber_Late %>%
  filter(!is.na(shannon))
pot_tuber_late$treatment <- factor(pot_tuber_late$treatment, levels = c("C", "G", "GP", "P"))

pot_tuber_late_test <- alphares_full %>%
  filter(samplecode == "Potato Tuber Late") %>%
  mutate(y.position = c(5, 5.25, 5.5, 5.5, 6.5, 6.75))

B <- ggboxplot(pot_tuber_late, x = "treatment", y = "shannon") + 
  stat_pvalue_manual(pot_tuber_late_test, label = "p.adj", hide.ns = T) +
  ylab("Shannon index") +
  xlab("Treatment")

oat_leaf_early <- indices$Oat_Leaf_Early
oat_leaf_early$treatment <- factor(oat_leaf_early$treatment, levels = c("C", "G", "GP", "P"))

oat_leaf_early_test <- alphares_full %>%
  filter(samplecode == "Oat Leaf Early") %>%
  mutate(y.position = c(4.5, 4.75, 5, 5.25, 5.5, 5.75))

C <- ggboxplot(oat_leaf_early, x = "treatment", y = "shannon") + 
  stat_pvalue_manual(oat_leaf_early_test, label = "p.adj", hide.ns = T) +
  ylab("Shannon index") +
  xlab("Treatment")

ggarrange(ggarrange(A, B,C, ncol = 3, labels = c("a", "b", "c"))) 

ggsave(path = "../output/", filename = "alpha_diversities.pdf", width = 10, height = 5, device='pdf', dpi=300)

##### full wilcoxon table #####

write_tsv(alphares_full, "../output/all_alpha_wilcoxon_tests.tsv")

##### Full alpha div results #####

indices <- list()
for(name in names(data_filt)){
  
  tmp_data <- as.data.frame(colData(data_filt[[name]])) %>%
    mutate(treatment = str_replace(treatment, "Control_No-phosphorus",
                                   "C")) %>%
    mutate(treatment = str_replace(treatment, "Control_Phosphorus",
                                   "P")) %>%
    mutate(treatment = str_replace(treatment, "Glyphosate_No-phosphorus",
                                   "G")) %>%
    mutate(treatment = str_replace(treatment, "Glyphosate_Phosphorus",
                                   "GP"))
  
  indices[[name]] <- tmp_data
  
}

all_shannon_indices <- bind_rows(indices, .id = "samplecode") %>%
  remove_rownames() %>%
  select(samplecode, treatment, replicate_plot = plot, Shannon = shannon) %>%
  mutate(samplecode = str_replace_all(samplecode, "_|-", " ")) %>%
  arrange(samplecode, treatment, replicate_plot, Shannon)

write_tsv(all_shannon_indices, "../output/Shannon_indices.tsv")

