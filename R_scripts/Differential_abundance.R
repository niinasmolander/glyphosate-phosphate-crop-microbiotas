library(tidyverse)
library(readxl)
library(phyloseq)
library(scuttle)
library(vegan)
library(mia)
library(ALDEx2)
library(ANCOMBC)
library(GUniFrac)
library(dacomp)
library(eBay)

source("Diff_abund_function.R")

rarefied_TSE <- read_rds("../input/rarefied_TSE.rds")

data_df <- read.csv("../input/barcode_sample_info.csv", quote="")
rownames(data_df) <- data_df$samplecode

#### Differential abundance analysis ####

# Run the DAA  with 5 different DAA estimators
res <- list()
for(sampletype in names(rarefied_TSE)){
  pl = str_split_i(sampletype, "_" ,1)
  pa = str_split_i(sampletype, "_" ,2)
  tp = str_split_i(sampletype, "_" ,3)

  print(sampletype)

  data_df_plant <- data_df %>% filter(plant == pl & part == pa & timepoint == tp)

  # for some reason mia used original counts when doing the subsetting by prevalent
  # hence replacing those with the rarefied version.

  assay(rarefied_TSE[[sampletype]], "counts") <- assay(rarefied_TSE[[sampletype]], "subsampled")


  altExp(rarefied_TSE[[sampletype]], "prevalent") <- subsetByPrevalent(rarefied_TSE[[sampletype]],
                                                                       assay.type = "subsampled",
                                                                       detection = 0,
                                                                       prevalence = 0.1)


  res[[sampletype]] <- five_DAA_estimators(altExp(rarefied_TSE[[sampletype]], "prevalent"), data_df_plant)

}


#### Combining ####

###### Aldex ######

aldex <- list()
for(name in names(res)){
  print(name)
  for(i in names(res[[name]]$aldex)){
    print(i)
    fullname <- paste(name, i, sep = "_")
    aldex[[fullname]] <- res[[name]]$aldex[[i]] %>%
      rownames_to_column(var = "taxon")
  }
}

aldex_all <- bind_rows(aldex, .id = "Comparison") %>%
  dplyr::select(Comparison, taxon, wi.eBH) %>%
  dplyr::rename("p-value" = "wi.eBH")

###### Ancombc2 ######

ancombc <- list()
for(name in names(res)){
  print(name)
  for(i in names(res[[name]]$ancombc2)){
    print(i)
    fullname <- paste(name, i, sep = "_")
    ancombc[[fullname]] <- res[[name]]$ancombc2[[i]]$res %>%
      dplyr::select(taxon, starts_with("q_treatment"))
    colnames(ancombc[[fullname]]) <- c("taxon", "p-value")
  }
}

ancombc_all <- bind_rows(ancombc, .id = "Comparison")

###### ZicoSeq ######

zicoseq <- list()
for(name in names(res)){
  print(name)
  for(i in names(res[[name]]$zicoseq)){
    print(i)
    fullname <- paste(name, i, sep = "_")
    zicoseq[[fullname]] <- as.data.frame(res[[name]]$zicoseq[[i]]$p.adj.fdr) %>%
      rownames_to_column(var = "taxon")
  }
}
zicoseq_all <- bind_rows(zicoseq, .id = "Comparison") %>%
  dplyr::rename("p-value" = "res[[name]]$zicoseq[[i]]$p.adj.fdr")

###### Dacomp ######

dacomp <- list()
for(name in names(res)){
  print(name)
  dacomp[[name]] <- bind_rows(res[[name]]$dacomp, .id = "Comp")
}

dacomp_all <- bind_rows(dacomp, .id = "Comparison") %>%
  mutate(Comparison = paste(Comparison, Comp, sep = "_")) %>%
  dplyr::select(-"Comp")

colnames(dacomp_all) <- c("Comparison", "taxon", "p-value")

##### eBay #####

ebay <- list()
for(name in names(res)){
  print(name)
  ebay[[name]] <- res[[name]]$eBay
}

ebay_all <- bind_rows(ebay, .id = "Comparison2") %>%
  mutate(Comparison = paste(Comparison2, Comparison, sep = "_")) %>%
  dplyr::select(-"Comparison2")


##### Combining #####
two_datas <- dplyr::full_join(aldex_all, ancombc_all, by = join_by(Comparison, taxon), suffix = c(".aldex", ".ancombc2"))
three_datas <- dplyr::full_join(two_datas, zicoseq_all, by = join_by(Comparison, taxon)) %>%
  dplyr::rename("p-value.zicoseq" = "p-value")
four_datas <- dplyr::full_join(three_datas, dacomp_all, by = join_by(Comparison, taxon)) %>%
  dplyr::rename("p-value.dacomp" = "p-value")

all_DAA_res <- dplyr::left_join(four_datas, ebay_all, by = join_by(Comparison, taxon))

#### Filtering and plots ####

all_DAA_res_filt <- all_DAA_res %>%
  dplyr::filter(`p-value.aldex` <= 0.05 | `p-value.ancombc2` <= 0.05 | `p-value.zicoseq` <= 0.05 | `p-value.dacomp` <= 0.05 | `p-value.eBay` <= 0.05) %>%
  mutate(Sampletype = paste(str_split_i(Comparison, "_", 1),
                            str_split_i(Comparison, "_", 2),
                            str_split_i(Comparison, "_", 3),
                            sep = "_"),
         Comparison = str_split_i(Comparison, "_", 4),
         Comparison = str_replace_all(Comparison, "AND", " vs. "))

all_DAA_res_filt$agreement <- rowSums(!(is.na(all_DAA_res_filt[,3:7])) & all_DAA_res_filt[,3:7] <= 0.05)

important <- all_DAA_res_filt %>%
  filter(agreement >= 3) %>%
  mutate(Comparison = str_replace_all(Comparison, "G vs. P", "P vs. G"))


###### Add info about the species ######

ASV_tax <- as.matrix(read.delim("../input/ASV_tax_species.tsv", row.names=1, quote=""))

taxas <- as.data.frame(ASV_tax) %>%
  rownames_to_column(var = "taxon") %>%
  dplyr::select(!(confidence:sequence)) %>%
  rownames_to_column("ASV_number")

important <- dplyr::left_join(important, taxas, by = join_by("taxon" == "taxon")) %>%
  mutate(Species = ifelse(Species_exact == "",
                          Species,
                          Species_exact)) %>%
  mutate(Level = ifelse(Family == "",
                        "Order",
                        ifelse(Genus == "",
                               "Family",
                               ifelse(Species == "",
                                      "Genus",
                                      "Species")))) %>%
  mutate(Taxa = ifelse(Level == "Order",
                       paste0(Level, ": ", Order, " ASV_", ASV_number),
                       ifelse(Level == "Family",
                              paste0(Level, ": ", Family, " ASV_", ASV_number),
                              ifelse(Level == "Genus",
                                     paste0(Level, ": ", Genus, " ASV_", ASV_number),
                                     paste0(Level, ": ", Genus, " ", Species, " ASV_", ASV_number))))) %>%
  mutate(Taxa = ifelse(Taxa == "Order: ",
                       taxon,
                       Taxa))


###### Direction ######

relabundances <- list()
for(name in names(rarefied_TSE)){
  relabundances[[name]] <- as.data.frame(assay(rarefied_TSE[[name]], "relabundance")) %>%
    rownames_to_column(var = "taxon") %>%
    pivot_longer(-1, names_to = "sample", values_to = "relabundance")
}

relabundances <- bind_rows(relabundances, .id = "sampletype")

relabundances <- dplyr::left_join(relabundances, dplyr::select(data_df, samplecode, treatment),
                                  by = join_by(sample == samplecode)) %>%
  mutate(treatment = str_replace(treatment, "Control_No-phosphorus",
                                 "C")) %>%
  mutate(treatment = str_replace(treatment, "Control_Phosphorus",
                                 "P")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_No-phosphorus",
                                 "G")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_Phosphorus",
                                 "GP"))

relabundances <- relabundances %>%
  group_by(sampletype, treatment, taxon) %>%
  mutate(median = median(relabundance),
         mean = mean(relabundance),
         n = n()) %>%
  ungroup() %>%
  mutate(treatment = str_replace(treatment, "Control_No-phosphorus",
                                 "C")) %>%
  mutate(treatment = str_replace(treatment, "Control_Phosphorus",
                                 "P")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_No-phosphorus",
                                 "G")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_Phosphorus",
                                 "GP"))

keepers <- important %>%
  mutate(first = str_split_i(Comparison, " vs. ", 1),
         second = str_split_i(Comparison, " vs. ", 2)) %>%
  dplyr::left_join(distinct(dplyr::select(relabundances, taxon, sampletype, treatment, median, mean)),
                   by = join_by(taxon, Sampletype == sampletype, first == treatment)) %>%
  dplyr::left_join(distinct(dplyr::select(relabundances, taxon, sampletype, treatment, median, mean)),
                   by = join_by(taxon, Sampletype == sampletype, second == treatment),
                   suffix = c("_first", "_second")) %>%
  mutate(direction_med = median_first - median_second,
         direction_mean = mean_first - mean_second) %>%
  filter(direction_mean >= 0.05 | direction_mean <= -0.05)

fortable <- important %>%
  mutate(first = str_split_i(Comparison, " vs. ", 1),
         second = str_split_i(Comparison, " vs. ", 2)) %>%
  dplyr::left_join(distinct(dplyr::select(relabundances, taxon, sampletype, treatment, median, mean, n)),
                   by = join_by(taxon, Sampletype == sampletype, first == treatment)) %>%
  dplyr::left_join(distinct(dplyr::select(relabundances, taxon, sampletype, treatment, median, mean, n)),
                   by = join_by(taxon, Sampletype == sampletype, second == treatment),
                   suffix = c("_first", "_second")) %>%
  mutate(direction_med = median_first - median_second,
         direction_mean = mean_first - mean_second) 

forplot <- relabundances %>%
  dplyr::semi_join(keepers, by = join_by(sampletype == Sampletype, taxon == taxon)) %>%
  dplyr::left_join(distinct(dplyr::select(keepers, taxon, Taxa)), by = join_by(taxon)) %>%
  mutate(sampletype = str_replace_all(sampletype, "_|-", " "),
         sampletype.f = factor(sampletype, levels = c("Potato Leaf Late",
                                                    "Faba bean Leaf Early",
                                                    "Faba bean Leaf Late"))) %>%
  mutate(Taxa1 = str_split_i(Taxa, " ASV", 1),
         Taxa2 = paste0("ASV", str_split_i(Taxa, " ASV", 2)))
  

test <- keepers %>%
  dplyr::select(taxon, Taxa, first, second, Sampletype) %>%
  mutate(y.value = c(0.5, 0.55, 0.6, 0.3, 0.35, 0.3, 0.4, 0.35, 0.45, 0.4, 0.5, 0.45, 0.55)) %>%
  rownames_to_column(var = "group") %>%
  pivot_longer(cols = c("first", "second")) %>%
  dplyr::select(-name) %>%
  dplyr::rename("sampletype" = "Sampletype") %>%
  mutate(sampletype = str_replace_all(sampletype, "_|-", " "),
         sampletype.f = factor(sampletype, levels = c("Potato Leaf Late",
                                                      "Faba bean Leaf Early",
                                                      "Faba bean Leaf Late"))) %>%
  mutate(Taxa1 = str_split_i(Taxa, " ASV", 1),
         Taxa2 = paste0("ASV", str_split_i(Taxa, " ASV", 2)))



potleaflate_asv9 <- ggplot() +
  geom_boxplot(data = filter(forplot, sampletype == "Potato Leaf Late" & str_detect(Taxa, "ASV_9")), aes(x = treatment, y = relabundance)) +
  geom_line(data = filter(test, sampletype == "Potato Leaf Late" & str_detect(Taxa, "ASV_9")), aes(x = value, y = y.value, group = group)) +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

potleaflate_asv3 <- ggplot() +
  geom_boxplot(data = filter(forplot, sampletype == "Potato Leaf Late" & str_detect(Taxa, "ASV_3")), aes(x = treatment, y = relabundance)) +
  geom_line(data = filter(test, sampletype == "Potato Leaf Late" & str_detect(Taxa, "ASV_3")), aes(x = value, y = y.value, group = group)) +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fabaleafearly_asv2 <- ggplot() +
  geom_boxplot(data = filter(forplot, sampletype == "Faba bean Leaf Early" & str_detect(Taxa, "ASV_2")), aes(x = treatment, y = relabundance)) +
  geom_line(data = filter(test, sampletype == "Faba bean Leaf Early" & str_detect(Taxa, "ASV_2")), aes(x = value, y = y.value, group = group)) +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fabaleaflate_asv3 <- ggplot() +
  geom_boxplot(data = filter(forplot, sampletype == "Faba bean Leaf Late" & str_detect(Taxa, "ASV_3")), aes(x = treatment, y = relabundance)) +
  geom_line(data = filter(test, sampletype == "Faba bean Leaf Late" & str_detect(Taxa, "ASV_3")), aes(x = value, y = y.value, group = group)) +
  ylab("Relative abundance") +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggpubr::ggarrange(potleaflate_asv9, potleaflate_asv3, fabaleafearly_asv2, fabaleaflate_asv3, labels = c(letters[1:4]))

ggsave(path = "../output/", filename = "DAA_top_species.pdf", width = 10, height = 8, device='pdf', dpi=300)

##### table #####

fortable <- fortable %>%
  dplyr::select(Sampletype, Taxa, Comparison, mean_first, mean_second, agreement, `p-value.aldex`, `p-value.ancombc2`, `p-value.zicoseq`, `p-value.dacomp`, `p-value.eBay`, n_first, n_second) %>%
  arrange(Sampletype, Taxa, Comparison) %>%
  mutate(mean_first = round(mean_first, digits = 3) * 100,
         mean_second = round(mean_second, digits = 3) * 100,
         Sampletype = str_replace_all(Sampletype, "-|_", " "))

write_tsv(fortable, "../output/DAA_table_means.tsv")
