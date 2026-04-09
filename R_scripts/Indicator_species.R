library(tidyverse)
library(phyloseq)
library(mia)
library(scuttle)
library(vegan)
library(ggpubr)
library(labdsv)

# for the rarefaction method please refer to the "CAP_raref.R" script
rarefied_TSE <- read_rds("../input/rarefied_TSE.rds")

data_df <- read.csv("../input/barcode_sample_info.csv", quote="")
rownames(data_df) <- data_df$samplecode

#### Multiple comparison loop of all sample types ####

indicator_results <- list()
res_df <- list()
for(name in names(rarefied_TSE)){

  print(name)

  test_df <- assay(rarefied_TSE[[name]], "relabundance") %>%
    t() %>%
    as.data.frame()

  test_df <- test_df[,colSums(test_df) > 0] %>%
    merge(data_df, ., by=0) %>%
    column_to_rownames(var = "Row.names")

  clust <- test_df$treatment

  indicator_results[[name]] <- withr::with_seed(112, indval(test_df[-(1:13)],clust, numitr=10000))

  res_df[[name]] <- cbind(indicator_results[[name]]$indval,
                          indicator_results[[name]]$pval)

  res_df[[name]] <- res_df[[name]] %>%
    dplyr::rename("p.value" = "indicator_results[[name]]$pval") %>%
    rownames_to_column(var = "ASV")

  res_df[[name]]$p.BH = p.adjust(res_df[[name]]$p.value, "BH")
  
}

all_res <- bind_rows(res_df, .id = "Sampletype")


#### Pair-wise comparisons loop ####

comparisons <- list(c("C", "G"), c("C", "P"), c("C", "GP"), c("P", "G"), c("P", "GP"), c("G", "GP"))

indicator_results2 <- list()
res_df2 <- list()
for(s_name in names(rarefied_TSE)){
  
  for(t in comparisons){
    
    name <- paste(s_name, t[1], t[2], sep = "_")
    
    print(name)
    
    test_data <- rarefied_TSE[[s_name]][,colData(rarefied_TSE[[s_name]])$treatment == t[1] | colData(rarefied_TSE[[s_name]])$treatment == t[2]]
    
    test_df <- assay(test_data, "relabundance") %>%
      t() %>%
      as.data.frame()
    
    test_df <- test_df[,colSums(test_df) > 0] %>%
      merge(data_df, ., by=0) %>%
      column_to_rownames(var = "Row.names")
    
    clust <- test_df$treatment
    
    indicator_results2[[name]] <- withr::with_seed(112, indval(test_df[-(1:13)],clust, numitr=10000))
    
    res_df2[[name]] <- cbind(indicator_results2[[name]]$indval,
                            indicator_results2[[name]]$pval)
    
    res_df2[[name]] <- res_df2[[name]] %>%
      dplyr::rename("p.value" = "indicator_results2[[name]]$pval") %>%
      rownames_to_column(var = "ASV")
    res_df2[[name]]$p.BH = p.adjust(res_df2[[name]]$p.value, "BH")
    
  }
}
 
all_res2 <- bind_rows(res_df2, .id = "Sampletype")

#### Agrochemical specific test ####

indicator_results3 <- list()
res_df3 <- list()
for(s_name in names(rarefied_TSE)){
  
  for(t in c("glyphosate", "phosphorus")){
    
    name <- paste(s_name, t, sep = "_")
    
    print(name)
    
    test_df <- assay(rarefied_TSE[[s_name]], "relabundance") %>%
      t() %>%
      as.data.frame()
    
    test_df <- test_df[,colSums(test_df) > 0] %>%
      merge(data_df, ., by=0) %>%
      column_to_rownames(var = "Row.names")
    
    clust <- test_df[,t]
    
    indicator_results3[[name]] <- withr::with_seed(112, indval(test_df[-(1:13)],clust, numitr=10000))
    
    res_df3[[name]] <- cbind(indicator_results3[[name]]$indval,
                             indicator_results3[[name]]$pval)
    
    res_df3[[name]] <- res_df3[[name]] %>%
      dplyr::rename("p.value" = "indicator_results3[[name]]$pval") %>%
      rownames_to_column(var = "ASV")
    res_df3[[name]]$p.BH = p.adjust(res_df3[[name]]$p.value, "BH")
    
  }
}

all_res3 <- bind_rows(res_df3, .id = "Sampletype")

#### Prepare the table ####

no_sep <- all_res %>%
  filter(p.BH <= 0.05)
sep <- all_res2 %>%
  filter(p.BH <= 0.05) %>%
  filter(!is.na(`Control_No-phosphorus`))
chem <- all_res3 %>%
  filter(p.BH <= 0.05) %>%
  pivot_longer(cols = c("Control", "Phosphorus", "No-phosphorus", "Glyphosate")) %>%
  dplyr::rename(Treatment = name,
                Indval = value) %>%
  filter(!is.na(Indval)) %>%
  group_by(Sampletype, ASV) %>%
  slice_max(Indval, n = 1) %>%
  ungroup()


##### ASV taxonomies #####

ASV_tax <- as.matrix(read.delim("../input/ASV_tax_species.tsv", row.names=1, quote=""))

taxas <- as.data.frame(ASV_tax) %>%
  rownames_to_column(var = "taxon") %>%
  dplyr::select(!(confidence:sequence)) %>%
  rownames_to_column("ASV_number")

chem_tax <- dplyr::left_join(chem, taxas, by = join_by(ASV == taxon)) %>%
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
                       paste0(Level, ": ", Order, " ASV", ASV_number),
                       ifelse(Level == "Family",
                              paste0(Level, ": ", Family, " ASV", ASV_number),
                              ifelse(Level == "Genus",
                                     paste0(Level, ": ", Genus, " ASV", ASV_number),
                                     paste0(Level, ": ", Genus, " ", Species, " ASV", ASV_number))))) %>%
  mutate(Taxa = ifelse(Taxa == "Order: ",
                       taxon,
                       Taxa)) %>%
  select(1:8, Taxa)


##### Table with Taxas #####

indval_table <- chem_tax %>% dplyr::select(Sampletype, Treatment, Taxa, Indval, p.BH) %>%
  mutate(Treatment = str_replace_all(Treatment, "Phosphorus", "Phosphate"),
         Treatment = str_replace_all(Treatment, "-phosphorus", "-phosphate"),
         Taxa = str_replace_all(Taxa, "ASV", "ASV_"),
         Sampletype = str_remove_all(Sampletype, "_phosphorus"),
         Sampletype = str_replace_all(Sampletype, "_|-", " ")) %>%
  group_by(Sampletype, Treatment) %>%
  arrange(desc(Indval), .by_group = T) %>%
  ungroup()


write_tsv(indval_table, "../output/Indval_values_table.tsv")

