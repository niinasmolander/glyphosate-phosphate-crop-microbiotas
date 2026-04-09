library(scuttle)
library(ggpubr)
library(car)
library(readxl)
library(tidyverse)

#### Prep the data ####
epsps_results <- read_excel("../input/epsps_results.xlsx", 
                            sheet = "EPSPS") %>%
  select(1:5) %>%
  dplyr::rename(OTU = `#OTU` ,status = ...5)


rarefied_TSE <- read_rds("../input/rarefied_TSE.rds")

samples <- read_tsv("../input/samples_info.tsv")

relab <- list()
for(name in names(rarefied_TSE)){
  relab[[name]] <- as.data.frame(assay(rarefied_TSE[[name]], "relabundance")) %>%
    rownames_to_column(var = "OTU") %>%
    relocate(OTU) %>%
    pivot_longer(-1)
  colnames(relab[[name]]) <- c("OTU", "sample", "relabundance")
}

relab <- bind_rows(relab, .id = "samplecode")

all_res <- dplyr::left_join(relab, epsps_results, by = "OTU") %>%
  dplyr::left_join(samples, by = join_by(sample == Sample)) %>%
  mutate(Treament = str_remove_all(treatment, "ontrol|lyphosate|hosphorus|\\+")) %>%
  group_by(plant, part, treatment, timepoint, plot, sample, status) %>%
  summarise(sum = sum(relabundance)) %>% ungroup() %>% group_by(sample) %>%
  mutate(freq = sum / sum(sum)) %>% ungroup() %>%
  mutate(samplecode = paste(plant, part, timepoint, sep = "_"))

comp <- list(c("Control", "Glyphosate"), c("Control", "Phosphorus"), c("Control", "Glyphosate+Phosphorus"))

groups <- unique(all_res$samplecode)

#### run tests ####
  
wilcox_res <- list()
levenes <- list()
for(i in groups){
  
  print(i)
  
  data <- all_res %>%
    dplyr::filter(samplecode == i) %>%
    filter(status == "S")
  
  for(c in comp){
    tmp_data <- data %>%
      filter(treatment %in% c)
    
    name <- paste(i, "S", sep = "_")
    name2 <- paste(name, paste(c, collapse = "_"), sep = "_")
    
    levenes[[name2]] <- rbind(levenes[[name2]], leveneTest(freq ~ as.factor(treatment), data = tmp_data))
    wilcox_res[[name]] <- rbind(wilcox_res[[name]], rstatix::wilcox_test(freq ~ treatment, data = tmp_data))
  }
  
  wilcox_res[[name]]$p.holm <- p.adjust(wilcox_res[[name]]$p, method = "holm")
  
}

#### tables ####

all_S <- bind_rows(wilcox_res, .id = "sample") %>%
  mutate(p.holm = round(p.holm, digits = 3))

# this is for checking the comparisons meet the assumptions
levenes <- bind_rows(levenes, .id = "sampletype")

write_tsv(all_S, "../output/EPSPS_sensitive_Wilcoxon_test.tsv")
