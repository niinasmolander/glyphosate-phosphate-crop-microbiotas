library(readxl)
library(tidyverse)
library(vegan)
library(dunn.test)
library(mia)

#### prep the data ####

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
  mutate(Treament = str_remove_all(treatment, "ontrol|lyphosate|hosphorus|\\+"))

##### plots #####

plot_data <- all_res %>%
  group_by(treatment, plant, part, timepoint, status) %>%
  summarise(sum = sum(relabundance)) %>% ungroup() %>% group_by(treatment, plant, part, timepoint) %>%
  mutate(freq = sum / sum(sum)) %>% ungroup() %>%
  mutate(treatment = str_remove_all(treatment, "ontrol|lyphosate|hosphorus|\\+")) %>%
  mutate(freq2 = ifelse(round(freq, digits = 2) == 0, "", round(freq, digits = 2)))

pot <- plot_data %>% filter(plant == "Potato") %>%
  ggplot(aes(x = treatment, y = freq, fill = status, label = freq2)) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5), size = 2)+
  facet_grid(. ~ part+timepoint) +
  ggpubr::theme_pubclean() +
  xlab("Treatment") +
  ylab("Proportion") +
  scale_fill_manual(values = c("orange","skyblue", "lightgray"))

faba <- plot_data %>% filter(plant == "Faba-bean") %>%
  ggplot(aes(x = treatment, y = freq, fill = status, label = freq2)) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5), size = 2)+
  facet_grid(. ~ part+timepoint) +
  ggpubr::theme_pubclean() +
  xlab("Treatment") +
  ylab("Proportion") +
  scale_fill_manual(values = c("orange","skyblue", "lightgray"))

oat <- plot_data %>% filter(plant == "Oat") %>%
  ggplot(aes(x = treatment, y = freq, fill = status, label = freq2)) +
  geom_bar(stat = "identity") +
  geom_text(position = position_stack(vjust = 0.5), size = 2)+
  facet_grid(. ~ part+timepoint) +
  ggpubr::theme_pubclean() +
  xlab("Treatment") +
  ylab("Proportion") +
  scale_fill_manual(values = c("orange","skyblue", "lightgray"))

ggpubr::ggarrange(pot, faba, oat, labels = c("a","b", "c"), common.legend = T,
          ncol = 1)

ggsave(path = "../output/", filename = "EPSPS_with_values.pdf", width = 8, height = 10, device='pdf', dpi=300)
