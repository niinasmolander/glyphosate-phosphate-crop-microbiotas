library(tidyverse)
library(mia)

# for the rarefaction method please refer to the "CAP_raref.R" script
raref <- read_rds("../input/rarefied_TSE.rds")

raref_counts <- list()
for(name in names(raref)){
   tmp <- as.data.frame(assay(raref[[name]], "subsampled"))
   raref_counts[[name]] <- as.data.frame(colSums(tmp)) %>%
     rownames_to_column(var = "sample")
   colnames(raref_counts[[name]]) <- c("sample", "raref_counts")
}

raref_counts <- bind_rows(raref_counts)

read_counts <- read_tsv("../input/overall_summary.tsv",
                        col_types = c("ccccccccc")) %>%
  select(sample, raw_n_reads = cutadapt_total_processed, filt_n_reads = nonchim) %>%
  mutate(raw_n_reads = as.numeric(str_remove_all(raw_n_reads, ",")),
         filt_n_reads = as.numeric(filt_n_reads)) %>%
  dplyr::left_join(raref_counts, by = "sample")

data_df <- read.csv("../input/barcode_sample_info.csv", quote="") %>%
  select(sample = samplecode, plant, part, timepoint, plot, glyphosate, phosphorus)

read_counts <- dplyr::left_join(read_counts, data_df, by = join_by(sample)) %>%
  mutate(glyphosate = str_sub(glyphosate, 1,1),
         phosphorus = str_sub(phosphorus, 1,1),
         treatment = paste0(glyphosate, phosphorus),
         treatment = str_remove_all(treatment, "N"),
         treatment = str_replace_all(treatment, "CP", "P")) %>%
  group_by(plant, part, timepoint, treatment) %>%
  mutate(n = n(),
         n_after_raref = sum(!is.na(raref_counts))) %>% ungroup()

finaltable <- read_counts %>%
  group_by(plant, part, timepoint, treatment) %>%
  summarise(raw_mean = mean(raw_n_reads),
            raw_median = median(raw_n_reads),
            raw_min = min(raw_n_reads),
            raw_max = max(raw_n_reads),
            raw_sd = sd(raw_n_reads),
            filt_mean = mean(filt_n_reads),
            filt_median = median(filt_n_reads),
            filt_min = min(filt_n_reads),
            filt_max = max(filt_n_reads),
            filt_sd = sd(filt_n_reads),
            rarefaction_depth = mean(raref_counts, na.rm = T),
            n = mean(n),
            n_after_raref = mean(n_after_raref)) %>% ungroup()

write_tsv(finaltable, "../output/read_counts.tsv")

