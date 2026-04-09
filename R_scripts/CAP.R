library(tidyverse)
library(phyloseq)
library(mia)
library(scuttle)
library(vegan)
library(ggpubr)

data_df <- read.csv("../input/barcode_sample_info.csv", quote="")

rownames(data_df) <- data_df$samplecode
data_df <- data_df %>%
  mutate(treatment = str_replace(treatment, "Control_No-phosphorus",
                                 "C")) %>%
  mutate(treatment = str_replace(treatment, "Control_Phosphorus",
                                 "P")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_No-phosphorus",
                                 "G")) %>%
  mutate(treatment = str_replace(treatment, "Glyphosate_Phosphorus",
                                 "GP")) %>%
  dplyr::rename("G" = "glyphosate",
                "P" = "phosphorus")

Sample_metadata <- data_df %>%
  dplyr::select(c("plant", "part", "treatment", "plot", "G", "P", "timepoint")) %>%
  dplyr::filter(!(str_detect(rownames(data_df), "_NA_")))
Sample_metadata$code <- paste0(Sample_metadata$plant,Sample_metadata$part,Sample_metadata$treatment,Sample_metadata$timepoint)

##### combining the asv abundance, taxonomy and sample datas #####

ASV_abundance <- as.matrix(read.delim("../input/ASV_table.tsv", row.names=1, quote=""))

ASV_tax <- as.matrix(read.delim("../input/ASV_tax_species.tsv", row.names=1, quote=""))

taxas <- as.data.frame(ASV_tax) %>%
  rownames_to_column(var = "taxon") %>%
  dplyr::select(!(confidence:sequence)) %>%
  rownames_to_column("ASV_number")

OTU <- phyloseq::otu_table(ASV_abundance, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(ASV_tax)
samples = phyloseq::sample_data(Sample_metadata)

alldata <- phyloseq(OTU, TAX, samples)

alldata_tse <- mia::convertFromPhyloseq(alldata)

# Remove cyanobacteria, NAs and eukaryota
alldata_tse <- alldata_tse[rowData(alldata_tse)$Domain != "Eukaryota" & rowData(alldata_tse)$Phylum != "Cyanobacteria" & !is.na(rowData(alldata_tse)$Phylum), ]

# Adding QC metrics
alldata_tse <- addPerCellQCMetrics(alldata_tse, use.altexps = TRUE)

groups <- as.data.frame(colData(alldata_tse)) %>%
  remove_rownames() %>%
  select(plant, part, timepoint) %>%
  distinct() %>%
  arrange(plant, part, timepoint)


##### using avgdist to calculate the bray-curtis dissimilarity #####

diss <- list()
betadisper <- list()
betadisper_anova <- list()
cap <- list()
cap_anova <- list()
data_filt <- list()
eig <- list()

withr::with_seed(112,
  for(i in 1:nrow(groups)){
  
  # get the plant, part, timepoint
  pl <- groups[i,1]
  pa <- groups[i,2]
  tp <- groups[i,3]

  name <- paste(pl, pa, tp, sep = "_")

  print(name)
  
  # filter tse based on plant, part and timepoint
  data_filt[[name]] <- alldata_tse[,colData(alldata_tse)$plant == pl & colData(alldata_tse)$part == pa & colData(alldata_tse)$timepoint == tp]

  # obtain the lowest read numbers for rarefying and median
  a <- as.numeric(sort(colSums(assay(data_filt[[name]], "counts")))[1])
  b <- as.numeric(sort(colSums(assay(data_filt[[name]], "counts")))[2])
  c <- as.numeric(sort(colSums(assay(data_filt[[name]], "counts")))[3])
  d <- median(sort(colSums(assay(data_filt[[name]], "counts"))))
  
  # community matrix
  data_tmp <- t(assay(data_filt[[name]], "counts"))
  
  if(a/d >= 0.15){
    # bray-curtis dissimilarity matrix
    diss[[name]] <- avgdist(data_tmp,
                            sample = a)
    # rarefied summarisedExperiment
    data_filt[[name]] <- rarefyAssay(data_filt[[name]],
                sample = a,
                name = "subsampled")
    
    # rarerified into relative abundances
    data_filt[[name]] <- transformAssay(data_filt[[name]],
                                        assay.type = "subsampled",
                                        method = "relabundance")
    
  } else if(b/d >= 0.15) {
    diss[[name]] <- avgdist(data_tmp,
                            sample = b)
    
    data_filt[[name]] <- rarefyAssay(data_filt[[name]],
                                     sample = b,
                                     name = "subsampled")
    
    data_filt[[name]] <- transformAssay(data_filt[[name]],
                                        assay.type = "subsampled",
                                        method = "relabundance")
    
  } else if(c/d >= 0.15){
    diss[[name]] <- avgdist(data_tmp,
                            sample = c)
    
    data_filt[[name]] <- rarefyAssay(data_filt[[name]],
                                     sample = c,
                                     name = "subsampled")
    
    data_filt[[name]] <- transformAssay(data_filt[[name]],
                                        assay.type = "subsampled",
                                        method = "relabundance")
    
  } else{
    diss[[name]] <- c
  }
  
  # rarefied abundances community matrix for capscale
  comm_tmp <- t(assay(data_filt[[name]], "relabundance"))
  
  # make sure meta order of samples same as the dissimilarity matrix
  meta <- data_df %>%
    filter(samplecode %in% rownames(as.matrix(diss[[name]]))) %>%
    arrange(factor(samplecode, levels = rownames(as.matrix(diss[[name]]))))
  
  # betadisper and anova
  betadisper[[name]] <- betadisper(diss[[name]], group = meta$treatment)
  betadisper_anova[[name]] <- anova(betadisper[[name]])
  
  # capscale with the average dissimilarity matrix. For species scores rarefied community matrix given (comm)
  cap[[name]] <- capscale(diss[[name]] ~ G*P,
           data = meta, comm = comm_tmp)
  cap_anova[[name]] <- anova.cca(cap[[name]], by = "terms") %>%
    rownames_to_column(var = "factor")
  
  eig[[name]] <- eigenvals(cap[[name]]) %>% summary()
  
}
)


for(name in names(eig)){
  eig[[name]] <- as.data.frame(eig[[name]]) %>%
    select(CAP2) %>%
    rownames_to_column(var = "items") %>%
    filter(items == "Cumulative Proportion")
}

eig <- bind_rows(eig, .id = "sampletype")

all_betadisp <- as.data.frame(bind_rows(betadisper_anova, .id = "samplecode")) %>%
  filter(!is.na(`Pr(>F)`)) %>%
  remove_rownames() %>%
  dplyr::rename("betadisper_p_value" = "Pr(>F)")

all_rda <- bind_rows(cap_anova, .id = "samplecode")

all_results <- dplyr::left_join(all_rda, select(all_betadisp, samplecode, betadisper_p_value), by = join_by(samplecode))

##### plots #####
# function may be more convenient, but here each produced separately as there was  
# a need to inspect the plots 

###### pot late leaf #####

baseplot <- plot(cap$Potato_Leaf_Late)
arrows <- as.data.frame(cap$Potato_Leaf_Late$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


pot_leaf_late <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$Potato_Leaf_Late) / sum(vegan::eigenvals(cap$Potato_Leaf_Late))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$Potato_Leaf_Late) / sum(vegan::eigenvals(cap$Potato_Leaf_Late))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_pot_leaf_late <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


###### pot late tuber #####

baseplot <- plot(cap$Potato_Tuber_Late)
arrows <- as.data.frame(cap$Potato_Tuber_Late$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


pot_tuber_late <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$Potato_Tuber_Late) / sum(vegan::eigenvals(cap$Potato_Tuber_Late))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$Potato_Tuber_Late) / sum(vegan::eigenvals(cap$Potato_Tuber_Late))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_pot_tuber_late <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


###### faba early leaf #####

baseplot <- plot(cap$`Faba-bean_Leaf_Early`)
arrows <- as.data.frame(cap$`Faba-bean_Leaf_Early`$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


faba_leaf_early <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$`Faba-bean_Leaf_Early`) / sum(vegan::eigenvals(cap$`Faba-bean_Leaf_Early`))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$`Faba-bean_Leaf_Early`) / sum(vegan::eigenvals(cap$`Faba-bean_Leaf_Early`))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_faba_leaf_early <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))



###### faba late leaf #####

baseplot <- plot(cap$`Faba-bean_Leaf_Late`)
arrows <- as.data.frame(cap$`Faba-bean_Leaf_Late`$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


faba_leaf_late <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$`Faba-bean_Leaf_Late`) / sum(vegan::eigenvals(cap$`Faba-bean_Leaf_Late`))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$`Faba-bean_Leaf_Late`) / sum(vegan::eigenvals(cap$`Faba-bean_Leaf_Late`))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_faba_leaf_late <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))



###### faba late root #####

baseplot <- plot(cap$`Faba-bean_Root_Late`)
arrows <- as.data.frame(cap$`Faba-bean_Root_Late`$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


faba_root_late <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$`Faba-bean_Root_Late`) / sum(vegan::eigenvals(cap$`Faba-bean_Root_Late`))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$`Faba-bean_Root_Late`) / sum(vegan::eigenvals(cap$`Faba-bean_Root_Late`))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_faba_root_late <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


###### oat early leaf #####

baseplot <- plot(cap$Oat_Leaf_Early)
arrows <- as.data.frame(cap$Oat_Leaf_Early$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


oat_leaf_early <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$Oat_Leaf_Early) / sum(vegan::eigenvals(cap$Oat_Leaf_Early))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$Oat_Leaf_Early) / sum(vegan::eigenvals(cap$Oat_Leaf_Early))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_oat_leaf_early <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


###### oat late root #####

baseplot <- plot(cap$Oat_Root_Late)
arrows <- as.data.frame(cap$Oat_Root_Late$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)


oat_root_late <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$Oat_Root_Late) / sum(vegan::eigenvals(cap$Oat_Root_Late))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$Oat_Root_Late) / sum(vegan::eigenvals(cap$Oat_Root_Late))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_oat_root_late <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


###### pot early root #####

baseplot <- plot(cap$Potato_Root_Early)
arrows <- as.data.frame(cap$Potato_Root_Early$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)

pot_root_early <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$Potato_Root_Early) / sum(vegan::eigenvals(cap$Potato_Root_Early))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$Potato_Root_Early) / sum(vegan::eigenvals(cap$Potato_Root_Early))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_pot_root_early <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


###### pot late root #####

baseplot <- plot(cap$Potato_Root_Late)
arrows <- as.data.frame(cap$Potato_Root_Late$CCA$biplot) %>%
  mutate(item = rownames(.)) %>%
  mutate(item = str_remove_all(item, pattern = "Glyphosate|Phosphorus"))

species <- data.frame(baseplot$species) %>%
  mutate(ASV_ID = rownames(.),
         dist = sqrt(CAP1^2 + CAP2^2)) %>%
  arrange(desc(dist)) %>%
  head(5) %>%
  dplyr::left_join(rownames_to_column(as.data.frame(ASV_tax), var = "ASV_ID")) %>%
  mutate(name = row.names(.))


species <- species %>%
  mutate(name2 = ifelse(Family == "",
                        paste0(name, ". ", Order),
                        ifelse(Species == "" & Genus == "",
                               paste0(name, ". ", Family),
                               ifelse(Species == "",
                                      paste0(name, ". ", Genus),
                                      paste0(name, ". ", Genus, " ", Species)))))

sites <- data.frame(baseplot$sites)
sites$samplecode <- rownames(sites)
sites <- dplyr::left_join(sites, data_df)
sites$plot <- as.character(sites$plot)

pot_root_late <- ggplot() +
  geom_point(data=sites, aes(x=CAP1, y=CAP2, color=treatment), alpha = 1.5, size = 2.5) +
  #geom_point(data=species, aes(x=CAP1, y=CAP2), size = 2) +
  geom_segment(data = arrows, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length=unit(.2, 'cm'))) +
  ggrepel::geom_text_repel(data = arrows,
                           aes(x= CAP1, y = CAP2, label = item),
                           size = 3,
                           color = "red",
                           max.overlaps = 50,
                           nudge_x = 0.02,
                           nudge_y = -0.005) +
  #stat_conf_ellipse(data=sites, aes(x=CAP1, y=CAP2, fill=plot), alpha = 0.1,
  #                  geom="polygon", bary=FALSE, level = 0.95,
  #                  inherit.aes = FALSE) +
  geom_text(data = species,
            aes(x= CAP1, y = CAP2, label = name),
            size = 3,
            color = "blue") +
  labs(x = paste0("CAP1 (", round((vegan::eigenvals(cap$Potato_Root_Late) / sum(vegan::eigenvals(cap$Potato_Root_Late))*100)[1], 2), "%)"),
       y = paste0("CAP2 (", round((vegan::eigenvals(cap$Potato_Root_Late) / sum(vegan::eigenvals(cap$Potato_Root_Late))*100)[2], 2), "%)")) +
  # scale_alpha_manual(values = c("1" = 1, "2" = 0.3)) +
  theme_light()

species_pot_root_late <- species %>%
  dplyr::left_join(., select(taxas, taxon, ASV_number), by = join_by(ASV_ID == taxon)) %>%
  mutate(name2 = paste0(name2, " ASV ", ASV_number))


#### make plots ####

plot1 <- ggarrange(pot_leaf_late, pot_root_early, pot_root_late, pot_tuber_late, labels = c("a","b", "c", "d"),  common.legend = T, nrow = 2,
                   ncol = 2)

ggsave(path = "../output/", filename = "CAP_plots_all1.pdf", plot = plot1,  width = 10, height = 10, device='pdf', dpi=300)

plot2 <- ggarrange(faba_leaf_early, faba_root_late, oat_leaf_early, oat_root_late, labels = c("a","b", "c", "d"),  common.legend = T, nrow = 2,
                   ncol = 2)
ggsave(path = "../output/", filename = "CAP_plots_all2.pdf", plot = plot2,  width = 10, height = 10, device='pdf', dpi=300)

#### save all_results ####

write_tsv(all_results, "../output/cap_anovas.tsv", na = "")


#### all species tables ####

species_pot_leaf_late <- species_pot_leaf_late %>%
  select(name, name2) %>%
  dplyr::rename("potato_leaf_late" = "name2")

species_pot_tuber_late <- species_pot_tuber_late %>%
  select(name, name2) %>%
  dplyr::rename("potato_tuber_late" = "name2")

species_pot_root_early <- species_pot_root_early %>%
  select(name, name2) %>%
  dplyr::rename("potato_root_early" = "name2")

species_pot_root_late <- species_pot_root_late %>%
  select(name, name2) %>%
  dplyr::rename("potato_root_late" = "name2")

species_faba_leaf_early <- species_faba_leaf_early %>%
  select(name, name2) %>%
  dplyr::rename("faba_leaf_early" = "name2")

species_faba_root_late <- species_faba_root_late %>%
  select(name, name2) %>%
  dplyr::rename("faba_root_late" = "name2")


species_oat_leaf_early <- species_oat_leaf_early %>%
  select(name, name2) %>%
  dplyr::rename("oat_leaf_early" = "name2")

species_oat_root_late <- species_oat_root_late %>%
  select(name, name2) %>%
  dplyr::rename("oat_root_late" = "name2")

all_tables <- dplyr::left_join(species_pot_leaf_late, species_pot_root_early) %>%
  dplyr::left_join(.,species_pot_root_late) %>%
  dplyr::left_join(.,species_pot_tuber_late) %>%
  dplyr::left_join(.,species_faba_leaf_early) %>%
  dplyr::left_join(.,species_faba_root_late) %>%
  dplyr::left_join(.,species_oat_leaf_early) %>%
  dplyr::left_join(.,species_oat_root_late)

final_table <- as.data.frame(t(all_tables))

final_table <- as.data.frame(apply(final_table, 2, function(x) substring(x, first = 4)))

final_table <- final_table %>% rownames_to_column(var = "samplecode") %>%
  filter(samplecode != "name")

colnames(final_table) <- c("samplecode",1:5)

write_tsv(final_table, "../output/CAP_plot_species.tsv")
