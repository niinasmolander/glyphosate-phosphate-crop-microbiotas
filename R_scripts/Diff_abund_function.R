five_DAA_estimators <- function(tsePrevalent_plant, data_df_plant){
  # Aldex
  
  tsePrevalent_plant_OTU <- tsePrevalent_plant
  
  data_split <- list()
  for(t in unique(colData(tsePrevalent_plant_OTU)$treatment)){
    for(tr in unique(colData(tsePrevalent_plant_OTU)$treatment)){
      if(t != tr & !(paste0(tr, "AND", t) %in% names(data_split))){
        name <- paste0(t, "AND", tr)
        data_split[[name]] <- tsePrevalent_plant_OTU[, colData(tsePrevalent_plant_OTU)$treatment == t | colData(tsePrevalent_plant_OTU)$treatment == tr]
      }
    }
  }
  
  data_split_df <- list()
  for(i in c(1:length(data_split))){
    name <- names(data_split[i])
    data_split_df[[name]] <- assay(data_split[[i]], "counts") %>%
      t() %>%
      as.data.frame() %>%
      merge(data_df_plant, ., by=0) %>%
      column_to_rownames(var = "Row.names") %>%
      dplyr::select(!(library:plot) & !(samplecode:timepoint)) %>%
      arrange(treatment) %>%
      t() %>%
      as.data.frame()
  }
  
  conditions <- list()
  for(i in c(1:length(data_split_df))){
    name <- names(data_split_df[i])
    conditions[[name]] <- as.character(data_split_df[[name]]["treatment",])
    data_split_df[[name]] <- data_split_df[[name]][!(row.names(data_split_df[[name]]) %in% "treatment"),]
    data_split_df[[name]] <- data_split_df[[name]] %>% mutate_all(as.numeric)
  }
  
  aldex_results <- list()
  for(i in names(data_split_df)){
    aldex_results[[i]] <- aldex(data_split_df[[i]], conditions[[i]], denom = "all",
                                test = "t", effect = TRUE, paired.test = F, verbose = T)
  }
  
  # Ancombc2
  
  tsePrevalent_plant_OTU <- tsePrevalent_plant
  
  data_split <- list()
  for(t in unique(colData(tsePrevalent_plant_OTU)$treatment)){
    for(tr in unique(colData(tsePrevalent_plant_OTU)$treatment)){
      if(t != tr & !(paste0(tr, "AND", t) %in% names(data_split))){
        name <- paste0(t, "AND", tr)
        data_split[[name]] <- tsePrevalent_plant_OTU[, colData(tsePrevalent_plant_OTU)$treatment == t | colData(tsePrevalent_plant_OTU)$treatment == tr]
      }
    }
  }
  
  
  ancombc2_results <- list()
  for(i in names(data_split)){
    print(i)
    ancombc2_results[[i]] <- ancombc2(data_split[[i]], 
                                      assay_name = "counts",
                                      fix_formula = "treatment",
                                      group = "treatment", verbose = TRUE,
                                      p_adj_method = "BH",
                                      prv_cut = 0)
  }
  
  
  # ZicoSeq
  
  tsePrevalent_plant_OTU <- tsePrevalent_plant
  
  data_split <- list()
  for(t in unique(colData(tsePrevalent_plant_OTU)$treatment)){
    for(tr in unique(colData(tsePrevalent_plant_OTU)$treatment)){
      if(t != tr & !(paste0(tr, "AND", t) %in% names(data_split))){
        name <- paste0(t, "AND", tr)
        data_split[[name]] <- tsePrevalent_plant_OTU[, colData(tsePrevalent_plant_OTU)$treatment == t | colData(tsePrevalent_plant_OTU)$treatment == tr]
      }
    }
  }
  
  data_split_df <- list()
  meta_split <- list()
  for(i in c(1:length(data_split))){
    name <- names(data_split[i])
    data_split_df[[name]] <- assay(data_split[[i]], "counts")
    data_split_df[[name]] <- data_split_df[[name]][base::rowSums(data_split_df[[name]]) != 0,]
    data_split_df[[name]] <- data_split_df[[name]][,base::colSums(data_split_df[[name]]) != 0]
    meta_split[[name]] <- data_df %>%
      filter(samplecode %in% colnames(data_split_df[[name]]))
    data_split_df[[name]] <- data_split_df[[name]] %>%
      as.data.frame(.) %>%
      relocate(base::rownames(meta_split[[name]])) %>%
      as.matrix(.)
  }
  
  
  ZicoSeq_results <- list()
  for(i in names(data_split_df)){
    print(i)
    ZicoSeq_results[[i]] <- ZicoSeq(meta.dat = meta_split[[i]], feature.dat = data_split_df[[i]],
                                    feature.dat.type = "count",
                                    grp.name = 'treatment', prev.filter = 0, is.winsor = FALSE,
                                    verbose = TRUE)
  }
  
  # Dacomp
  
  tsePrevalent_plant_OTU <- tsePrevalent_plant
  
  data_split <- list()
  for(t in unique(colData(tsePrevalent_plant_OTU)$treatment)){
    for(tr in unique(colData(tsePrevalent_plant_OTU)$treatment)){
      if(t != tr & !(paste0(tr, "AND", t) %in% names(data_split))){
        name <- paste0(t, "AND", tr)
        data_split[[name]] <- tsePrevalent_plant_OTU[, colData(tsePrevalent_plant_OTU)$treatment == t | colData(tsePrevalent_plant_OTU)$treatment == tr]
      }
    }
  }
  
  data_split_df <- list()
  meta_split <- list()
  dsr <- list()
  for(i in c(1:length(data_split))){
    name <- names(data_split[i])
    data_split_df[[name]] <- assay(data_split[[i]], "counts")
    data_split_df[[name]] <- data_split_df[[name]][base::rowSums(data_split_df[[name]]) != 0,]
    data_split_df[[name]] <- data_split_df[[name]][,base::colSums(data_split_df[[name]]) != 0] %>%
      t()
    
    meta_split[[name]] <- data_df %>%
      filter(samplecode %in% rownames(data_split_df[[name]])) %>%
      dplyr::select(samplecode, treatment)
    meta_split[[name]] <- meta_split[[name]][base::order(base::match(meta_split[[name]]$samplecode, base::rownames(data_split_df[[name]]))),]
    meta_split[[name]] <- base::as.vector(meta_split[[name]]$treatment)
    
    dsr[[name]] <- dacomp.select_references(data_split_df[[name]])
  }
  
  Dacomp_results <- list()
  for(i in names(data_split_df)){
    Dacomp_results[[i]] <- dacomp.test(data_split_df[[i]], meta_split[[i]], dsr[[i]],
                                       test = DACOMP.TEST.NAME.WILCOXON, disable_DSFDR = TRUE,
                                       verbose = TRUE)
  }
  
  
  dacomp_all <- list()
  for(i in names(Dacomp_results)){
    dacomp_all[[i]] <- data_split_df[[i]] %>% t() %>% base::rownames() %>%
      bind_cols(., Dacomp_results[[i]]$p.values.test.adjusted)
  }
  
  # eBay
  
  tsePrevalent_plant_OTU <- tsePrevalent_plant
  
  # split
  data_split <- list()
  for(t in unique(colData(tsePrevalent_plant_OTU)$treatment)){
    for(tr in unique(colData(tsePrevalent_plant_OTU)$treatment)){
      if(t != tr & !(paste0(tr, "AND", t) %in% names(data_split))){
        name <- paste0(t, "AND", tr)
        data_split[[name]] <- tsePrevalent_plant_OTU[, colData(tsePrevalent_plant_OTU)$treatment == t | colData(tsePrevalent_plant_OTU)$treatment == tr]
      }
    }
  }
  
  #create an OTU table with rows as samples and columns as taxa
  
  data_split_df <- list()
  for(i in c(1:length(data_split))){
    name <- names(data_split[i])
    data_split_df[[name]] <- assay(data_split[[i]], "counts")
    data_split_df[[name]] <- data_split_df[[name]][,colSums(data_split_df[[name]]) != 0]
    data_split_df[[name]] <- data_split_df[[name]][rowSums(data_split_df[[name]]) != 0,] %>%
      t() %>%
      as.data.frame() %>%
      merge(data_df_plant, ., by=0) %>%
      column_to_rownames(var = "Row.names") %>%
      dplyr::select(!(library:plot) & !(samplecode:timepoint)) %>%
      arrange(treatment) %>%
      as.data.frame()
  }
  
  # treatment conditions
  conditions <- list()
  for(i in c(1:length(data_split_df))){
    name <- names(data_split_df[i])
    conditions[[name]] <- dplyr::select(data_split_df[[name]], treatment)
    swap <- cbind(base::unique(conditions[[name]]), c(0,1))
    rownames(swap) <- c(0,1)
    colnames(swap) <- c("treatment", "condition")
    conditions[[name]] <- dplyr::left_join(conditions[[name]], swap)
    
    data_split_df[[name]] <- data_split_df[[name]] %>%
      dplyr::select(-c("treatment"))
  }
  
  
  eBay_results <- list()
  for(name in names(data_split_df)){
    eBay_results[[name]] <- eBay(otu.data=data_split_df[[name]],
                                 group=conditions[[name]]$condition,
                                 test.method="wilcoxon",
                                 cutf=0.05,
                                 adj.m="BH")
  }
  
  
  eBay_pvalues <- list()
  for(name in names(eBay_results)){
    eBay_pvalues[[name]] <- as.data.frame(eBay_results[[name]]$final.p) %>%
      rownames_to_column(var = "taxon")
  }
  
  eBay_pvalues <- bind_rows(eBay_pvalues, .id = "Comparison")
  colnames(eBay_pvalues) <- c("Comparison", "taxon", "p-value.eBay")
  
  
  return(lst("aldex" = aldex_results, 
             "ancombc2" = ancombc2_results,
             "zicoseq" = ZicoSeq_results,
             "dacomp" = dacomp_all,
             "eBay" = eBay_pvalues))
  
}
