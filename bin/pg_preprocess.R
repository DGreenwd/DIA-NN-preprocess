require(data.table)
require(stringr)
require(dplyr)
require(tidyr)
require(purrr)
require(tibble)

read_pg_matrix = function(filename,
                          col_suffix_trim = ".d",
                          col_prefix_trim = "E:\\\\Long Covid\\\\",
                          expected_cols = c("Protein.Group","Protein.Names","Genes","First.Protein.Description")){
  raw = data.table::fread(filename,check.names = F)
  if(!all(expected_cols%in%colnames(raw))){
    print(paste0("Warning, expected columns absent from input: ",setdiff(expected_cols,colnames(raw))))
  }
  
  
  if(!is.null(col_suffix_trim)){
    colnames(raw) = stringr::str_remove(colnames(raw),paste0(col_suffix_trim,"$"))
  }
  
  if(!is.null(col_prefix_trim)){
    colnames(raw) = str_replace(colnames(raw),paste0("^",col_prefix_trim),"")
  }
  
  return(raw)
  
}

create_protein_dictionary <-
  function(prot_raw_list, pg_cols) {
    protein_dictionary <-
      bind_rows(map(prot_raw_list, ~ select(.x, all_of(pg_cols)))) %>%
      distinct()
    
    # Genes contains a gene per protein in group, extract first gene
    protein_dictionary <- protein_dictionary %>%
      mutate(Genes_first = str_split_fixed(Genes, ";", 2)[, 1])
    
    # Map gene names (first in group) to ENTREZID
    protein_dictionary <-
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        unique(protein_dictionary$Genes_first),
        'ENTREZID',
        'SYMBOL'
      ) %>%
      bind_rows(., id = "Genes") %>%
      pivot_longer(
        .,
        cols = dplyr::everything(),
        names_to = "Genes_first",
        values_to = "ENTREZID"
      ) %>%
      left_join(protein_dictionary, ., by = "Genes_first")
    
    if (protein_dictionary %>% filter(is.na(ENTREZID)) %>% nrow() > 0) {
      print(
        paste0(
          "Number of protein groups not mapped to an ENTREZID by GENE SYMBOL: ",
          protein_dictionary %>%
            filter(is.na(ENTREZID)) %>%
            nrow()
        )
      )
      
      print("Mapping remaining based on protein UNIPROT")
      
      
      # Map id for the first gene associated with each protein to ENTREZID
      protein_dictionary <-
        AnnotationDbi::mapIds(
          org.Hs.eg.db::org.Hs.eg.db,
          unique(
            protein_dictionary %>% filter(is.na(ENTREZID)) %>%  pull(Protein.Group)
          ),
          'ENTREZID',
          'UNIPROT'
        ) %>%
        bind_rows(.) %>%
        pivot_longer(
          .,
          cols = dplyr::everything(),
          names_to = "Protein.Group",
          values_to = "ENTREZID"
        ) %>%
        filter(!is.na(ENTREZID)) %>%
        rows_update(protein_dictionary, ., by = "Protein.Group")
      
      if (protein_dictionary %>% filter(is.na(ENTREZID)) %>% nrow() > 0) {
        print(
          paste0(
            "Number of protein groups not mapped to an ENTREZID by GENE SYMBOL or UNIPROT: ",
            protein_dictionary %>%
              filter(is.na(ENTREZID)) %>%
              nrow()
          )
        )
      }
      
    }
    return(protein_dictionary)
  }


extract_colData <-
  function(prot_raw = NULL,
           meta_batch = NULL,
           pg_cols = NULL,
           sample_string_remove = "^RSM2") {
    
    raw_cols = colnames(prot_raw)
    # Extract samples and batch
    raw_cols = setdiff(raw_cols, pg_cols)
    sample_barcode = stringr::str_split_fixed(raw_cols, "_", n = 2)[, 1]
    batch = stringr::str_split_fixed(raw_cols, "_", n = 2)[, 2]
    batch = stringr::str_split_fixed(batch, "-", n = 2)[, 1]
    
    # Manual error fix: one sample uses an RSM2 prefix which causes issues later on
    if (!is.null(sample_string_remove)) {
      sample_barcode = gsub(sample_string_remove, "RSM", sample_barcode)
    }
    
    if (!is.null(prefix)) {
      sample_barcode_core =  stringr::str_remove(sample_barcode, pattern = prefix)
    }
    if (!is.null(suffix)) {
      sample_barcode_core =  stringr::str_remove(sample_barcode_core, pattern = suffix)
    }
    
    sample_data <-       
      data.frame(
        filename = raw_cols,
        sample_barcode = sample_barcode,
        sample_barcode_core = sample_barcode_core,
        batch = batch,
        meta_batch = meta_batch
      )
    
    ### Count NA and non-NA values 
    sample_data <- 
      prot_raw %>%
      select(-c(pg_cols)) %>%
      summarise(across(everything(), ~ sum(!is.na(.)))) %>%
      pivot_longer(everything(), names_to = "filename", values_to = "non_na_count") %>%
      left_join(sample_data, ., by = "filename")
    
    
    sample_data <- 
      prot_raw %>%
      select(-c(pg_cols)) %>%
      summarise(across(everything(), ~ sum(is.na(.)))) %>%
      pivot_longer(everything(), names_to = "filename", values_to = "na_count") %>%
      left_join(sample_data, ., by = "filename")
    
    return(sample_data)
    
  }

pg_non_na_by_group <- 
  function(prot_raw,colData,pg_cols){
    
    prot_raw %>%
      
      ## Pivot longer - value per filename 
      pivot_longer(-pg_cols, names_to = "filename", values_to = "value") %>%
      
      ## Join sample info - filename and group (i.e. study or experimental group)
      left_join(colData %>% select(filename,group), by = "filename") %>%
      
      ## Explicit missing values 
      # Expand each data set so it has all proteins groups for all samples (not just those in its metabatch )
      complete(filename, nesting(Protein.Group = protein_dictionary %>% pull(Protein.Group))) %>% 
      
      ## Count non-missing values and total
      group_by(Protein.Group, group) %>%
      summarise(non_na_count = sum(!is.na(value)),
                total_count = n(),
                .groups = "drop") %>%
      
      ## Proportion of total
      mutate(non_na_proportion = non_na_count / total_count) %>% 
      
      ## Pivot wider
      pivot_wider(names_from = group, values_from = non_na_proportion, id_cols = Protein.Group,
                  names_glue = "{group}_{.value}")
  }

library(tibble)
manual_impute <- function(mat, scale = 0.3, shift = 1.8, seed = NULL) {
  # https://rdrr.io/bioc/DEP/src/R/functions.R#sym-manual_impute
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  if(is.integer(scale)) scale <- is.numeric(scale)
  if(is.integer(shift)) shift <- is.numeric(shift)
  
  
  # Get descriptive parameters of the current sample distributions
  stat <- mat %>%
    data.frame(check.names = F) %>%
    rownames_to_column() %>%
    gather(samples, value, -rowname) %>%
    filter(!is.na(value))  %>%
    group_by(samples) %>%
    summarise(mean = mean(value),
              median = median(value),
              sd = sd(value),
              n = n(),
              infin = nrow(mat) - n)
  # Impute missing values by random draws from a distribution
  # which is left-shifted by parameter 'shift' * sd and scaled by parameter 'scale' * sd.
  for (a in seq_len(nrow(stat))) {
    mat[is.na(mat[, stat$samples[a]]), stat$samples[a]] <-
      rnorm(stat$infin[a],
            mean = stat$median[a] - shift * stat$sd[a],
            sd = stat$sd[a] * scale)
  }
  return(mat)
}


