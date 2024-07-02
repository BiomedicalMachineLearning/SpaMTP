
#' This the function used to compute the exact fisher test for over-represntation based pathway analysis
#'
#' @param seurat A seurat object contains spatial metabolomics/transcriptolomics features or both.
#' @param polarity The polarity of the MALDI experiment
#' @param ... The arguments pass to FisherexactTest
#'
#' @return A dataframe with the relevant pathway information
#' @export
#'
pathway_analysis = function(seurat,
                            polarity,
                            ...){
  met_analytes = row.names(seurat@assays[["SPM"]]@features)
  rna_analytes = row.names(seurat@assays[["SPT"]]@features)
  analytes =  list(mz = met_analytes,
                   genes = paste0("gene_symbol:",toupper(rna_analytes)))
  result = FisherexactTest(analytes,
                           polarity = polarity,...)
  return(result)
}

#' Calculates Significant Metabolic Pathways using a Fisher Exact Test
#'
#' @param Analyte A list of analytes with 3 elements, namely "mz", "genes" and "metabolites", each is comprised of the corresponding labels, for metabolites,
#' Supported metabolites format including, X stands for upper case of the cooresponding ID in each database: "hmdb:HMDBX", "chebi:X", "pubchem:X","wikidata:X" ,"kegg:X" ,"CAS:X","lipidbank:X","chemspider:X","	LIPIDMAPS:X"
#' Supported gene data format including: "entrez:X", "gene_symbol:X", "uniprot:X", "ensembl:X", "hmdb:HMDBPX"
#' Supported mz format: any string or numeric vector contains the m/z
#' @param analyte_type = "metabolites" or "gene" or "mz", or a vector contains any combinations of them
#' @param max_path_size The max number of metabolites in a specific pathway
#' @param min_path_size The min number of metabolites in a specific pathway
#' @param alternative The hypothesis of the fisher exact test
#' @param polarity The polarity of the MALDI test
#' @param pathway_all_info Whether to included all genes/metabolites screened in the return
#' @param pval_cutoff The cut off of raw p value to retain the pathways

#' @return a dataframe with the relevant pathway information

FisherexactTest = function (Analyte,
                            analyte_type = c("mz","genes"),
                            polarity = NULL,
                            ppm_error = 10,
                            max_path_size = 500,
                            alternative = "greater",
                            min_path_size = 5,
                            pathway_all_info = F,
                            pval_cutoff = 0.05)
{
  require(dplyr)
  require(stringr)

  if((!"mz" %in% analyte_type) & (!"metabolites" %in% analyte_type) & (!"genes" %in% analyte_type)){
    stop(
      "analyte_type was not specified correctly.  Please specify one of the following options: metabolites, genes"
    )
  }
  now <- proc.time()
  print("Fisher Testing ......")
  pathwayRampId <- rampId <- c()
  # Get the RaMP ids for metabolites/genes
  convert_to_rows <- function(row,
                              pattern) {
    identifiers = data.frame()
    for (i in which(grepl(row,
                          pattern = pattern))) {
      temp = unlist(strsplit(unlist(row[i]), pattern))
      identifiers <- rbind(identifiers,
                           temp)
    }
    identifiers = t(identifiers)
    colnames(identifiers) = colnames(row)[which(grepl(row,
                                                      pattern = pattern) ==
                                                  T)]
    return(cbind(row[which(grepl(row,
                                 pattern = pattern) == F)][rep(1, times = nrow(identifiers)), ],
                 identifiers))
  }

  print("Loading files ......")
  analytehaspathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/analytehaspathway.rds"))
  pathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/pathway.rds"))
  source = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/source.rds"))
  chem_props = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/chem_props.rds"))
  analyte = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/analyte.rds"))
  pathway = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/pathway.rds"))

  print("Loading files finished!")

  if ( "metabolites" %in% analyte_type) {
    analytes_met = Analyte[["metabolites"]]
    source_met = source[which(grepl(source$rampId, pattern = "RAMP_C") == T),]
    analytehaspathway_met = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_C") == T),]
    analyte_met = analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),]
    adduct_file = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/adduct_file.rds"))
  }
  if ("genes" %in% analyte_type) {
    analytes_rna = Analyte[["genes"]]
    source_rna = source[which(grepl(source$rampId, pattern = "RAMP_G") == T),]
    analytehaspathway_rna = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_G") == T),]
    analyte_rna = analyte[which(grepl(analyte$rampId, pattern = "RAMP_G") == T),]
  }

  if ("mz" %in% analyte_type) {
    analytes_mz = Analyte[["mz"]]
    adduct_file = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/adduct_file.rds"))
    load(paste0(dirname(system.file(package = "SpaMTP")), "/data/Chebi_db.rda"))
    load(paste0(dirname(system.file(package = "SpaMTP")), "/data/Lipidmaps_db.rda"))
    load(paste0(dirname(system.file(package = "SpaMTP")), "/data/HMDB_db.rda"))
    # Since this file was tested in positive ion mode
    db = rbind(Chebi_db,
               Lipidmaps_db,
               HMDB_db)
    rm(Chebi_db)
    rm(Lipidmaps_db)
    rm(HMDB_db)
    gc()
    if (polarity == "Positive") {
      test_add_pos <- adduct_file$adduct_name[which(adduct_file$charge > 0)]
      # Using Chris' pipeline for annotation
      # 1) Filter DB by adduct.
      db_1 <- db_adduct_filter(db, test_add_pos, polarity = "pos")
    } else if (polarity == "Negative") {
      test_add_neg <- adduct_file$adduct_name[which(adduct_file$charge < 0)]
      # Using Chris' pipeline for annotation
      # 1) Filter DB by adduct.
      db_1 <- db_adduct_filter(db, test_add_neg, polarity = "neg")
    } else if (polarity == "Neutral") {
      # Using Chris' pipeline for annotation
      # 1) Filter DB by adduct.
      db_1 <- db %>% mutate("M" = `M-H ` + 1.007276)
    }else{
      stop("Please enter correct polarity from: 'Positive', 'Negative', 'Neutral'")
    }

    # 2) only select natural elements
    db_2 <- formula_filter(db_1)

    # 3) search db against mz df return results
    print("search db against mz df return results")
    input_mz = data.frame(cbind(
      row_id = 1:length(analytes_mz),
      mz = as.numeric(str_extract(analytes_mz, pattern = "\\d+\\.?\\d*"))
    ))
    ppm_error <- ppm_error
    db_3 <- proc_db(data.frame(input_mz), db_2, ppm_error)
    # Expand the isomer entries
    print("Expanding database to extract all potential metabolites")
    db_3list = pblapply(1:nrow(db_3), function(i){
      if (any(grepl(db_3[i, ], pattern = ";"))) {
        # Only take the first row
        ids = unlist(str_split(db_3[i, ]$Isomers, pattern = ";"))
        ids[which(grepl(ids, pattern = "HMDB"))] = paste0("HMDB:",ids[which(grepl(ids, pattern = "HMDB"))])
        ids[which(grepl(ids, pattern = "LM"))] = paste0("LIPIDMAPS:",ids[which(grepl(ids, pattern = "LM"))])
        return()
      } else{
        return(db_3[i,]$Isomers)
      }
    })
    expand_db3 = do.call(c,db_3list)
    analytes_mz = sub(" ", "", expand_db3)
    source_mz = source[which(grepl(source$rampId, pattern = "RAMP_C") == T),]
    analytehaspathway_mz = analytehaspathway[which(grepl(analytehaspathway$rampId, pattern = "RAMP_C") == T),]
  }

  print("Parsing the information of given analytes class")
  analyte_new = analytehaspathway_new = source_new =  data.frame()
  analytes_new = c()
  if("mz" %in% analyte_type){
    analyte_new = rbind(analyte_new,
                        analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),])
    analytehaspathway_new = rbind(analytehaspathway_new,
                                  analytehaspathway_mz)
    source_new = rbind(source_new,
                       source_mz)
    analytes_new = c(analytes_new,
                     analytes_mz)
  }

  if("metabolites" %in% analyte_type){
    analyte_new = rbind(analyte_new,
                        analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),])
    analytehaspathway_new = rbind(analytehaspathway_new,
                                  analytehaspathway_met)
    source_new = rbind(source_new,
                       source_met)
    analytes_new = c(analytes_new,
                     analytes_met)
  }

  if("genes" %in% analyte_type){
    analyte_new = rbind(analyte_new,
                        analyte[which(grepl(analyte$rampId, pattern = "RAMP_C") == T),])

    analytehaspathway_new = rbind(analytehaspathway_new,
                                  analytehaspathway_rna)
    source_new = rbind(source_new,
                       source_rna)
    analytes_new = c(analytes_new,
                     analytes_rna)
  }

  analytehaspathway_new = unique(analytehaspathway_new)
  source_new = unique(source_new)
  analytes_new = unique(analytes_new)

  # Pathway enrichment

  ############ Metabolites pathway analysis ##############
  print("Begin metabolic pathway analysis ......")
  analytes_rampids = c()
  unique_analytes = unique(analytes_new)
  analytes_rampids = unique(source_new$rampId[which(tolower(source_new$sourceId) %in% tolower(unique_analytes))])
  # for(k in 1:length(unique_analytes)){
  #   pattern = unique_analytes[k]
  #   analytes_rampids = c(analytes_rampids,
  #                        unique(source$rampId[which(grepl(source$sourceId, pattern = pattern,
  #                                                         ignore.case = T))]))
  # }
  # analytes_rampids = unique(na.omit(analytes_rampids))
  # (1) Get candidate pathways
  # Get all analytes and number of analytes within a specific pathway
  pathway_db = readRDS(paste0(dirname(system.file(package = "SpaMTP")), "/data/pathway.rds"))
  source_non_duplicated = source_new[which(!duplicated(source_new$rampId)&(source_new$rampId %in% analytes_rampids)),]
  # rampid = the subset of the database with our query data
  pathway_rampids = analytehaspathway_new[which(analytehaspathway_new$rampId %in% analytes_rampids),]
  pathway_rampids_count = pathway_rampids %>% group_by(pathwayRampId) %>% dplyr::mutate(analytes_in_pathways  = n())
  # analytespathway_new. = the subset of the database with all pathways
  analytehaspathway_full = analytehaspathway_new %>%
    group_by(pathwayRampId) %>% dplyr::mutate(total_in_pathways = n())
  analytehaspathway_full =analytehaspathway_full[which(analytehaspathway_full$total_in_pathways>= min_path_size & analytehaspathway_full$total_in_pathways <= max_path_size),]

  # Generate a dataframe contains: the list of metabolites IDs, the list of metabolites names, the number of elements in pathway, the number of elements in our dataset, for each pathway
  if(pathway_all_info == T){
    enrichment_df = pblapply(1:length(unique(pathway_rampids_count$pathwayRampId)), function(x){
      #ananlytes_id_df = analytehaspathway_new[which(analytehaspathway_new$pathwayRampId == unique(pathway_rampids$pathwayRampId)[x]),]
      pathway_id = unique(pathway_rampids_count$pathwayRampId)[x]
      pathway_info = pathway_db[which(pathway_db$pathwayRampId == pathway_id),]
      # get rampids associated with the pathway
      full_list = analytehaspathway_full[which(analytehaspathway_full$pathwayRampId == pathway_id),]
      screened_List = pathway_rampids_count[which(pathway_rampids_count$pathwayRampId == pathway_id),]
      # Get screened index
      source_index_met = which(source_non_duplicated$rampId %in% screened_List$rampId[which(grepl(screened_List$rampId,
                                                                                                  pattern = "RAMP_C_"))])
      source_index_gene = which(source_non_duplicated$rampId %in% screened_List$rampId[which(grepl(screened_List$rampId,
                                                                                                   pattern = "RAMP_G_"))])
      #met
      ananlytes_name_list_met = list(source_non_duplicated$commonName[source_index_met])
      ananlytes_id_list_met = list(source_non_duplicated$sourceId[source_index_met])
      #gene
      ananlytes_name_list_gene = list(source_non_duplicated$commonName[source_index_gene])
      ananlytes_id_list_gene = list(source_non_duplicated$sourceId[source_index_gene])

      analytes_in_pathways = screened_List$count[1]
      total_in_pathways = full_list$count[1]
      return_df = data.frame(pathway_name = pathway_info$pathwayName,
                             pathway_id = pathway_info$sourceId) %>% mutate(metabolite_name_list=ananlytes_name_list_met) %>%
        mutate(metabolite_id_list= ananlytes_id_list_met) %>% mutate(total_in_pathways = total_in_pathways) %>%
        mutate(gene_name_list = ananlytes_name_list_gene) %>% mutate(gene_id_list = ananlytes_id_list_gene) %>%
        mutate(analytes_in_pathways = analytes_in_pathways)
      return(return_df)
    })
  }else{
    print("Merging datasets")
    enrichment_df = merge(pathway_rampids_count[which(!duplicated(pathway_rampids_count$pathwayRampId)),], analytehaspathway_full[which(!duplicated(analytehaspathway_full$pathwayRampId)),],
                          by = "pathwayRampId")
    #colnames(enrichment_df)[which(colnames(enrichment_df) == "count")] = "total_in_pathways"
    enrichment_df = merge(enrichment_df, pathway_db, by = "pathwayRampId")
    enrichment_df = enrichment_df %>% filter(!duplicated(pathwayName))
  }

  # abbb = analytehaspathway_new[which(analytehaspathway_new$rampId %in% analytes_rampids)[1:100],] %>% rowwise() %>%
  # dplyr::mutate(ananlytes_id_list = list(analytehaspathway_new$rampId[which(analytehaspathway_new$pathwayRampId == pathwayRampId)])) %>%
  # dplyr::count(pathwayRampId,
  #              ananlytes_id_list,
  #              sort = T,
  #              name = "analytes_in_pathways")
  print("Running test")
  # (2) Conduct pathway enrichment
  total_analytes = length(unique(pathway_rampids_count$rampId))
  total_in_selected_pathways = length(unique(analytehaspathway_full$rampId))
  print("Calculating p value......")

  enrichment_df = enrichment_df %>% rowwise() %>% mutate(p_val = fisher.test(matrix(
    c(
      analytes_in_pathways,
      # Detected metabolites in pathway
      total_analytes -
        analytes_in_pathways,
      # Detected metabolites not in pathway
      total_in_pathways -
        analytes_in_pathways,
      # Pathway elements not detected
      total_in_selected_pathways -
        total_in_pathways - total_analytes + analytes_in_pathways
    ),
    2,
    2
  ),
  alternative = alternative)$p.value)
  enrichment_df = cbind(enrichment_df,
                        fdr = p.adjust(enrichment_df$p_val, method = "fdr")) %>% mutate(background_analytes_number = total_in_selected_pathways)
  print("P value obtained")
  # (5) Append pathway information to the original df
  gc()
  # (6) Append metabolites information to the original df
  # Paste back the original Ids
  # (7) Reduce the dataframe with respected to the User input pathway size

  print("Done")
  rm(chem_props)
  gc()
  #return(enrichment_df_with_both_info %>% select(-c(
  #  pathwayRampId,
  #  ananlytes_id_list,
  #  screened_analytes
  #)))
  return =enrichment_df %>% select(-c(pathwayRampId,rampId.y, pathwaySource.y)) %>% filter(p_val<=pval_cutoff)
  return(return)
}
