############################################
## load libraries
library(tidyverse)  ##needed for general table operations
library(jsonlite)  ##needed for HGNC requests
library(rvest)    ##needed for scraping
library(httr)    ##needed for scraping
library(config) # needed for config loading
library("R.utils")  ## gzip downloaded and result files
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/G02_HPO/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
source("../functions/hpo-functions.R", local = TRUE)
source("../functions/file-functions.R", local = TRUE)
############################################


############################################
## get all children of term neoplasm HP:0002664 and annotating them with name and definition.

current_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

# walk through the ontology tree and add all unique terms descending from
# HP:0002664 (neoplasm)

if (check_file_age("hpo_list_neoplasm", "data/", 1)) {
  hpo_list_neoplasm <- read_csv(get_newest_file("hpo_list_neoplasm", "data"))
} else {
  all_hpo_children_list_neoplasm <- HPO_all_children_from_term("HP:0002664")

  # transform hte list into a tibble
  hpo_list_neoplasm <- all_hpo_children_list_neoplasm %>%
    unlist() %>%
    tibble(`term` = .) %>%
    unique() %>%
    mutate(query_date = current_date)

  write_csv(hpo_list_neoplasm,
    file = paste0("data/hpo_list_neoplasm.",
      current_date,
      ".csv"),
    na = "NULL")

  gzip(paste0("data/hpo_list_neoplasm.", current_date, ".csv"),
    overwrite = TRUE)
}
############################################


############################################
## download all required database sources from HPO and OMIM
# we load and use the results of previous walks through the ontology tree if not older then 1 month

if (check_file_age("phenotype", "data/", 1)) {
  phenotype_hpoa_filename <- get_newest_file("phenotype", "data/")
} else {
  # TODO: compute date only once or somehow in config
  file_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

  # disease ontology annotations from HPO
  # TODO: this should be a config variable
  phenotype_hpoa_url <- "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"

  phenotype_hpoa_filename <- paste0("data/phenotype.",
    file_date,
    ".hpoa")

  download.file(phenotype_hpoa_url, phenotype_hpoa_filename, mode = "wb")

  gzip(phenotype_hpoa_filename,
    overwrite = TRUE)
}

if (check_file_age("omim_genemap2", "data/", 1)) {
  omim_genemap2_filename <- get_newest_file("omim_genemap2", "data/")
} else {
  # OMIM links to genemap2 file needs to be set in config and applied for at
  # https://www.omim.org/downloads
  omim_genemap2_url <- config_vars_proj$omim_genemap2_url

  omim_genemap2_filename <- paste0("data/omim_genemap2.",
    file_date,
    ".txt")

  download.file(omim_genemap2_url, omim_genemap2_filename, mode = "wb")

  gzip(omim_genemap2_filename,
    overwrite = TRUE)
}
############################################


############################################
## search OMIM and Orphanet for HPO terms and filter
phenotype_hpoa <- read_delim(phenotype_hpoa_filename,
    delim = "\t",
  escape_double = FALSE,
    trim_ws = TRUE,
  skip = 4)

omim_genemap2 <- read_delim(omim_genemap2_filename, "\t",
    escape_double = FALSE,
    col_names = FALSE,
    comment = "#",
    trim_ws = TRUE) %>%
  select(Chromosome = X1,
    Genomic_Position_Start = X2,
    Genomic_Position_End = X3,
    Cyto_Location = X4,
    Computed_Cyto_Location = X5,
    MIM_Number = X6,
    Gene_Symbols = X7,
    Gene_Name = X8,
    Approved_Symbol = X9,
    Entrez_Gene_ID = X10,
    Ensembl_Gene_ID = X11,
    Comments = X12,
    Phenotypes = X13,
    Mouse_Gene_Symbol_ID = X14) %>%
  select(Approved_Symbol, Phenotypes) %>%
  separate_rows(Phenotypes, sep = "; ") %>%
  separate(Phenotypes, c("disease_ontology_name", "hpo_mode_of_inheritance_term_name"), "\\), (?!.+\\))") %>%
  separate(disease_ontology_name, c("disease_ontology_name", "Mapping_key"), "\\((?!.+\\()") %>%
  mutate(Mapping_key = str_replace_all(Mapping_key, "\\)", "")) %>%
  separate(disease_ontology_name, c("disease_ontology_name", "MIM_Number"), ", (?=[0-9][0-9][0-9][0-9][0-9][0-9])") %>%
  mutate(Mapping_key = str_replace_all(Mapping_key, " ", "")) %>%
  mutate(MIM_Number = str_replace_all(MIM_Number, " ", "")) %>%
  filter(!is.na(MIM_Number))  %>%
  filter(!is.na(Approved_Symbol))  %>%
  mutate(disease_ontology_id = paste0("OMIM:", MIM_Number)) %>%
  separate_rows(hpo_mode_of_inheritance_term_name, sep = ", ") %>%
  mutate(hpo_mode_of_inheritance_term_name = str_replace_all(hpo_mode_of_inheritance_term_name, "\\?", "")) %>%
  select(-MIM_Number) %>%
  unique() %>%
  mutate(hpo_mode_of_inheritance_term_name = case_when(hpo_mode_of_inheritance_term_name == "Autosomal dominant" ~ "Autosomal dominant inheritance",
    hpo_mode_of_inheritance_term_name == "Autosomal recessive" ~ "Autosomal recessive inheritance",
    hpo_mode_of_inheritance_term_name == "Digenic dominant" ~ "Digenic inheritance",
    hpo_mode_of_inheritance_term_name == "Digenic recessive" ~ "Digenic inheritance",
    hpo_mode_of_inheritance_term_name == "Isolated cases" ~ "Sporadic",
    hpo_mode_of_inheritance_term_name == "Mitochondrial" ~ "Mitochondrial inheritance",
    hpo_mode_of_inheritance_term_name == "Multifactorial" ~ "Multifactorial inheritance",
    hpo_mode_of_inheritance_term_name == "Pseudoautosomal dominant" ~ "X-linked dominant inheritance",
    hpo_mode_of_inheritance_term_name == "Pseudoautosomal recessive" ~ "X-linked recessive inheritance",
    hpo_mode_of_inheritance_term_name == "Somatic mosaicism" ~ "Somatic mosaicism",
    hpo_mode_of_inheritance_term_name == "Somatic mutation" ~ "Somatic mutation",
    hpo_mode_of_inheritance_term_name == "X-linked" ~ "X-linked inheritance",
    hpo_mode_of_inheritance_term_name == "X-linked dominant" ~ "X-linked dominant inheritance",
    hpo_mode_of_inheritance_term_name == "X-linked recessive" ~ "X-linked recessive inheritance",
    hpo_mode_of_inheritance_term_name == "Y-linked" ~ "Y-linked inheritance"))

phenotype_hpoa_filter <- phenotype_hpoa %>%
   filter(hpo_id %in% hpo_list_neoplasm$term) %>%
   select(database_id, hpo_id) %>%
   unique() %>%
  group_by(database_id) %>%
  summarise(neoplasm_hpo_terms = paste(hpo_id, collapse = " | "),
    .groups = "keep") %>%
  ungroup()

hpo_gene_list <- phenotype_hpoa_filter %>%
  left_join(omim_genemap2, by = c("database_id" = "disease_ontology_id")) %>%
  filter(!is.na(Approved_Symbol)) %>%
  mutate(database_and_hpo_id = paste0(database_id, " (", neoplasm_hpo_terms, ")")) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Approved_Symbol)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(-disease_ontology_name, -Mapping_key, -hpo_mode_of_inheritance_term_name) %>%
  group_by(approved_symbol) %>%
  summarise(Approved_Symbol = paste(unique(Approved_Symbol), collapse = " | "),
    hgnc_id = paste(unique(hgnc_id), collapse = " | "),
    neoplasm_hpo_terms = paste(neoplasm_hpo_terms, collapse = "; "),
    database_id = paste(database_id, collapse = "; "),
    database_and_hpo_id = paste(database_and_hpo_id, collapse = "; "),
    source_count = n(),
    .groups = "keep") %>%
  ungroup() %>%
  mutate(at_least_one_database = (source_count > 0)) %>%
  select(approved_symbol,
    hgnc_id,
    gene_name_reported = Approved_Symbol,
    source = database_and_hpo_id,
    source_count,
    source_evidence = at_least_one_database)

# TODO: normalize source_evidence to 0/1 as percentiles
# TODO: write a function for this normalization step

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(hpo_gene_list,
  file = paste0("results/G02_HPO_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/G02_HPO_genes.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
