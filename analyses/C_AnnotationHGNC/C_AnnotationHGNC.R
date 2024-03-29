############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(biomaRt)    ## needed to get gene coordinates
library(STRINGdb)  ## needed to compute StringDB identifiers
library("R.utils")  ## gzip downloaded and result files
library(config)     ## needed for config loading
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/C_AnnotationHGNC/"

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
source("../functions/ensembl-functions.R", local = TRUE)
source("../functions/gnomad-functions.R", local = TRUE)
############################################


############################################
## download HGNC file
file_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

# define link and file name
# TODO: this should be a config variable
# TODO: add check if files already exist
hgnc_link <-
  "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt"

hgnc_file <- "data/non_alt_loci_set.txt"

# download and gzip file to save space
download.file(hgnc_link, hgnc_file, mode = "wb")
gzip(hgnc_file, overwrite = TRUE)
############################################


############################################
## load the downloaded HGNC file
# TODO: specify column specifications to suppress warnings
non_alt_loci_set <- read_delim(paste0(hgnc_file, ".gz"),
    "\t",
    col_names = TRUE) %>%
  mutate(update_date = file_date)
############################################


############################################
## load STRINGdb database
string_db <- STRINGdb$new(version = "11.5",
  species = 9606,
  score_threshold = 200,
  input_directory = "data/")
############################################


############################################
## map gene symbols to StringDB identifiers
non_alt_loci_set_table <- non_alt_loci_set %>%
  dplyr::select(symbol) %>%
  unique()

non_alt_loci_set_df <- non_alt_loci_set_table %>%
    as.data.frame()

non_alt_loci_set_mapped <- string_db$map(non_alt_loci_set_df, "symbol")
non_alt_loci_set_mapped_tibble <- as_tibble(non_alt_loci_set_mapped) %>%
  filter(!is.na(STRING_id)) %>%
  group_by(symbol) %>%
  summarise(STRING_id = str_c(STRING_id, collapse = ";")) %>%
  ungroup %>%
  unique()

## join with String identifiers
non_alt_loci_set_string <- non_alt_loci_set %>%
  left_join(non_alt_loci_set_mapped_tibble, by = "symbol")
############################################


############################################
## add gene coordinates from ensembl
# TODO: fix warning "! Ensembl will soon enforce the use of https. Ensure the 'host' argument includes https://""
non_alt_loci_set_coordinates <- non_alt_loci_set_string %>%
  mutate(hg19_coordinates_from_ensembl =
    gene_coordinates_from_ensembl(ensembl_gene_id)) %>%
  mutate(hg19_coordinates_from_symbol =
    gene_coordinates_from_symbol(symbol)) %>%
  mutate(hg38_coordinates_from_ensembl =
    gene_coordinates_from_ensembl(ensembl_gene_id, reference = "hg38")) %>%
  mutate(hg38_coordinates_from_symbol =
    gene_coordinates_from_symbol(symbol, reference = "hg38")) %>%
  mutate(bed_hg19 =
    case_when(
      !is.na(hg19_coordinates_from_ensembl$bed_format) ~
        hg19_coordinates_from_ensembl$bed_format,
      is.na(hg19_coordinates_from_ensembl$bed_format) ~
        hg19_coordinates_from_symbol$bed_format,
    )
  ) %>%
  mutate(bed_hg38 =
    case_when(
      !is.na(hg38_coordinates_from_ensembl$bed_format) ~
        hg38_coordinates_from_ensembl$bed_format,
      is.na(hg38_coordinates_from_ensembl$bed_format) ~
        hg38_coordinates_from_symbol$bed_format,
    )
  ) %>%
  dplyr::select(-hg19_coordinates_from_ensembl,
    -hg19_coordinates_from_symbol,
    -hg38_coordinates_from_ensembl,
    -hg38_coordinates_from_symbol)
############################################


############################################
## annotate genomic and CDS size
# TODO: add mane transcript length (currently only cds length)
non_alt_loci_set_coordinates_length <- non_alt_loci_set_coordinates %>%
  mutate(hg19_genomic_size_from_symbol = gene_length_from_symbol(symbol, reference = "hg19")) %>%
  mutate(hg19_genomic_size_from_ensembl = gene_length_from_ensembl(ensembl_gene_id, reference = "hg19")) %>%
  mutate(hg19_cds_size_from_mane = transcript_length_from_ensembl(mane_select, reference = "hg19")) %>%
  mutate(hg38_genomic_size_from_symbol = gene_length_from_symbol(symbol, reference = "hg38")) %>%
  mutate(hg38_genomic_size_from_ensembl = gene_length_from_ensembl(ensembl_gene_id, reference = "hg38")) %>%
  mutate(hg38_cds_size_from_mane = transcript_length_from_ensembl(mane_select, reference = "hg38")) %>%
  mutate(hg19_genomic_size =
    case_when(
      !is.na(hg19_genomic_size_from_symbol$genomic_length) ~
        hg19_genomic_size_from_symbol$genomic_length,
      is.na(hg19_genomic_size_from_symbol$genomic_length) ~
        hg19_genomic_size_from_ensembl$genomic_length,
    )
  ) %>%
  mutate(hg38_genomic_size =
    case_when(
      !is.na(hg38_genomic_size_from_symbol$genomic_length) ~
        hg38_genomic_size_from_symbol$genomic_length,
      is.na(hg38_genomic_size_from_symbol$genomic_length) ~
        hg38_genomic_size_from_ensembl$genomic_length,
    )
  ) %>%
  mutate(hg19_cds_size_mane = hg19_cds_size_from_mane$transcript_length) %>%
  mutate(hg38_cds_size_mane = hg38_cds_size_from_mane$transcript_length) %>%
  dplyr::select(-hg19_genomic_size_from_symbol,
    -hg19_genomic_size_from_ensembl,
    -hg19_cds_size_from_mane,
    -hg38_genomic_size_from_symbol,
    -hg38_genomic_size_from_ensembl,
    -hg38_cds_size_from_mane)
############################################


############################################
## export table as csv with date of creation
creation_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

write_csv(non_alt_loci_set_coordinates_length,
  file = paste0("results/non_alt_loci_set_coordinates_length.", creation_date, ".csv"))

gzip(paste0("results/non_alt_loci_set_coordinates_length.", creation_date, ".csv"),
  overwrite = TRUE)
############################################