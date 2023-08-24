############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library("R.utils")  ## gzip downloaded and result files
library(config)     ## needed for config loading
library(biomaRt)    ## needed to get gene coordinates
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/V04_DeepIntronicClinVar/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################


############################################
# load global functions
# hgnc functions
source("../functions/file-functions.R", local = TRUE)
source("../functions/ensembl-functions.R", local = TRUE)
source("../functions/gnomad-functions.R", local = TRUE)
############################################


############################################
## download all required database sources from clinvar

current_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

# clinvar VCF file download
if (check_file_age("clinvar", "data/downloads/", 1)) {
  clinvar_vcf_hg19_filename <- get_newest_file("clinvar", "data/downloads/")
} else {
  # clinvar VCF file links to genemap2 file needs to be set in config
  clinvar_vcf_hg19_url <- config_vars_proj$clinvar_vcf_url

  clinvar_vcf_hg19_filename <- paste0("data/downloads/clinvar.",
    current_date,
    ".vcf.gz")

  download.file(clinvar_vcf_hg19_url, clinvar_vcf_hg19_filename, mode = "wb")
}
############################################


############################################
# load or process clinvar VCF file into table
# and then safe/load the result as a table for future use
if (check_file_age("clinvar_table", "results/", 1)) {
  clinvar_table_filename <- get_newest_file("clinvar_table", "results/")

  clinvar_table <- read_csv(clinvar_table_filename)
} else {
  ## load the downloaded clinvar VCF file
  clinvar_vcf <- read_delim(clinvar_vcf_hg19_filename,
      delim = "\t", escape_double = FALSE,
      trim_ws = TRUE,
      comment = "##")

  # reformat the vcf into a table
  clinvar_table <- clinvar_vcf %>% 
    separate_rows(INFO, sep = ";") %>%
    separate(INFO, into = c("key", "value"), sep = "=") %>%
    spread(key, value) %>%
    separate_rows(GENEINFO, sep = "\\|") %>%
    separate_rows(CLNSIG, sep = "/")

  # write the table to a csv file
  write_csv(clinvar_table,
    file = paste0("results/clinvar_table.", current_date, ".csv"))

  gzip(paste0("results/clinvar_table.", current_date, ".csv"),
    overwrite = TRUE)
}
############################################


############################################
# load the merged gene table and extract gene names
# compute exon coordinates using ensembl functions
germline_merged_filename <- get_newest_file("CancerPanel_Germline_MergeAnalysesSources", "../D_MergeAnalysesSources/results/")

germline_merged <- read_csv(germline_merged_filename)

germline_merged_exons <- exon_coordinates_from_symbol(germline_merged$approved_symbol, padding = 20) %>%
  dplyr::select(symbol = external_gene_name, chromosome = chromosome_name, start = padded_start, end = padded_end)
############################################


############################################
# extract variants from clinvar table found in germline_merged
# filter for pathogenic and likely pathogenic variants
clinvar_table_variants <- clinvar_table %>%
  separate(GENEINFO, into = c("symbol", "Gene_ID"), sep = ":") %>%
  filter(CLNSIG %in% c("Pathogenic", "Likely_pathogenic")) %>%
  filter(symbol %in% germline_merged$approved_symbol) %>%
  dplyr::select(symbol, chromosome = `#CHROM`, start = POS, end = POS, ref = REF, alt = ALT, clinsig = CLNSIG, MC) %>%
  unique()
############################################




clinvar_variants_outside_target <- genome_anti_join(clinvar_table_variants, germline_merged_exons, by = c("chromosome", "start", "end"))











############################################
## export table as csv with date of creation
creation_date <- strftime(as.POSIXlt(Sys.time(),
    "UTC", "%Y-%m-%dT%H:%M:%S"),
  "%Y-%m-%d")

write_csv(non_alt_loci_set_coordinates_gnomad_omim_gencc_clinvar,
  file = paste0("results/non_alt_loci_set_coordinates.", creation_date, ".csv"))

gzip(paste0("results/non_alt_loci_set_coordinates.", creation_date, ".csv"),
  overwrite = TRUE)
############################################