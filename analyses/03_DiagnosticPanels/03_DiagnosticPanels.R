############################################
## load libraries
library(readr)
library(readxl)
library(tidyverse)
library(rvest)
library(jsonlite)
library(curl)
library(httr)
library(config) # needed for config loading
library(webdriver) # needed for headless browsing
library("R.utils")  ## gzip downloaded and result files
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/03_DiagnosticPanels/"

## read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = "default")
config_vars_path <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_path$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
source("../functions/phantomjs-functions.R", local = TRUE)
source("../functions/blueprintgenetics-functions.R", local = TRUE)
############################################


############################################
## download web urls

# load the list of sources
diagnostic_panels_list <- read_excel("data/cancer_diagnostic_panels_list.xlsx") %>%
  filter(use == "yes")

# download using phantomJS
diagnostic_panels <- diagnostic_panels_list %>%
  rowwise() %>%
  mutate(filename_download = download_url_by_phantomjs(diagnostic_panel_source,
    diagnostic_panel_name, type)) %>%
  ungroup()
############################################


############################################
## compute standardized gene list per panel

## 1) blueprintgenetics_hereditary_cancer
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "blueprintgenetics_hereditary_cancer"))$diagnostic_panel_source

# generate list of blueprintgenetics sub panels
blueprintgenetics_panel_list <- diagnostic_panels %>%
  filter(diagnostic_panel_name == "blueprintgenetics_hereditary_cancer") %>%
  select(diagnostic_panel_name, subpanel_name, subpanel_source) %>%
  separate_rows(c(subpanel_name, subpanel_source), sep = ", ", convert = TRUE) %>%
  mutate(subpanel_name = paste0(diagnostic_panel_name, "_", subpanel_name)) %>%
  mutate(type = "html")

# download using phantomJS
blueprintgenetics_panel_panels <- blueprintgenetics_panel_list %>%
  rowwise() %>%
  mutate(filename_download = download_url_by_phantomjs(subpanel_source, subpanel_name, type)) %>%
  ungroup()

# read all blueprintgenetics nephrology sub panels
blueprintgenetics_panel_panels_genes <- blueprintgenetics_panel_panels %>%
  rowwise() %>%
  mutate(gene_list = list(blueprintgenetics_panel_extract_genes(filename_download))) %>%
  ungroup()

blueprintgenetics_hereditary_cancer_genes <- blueprintgenetics_panel_panels_genes$gene_list %>%
  unlist() %>%
  tibble(`Genes` = .) %>%
  mutate(Panel = "blueprintgenetics_hereditary_cancer") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#\\.]")) %>%
  unique() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 2) centogene_solid_tumor_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "centogene_solid_tumor_panel"))$filename_download

centogene_solid_tumor_panel <- read_html(url)

centogene_solid_tumor_panel_genes_nodes <- centogene_solid_tumor_panel %>%
  html_nodes(xpath = '//*[@id="t3m-Modal--74420"]//table[1]') %>%
  html_table()

centogene_solid_tumor_panel_genes <- centogene_solid_tumor_panel_genes_nodes[[1]] %>%
  select(Genes = Gene) %>%
  unique() %>%
  mutate(Panel = "centogene_solid_tumor_panel") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 3) fulgentgenetics_comprehensivecancer_full
url <- (diagnostic_panels %>%
      filter(diagnostic_panel_name == "fulgentgenetics_comprehensivecancer_full"))$filename_download

fulgentgenetics_comprehensivecancer_full <- read_html(url)

fulgentgenetics_comprehensivecancer_full_genes <- fulgentgenetics_comprehensivecancer_full %>%
    html_nodes(xpath = '//meta[contains(@content,"AIP")]') %>%
    html_attr("content") %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, convert = TRUE) %>%
  mutate(Panel = "fulgentgenetics_comprehensivecancer_full") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
  filter(!(Genes %in% c("Full", "Comprehensive", "Cancer", "Panel", "COMPREHENSIVE"))) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


#TODO: fix this as it is not working
## 4) cegat_tumor_syndromes_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "cegat_tumor_syndromes_panel"))$filename_download

cegat_tumor_syndromes_panel <- read_html(url)

cegat_tumor_syndromes_panel_genes <- cegat_tumor_syndromes_panel %>%
    html_nodes(xpath = '//h2[contains(text(),"Gene Directory")]//following::em') %>%
  html_text() %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, convert = TRUE) %>%
  mutate(Panel = "cegat_tumor_syndromes_panel") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 5) invitae_tumor_syndromes_panel
url <- (diagnostic_panels %>%
      filter(diagnostic_panel_name == "invitae_tumor_syndromes_panel"))$filename_download

invitae_tumor_syndromes_panel <- read_html(url)

invitae_tumor_syndromes_panel_genes <- invitae_tumor_syndromes_panel %>%
    html_nodes(xpath = '//meta[contains(@content,"AIP")]') %>%
    html_attr("content") %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, convert = TRUE) %>%
  mutate(Panel = "invitae_tumor_syndromes_panel") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 6) mgz_erbliche_tumorerkrankungen_umfassende_diagnostik
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "mgz_erbliche_tumorerkrankungen_umfassende_diagnostik"))$filename_download

mgz_erbliche_tumorerkrankungen_umfassende_diagnostik <- read_html(url)

mgz_erbliche_tumorerkrankungen_umfassende_diagnostik_genes <- mgz_erbliche_tumorerkrankungen_umfassende_diagnostik %>%
    html_nodes(xpath = '//div[contains(@class,"panel_gene")]') %>%
  html_text() %>%
  str_remove_all(., "[\\n\\r ]+") %>%
  tibble(`Genes` = .) %>%
  mutate(Panel = "mgz_erbliche_tumorerkrankungen_umfassende_diagnostik") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
  filter(!str_detect(Genes, "\\(")) %>%
  filter(Genes != "") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel) %>%
  unique()


## 7) uchicago_comprehensive_hereditary_cancer_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "uchicago_comprehensive_hereditary_cancer_panel"))$filename_download

uchicago_comprehensive_hereditary_cancer_panel <- read_html(url)

uchicago_comprehensive_hereditary_cancer_panel_genes <- uchicago_comprehensive_hereditary_cancer_panel %>%
    html_nodes(xpath = '//a[contains(@href,"/gene/")]') %>%
  html_text() %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, convert = TRUE) %>%
  mutate(Panel = "uchicago_comprehensive_hereditary_cancer_panel") %>%
  mutate(Source = url) %>%
  filter(Genes != "") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 8) preventiongenetics_cancer_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "preventiongenetics_cancer_panel"))$filename_download

preventiongenetics_cancer_panel <- read_html(url)

preventiongenetics_cancer_panel_genes <- preventiongenetics_cancer_panel %>%
  html_nodes(xpath = '//*[@id="genes-div"]//table[1]') %>%
  html_table()

preventiongenetics_cancer_panel_genes <- preventiongenetics_cancer_panel_genes[[1]] %>%
  select(Genes = `Official Gene Symbol`) %>%
  unique() %>%
  mutate(Panel = "preventiongenetics_cancer_panel") %>%
  mutate(Source = url) %>%
  mutate(Genes = str_remove_all(Genes, "[\\*\\#]")) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


#TODO: fix this as it is not working
## 9) mayocliniclabs_xcp_hereditary_expanded_cancer_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "mayocliniclabs_xcp_hereditary_expanded_cancer_panel"))$filename_download

mayocliniclabs_xcp_hereditary_expanded_cancer_panel <- read_html(url)

mayocliniclabs_xcp_hereditary_expanded_cancer_panel_genes <- mayocliniclabs_xcp_hereditary_expanded_cancer_panel %>%
    html_nodes(xpath = '//span[contains(text(),"Genes analyzed:")]') %>%
  html_text() %>%
  str_remove_all("Genes analyzed: ") %>%
  str_remove_all(" \\(including promoters 1A \\& 1B\\)") %>%
  str_remove_all(" \\(upstream enhancer region duplication only\\)") %>%
  str_remove_all(" \\(c\\.952G>A p\\.E318K variant only\\)") %>%
  str_remove_all(" \\(including promoter\\)") %>%
  str_remove_all(" \\(copy number variants only\\)") %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, convert = TRUE) %>%
  mutate(Panel = "mayocliniclabs_xcp_hereditary_expanded_cancer_panel") %>%
  mutate(Source = url) %>%
  filter(Genes != "and") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 10) myriadgenetics_myrisk_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "myriadgenetics_myrisk_panel"))$filename_download

myriadgenetics_myrisk_panel <- read_html(url)

myriadgenetics_myrisk_panel_genes <- myriadgenetics_myrisk_panel %>%
    html_nodes(xpath = '//td[contains(@class,"gene")]') %>%
  html_text() %>%
  str_remove_all(., "[\\n\\r ]+") %>%
  str_remove_all("Monoallelic") %>%
  str_remove_all("Biallelic") %>%
  str_remove_all("\\(p16INK4a\\)") %>%
  str_remove_all("\\(p14ARF\\)") %>%
  tibble(`Genes` = .) %>%
  mutate(Panel = "myriadgenetics_myrisk_panel") %>%
  mutate(Source = url) %>%
  unique() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel) %>%
  unique()


## 11) genedx_oncogenedx_custom_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "genedx_oncogenedx_custom_panel"))$filename_download

genedx_oncogenedx_custom_panel <- read_html(url)

genedx_oncogenedx_custom_panel_genes <- genedx_oncogenedx_custom_panel %>%
    html_nodes(xpath = '//span[contains(@class,"gene-name")]') %>%
  html_text() %>%
  str_remove_all(" \\(GREM1\\)") %>%
  tibble(`Genes` = .) %>%
  mutate(Panel = "genedx_oncogenedx_custom_panel") %>%
  mutate(Source = url) %>%
  unique() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel) %>%
  unique()


## 12) arupconsult_hereditary_cancer_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "arupconsult_hereditary_cancer_panel"))$filename_download

arupconsult_hereditary_cancer_panel <- read_html(url)

arupconsult_hereditary_cancer_panel_genes <- arupconsult_hereditary_cancer_panel %>%
  html_nodes(xpath = '//*[@id="genes-tested1"]//following::table') %>%
  html_table()

arupconsult_hereditary_cancer_panel_genes <- arupconsult_hereditary_cancer_panel_genes[[1]] %>%
  select(Genes = Gene) %>%
  unique() %>%
  mutate(Panel = "arupconsult_hereditary_cancer_panel") %>%
  mutate(Source = url) %>%
  filter(!str_detect(Genes, "Association")) %>%
  mutate(Genes = str_remove_all(Genes, "\\(Exon 9 deletion\\/duplications only\\)")) %>%
  mutate(Genes = str_trim(Genes)) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 13) cancergeneticslab_hereditary_cancer_panel
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "cancergeneticslab_hereditary_cancer_panel"))$filename_download

cancergeneticslab_hereditary_cancer_panel <- read_html(url)

cancergeneticslab_hereditary_cancer_panel_genes <- cancergeneticslab_hereditary_cancer_panel %>%
    html_nodes(xpath = '//h2[contains(text(),"GENES TARGETED")]//following::p') %>%
  html_text() %>%
  tibble(`Genes` = .) %>%
  filter(!str_detect(Genes, "Single nucleotide variants")) %>%
  filter(!str_detect(Genes, "Variants not presumed by nature")) %>%
  filter(!str_detect(Genes, "Genomic DNA is subjected to FFPE repair")) %>%
  filter(!str_detect(Genes, "Variants are validated")) %>%
  mutate(Genes = str_remove_all(Genes, "Entire coding region\\:  ")) %>%
  mutate(Genes = str_remove_all(Genes, "Partial Genes\\: ")) %>%
  mutate(Genes = str_remove_all(Genes, "Copy Number Variants \\(not applicable to FFPE specimens\\)\\: Copy number variants are called in all of the above genes as well as:")) %>%
  separate_rows(., Genes, sep = "; ", convert = TRUE) %>%
  mutate (Genes = str_remove_all(Genes, "_.+")) %>%
  filter(Genes != "") %>%
  mutate(Genes = str_trim(Genes)) %>%
  mutate(Panel = "cancergeneticslab_hereditary_cancer_panel") %>%
  mutate(Source = url) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 14) neogenomics_full_comprehensive_cancer_panel_germline
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "neogenomics_full_comprehensive_cancer_panel_germline"))$filename_download

neogenomics_full_comprehensive_cancer_panel_germline <- read_html(url)

neogenomics_full_comprehensive_cancer_panel_germline_genes <- neogenomics_full_comprehensive_cancer_panel_germline %>%
    html_nodes(xpath = '//p[contains(text(),"AIP")]') %>%
  html_text() %>%
  str_remove_all("Gene list: ") %>%
  str_remove_all("\\( 127 genes \\)") %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, sep = ",", convert = TRUE) %>%
  mutate(Genes = str_trim(Genes)) %>%
  mutate(Panel = "neogenomics_full_comprehensive_cancer_panel_germline") %>%
  mutate(Source = url) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)


## 15) natera_hereditary_cancer_test_comprehensive
url <- (diagnostic_panels %>%
  filter(diagnostic_panel_name == "natera_hereditary_cancer_test_comprehensive"))$filename_download

natera_hereditary_cancer_test_comprehensive <- read_html(url)

natera_hereditary_cancer_test_comprehensive_genes <- natera_hereditary_cancer_test_comprehensive %>%
    html_nodes(xpath = '//strong[contains(text(),"Genes:")]/following-sibling::span|//strong[contains(text(),"Genes:")]/following-sibling::text()') %>%
  html_text() %>%
  tibble(`Genes` = .) %>%
  separate_rows(., Genes, sep = ", ", convert = TRUE) %>%
  mutate(Genes = str_trim(Genes)) %>%
  mutate(Genes = str_remove_all(Genes, ",")) %>%
  unique() %>%
  mutate(Panel = "natera_hereditary_cancer_test_comprehensive") %>%
  mutate(Source = url) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = Genes, source = Source, panel = Panel)

############################################


############################################
## bind all tables and summarize
## compute count of panels a gene is reported in
## compute diagnostic panel evidence as genes reported in at least 2 panel sources
all_diagnostic_panels_genes <- bind_rows(blueprintgenetics_hereditary_cancer_genes,
  centogene_solid_tumor_panel_genes,
  fulgentgenetics_comprehensivecancer_full_genes,
  cegat_tumor_syndromes_panel_genes,
  invitae_tumor_syndromes_panel_genes,
  mgz_erbliche_tumorerkrankungen_umfassende_diagnostik_genes,
  uchicago_comprehensive_hereditary_cancer_panel_genes,
  preventiongenetics_cancer_panel_genes,
  mayocliniclabs_xcp_hereditary_expanded_cancer_panel_genes,
  myriadgenetics_myrisk_panel_genes,
  genedx_oncogenedx_custom_panel_genes,
  arupconsult_hereditary_cancer_panel_genes,
  cancergeneticslab_hereditary_cancer_panel_genes,
  neogenomics_full_comprehensive_cancer_panel_germline_genes,
  natera_hereditary_cancer_test_comprehensive_genes) %>%
  group_by(approved_symbol) %>%
  summarise(panel_diagnostic = paste(unique(panel), collapse = "; "),
    panel_diagnostic_count = n(),
    gene_name_reported = paste(unique(gene_name_reported), collapse = " | ")) %>%
  ungroup() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(approved_symbol)) %>%
  mutate(at_least_two_panels = (panel_diagnostic_count > 1))

all_diagnostic_panels_genes_format <- all_diagnostic_panels_genes %>%
  select(approved_symbol, hgnc_id, gene_name_reported, source = panel_diagnostic, source_count = panel_diagnostic_count, source_evidence = at_least_two_panels)

############################################


############################################
## save results
# TODO: gzip csv result files
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(all_diagnostic_panels_genes_format,
  file = paste0("results/03_DiagnosticPanels_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

write_csv(diagnostic_panels,
  file = paste0("results/03_DiagnosticPanels_list.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/03_DiagnosticPanels_list.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
