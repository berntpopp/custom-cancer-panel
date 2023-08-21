#### This file holds analyses functions for the ensembl database name using biomart

#' Retrieve gene coordinates in BED format from gene symbols
#'
#' This function retrieves the gene coordinates in BED format for the given gene
#' symbols. The coordinates are obtained from the specified reference genome.
#'
#' @param gene_symbols A vector or tibble containing the gene symbols.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the gene symbols and their corresponding coordinates in BED format.
#'
#' @examples
#' gene_symbols <- c("ARID1B ", "GRIN2B", "NAA10")
#' gene_coordinates_from_symbol(gene_symbols, reference = "hg19")
#'
#' @export
gene_coordinates_from_symbol <- function(gene_symbols, reference = "hg19") {
  gene_symbol_list <- as_tibble(gene_symbols) %>%
    dplyr::select(hgnc_symbol = value)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position")
  filters <- c("hgnc_symbol")

  values <- list(hgnc_symbol = gene_symbol_list$hgnc_symbol)

  gene_coordinates_hg19 <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
    group_by(hgnc_symbol) %>%
    summarise(hgnc_symbol = max(hgnc_symbol), chromosome_name = max(chromosome_name), start_position = max(start_position), end_position = max(end_position)) %>%
    mutate(bed_format = paste0("chr", chromosome_name, ":", start_position, "-", end_position)) %>%
    dplyr::select(hgnc_symbol, bed_format)

  gene_symbol_list_return <- gene_symbol_list %>%
  left_join(gene_coordinates_hg19, by = ("hgnc_symbol"))

  return(gene_symbol_list_return)
}


#' Retrieve gene coordinates in BED format from Ensembl IDs
#'
#' This function retrieves the gene coordinates in BED format for the given Ensembl
#' gene IDs. The coordinates are obtained from the specified reference genome.
#'
#' @param ensembl_id A vector or tibble containing the Ensembl gene IDs.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the Ensembl gene IDs and their corresponding coordinates in BED format.
#'
#' @examples
#' ensembl_id <- c("ENSG00000049618", "ENSG00000273079", "ENSG00000102030")
#' gene_coordinates_from_ensembl(ensembl_id, reference = "hg19")
#'
#' @export
gene_coordinates_from_ensembl <- function(ensembl_id, reference = "hg19") {
  ensembl_id_list <- as_tibble(ensembl_id) %>%
    dplyr::select(ensembl_gene_id = value)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
  filters <- c("ensembl_gene_id")

  values <- list(ensembl_gene_id = ensembl_id_list$ensembl_gene_id)

  gene_coordinates_hg19 <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
    group_by(ensembl_gene_id) %>%
    summarise(ensembl_gene_id = max(ensembl_gene_id), chromosome_name = max(chromosome_name), start_position = max(start_position), end_position = max(end_position)) %>%
    mutate(bed_format = paste0("chr", chromosome_name, ":", start_position, "-", end_position)) %>%
    dplyr::select(ensembl_gene_id, bed_format)

  ensembl_id_list_return <- ensembl_id_list %>%
  left_join(gene_coordinates_hg19, by = ("ensembl_gene_id"))

  return(ensembl_id_list_return)
}


#' Retrieve genomic length from gene symbols
#'
#' This function retrieves the genomic length for the given gene
#' symbols. The length is obtained by subtracting start_position from end_position.
#'
#' @param gene_symbols A vector or tibble containing the gene symbols.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the gene symbols and their corresponding genomic lengths.
#'
#' @examples
#' gene_symbols <- c("ARID1B ", "GRIn2B", "NAA10")
#' gene_length_from_symbol(gene_symbols, reference = "hg19")
#'
#' @export
gene_length_from_symbol <- function(gene_symbols, reference = "hg19") {
  gene_symbol_list <- as_tibble(gene_symbols) %>%
    dplyr::select(hgnc_symbol = value)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("hgnc_symbol", "start_position", "end_position")
  filters <- c("hgnc_symbol")

  values <- list(hgnc_symbol = gene_symbol_list$hgnc_symbol)

  gene_coordinates <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
    group_by(hgnc_symbol) %>%
    summarise(hgnc_symbol = max(hgnc_symbol), start_position = max(start_position), end_position = max(end_position)) %>%
    mutate(genomic_length = end_position - start_position) %>%
    dplyr::select(hgnc_symbol, genomic_length)

  gene_symbol_list_return <- gene_symbol_list %>%
  left_join(gene_coordinates, by = ("hgnc_symbol"))

  return(gene_symbol_list_return)
}


#' Retrieve genomic length from Ensembl IDs
#'
#' This function retrieves the genomic length for the given Ensembl
#' gene IDs. The length is obtained by subtracting start_position from end_position.
#'
#' @param ensembl_id A vector or tibble containing the Ensembl gene IDs.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the Ensembl gene IDs and their corresponding genomic lengths.
#'
#' @examples
#' ensembl_id <- c("ENSG00000049618", "ENSG00000273079", "ENSG00000102030")
#' gene_length_from_ensembl(ensembl_id, reference = "hg19")
#'
#' @export
gene_length_from_ensembl <- function(ensembl_id, reference = "hg19") {
  ensembl_id_list <- as_tibble(ensembl_id) %>%
    dplyr::select(ensembl_gene_id = value)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("ensembl_gene_id", "start_position", "end_position")
  filters <- c("ensembl_gene_id")

  values <- list(ensembl_gene_id = ensembl_id_list$ensembl_gene_id)

  gene_coordinates <- getBM(attributes=attributes, filters=filters, values=values, mart=mart) %>%
    group_by(ensembl_gene_id) %>%
    summarise(ensembl_gene_id = max(ensembl_gene_id), start_position = max(start_position), end_position = max(end_position)) %>%
    mutate(genomic_length = end_position - start_position) %>%
    dplyr::select(ensembl_gene_id, genomic_length)

  ensembl_id_list_return <- ensembl_id_list %>%
  left_join(gene_coordinates, by = ("ensembl_gene_id"))

  return(ensembl_id_list_return)
}


#' Retrieve transcript length from Ensembl transcript IDs
#'
#' This function retrieves the transcript length for the given Ensembl
#' transcript IDs. The length is obtained from the specified reference genome.
#'
#' @param transcript_ids A vector or tibble containing the Ensembl transcript IDs.
#' @param reference The reference genome to use (default: "hg19").
#'
#' @return A tibble with the Ensembl transcript IDs and their corresponding lengths.
#'
#' @examples
#' transcript_ids <- c("ENST00000636930.2|NM_001374828.1", "ENST00000609686.4|NM_000834.5", "ENST00000464845.6|NM_003491.4")
#' transcript_length_from_ensembl(transcript_ids, reference = "hg19")
#'
#' @export
transcript_length_from_ensembl <- function(transcript_ids, reference = "hg19") {
  transcript_id_list <- as_tibble(transcript_ids) %>%
    dplyr::mutate(transcript_id = str_extract(value, "ENST\\d+")) %>%
    dplyr::select(ensembl_transcript_id = transcript_id)

  # define mart
  mart_hg19 <- useMart("ensembl", host="grch37.ensembl.org")
  mart_hg19 <- useDataset("hsapiens_gene_ensembl", mart_hg19)

  mart_hg38 <- useMart("ensembl", host="ensembl.org")
  mart_hg38 <- useDataset("hsapiens_gene_ensembl", mart_hg38)

  if (reference == "hg19") {
    mart <- useMart("ensembl", host = "grch37.ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg19)
  } else {
    mart <- useMart("ensembl", host = "ensembl.org")
    mart <- useDataset("hsapiens_gene_ensembl", mart_hg38)
  }

  attributes <- c("ensembl_transcript_id", "transcript_length")
  filters <- c("ensembl_transcript_id")

  values <- list(ensembl_transcript_id = transcript_id_list$ensembl_transcript_id)

  transcript_lengths <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)

  transcript_id_list_return <- transcript_id_list %>%
  left_join(transcript_lengths, by = ("ensembl_transcript_id"))

  return(transcript_id_list_return)
}


#' Retrieve SNP positions from rsIDs using Ensembl
#'
#' This function retrieves the position and chromosome information for the given SNP rsIDs.
#' The information is obtained from the specified reference genome via the Ensembl API.
#'
#' @param rs_ids A vector or tibble containing the SNP rsIDs.
#' @param reference The reference genome to use. Valid values are "hg38" and "hg19" (default: "hg38").
#'
#' @return A tibble with the SNP rsIDs and their corresponding chromosome, chrom_start, and chrom_end.
#'
#' @examples
#' rs_ids <- c("rs1234567", "rs9876543", "rs11223344")
#' snp_position_from_rs(rs_ids, reference = "hg38")
#'
#' @export
snp_position_from_rs <- function(rs_ids, reference = "hg38") {
  # transform rs_ids to tibble
  rs_ids_list <- as_tibble(rs_ids) %>%
    dplyr::select(snp_id = value)

  # define mart
  if (reference == "hg19") {
    mart <- useEnsembl("snp", dataset = "hsapiens_snp", GRCh = "37")
  } else {
    mart <- useEnsembl("snp", dataset = "hsapiens_snp")
  }

  # get snp position and chromosome
  snps <- getBM(attributes=c("refsnp_id",
                    "chr_name",
                    "chrom_start",
                    "chrom_end"),
      filters = "snp_filter", values = rs_ids_list$snp_id, mart = mart, uniqueRows = TRUE) %>%
      tibble() %>%
      filter(!str_detect(chr_name, pattern = "^H"))

  # join back to rs_ids_list
  snps_return <- left_join(rs_ids_list, snps, by = c("snp_id" = "refsnp_id")) %>%
    dplyr::select(chr = chr_name, chrom_start, chrom_end)

  return(snps_return)
}