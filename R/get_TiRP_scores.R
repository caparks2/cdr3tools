#' Calculate TiRP Scores
#'
#' A scoring system to quantify TCR-intrinsic regulatory potential (TiRP)
#'
#' A special note should be made to consider that this function, as written by
#'   the authors, evaluates only sequences that are between 12 and 17 amino
#'   acids in length and that evaluation is performed on sequences as they are,
#'   in their linear, native states without first being adjusted for length by
#'   IMGT unique numbering rules for JUNCTIONS. It was also discovered that the
#'   original version of this function as written by the authors would fail to
#'   produce "vgene scores", resulting in incorrect overall scores, if allele
#'   identifiers were part of the V gene names. This was corrected in this
#'   implementation of the calculation.
#'
#' @param .data A data frame (and/or a tibble). This is the input TCR sequence
#'   data that will be used to calculate TiRP scores. The format must conform to
#'   the following:
#'   \describe{
#'   \item{First Column}{A character vector of TRBV gene names}
#'   \item{Second Column}{A character vector of IMGT JUNCTION sequences (CDR3
#'     amino acid sequences from 2nd-CYS p104 to J-TRP/J-PHE p118)}
#'   }
#'   Column names are not important.
#'
#'   `.data` can alternativeley be a list containing multiple of such data
#'     frames.
#' @param .details A logical. `FALSE` (the default) returns results as numeric
#'   vectors of TiRP scores. `TRUE` returns results as data frames containing
#'   all evaluated V genes and their CDR3 sequences, along with all
#'   intermediary scores plus the final TiRP scores.
#' @returns A numeric vector of TiRP scores or a data frame of detailed results
#'   for each evaluated unique sequence.
#' @examples
#' # Confirm the proper format for the input data
#' head(tcr_seqs)
#'
#' get_TiRP_scores(head(tcr_seqs))
#' @author Kaitlyn A. Lagattuta, Joyce B. Kang, Aparna Nathan (Center for Data
#'   Sciences, Brigham and Women’s Hospital, Boston, MA, USA)
#'   Christopher Parks (rewritten and improved for use in this package.)
#' @references
#' Lagattuta KA, Kang JB, Nathan A, Pauken KE, Jonsson AH, Rao DA, Sharpe AH,
#'   Ishigaki K, Raychaudhuri S. Repertoire analyses reveal T cell antigen
#'   receptor sequence features that influence T cell fate. Nat Immunol. 2022
#'   Mar;23(3):446-457. doi: 10.1038/s41590-022-01129-x. Epub 2022 Feb 17. PMID:
#'   35177831; PMCID: PMC8904286.
#'
#' https://www.nature.com/articles/s41590-022-01129-x
#'
#' https://github.com/immunogenomics/TiRP
#' @export get_TiRP_scores
get_TiRP_scores <- function(.data, .details = FALSE) {

  if (inherits(.data, "list")) {
    lapply(.data, get_TiRP_scores_internal, .details = .details)
  } else if (inherits(.data, "data.frame")) {
    get_TiRP_scores_internal(.data, .details = .details)
  } else {
    rlang::abort(
      paste(
        ".data must either be a single data frame (with V genes in the first",
        "column and CDR3 amino acid sequences in the second column) or a list",
        "containing multiple of such data frames."
      )
    )
  }
}

get_TiRP_scores_internal <- function(.data, .details) {
  data <- as.data.frame(.data)

  # v_gene_col <- data[, apply(data, 2, function(x) {
  #   x <- x[!is.na(x)]
  #   sum(grepl("^TRBV[0-9\\*-]+$|^TCRBV[0-9\\*-]+$", x, ignore.case = TRUE)) > 0
  # }), drop = TRUE]
  #
  # cdr3_col <- data[, apply(data, 2, function(x) {
  #   x <- x[!is.na(x)]
  #   not_v_genes <- !grepl("^TRBV[0-9\\*-]+$|^TCRBV[0-9\\*-]+$", x, ignore.case = TRUE)
  #   amino_acids <- grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", x, ignore.case = TRUE)
  #   not_nucleotides <- grepl("[^ACTG]", x, ignore.case = TRUE)
  #   len <- c(length(not_v_genes) > 0, length(amino_acids) > 0, length(not_nucleotides) > 0)
  #   if (!all(len)) return(FALSE)
  #   sum(c(not_v_genes, amino_acids, not_nucleotides)) == length(c(not_v_genes, amino_acids, not_nucleotides))
  # }), drop = TRUE]
  #
  # data <- data.frame(v_gene = v_gene_col, cdr3 = cdr3_col)
  # data <- data[!is.na(data$v_gene), ]
  # data <- data[!is.na(data$cdr3), ]

  weights <- cdr3tools::TiRP_weights

  # suppressPackageStartupMessages({
  #   library(stringr)
  #   library(dplyr)
  # })

  ######################## defining TCR features ###########################

  ## standardizing TRBV gene names

  # ex <- data[, 1][!(grepl("nresolved", data[, 1]))][1]
  # # if (substr(data[1,1], 1, 4) =="TRBV"){
  # if (!(grepl("-0", ex))) {
  #   data$vgene <- reformat_vgene_cp_modified(data[, 1])
  # } else {
  #   data$vgene <- as.character(data[, 1])
  # }

  # if (any(stringr::str_detect(data[, 1], "TCRBV|\\*[0-9][0-9]"))) {
    data$vgene <- fix_vdj_genes2(data[, 1])
  # } else {
  #   data$vgene <- as.character(data[, 1])
  # }


  ## CDR3 sequence
  data$cdr3 <- as.character(data[, 2])
  data <- data[!(is.na(data$cdr3)), ]

  ## CDR3 length
  data$length <- sapply(data$cdr3, function(x) nchar(x))
  data <- data[data$length >= 12 & data$length <= 17, ]

  ## Jmotif
  jmotifs <- gsub("Jmotif", "", weights$feat[grepl("Jmotif", weights$feat)])
  data$Jmotif <- sapply(data$cdr3, function(x) substr(x, nchar(x) - 4, nchar(x)))
  data$Jmotif[!(data$Jmotif %in% jmotifs)] <- "other"

  ## amino acid composition in the CDR3 middle region
  data$cdr3MR <- sapply(data$cdr3, function(x) substr(x, 5, nchar(x) - 6))
  aminos <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  ref_stats <- cdr3tools::heldout_means_sds
  for (i in 1:length(aminos)) {
    new <- sapply(data$cdr3MR, function(x) stringr::str_count(x, aminos[i]) / nchar(x))
    mn <- ref_stats$mean[ref_stats$amino == aminos[i]]
    sd <- ref_stats$sd[ref_stats$amino == aminos[i]]
    new <- (new - mn) / sd
    data <- cbind(data, new)
    colnames(data)[ncol(data)] <- paste("perc_mid", aminos[i], sep = "_")
  }


  ######################## scoring TCR features ###########################

  ## including only terms that were Bonferroni-significant in meta-analysis:
  touse <- weights[weights$metap < (0.05 / nrow(weights)), ]

  perc_terms <- touse$feat[grepl("perc_mid", touse$feat)]
  data$perc_score <- sapply(1:nrow(data), function(x) sum(data[x, perc_terms] * touse[perc_terms, c("metabeta")]))

  pos_terms <- list()
  pos_terms[[1]] <- c("12_5", "12_6")
  pos_terms[[2]] <- c("13_5", "13_6", "13_7")
  pos_terms[[3]] <- c("14_5", "14_6", "14_7", "14_-7")
  pos_terms[[4]] <- c("15_5", "15_6", "15_7", "15_8", "15_-7")
  pos_terms[[5]] <- c("16_5", "16_6", "16_7", "16_8", "16_-8", "16_-7")
  pos_terms[[6]] <- c("17_5", "17_6", "17_7", "17_8", "17_9", "17_-8", "17_-7")
  data$pos_score <- sapply(data$cdr3MR, function(x) get_pos_score(x, touse, pos_terms))

  data$feat <- paste("vgene", data$vgene, sep = "")
  data <- dplyr::left_join(data, touse[, c("feat", "metabeta")], by = "feat")
  colnames(data)[ncol(data)] <- "vgene_score"
  data$vgene_score <- sapply(data$vgene_score, function(x) ifelse(is.na(x), 0, x))
  # if (nrow(data[data$vgene_score != 0, ]) == 0) {
  #   stop("Unidentifiable TRBV genes.")
  # }

  data$feat <- paste("p107", substr(data$cdr3, 4, 4), sep = "")
  data <- dplyr::left_join(data, touse[, c("feat", "metabeta")], by = "feat")
  colnames(data)[ncol(data)] <- "p107_score"
  data$p107_score <- sapply(data$p107_score, function(x) ifelse(is.na(x), 0, x))

  data$feat <- paste("Jmotif", data$Jmotif, sep = "")
  data <- dplyr::left_join(data, touse[, c("feat", "metabeta")], by = "feat")
  colnames(data)[ncol(data)] <- "Jmotif_score"
  data$Jmotif_score <- sapply(data$Jmotif_score, function(x) ifelse(is.na(x), 0, x))

  data$feat <- sapply(data$cdr3, function(x) substr(x, nchar(x) - 5, nchar(x) - 5))
  data$feat <- paste("p113", data$feat, sep = "")
  data <- dplyr::left_join(data, touse[, c("feat", "metabeta")], by = "feat")
  colnames(data)[ncol(data)] <- "p113_score"
  data$"p113_score" <- sapply(data$"p113_score", function(x) ifelse(is.na(x), 0, x))

  data$feat <- paste("length", data$length, sep = "")
  data <- dplyr::left_join(data, touse[, c("feat", "metabeta")], by = "feat")
  colnames(data)[ncol(data)] <- "length_score"
  data$length_score <- sapply(data$length_score, function(x) ifelse(is.na(x), 0, x))
  data$feat <- NULL

  ######################## summation ###########################
  data$total_score <- data$vgene_score + data$Jmotif_score + data$p107_score + data$p113_score + data$pos_score + data$perc_score + data$length_score
  ## scaling by the mean and standard deviation of originally held-out data to standardize the TiRP scale
  data$vTiRP <- (data$vgene_score + data$p107_score + 0.1459054) / 0.2364
  data$mTiRP <- (data$perc_score + data$pos_score + data$length_score - 0.03846178) / 0.2364
  data$jTiRP <- (data$Jmotif_score + data$p113_score - 0.07454362) / 0.2364
  data$TiRP <- (data$total_score + 0.0329) / 0.2364
  data <- data[, !(grepl("perc_mid", colnames(data)))]

  data <- tibble::as_tibble(data)

  if (!.details) {
    data <- data$TiRP
  }

  return(data)
}

reformat_vgene_cp_modified <- function(vg) {
  vg <- as.character(vg)
  vgene <- gsub("TRBV", "TCRBV", vg)
  vgene <- gsub("\\*0[[:digit:]]$", "", vgene)
  info <- strsplit(vgene, "-")
  res <- vapply(X = info, FUN.VALUE = character(1), FUN = function(.info) {
    family <- .info[1]
    member <- ifelse(length(.info) == 2, .info[2], "01")
    family <- ifelse(nchar(family) == 6, paste("TCRBV0", substr(family, 6, 6), sep = ""), family)
    member <- ifelse(substr(member, 1, 1) == "0", member, paste("0", member, sep = ""))
    paste(family, member, sep = "-")
  })
  return(res)
}

fix_vdj_genes2 <- function(.x) {

  res <- stringr::str_replace_all({{ .x }}, ",", ", ")
  res <- stringr::str_replace_all(res, "-([0])([0-9])", "-\\2")
  res <- stringr::str_replace_all(res, "([VDJ])([0])([0-9])", "\\1\\3")
  res <- stringr::str_replace_all(res, "TCR", "TR")
  res <- stringr::str_replace_all(res, "\\*[0-9][0-9]", "")

  return(res)
}

get_perc_score <- function(percents, touse, perc_term) {
  df <- data.frame(perc_term, percents)
  df <- df[df$perc_term %in% touse$feat, ]
  df <- dplyr::left_join(df, touse, by = c("perc_term" = "feat"))
  df$prod <- df$percents * df$metabeta
  return(sum(df$prod))
}

get_pos_score <- function(x, touse, pos_terms) {
  terms <- pos_terms[[nchar(x) - 1]]
  terms <- sapply(1:length(terms), function(y) paste(terms[y], substr(x, y, y), sep = "_"))
  return(sum(touse$metabeta[touse$feat %in% terms]))
}
