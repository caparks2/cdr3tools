#' Calculate TiRP Scores
#'
#' A scoring system to quantify TCR-intrinsic regulatory potential (TiRP)
#'
#' A special note should be made to consider that this function, as written by
#'   the authors, evaluates only sequences that are between 12 and 17 amino
#'   acids in length and that evaluation is performed on sequences as they are,
#'   in their linear, native states without first being adjusted for length by
#'   IMGT unique numbering rules for JUNCTIONS. It was also discovered that the
#'   original version of this function as would fail to produce "vgene scores",
#'   resulting in incorrect overall scores, if allele identifiers were part of
#'   the V gene names. This was corrected in this implementation of the
#'   calculation. Additionally, the function has been modified to use IMGT
#'   conventions in naming of V genes and will check this for you. Further
#'   changes have been added so that errors result in informative messages, the
#'   input can now accept multiple data files, the output can be simple or
#'   detailed depending on user's choice, and a much speedier overall execution
#'   time has been achieved.
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
#' @author Kaitlyn A. Lagattuta
#' @author Joyce B. Kang
#' @author Aparna Nathan
#' @author Center for Data Sciences, Brigham and Women’s Hospital, Boston, MA, USA
#' @author Christopher Parks (I rewrote and improved for use in this package.)
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
#' @export
get_TiRP_scores <- function(.data, .details = FALSE) {

  data <- .data
  details <- .details

  if (!(inherits(data, "list") || inherits(data, "data.frame"))) {
    rlang::abort(
      paste(
        ".data must either be a single data frame (with V genes in the first",
        "column and CDR3 amino acid sequences in the second column) or a list",
        "containing multiple of such data frames."
      )
    )
  }

  if (inherits(data, "list")) {
    result <- lapply(data, function(data) get_TiRP_scores_internal(data, details))
    return(result)
  }

  if (inherits(data, "data.frame")) {
    result <- get_TiRP_scores_internal(data, details)
    return(result)
  }

}

get_TiRP_scores_internal <- function(data, details) {

  data <- as.data.frame(data)

  if (!all(grepl("TRBV|TCRBV", data[!is.na(data[, 1]), 1]))) {
    rlang::abort(
      paste(
        "Column 1 must be a character vector of TRBV genes."
      )
    )
  }

  if (!all(grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", data[!is.na(data[, 2]), 2]))) {
    rlang::abort(
      paste(
        "Column 2 must be a character vector of IMGT JUNCTION sequences.",
        "The only characters allowed are in the set [ACDEFGHIKLMNPQRSTVWY]."
      )
    )
  }

  data <- data[!(is.na(data[, 1]) | is.na(data[, 2])), ]

  weights <- cdr3tools::TiRP_weights

  # defining TCR features

  # standardizing TRBV gene names
  check_v_genes <- sub("vgene", "", weights[grep("vgene", weights$feat), "feat"])

  if (!(sum(data[, 1] %in% check_v_genes) == length(data[, 1]))) {
    data$vgene <- cdr3tools::imgt_format_gene_names(data[, 1])
  } else {
    data$vgene <- as.character(data[, 1])
  }

  # CDR3 sequence
  data$cdr3 <- as.character(data[, 2])

  # CDR3 length
  data$length <- nchar(data$cdr3)
  data <- data[data$length >= 12 & data$length <= 17, ]

  if (nrow(data) == 0) {
    rlang::abort("Only sequences of lengths from 12 to 17 amino acids are analyzed. No such sequences were detected.")
  }

  # Jmotif
  jmotifs <- gsub("Jmotif", "", weights$feat[grepl("Jmotif", weights$feat)])
  data$Jmotif <- substr(data$cdr3, nchar(data$cdr3) - 4, nchar(data$cdr3))
  data$Jmotif[!(data$Jmotif %in% jmotifs)] <- "other"

  # amino acid composition in the CDR3 middle region
  data$cdr3MR <- substr(data$cdr3, 5, nchar(data$cdr3) - 6)
  aminos <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
              "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  ref_stats <- cdr3tools::heldout_means_sds

  # For loops are very slow in R. Rewriting using a vectorised approach.
  # for (i in 1:length(aminos)) {
  #   new <- sapply(data$cdr3MR, function(x) stringr::str_count(x, aminos[i]) / nchar(x))
  #   mn <- ref_stats$mean[ref_stats$amino == aminos[i]]
  #   sd <- ref_stats$sd[ref_stats$amino == aminos[i]]
  #   new <- (new - mn) / sd
  #   data <- cbind(data, new)
  #   colnames(data)[ncol(data)] <- paste("perc_mid", aminos[i], sep = "_")
  # }

  perc_mid_cols <- do.call(cbind, lapply(aminos, function(aminos) {
    new <- stringr::str_count(data$cdr3MR, aminos) / nchar(data$cdr3MR)
    mn <- ref_stats$mean[ref_stats$amino == aminos]
    sd <- ref_stats$sd[ref_stats$amino == aminos]
    new <- (new - mn) / sd
    colname <- paste("perc_mid", aminos, sep = "_")
    cols <- data.frame(new)
    names(cols) <- colname
    cols
  }))

  data <- cbind(data, perc_mid_cols)


  # scoring TCR features

  # including only terms that were Bonferroni-significant in meta-analysis:
  touse <- weights[weights$metap < (0.05 / nrow(weights)), ]

  perc_terms <- touse$feat[grepl("perc_mid", touse$feat)]
  # data$perc_score <- sapply(1:nrow(data), function(x) sum(data[x, perc_terms] * touse[perc_terms, c("metabeta")]))
  perc_score <- do.call(rbind, as.list(apply(data[, perc_terms], 1, function(rows) {
    rows * touse[perc_terms, c("metabeta")]
  }, simplify = FALSE)))
  data$perc_score <- rowSums(perc_score)

  pos_terms <- list()
  pos_terms[[1]] <- c("12_5", "12_6")
  pos_terms[[2]] <- c("13_5", "13_6", "13_7")
  pos_terms[[3]] <- c("14_5", "14_6", "14_7", "14_-7")
  pos_terms[[4]] <- c("15_5", "15_6", "15_7", "15_8", "15_-7")
  pos_terms[[5]] <- c("16_5", "16_6", "16_7", "16_8", "16_-8", "16_-7")
  pos_terms[[6]] <- c("17_5", "17_6", "17_7", "17_8", "17_9", "17_-8", "17_-7")
  data$pos_score <- get_pos_score_CP(data$cdr3MR, touse, pos_terms)

  data$feat <- paste("vgene", data$vgene, sep = "")
  data <- dplyr::left_join(
    dtplyr::lazy_dt(data),
    touse[, c("feat", "metabeta")],
    by = "feat"
  )
  data <- as.data.frame(data)
  colnames(data)[ncol(data)] <- "vgene_score"
  data$vgene_score[is.na(data$vgene_score)] <- 0
  # data$vgene_score <- sapply(data$vgene_score, function(x) ifelse(is.na(x), 0, x))
  if (nrow(data[data$vgene_score != 0, ]) == 0) {
    rlang::warn(
      paste(
        "TRBV genes could not be matched to internal data. Make sure TRBV",
        "gene names conform to IMGT conventions AND that they lack allele",
        "identifiers (ex. take out the trailing *01). Scores are returned",
        "without including V gene scores and may not be informative."
      )
    )
  }

  data$feat <- paste("p107", substr(data$cdr3, 4, 4), sep = "")
  LHS <- dtplyr::lazy_dt(data, key_by = "feat")
  RHS <- dtplyr::lazy_dt(touse[, c("feat", "metabeta")], key_by = "feat")
  data <- dplyr::left_join(LHS, RHS, by = "feat")
  # data <- dplyr::left_join(
  #   dtplyr::lazy_dt(data),
  #   touse[, c("feat", "metabeta")],
  #   by = "feat"
  # )
  data <- as.data.frame(data)
  colnames(data)[ncol(data)] <- "p107_score"
  data$p107_score[is.na(data$p107_score)] <- 0
  # data$p107_score <- sapply(data$p107_score, function(x) ifelse(is.na(x), 0, x))

  data$feat <- paste("Jmotif", data$Jmotif, sep = "")
  LHS <- dtplyr::lazy_dt(data, key_by = "feat")
  RHS <- dtplyr::lazy_dt(touse[, c("feat", "metabeta")], key_by = "feat")
  data <- dplyr::left_join(LHS, RHS, by = "feat")
  # data <- dplyr::left_join(
  #   dtplyr::lazy_dt(data),
  #   touse[, c("feat", "metabeta")],
  #   by = "feat"
  # )
  data <- as.data.frame(data)
  colnames(data)[ncol(data)] <- "Jmotif_score"
  data$Jmotif_score[is.na(data$Jmotif_score)] <- 0
  # data$Jmotif_score <- sapply(data$Jmotif_score, function(x) ifelse(is.na(x), 0, x))

  # data$feat <- sapply(data$cdr3, function(x) substr(x, nchar(x) - 5, nchar(x) - 5), USE.NAMES = F)
  data$feat <- paste("p113", substr(data$cdr3, nchar(data$cdr3) - 5, nchar(data$cdr3) - 5), sep = "")
  LHS <- dtplyr::lazy_dt(data, key_by = "feat")
  RHS <- dtplyr::lazy_dt(touse[, c("feat", "metabeta")], key_by = "feat")
  data <- dplyr::left_join(LHS, RHS, by = "feat")
  # data <- dplyr::left_join(
  #   dtplyr::lazy_dt(data),
  #   touse[, c("feat", "metabeta")],
  #   by = "feat"
  # )
  data <- as.data.frame(data)
  colnames(data)[ncol(data)] <- "p113_score"
  data$p113_score[is.na(data$p113_score)] <- 0
  # data$"p113_score" <- sapply(data$"p113_score", function(x) ifelse(is.na(x), 0, x))

  data$feat <- paste("length", data$length, sep = "")
  LHS <- dtplyr::lazy_dt(data, key_by = "feat")
  RHS <- dtplyr::lazy_dt(touse[, c("feat", "metabeta")], key_by = "feat")
  data <- dplyr::left_join(LHS, RHS, by = "feat")
  # data <- dplyr::left_join(
  #   dtplyr::lazy_dt(data),
  #   touse[, c("feat", "metabeta")],
  #   by = "feat"
  # )
  data <- as.data.frame(data)
  colnames(data)[ncol(data)] <- "length_score"
  data$length_score[is.na(data$length_score)] <- 0
  # data$length_score <- sapply(data$length_score, function(x) ifelse(is.na(x), 0, x))
  # data$feat <- NULL
  # data$feat <- NA_character_

  data <- data[, -which(names(data) == "feat")]

  # Summation

  # data$total_score <- data$vgene_score + data$Jmotif_score + data$p107_score + data$p113_score + data$pos_score + data$perc_score + data$length_score
  data$total_score <- with(data,
    vgene_score + Jmotif_score + p107_score + p113_score + pos_score + perc_score + length_score
  )

  # scaling by the mean and standard deviation of originally held-out data to standardize the TiRP scale
  data$vTiRP <- (data$vgene_score + data$p107_score + 0.1459054) / 0.2364
  data$mTiRP <- (data$perc_score + data$pos_score + data$length_score - 0.03846178) / 0.2364
  data$jTiRP <- (data$Jmotif_score + data$p113_score - 0.07454362) / 0.2364
  data$TiRP <- (data$total_score + 0.0329) / 0.2364
  data <- data[, !(grepl("perc_mid", colnames(data)))]

  data <- tibble::as_tibble(data)

  if (!details) {
    data <- data$TiRP
  }

  return(data)
}

get_pos_score_CP <- function(data, touse, pos_terms) {
  term_indices <- pos_terms[nchar(data) - 1]
  residues <- strsplit(data, "")
  terms <- mapply(function(.x, .y) paste(.x, .y, sep = "_"), .x = term_indices, .y = residues, SIMPLIFY = FALSE)
  vapply(terms, function(terms) sum(touse$metabeta[touse$feat %in% terms]), numeric(1))
}

# get_perc_score <- function(percents, touse, perc_term) {
#   df <- data.frame(perc_term, percents)
#   df <- df[df$perc_term %in% touse$feat, ]
#   df <- dplyr::left_join(df, touse, by = c("perc_term" = "feat"))
#   df$prod <- df$percents * df$metabeta
#   return(sum(df$prod))
# }

# get_pos_score <- function(x, touse, pos_terms) {
#   terms <- pos_terms[[nchar(x) - 1]]
#   terms <- sapply(1:length(terms), function(y) paste(terms[y], substr(x, y, y), sep = "_"))
#   return(sum(touse$metabeta[touse$feat %in% terms]))
# }
