#' Remove Contaminant Clones
#'
#' Remove contaminant TCR sequences from Immunoseq or Immunarch repertoire data.
#'
#' Accepts amino acid sequences and nucleotide sequences corresponding to
#'   Immunoseq export v1 CDR3_amino_acid and CDR3_rearrangement/rearrangement
#'   columns and Immunarch CDR3.aa and CDR3.nt/Sequence columns.
#'
#' @param data A list or a single data frame as returned by
#'   [cdr3tools::read_immunoseq()] or [cdr3tools::format_immunarch].
#' @param seqs Character vector of CDR3 nucleotide and amino acid sequences.
#' @returns A list or a single data frame of repertoire data in the same format
#'   as `data`, but with contaminant sequence rows removed.
#' @author Christopher Parks
#' @examples
#' remove_contaminants(
#'   cfselo_seqs,
#'   c("CASRDRGSYGYTF", "TGTGCCAGCAGGGCGACTAGCGGGAGGGCCTACGAGCAGTACTTC")
#' )
#' @family Immunoseq
#' @export
remove_contaminants <- function(data, seqs) {

  if (!inherits(data, "list") && inherits(data, "data.frame")) {
    data <- list(data)
  }

  if (!inherits(seqs, "character")) {
    rlang::abort("... must contain only nucleotide and amino acid sequences")
  }

  nt <- seqs[grepl("^[ACTG]+$", seqs, ignore.case = TRUE)]
  nt <- paste(nt, collapse = "|")

  aa <- seqs[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", seqs, ignore.case = TRUE) &
               !grepl("^[ACTG]+$", seqs, ignore.case = TRUE)]
  aa <- paste(aa, collapse = "|")

  data <- lapply(data, function(data) {

      if (attr(data, "repertoire_data_format") == "cdr3tools_immunoseq_v1") {

        data <- data[!grepl(aa, data$cdr3_amino_acid, ignore.case = TRUE), ]
        data <- data[!grepl(nt, data$rearrangement, ignore.case = TRUE), ]
        attr(data, "repertoire_data_format") <- "cdr3tools_immunoseq_v1"
        data

      } else if (attr(data, "repertoire_data_format") == "cdr3tools_immunarch") {

        data <- data[!grepl(aa, data$CDR3.aa, ignore.case = TRUE), ]
        data <- data[!grepl(nt, data$Sequence, ignore.case = TRUE), ]
        attr(data, "repertoire_data_format") <- "cdr3tools_immunarch"
        data

      }

    }
  )

  if (length(data) == 1) {
    data <- data[[1]]
  }

  return(data)
}
