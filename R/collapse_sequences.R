#' Collapse Unique Sequences
#'
#' Collapse repertoires using various combinations of unique sequence
#'   definitions.
#'
#' Rearrangement sequences are grouped according to the user selected unique
#'   sequence definition and returned as list columns. The returned data will
#'   contain one row per unique sequence. A column `duplicates` is added to
#'   show the number of rearrangement sequences that match to a single
#'   unique sequence. This is useful for minimizing PCR and sequencing errors that
#'   result in artificially inflated unique sequence counts.
#'
#' @param data A list or a single data frame as returned by
#'   [cdr3tools::read_immunoseq].
#' @param .group_by Character. The combination of T cell receptor sequences and
#'   V and J genes to use in order to group unique sequences. If `NULL`, the
#'   default is to group unique sequences by CDR3 nucleotide sequence plus the
#'   most resolved V and J gene calls. Allowed inputs to `.group_by` are:
#'     \describe{
#'       \item{"nt"}{CDR3 nucleotide sequence}
#'       \item{"ntvj"}{*Default*. CDR3 nucleotide sequence plus most resolved V and J gene calls}
#'       \item{"aa"}{CDR3 amino acid sequence}
#'       \item{"aavj"}{CDR3 amino acid sequence plus most resolved V and J gene calls}
#'     }
#' @returns A list or a single data frame of Immunoseq data with
#' @examples
#' collapse_sequences(cfselo_seqs)
#' @author Christopher Parks
#' @export
collapse_sequences <- function(data, .group_by = NULL) {
  if (!inherits(data, "list") && inherits(data, "data.frame")) {
    data <- list(data)
  }

  if (is.null(.group_by)) .group_by <- "ntvj"

  .group_by <- rlang::arg_match(.group_by, c("nt", "ntvj", "aa", "aavj"))

  .group_by <- switch(.group_by,
    "nt"   = c("cdr3_rearrangement"),
    "ntvj" = c("cdr3_rearrangement", "v_resolved", "j_resolved"),
    "aa"   = c("cdr3_amino_acid"),
    "aavj" = c("cdr3_amino_acid", "v_resolved", "j_resolved")
  )

  data <- purrr::map(data, ~ .x %>%
    dplyr::group_by(
      dplyr::across(
        c(tidyselect::all_of(.group_by), sample_name, frame_type)
      )
    ) %>%
    dplyr::summarise(
      duplicates = dplyr::n(),
      templates = sum(templates),
      seq_reads = sum(seq_reads),
      cdr3_rearrangement = cdr3_rearrangement[1],
      cdr3_amino_acid = cdr3_amino_acid[1],
      v_resolved = v_resolved[1],
      d_resolved = d_resolved[1],
      j_resolved = j_resolved[1],
      n1_index = n1_index[1],
      v_index = v_index[1],
      d_index = d_index[1],
      n2_index = n2_index[1],
      j_index = j_index[1],
      n1_insertions = n1_insertions[1],
      n2_insertions = n2_insertions[1],
      rearrangement = list(rearrangement),
      rearrangement_trunc = list(rearrangement_trunc),
      .groups = "drop"
    ) %>%
    dplyr::mutate(frequency = templates / sum(templates)) %>%
    dplyr::arrange(dplyr::desc(frequency)) %>%
    dplyr::select(
      sample_name, frame_type, duplicates, templates, frequency, seq_reads,
      cdr3_rearrangement, cdr3_amino_acid, v_resolved, d_resolved,
      j_resolved, n1_index, v_index, d_index, n2_index, j_index,
      n1_insertions, n2_insertions, rearrangement, rearrangement_trunc
    ))

  if (length(data) == 1) {
    return(data[[1]])
  } else {
    return(data)
  }
}
