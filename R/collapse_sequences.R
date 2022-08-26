#' Collapse Unique Sequences
#'
#' Collapse repertoires using various combinations of unique sequence
#'   definitions.
#'
#' Rearrangement sequences are grouped according to the user selected unique
#'   sequence definition and returned as list columns. The returned data will
#'   contain one row per unique sequence. A column `duplicates` is added to
#'   show the number of rearrangement sequences that match to a single
#'   unique sequence. This function is useful for minimizing PCR and sequencing
#'   errors that result in artificially inflated unique sequence counts.
#'
#' @param data A list or a single data frame as returned by
#'   [cdr3tools::read_immunoseq()].
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
#' @family Immunoseq
#' @export
collapse_sequences <- function(data, .group_by = NULL) {

  if (!inherits(data, "list") && inherits(data, "data.frame")) {
    data <- list(data)
  }

  if (is.null(.group_by)) {
    .group_by <- "ntvj"
  }

  .group_by <- rlang::arg_match(.group_by, c("nt", "ntvj", "aa", "aavj"))

  data <- purrr::map(data, function(data) {

    if (attr(data, "repertoire_data_format") == "cdr3tools_immunoseq_v1") {

      .group_by <- switch(
        .group_by,
        "nt"   = c("cdr3_rearrangement"),
        "ntvj" = c("cdr3_rearrangement", "v_resolved", "j_resolved"),
        "aa"   = c("cdr3_amino_acid"),
        "aavj" = c("cdr3_amino_acid", "v_resolved", "j_resolved")
      )

      data <- data %>%

        dtplyr::lazy_dt(key_by = c(.group_by, "sample_name", "frame_type")) %>%

        dplyr::group_by(
          dplyr::across(
            c(
              tidyselect::all_of({{ .group_by }}),
              .data$sample_name,
              .data$frame_type
            )
          )
        ) %>%

        dplyr::summarise(
          duplicates = dplyr::n(),
          templates = sum(.data$templates),
          seq_reads = sum(.data$seq_reads),
          cdr3_rearrangement = .data$cdr3_rearrangement[1],
          cdr3_amino_acid = .data$cdr3_amino_acid[1],
          v_resolved = .data$v_resolved[1],
          d_resolved = .data$d_resolved[1],
          j_resolved = .data$j_resolved[1],
          n1_index = .data$n1_index[1],
          v_index = .data$v_index[1],
          d_index = .data$d_index[1],
          n2_index = .data$n2_index[1],
          j_index = .data$j_index[1],
          n1_insertions = .data$n1_insertions[1],
          n2_insertions = .data$n2_insertions[1],
          rearrangement = list(.data$rearrangement),
          rearrangement_trunc = list(.data$rearrangement_trunc),
          .groups = "drop"
        ) %>%

        dplyr::mutate(frequency = .data$templates / sum(.data$templates)) %>%

        dplyr::arrange(-.data$frequency) %>%

        dplyr::select(
          .data$sample_name, .data$frame_type, .data$duplicates, .data$templates,
          .data$frequency, .data$seq_reads, .data$cdr3_rearrangement,
          .data$cdr3_amino_acid, .data$v_resolved, .data$d_resolved,
          .data$j_resolved, .data$n1_index, .data$v_index, .data$d_index,
          .data$n2_index, .data$j_index, .data$n1_insertions, .data$n2_insertions,
          .data$rearrangement, .data$rearrangement_trunc
        ) %>%

        tibble::as_tibble()

      attr(data, "repertoire_data_format") <- "cdr3tools_immunoseq_v1"
      data

    } else if (attr(data, "repertoire_data_format") == "cdr3tools_immunarch") {

      .group_by <- switch(
        .group_by,
        "nt"   = c("CDR3.nt"),
        "ntvj" = c("CDR3.nt", "V.name", "J.name"),
        "aa"   = c("CDR3.aa"),
        "aavj" = c("CDR3.aa", "V.name", "J.name")
      )

      data <- data %>%

        dtplyr::lazy_dt(key_by = .group_by) %>%

        dplyr::group_by(
          dplyr::across(
            tidyselect::all_of({{ .group_by }})
          )
        ) %>%

        dplyr::summarise(
          Duplicates = dplyr::n(),
          Clones = sum(.data$Clones),
          CDR3.nt = .data$CDR3.nt[1],
          CDR3.aa = .data$CDR3.aa[1],
          V.name = .data$V.name[1],
          D.name = .data$D.name[1],
          J.name = .data$J.name[1],
          V.end = .data$V.end[1],
          D.start = .data$D.start[1],
          D.end = .data$D.end[1],
          J.start = .data$J.start[1],
          VJ.ins = .data$VJ.ins[1],
          VD.ins = .data$VD.ins[1],
          DJ.ins = .data$DJ.ins[1],
          Sequence = list(.data$Sequence),
          .groups = "drop"
        ) %>%

        dplyr::mutate(Proportion = .data$Clones / sum(.data$Clones)) %>%

        dplyr::arrange(-.data$Proportion) %>%

        dplyr::select(
          .data$Clones, .data$Proportion, .data$Duplicates, .data$CDR3.nt,
          .data$CDR3.aa, .data$V.name, .data$D.name, .data$J.name, .data$V.end,
          .data$D.start, .data$D.end, .data$J.start, .data$VJ.ins, .data$VD.ins,
          .data$DJ.ins, .data$Sequence
        ) %>%

        tibble::as_tibble()

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
