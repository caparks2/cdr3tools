#' Collapse Unique Sequences
#'
#' Collapse repertoires using various combinations of unique sequence
#'   definitions
#' @param .group_by Character. The combination of T cell receptor sequences and
#'   V and J genes to use in order to group unique sequences. If `NULL`, the
#'   default is to group unique sequences by CDR3 nucleotide sequence plus the
#'   most resolved V and J gene calls. Allowed inputs to `.group_by` are:
#'     \describe{
#'       \item{"nt"}{CDR3 nucleotide sequence}
#'       \item{"ntvj"}{*Default*. CDR3 nucleotide sequence plus most resolved V and J gene calls}
#'       \item{"aa"}{CDR3 amino acid sequence}
#'       \item{"aavj"}{CDR3 amino acid sequence plus most resolved V and J gene calls}
#'       \item{"rearrangement"}{full nucleotide sequence (not suitable for comparisons across platforms)}
#'       \item{"extended_rearrangement"}{extended full nucleotide sequence (not suitable for comparisons across platforms)}
#'       \item{"rearrangement_trunc"}{truncated full nucleotide sequence}
#'     }
#' @param .use_reads Logical. `FALSE` (the default) indicates that repertoire
#'   frequencies should be calculated using templates. `TRUE` indicates that
#'   repertoire frequencies should be calculated using reads. Immunoseq files
#'   are sometimes exported without reads information if templates are present
#'   in the data. In this case `` will throw an
#'   error if `.use_reads` = TRUE.
#' @name collapse_sequences
NULL
# if (is.null(.use_reads)) {
#   .use_reads <- FALSE
# }
# if (.use_reads) {
#   purrr::map(data, function(data) {
#     if (all(is.na(data$seq_reads))) {
#       rlang::abort("data$seq_reads has no values. Set .use_reads = FALSE to use templates.")
#     }
#   })
# }
#
# if (!is.null(.group_by)) {
#
#   .group_by <- rlang::arg_match(.group_by,
#                                 c(
#                                   "nt", "ntvj", "aa", "aavj",
#                                   "rearrangement",
#                                   "rearrangement_trunc",
#                                   "extended_rearrangement"
#                                 )
#   )
#   .group_by <- switch(.group_by,
#                       nt                     = "cdr3_rearrangement",
#                       ntvj                   = c("cdr3_rearrangement", "v_resolved", "j_resolved"),
#                       aa                     = "cdr3_amino_acid",
#                       aavj                   = c("cdr3_amino_acid", "v_resolved", "j_resolved"),
#                       rearrangement          = "rearrangement",
#                       rearrangement_trunc    = "rearrangement_trunc",
#                       extended_rearrangement = "extended_rearrangement"
#   )
#   data <- purrr::map(data, function(data) {
#     grouped_data <- data %>%
#       dplyr::group_by(dplyr::across(tidyselect::all_of({{ .group_by }}))) %>%
#       dplyr::summarise(templates = sum(.data$templates),
#                        seq_reads = sum(.data$seq_reads)) %>%
#       dplyr::ungroup()
#     data <- data %>%
#       dplyr::select(-templates, -seq_reads, -rearrangement, -rearrangement_trunc, -extended_rearrangement) %>%
#       dplyr::inner_join(grouped_data, by = {{ .group_by }})
#     data
#   }
#   )
# }
#
# return(data)
#
# if (.use_reads) {
#   data <- purrr::map(data, function(x) {
#     x %>%
#       # add a seq_read frequency column
#       dplyr::mutate(
#         frequency = .data$seq_reads / sum(.data$seq_reads),
#         .after = .data$seq_reads
#       ) %>%
#       # sort by frequency
#       dplyr::arrange(dplyr::desc(.data$frequency))
#   }
#   )
# } else {
#   data <- purrr::map(data, function(x) {
#     x %>%
#       # add a template frequency column
#       dplyr::mutate(
#         frequency = .data$templates / sum(.data$templates),
#         .after = .data$templates
#       ) %>%
#       # sort by frequency
#       dplyr::arrange(dplyr::desc(.data$frequency))
#   }
#   )
# }
