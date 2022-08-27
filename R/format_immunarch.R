#' Format Imunoseq Data for use with Immunarch
#'
#' Prepare Immunoseq data according to Immunarch data format
#'
#' Works on Immunoseq data as returned by the [cdr3tools::read_immunoseq]
#'   function. If unique sequences were collapsed using
#'   [cdr3tools::collapse_sequences] the resulting column `sequence` will be a
#'   list column containing the rearrangment sequences from each of the
#'   duplicate unique sequences.
#'
#'   __*A special note about list columns:*__ If formatting data as read by
#'   [cdr3tools::read_immunoseq()] that has had unique sequences collapsed
#'   using [cdr3tools::collapse_sequences()], there will be list columns in the
#'   output. These list columns contain the duplicated (but error containing)
#'   full length sequences. [cdr3tools::format_immunarch()] will select the
#'   first sequence from each list column row as the representative sequence
#'   in order to dissolve the list columns. This ensures compatibility with
#'   Immunarch functions.
#'
#' @param data A list or a single data frame as returned by
#'   [cdr3tools::read_immunoseq()].
#' @returns A list or a single data frame of Immunarch formatted data.
#' @author Christopher Parks
#' @references \url{https://immunarch.com/reference/immunr_data_format.html}
#' @examples
#' format_immunarch(cfselo_seqs)
#' @family Immunoseq
#' @export
format_immunarch <- function(data) {

  # add data to a list if it's not already
  if (!inherits(data, "list") && inherits(data, "data.frame")) {
    data <- list(data)
  }

  if (!inherits(data, "list")) {
    rlang::abort("data must be a list or single data frame as returned by cdr3tools::read_immunoseq.")
  }

  # detect list columns in data
  list_cols <- any(
    vapply(
      data,
      function(data) {
        any(
          vapply(
            data,
            is.list,
            logical(1),
            USE.NAMES = FALSE
          )
        )
      },
      logical(1),
      USE.NAMES = FALSE
    )
  )

  # unnest list columns by selecting the first element in each list_column
  #   and converting to character vector.
  if (list_cols) {
    data <- purrr::map(data, function(data) {
      data <- data %>%
        dplyr::mutate(dplyr::across(
          .cols = where(is.list),
          .fns = ~ purrr::map_chr(.x, dplyr::first))
          )
    }
    )
  }

  # Format immunoseq files to immunarch format
  data <- purrr::map(data, function(data) {
    data <- data %>%

      # define CDR3 positions the immunarch way
      dplyr::mutate(
        V.end = .data$n1_index - .data$v_index,
        D.start = .data$d_index - .data$v_index,
        D.end = .data$n2_index - .data$v_index,
        J.start = .data$j_index - .data$v_index
      ) %>%

      # where a region isn't defined, replace the value with NA
      #   (I validated this with IMGT V-Quest)
      dplyr::mutate(
        dplyr::across(
          .cols = c(.data$V.end, .data$D.start, .data$D.end, .data$J.start),
          .fns = ~ dplyr::case_when(
            .x < 0 ~ NA_real_,
            TRUE   ~ .x
          )
        )
      ) %>%

      # select the columns needed for immunarch format. Rename then as we go.
      dplyr::select(
        Clones = .data$templates,
        Duplicates = tidyselect::any_of("duplicates"),
        CDR3.nt = .data$cdr3_rearrangement,
        CDR3.aa = .data$cdr3_amino_acid,
        V.name = .data$v_resolved,
        D.name = .data$d_resolved,
        J.name = .data$j_resolved,
        .data$V.end,
        .data$D.start,
        .data$D.end,
        .data$J.start,
        VD.ins = .data$n1_insertions,
        DJ.ins = .data$n2_insertions,
        Sequence = .data$rearrangement
      ) %>%

      # calculate proportion
      dplyr::mutate(Proportion = .data$Clones / sum(.data$Clones), .after = .data$Clones) %>%

      # set V-J insertions to NA by definition since these are TRB sequences
      dplyr::mutate(VJ.ins = NA, .after = .data$J.start) %>%

      # convert immunoseq VDJ gene names to IMGT gene names
      dplyr::mutate(V.name = cdr3tools::imgt_format_gene_names(.data$V.name)) %>%
      dplyr::mutate(D.name = cdr3tools::imgt_format_gene_names(.data$D.name)) %>%
      dplyr::mutate(J.name = cdr3tools::imgt_format_gene_names(.data$J.name)) %>%

      # sort the data according to abundance
      dplyr::arrange(dplyr::desc(.data$Proportion))

    attr(data, "repertoire_data_format") <- "cdr3tools_immunarch"
    data
  }
  )

  if (length(data) == 1) {
    data <- data[[1]]
  }

  return(data)

}
