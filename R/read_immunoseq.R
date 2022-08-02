#' Read Adaptive Biotechnologies Immunoseq Files
#'
#' A function to read Immunoseq v1 exported .tsv files.
#'
#' This function reads Adaptive Biotechnologies Immunoseq files into R. Only
#'   necessary columns from the Immunoseq files are read, with the others being
#'   dropped before being read. This was done to preserve RAM. Given a directory
#'   containing multiple files, the function reads them all into a named list.
#'   If only a path to one file is given, the function reads the file to a data
#'   frame. At this time, only Immunoseq v1 .tsv file exports are supported.
#'
#' @param .path A character vector (string). The full or relative (for example,
#'  if the working directory was set with `setwd()`) file path to a directory
#'  containing multiple Immunoseq v1 .tsv file exports. Alternatively, the full
#'  or relative file path to a single Immunoseq v1 .tsv file export.
#' @param .functional A boolean. `TRUE` (the default) filters sequences for
#'   in-frame, functional sequences. `FALSE` leaves non-coding, out-of-frame
#'   sequences in.
#'
#' @author Christopher Parks
#'
#' @returns if `.path` is a directory containing multiple files, a list of same
#'   length as the directory. if `.path` is a single file, a data frame.
#' @export
read_immunoseq <- function(.path, .functional = TRUE) {

  if (!dir.exists(.path)) {
    files <- .path
    file_names <- basename(.path)
    file_names <- stringr::str_remove_all(file_names, stringr::regex("\\.tsv", ignore_case = TRUE))
    file_names <- stringr::str_remove_all(file_names, stringr::regex("_TCR[AB]$", ignore_case = TRUE))
    names(files) <- file_names
  } else if (dir.exists(.path)) {
    files <- list.files(.path, pattern = "\\.tsv", full.names = TRUE, ignore.case = TRUE)
    file_names <- list.files(.path, pattern = "\\.tsv", full.names = FALSE, ignore.case = TRUE)
    file_names <- stringr::str_remove_all(file_names, stringr::regex("\\.tsv", ignore_case = TRUE))
    file_names <- stringr::str_remove_all(file_names, stringr::regex("_TCR[AB]$", ignore_case = TRUE))
    names(files) <- file_names
  }

  # check for immunoseq version
  if (!check_immunoseq_version(files)) {
    rlang::abort(
      paste(
        ".path must be a path to a directory of immunoseq v1 exported .tsv",
        "text files or a path to a single immunoseq v1 exported .tsv text file"
      )
    )
  }

  # read in immunoseq files while selecting needed columns and dropping the rest.
  data <- purrr::map(
    files, function(x) {
      readr::read_tsv(
        x,
        show_col_types = FALSE,
        progress = FALSE,
        col_select = c(
          "sample_name", "frame_type", "templates", "seq_reads", "cdr3_amino_acid",
          "cdr3_rearrangement", "v_resolved", "d_resolved", "j_resolved",
          "n1_index", "v_index", "d_index", "n2_index", "j_index",
          "n1_insertions", "n2_insertions", "rearrangement"
        )
      )
    }
  )

  # filter for in-frame, functional sequences.
  if (.functional) {
    data <- purrr::map(data, function(x) {
      x %>% dplyr::filter(.data$frame_type == "In")
    })
  }

  # correct VDJ gene names to IMGT conventions.
  data <- purrr::map(
    data, function(x) {
      x %>%
        dplyr::mutate(
          v_resolved = cdr3tools::imgt_format_gene_names(.data$v_resolved),
          d_resolved = cdr3tools::imgt_format_gene_names(.data$d_resolved),
          j_resolved = cdr3tools::imgt_format_gene_names(.data$j_resolved)
        ) %>%
        # add a template frequency column
        dplyr::mutate(
          frequency = .data$templates / sum(.data$templates),
          .after = .data$templates
        ) %>%
        # sort by frequency
        dplyr::arrange(dplyr::desc(.data$frequency))
    }
  )

  # if (length(data) == 1) {
  #   data <- data[[1]]
  # }

  return(data)
}

check_immunoseq_version <- function(files) {

  if (!any(stringr::str_detect(files, "\\.tsv"))) {
    rlang::abort(
      paste(
        ".path must be a path to a directory of immunoseq v1 exported .tsv",
        "text files or a path to a single immunoseq v1 exported .tsv text file"
      )
    )
  }

  version_check <- lapply(
    files, function(x) {

      con <- tryCatch(
        { file(x, "r") },
        error = function(c) return(FALSE),
        warning = function(c) return(FALSE)
      )

      if (is.logical(con)) {
        if (!con) {
          return(FALSE)
        }
      }

      header <- readLines(con = con, n = 1)
      close(con)
      header <- strsplit(header, "\t")[[1]]
      matches <- header %in% cdr3tools::immunoseq_export_v1_col_names
      all(matches)
    }
  )

  return(all(unlist(version_check)))
}
