#' Sample Repertoires
#'
#' Sample repertoires to common template, read, or copy counts.
#' @name sample_repertoires
NULL
# downsample <- function(.data, .col = NULL, .n = NULL, .counts = TRUE) {
#   templates <- data[[col]]
#   indices <- seq_along(templates)
#   # templates <- data$templates
#   expanded <- rep.int(indices, templates)
#   samples <- sample(expanded, n)
#   collapsed <- as.data.frame(table(samples))
#   sampled_indices <- as.numeric(collapsed[[1]])
#   sampled_templates <- as.numeric(collapsed[[2]])
#   # result <- data[sampled_indices, 'rearrangement']
#   # result['templates'] <- sampled_templates
#   result <- data.frame(Clones = sampled_templates)
#   result['Proportion'] <- result$Clones/sum(result$Clones)
#   result <- result[order(result$Clones, decreasing = TRUE), ]
#   tibble::as_tibble(result)
# }
#
# repertoire_downsample <- function(.data, .col = NULL, .n = NULL, .counts = TRUE) {
#
#   if (!(inherits(.data, "list") || inherits(.data, "data.frame")))
#     rlang::abort(
#       paste(
#         ".data must be a data frame as returned by cdr3tools::read_immunoseq",
#         "cdr3tools::format_immunarch."
#       )
#     )
#
#   if (inherits(.data, "list")) {
#     result <- lapply(.data, function(.data) {
#       .data <- as.data.frame(.data)
#       result <- repertoire_downsample_internal(.data, .col, .n, .counts)
#       tibble::as_tibble(result)
#       }
#     )
#     return(result)
#   }
#
#   if (inherits(.data, "data.frame")) {
#     result <- repertoire_downsample_internal(.data, .col, .n, .counts)
#     result <- tibble::as_tibble(result)
#     return(result)
#   }
# }
#
# repertoire_downsample_internal <- function(.data, .col, .n, .counts) {
#
#   .data <- as.data.frame(.data)
#
#   if (nrow(.data) == 0)
#     .data <- .data[NA, ]
#
#   .
#
#   if (!any(grepl("templates|reads|seq_reads|copies|clones|counts|frequency|proportion", names(.data),
#                  ignore.case = TRUE), na.rm = TRUE)) {
#     rlang::abort(
#       paste0(
#         ".data must be a single data frame containing a column of sequence ",
#         "counts or freuqncies (named `templates`, `frequency`, `seq_reads`, `reads`, `copies`,",
#         "`clones`, `proportion`, or `counts`)",
#         ", not case sensitive, or a list containing multiple of such data frames."
#       )
#     )
#   }
#
#   if (is.null(.col)) {
#     .col <- grep("proportion|frequency", names(.data), ignore.case = TRUE, value = TRUE)
#   }
#   if (!any(grepl("templates|reads|seq_reads|copies|clones|counts", .col,
#                  ignore.case = TRUE))) {
#     rlang::abort('`.col` must be one of "templates", "reads", "seq_reads", "copies", "clones", or "counts". Not case sensitive.')
#   }
#   counts_col <- col
#   if (length(counts_col) > 1) {
#     rlang::abort("There must be only one column containing counts.")
#   }
#
#   .data <- .data[[counts_col]]
#
#   if (length(.data) == 0) {
#     .data <- NA_real_
#   }
#
# }
