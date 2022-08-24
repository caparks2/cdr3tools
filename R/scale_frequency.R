#' Define Frequency Limits for Plot Scales
#'
#' Given a vector a frequencies, defines the minimum and maximum limits and
#'   rounds (down for the minimum and up for the maximum).
#'
#' Particularly useful in the context of ggplot2. For example:
#'   `scale_y_continuous(limits = c(scale_frequency_min(x), scale_frequency_max(x))`
#'
#' @param x A numeric vector of frequencies. Must be in the set 0-1.
#' @returns A numeric of length = 1 representing the minimum or maximum, rounded
#'   either up or down.
#' @examples
#' scale_frequency_min(cfselo_seqs$frequency)
#'
#' scale_frequency_max(cfselo_seqs$frequency)
#' @author Christopher Parks
#' @export
scale_frequency_min <- function(x) {
  x <- min({{ x }}, na.rm = TRUE)
  x <- as.character(x)
  if (grepl("e", x)) {
    expnt <- sub("(.*\\.[[:digit:]]+e)(.*)", "\\2", x)
    expnt <- as.numeric(expnt)
    x <- as.numeric(x)
    x <- x * `^`(10, abs(expnt))
    x <- floor(x)
    x * `^`(10, expnt)
  } else {
    zeros <- 1 + nchar(sub("(^0\\.)([0]+)([[:digit:]]+$)", "\\2", x))
    x <- as.numeric(x)
    x <- x * `^`(10, zeros)
    x <- floor(x)
    x * `^`(10, -1 * zeros)
  }
}

#' @rdname scale_frequency_min
#' @export
scale_frequency_max <- function(x) {
  x <- max({{ x }}, na.rm = TRUE)
  x <- as.character(x)
  if (grepl("e", x)) {
    expnt <- sub("(.*\\.[[:digit:]]+e)(.*)", "\\2", x)
    expnt <- as.numeric(expnt)
    x <- as.numeric(x)
    x <- x * `^`(10, abs(expnt))
    x <- ceiling(x)
    x * `^`(10, expnt)
  } else {
    zeros <- 1 + nchar(sub("(^0\\.)([0]+)([[:digit:]]+$)", "\\2", x))
    x <- as.numeric(x)
    x <- x * `^`(10, zeros)
    x <- ceiling(x)
    x * `^`(10, -zeros)
  }
}
