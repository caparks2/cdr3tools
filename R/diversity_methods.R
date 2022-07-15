#' Diversity Calculation Methods
#'
#' For use with the [cdr3tools::repertoire_diversity()] function.
#'
#' These are the various methods of diversity calculation that can be used to
#'   quantify diversity of T cell receptor sequence repertoires. These methods
#'   must be passed as individual character strings to the `.method` argument in
#'   the [cdr3tools::repertoire_diversity()] function.
#' @export
diversity_methods <- function() {
  c("shannons.clonality", "shannons.entropy", "shannons.diversity",
    "gini.simpson", "simpsons.clonality", "simpsons.index",
    "simpsons.dominance", "simpsons.equitability", "r20", "slope")
}
#' @seealso [cdr3tools::repertoire_diversity()]
#' @family Repertoire Diversity
