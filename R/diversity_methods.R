#' Diversity Calculation Methods
#'
#' For use with the [cdr3tools::repertoire_diversity()] function.
#'
#' These are the various methods of diversity calculation that can be used to
#'   quantify diversity of T cell receptor sequence repertoires. These methods
#'   must be passed as individual character strings to the `.method` argument in
#'   the [cdr3tools::repertoire_diversity()] function. Available options to
#'   `.method` are :
#'     * \code{"shannons.clonality"}
#'     * \code{"shannons.entropy"}
#'     * \code{"shannons.diversity"}
#'     * \code{"gini.simpson"}
#'     * \code{"simpsons.clonality"}
#'     * \code{"simpsons.index"}
#'     * \code{"simpsons.dominance"}
#'     * \code{"simpsons.equitability"}
#'     * \code{"r20"}
#'     * \code{"slope"}
#' @seealso [cdr3tools::repertoire_diversity()]
#' @family Repertoire Diversity
#' @export
diversity_methods <- function() {
  c("shannons.clonality", "shannons.entropy", "shannons.diversity",
    "gini.simpson", "simpsons.clonality", "simpsons.index",
    "simpsons.dominance", "simpsons.equitability", "r20", "slope")
}
