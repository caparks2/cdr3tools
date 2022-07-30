#' Compute the hydrophobicity index of a protein sequence
#'
#' This function calculates the GRAVY hydrophobicity index of an amino acids
#'   sequence using one of the 38 scales from different sources.
#'
#' The hydrophobicity is an important stabilization force in protein folding;
#'   this force changes depending on the solvent in which the protein is found.
#'   The hydrophobicity index is calculated adding the hydrophobicity of
#'   individual amino acids and dividing this value by the length of the
#'   sequence.
#'
#' @param seq An amino-acids sequence
#' @param scale A character string specifying the hydophobicity scale to be
#'   used; must be one of \code{"Aboderin"}, \code{"AbrahamLeo"}, \code{"Argos"},
#'   \code{"BlackMould"}, \code{"BullBreese"}, \code{"Casari"}, \code{"Chothia"},
#'   \code{"Cid"}, \code{"Cowan3.4"}, \code{"Cowan7.5"}, \code{"Eisenberg"},
#'   \code{"Engelman"}, \code{"Fasman"}, \code{"Fauchere"}, \code{"Goldsack"},
#'   \code{"Guy"}, \code{"HoppWoods"}, \code{"Janin"}, \code{"Jones"},
#'   \code{"Juretic"}, \code{"Kidera"}, \code{"Kuhn"}, \code{"KyteDoolittle"},
#'   \code{"Levitt"}, \code{"Manavalan"}, \code{"Miyazawa"}, \code{"Parker"},
#'   \code{"Ponnuswamy"}, \code{"Prabhakaran"}, \code{"Rao"}, \code{"Rose"},
#'   \code{"Roseman"}, \code{"Sweet"}, \code{"Tanford"}, \code{"Welling"},
#'   \code{"Wilson"}, \code{"Wolfenden"}, \code{"Zimmerman"},
#'   \code{"interfaceScale_pH8"}, \code{"interfaceScale_pH2"},
#'   \code{"octanolScale_pH8"}, \code{"octanolScale_pH2"}, \code{"oiScale_pH8"},
#'   \code{"oiScale_pH2"}, or \code{"Wimley"}.
#' @returns The computed GRAVY index for a given amino-acid sequence
#' @examples
#' hydrophobicity(seq = "QWGRRCCGWGPGRRYCVRWC", scale = "Wimley")
#' @source \url{https://github.com/dosorio/Peptides}
#' @export
hydrophobicity <- function(seq, scale = "Wimley") {
  # Loading hydrophobicity scales
  Hydrophobicity <- cdr3tools::AAdata[["Hydrophobicity"]]
  # Setting the hydrophobicity scale
  scale <- match.arg(scale, names(Hydrophobicity))
  # Split sequence by aminoacids
  seq <- aaCheck(seq)
  # Sum the hydrophobicity of each amino acid and divide them between the sequence length
  # Return the GRAVY value
  h <-
    lapply(seq, function(seq) {
      (sum(Hydrophobicity[[scale]][seq], na.rm = TRUE) / length(seq))
    })
  return(unlist(h))
}

aaCheck <- function(seq){
  if(!any(lengths(seq) > 1)){
    seq <- toupper(seq)
    seq <- gsub(pattern = "[[:space:]]+",replacement = "",x = seq)
    seq <- strsplit(x = seq,split = "")
  } else {
    seq <- lapply(seq,toupper)
  }
  check <- unlist(lapply(seq,function(sequence){
    !all(sequence%in%c("A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N" ,"P" ,"Q" ,"R" ,"S" ,"T" ,"V" ,"W" ,"Y", "-"))
  }))
  if(sum(check) > 0){
    sapply(which(check == TRUE),function(sequence){warning(paste0("Sequence ",sequence," has unrecognized amino acid types. Output value might be wrong calculated"),call. = FALSE)})
  }
  return(seq)
}
