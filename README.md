
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cdr3tools

<!-- badges: start -->
<!-- badges: end -->

Adjust IMGT JUNCTION sequences for length by inserting gaps according to
the IMGT unique numbering rules. Vectorised for JUNCTION sequences and
their gap-inserted, length-adjusted replacements.

## Installation

You can install cdr3tools like so:

``` r
if(!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools")

devtools::install_github("caparks2/cdr3")
```

Example useage:

``` r
cdr3_seqs <- c(
   "CASSF",
   "CASSGEKLFF",
   "CASSKPDRGIYGYTF",
   "TGPLHF",
   "CASSQETRYDFLTIDTGGKKKNTEAFF"
)

add_imgt_junction_gaps(cdr3_seqs)
add_imgt_junction_gaps(cdr3_seqs, remove_non_canonicals = TRUE)
```
