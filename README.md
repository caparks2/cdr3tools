
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cdr3tools

An R package containing helpful functions for the analysis of TCRseq
data, mainly as performed in the Sykes Lab at the Columbia Center for
Translational Immunology, at the Columbia University Medical Center,
NYC, NY, USA.

The alloreactive repertoire analysis and repertoire diversity core
functionality of this package are based on the functions published in
[Obradovic et. al., 2021](https://pubmed.ncbi.nlm.nih.gov/35291378/).
The functions available in this package produce the same results, but
are written to process multiple samples quickly rather than one at a
time.

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install cdr3tools like so:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("caparks2/cdr3tools")

library(cdr3tools)
```

## Example useage

### Reading repertoire data

Read many Adaptive Biotechnologies Immunoseq files into R

``` r
?read_immunoseq()
```

### Repertoire data utilities

Format for use with the [immunarch](https://immunarch.com/) package

``` r
?format_immunarch()
```

Collapse unique sequences together, adding their reads (or template
counts) together, while reducing sequencing/PCR error.

``` r
?collapse_sequences()
```

Remove contaminant sequences from repertoire data

``` r
?remove_contaminants()
```

### Alloreactive TCR tools

Define the unique sequences that are alloreactive

``` r
?get_alloreactives()
```

### Repertoire diversity measures

Calculate repertoire diversity using several different methods for
multiple repertoire files

``` r
?repertoire_diversity()
```

Calculate Jensen Shannon Divergence (or Distance) for a reads (or
template counts) matrix of two repertoires. R vector recycling rules are
different and can confound JSD calculations. This function deals with R
recycling gracefully.

``` r
?jensen_shannon()
```

### IMGT tools

Fetch VDJ gene reference sequences from the IMGT reference sequence
database.

``` r
?imgt_get_ref_seqs()
```

Align CDR3 amino acid sequences according to IMGT unique numbering
rules.

``` r
cdr3_seqs <- c(
   "CASSF",
   "CASSGEKLFF",
   "CASSKPDRGIYGYTF",
   "TGPLHF",
   "CASSQETRYDFLTIDTGGKKKNTEAFF"
)

imgt_align_junctions(cdr3_seqs)
imgt_align_junctions(cdr3_seqs, remove_non_canonicals = TRUE)
```

Simplify IMGT reference sequence FASTA headers

``` r
?imgt_simple_headers()
```

Extract sub sequences from CDR3 sequences using IMGT unique numbering.

``` r
?imgt_extract()
```

Re-format VDJ gene names to conform to the IMGT standard gene names

``` r
?imgt_format_gene_names()
```

### CDR3 sequence characterization

Calculate hydrophobicity based on the [Wimley and
White](https://pubmed.ncbi.nlm.nih.gov/8836100/) hydrophobicity scale.
Ported from the [peptides](https://github.com/dosorio/Peptides) R
package.

``` r
?hydrophobicity
```

Calculate TCR-intrinsic regulatory potential scores based on [Lagattuta
et. al., 2022](https://pubmed.ncbi.nlm.nih.gov/35177831/).

``` r
?get_TiRP_scores()
```

### Miscellaneaous functions

Make pretty scale labels for frequencies in `ggplot2`.

``` r
?scale_frequency()
```
