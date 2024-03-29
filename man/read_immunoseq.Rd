% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_immunoseq.R
\name{read_immunoseq}
\alias{read_immunoseq}
\title{Read Adaptive Biotechnologies Immunoseq Files}
\usage{
read_immunoseq(.path, .functional = NULL, .format_IMGT = NULL)
}
\arguments{
\item{.path}{A character vector (string). The full or relative (for example,
if the working directory was set with \code{setwd()}) file path to a directory
containing multiple Immunoseq v1 .tsv file exports. Alternatively, the full
or relative file path to a single Immunoseq v1 .tsv file export.}

\item{.functional}{A logical. \code{TRUE} (the default) filters sequences for
in-frame, functional sequences. \code{FALSE} leaves non-coding, out-of-frame
sequences in.}

\item{.format_IMGT}{A logical. \code{TRUE} (the default) indicates that VDJ gene
names will be reformatted according to IMGT gene name conventions. \code{FALSE}
leaves them as is.}
}
\value{
if \code{.path} is a directory containing multiple files, a list of same
length as the directory. if \code{.path} is a single file, a data frame.
}
\description{
A function to read Immunoseq v1 exported .tsv files.
}
\details{
This function reads Adaptive Biotechnologies Immunoseq files into R. Only
necessary columns from the Immunoseq files are read, with the others being
dropped before being read. This was done to preserve RAM. Given a directory
containing multiple files, the function reads them all into a named list.
If only a path to one file is given, the function reads the file to a data
frame. At this time, only Immunoseq v1 .tsv file exports are supported.
}
\seealso{
Other Immunoseq: 
\code{\link{collapse_sequences}()},
\code{\link{diversity_methods}()},
\code{\link{format_immunarch}()},
\code{\link{get_alloreactives}()},
\code{\link{remove_contaminants}()},
\code{\link{repertoire_diversity}()}
}
\author{
Christopher Parks
}
\concept{Immunoseq}
