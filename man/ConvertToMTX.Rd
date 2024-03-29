% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{ConvertToMTX}
\alias{ConvertToMTX}
\title{ConvertToMTX}
\usage{
ConvertToMTX(X)
}
\arguments{
\item{X}{A character variable. Specifies the name of the file that contains the raw counts data (should be in .csv, .txt, or .tsv format). Make sure it is in the genes along rows and cells along columns format.}
}
\value{
Generates an .mtx.gz file (and features and barcodes files .tsv files) in the working directory that contains the input file.
}
\description{
This function converts .csv or .txt counts files to .mtx format. Prepares the corresponding features.tsv and barcodes.tsv files as well.
}
\examples{
\dontrun{
ConvertToMTX(X = "10X_PBMC3k.csv")
}
}
