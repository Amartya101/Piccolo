% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{SelectFeatures}
\alias{SelectFeatures}
\title{SelectFeatures function}
\usage{
SelectFeatures(
  PiccoloList,
  NoOfHVG = NULL,
  Batch = NULL,
  MinPercNonZeroCells = 0.5,
  Reference = NULL,
  Out = F,
  verbose = T
)
}
\arguments{
\item{PiccoloList}{A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function and processed through the \link[Piccolo]{FilterCells}}

\item{NoOfHVG}{A numeric (integer) variable. The number of variable features to be shortlisted. If left NULL (default), will shortlist based on the threshold of ReferenceLevel.}

\item{Batch}{An optional character vector. Specifies the batch labels for the cells. The order of the batch labels should match the order of the cells in the barcodes file.}

\item{MinPercNonZeroCells}{A numeric variable. Specifies the minimum percentage of cells that must have non-zero counts for each gene in the data set. Default is 0.5 (\%).}

\item{Reference}{A numeric variable (value should be greater than 0 but less than 1). Specifies the reference level against which features are identified as variable. Default is the 10th quantile for each bin (Reference = 0.1).}

\item{Out}{A logical variable. Specifies whether to return output files (.csv) with the HVGs and the stable genes (if set to T), or not (if set to F). Default is FALSE.}

\item{verbose}{A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).}
}
\value{
An updated PiccoloList object containing data frames with the HVGs and the stable genes.
}
\description{
This function performs feature selection by identifying highly variable genes and stable genes.
}
\examples{
\dontrun{
pbmc3k <- SelectFeatures(PiccoloList = pbmc3k)
pbmc3k <- SelectFeatures(PiccoloList = pbmc3k,
NoOfHVG = 3000, Out = T)
}
}
