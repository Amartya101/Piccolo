% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{GenePairsSingleCell}
\alias{GenePairsSingleCell}
\title{GenePairsSingleCell function}
\usage{
GenePairsSingleCell(
  PiccoloList,
  PercSetSize = NULL,
  JcdInd,
  Stop = NULL,
  Out = F
)
}
\arguments{
\item{PiccoloList}{A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function}

\item{PercSetSize}{A numeric variable. Specifies the size of the top percentile set (in percent). Default value is 5.}

\item{JcdInd}{A numeric variable. Specifies the minimum extent of the overlap between percentile sets of any given gene-pair.}

\item{Stop}{A numeric variable. Specifies the stop point for the sliding window. For use as TuBA, leave this NULL (default).}

\item{Out}{A logical variable. If set to T (default is F), will prepare .csv files containing the information about the gene-pairs, cells in percentile sets etc.}
}
\value{
A list object containing the gene-pairs data frame, the genes percentile sets data frame, and a data frame with cell IDs.
}
\description{
This function identifies the gene-pairs with significant overlaps between their percentile sets
}
\examples{
\dontrun{
pbmc3k <- GenePairsSingleCell(PiccoloList = pbmc3k,PercSetSize = 5,JcdInd = 0.3)
}
}
