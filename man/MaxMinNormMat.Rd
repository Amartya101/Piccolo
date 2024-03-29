% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{MaxMinNormMat}
\alias{MaxMinNormMat}
\title{Max-Min Normalization Function}
\usage{
MaxMinNormMat(PiccoloList, Out = F)
}
\arguments{
\item{PiccoloList}{A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function}

\item{Out}{A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T). Default is F.}
}
\value{
A numeric matrix containing the normalized values.
}
\description{
This function will normalize the standardized values to normalized values in the range 0-1
}
\examples{
\dontrun{
NormMat <- MaxMinNormMat(PiccoloList = pbmc3k,
Out = F)
}
}
