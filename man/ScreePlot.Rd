% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{ScreePlot}
\alias{ScreePlot}
\title{ScreePlot function}
\usage{
ScreePlot(PiccoloList, MaxPC = NULL, cex = 1.25)
}
\arguments{
\item{PiccoloList}{A PiccoloList object after using the \link[Piccolo]{ComputePC} function.}

\item{MaxPC}{An integer variable. Specifies the maximum number of PCs to use for the Scree plot.}

\item{cex}{A numeric variable. Specifies the size of the elements in the plot.}
}
\value{
A Scree plot showing the percentage of the overall variance explained by each PC.
}
\description{
This function creates the Scree plot to show the percentage of total variation explained by each principal component (PC).
}
\examples{
\dontrun{
ScreePlot(PiccoloList = pbmc3k)
ScreePlot(PiccoloList = pbmc3k, MaxPC = 15)
}
}
