% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{LabeltSNE}
\alias{LabeltSNE}
\title{tSNE Plot With Cell Labels}
\usage{
LabeltSNE(
  PiccoloList,
  Labels,
  Levels = NULL,
  Alpha = 0.7,
  Size = 1.4,
  BaseSize = 24,
  Title = "Piccolo",
  LegendPosition = "right"
)
}
\arguments{
\item{PiccoloList}{A list object. Piccolo list object obtained after applying the \link[Piccolo]{tSNEcoords} function.}

\item{Labels}{A character variable. Specifies the cell-type/group/batch labels for the cells in the same order that they are listed in the barcodes file.}

\item{Levels}{A character variable. Specifies the order in which the cluster labels should be listed in the tSNE plot.}

\item{Alpha}{A numeric variable. Specifies the transparency of the dots in the tSNE plot. Smaller values lead to greater transparency. Default is 0.7.}

\item{Size}{A numeric variable. Specifies the size of the dots in the tSNE plot. Default is 1.4.}

\item{BaseSize}{A numeric variable. Specifies the base size of the text elements in the tSNE plot. Default is 28.}

\item{Title}{A character variable. Specifies the title of the tSNE plot.}

\item{LegendPosition}{A character variable. Specifies the position in the plot where the legend should be placed. Default is "right".}
}
\value{
A tSNE plot with the cells colored according to the cell labels provided by the user.
}
\description{
This function generates the tSNE plot with the cells labeled by the cluster labels identitied by Leiden clustering.
}
\examples{
\dontrun{
LabeltSNE(PiccoloList = pbmc3k)
LabeltSNE(PiccoloList = pbmc3k, Title = "PBMC3k")
}
}
