% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{UMAPcoords}
\alias{UMAPcoords}
\title{UMAP coordinates}
\usage{
UMAPcoords(PiccoloList, NoOfPC = NULL, Out = F)
}
\arguments{
\item{PiccoloList}{A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} and the \link[Piccolo]{ComputePC} functions.}

\item{NoOfPC}{An integer variable. Specifies the number of PCs to use for the UMAP. The default is NULL, which corresponds to the use of all the shortlisted principal components.}

\item{Out}{A logical variable. Specifies whether to return an output file (.csv) with the UJMAP coordinates (if set to T), or not (if set to F). Default is F.}
}
\value{
An updated PiccoloList with a data frame containing the UMAP coordinates of the cells.
}
\description{
This function generates the 2-dimensional coordinates for the cells using UMAP.
}
\examples{
\dontrun{
pbmc3k <- UMAPcoords(PiccoloList = pbmc3k,
Out = T)
}
}
