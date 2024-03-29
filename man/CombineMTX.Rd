% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Piccolo.R
\name{CombineMTX}
\alias{CombineMTX}
\title{CombineMTX}
\usage{
CombineMTX(
  Path,
  MinFeaturesPerCell = 10,
  MT.Perc = 100,
  RP.Perc = 100,
  verbose = T
)
}
\arguments{
\item{Path}{A character variable. Specifies the full path of the directory that contains the folders with the .mtx or .mtx.gz files that need to be combined into one file.}

\item{MinFeaturesPerCell}{An integer variable. The minimum number of features with non-zero counts that any given cell is expected to contain. Cells with fewer than the specified number will be filtered out (Default is 10).}

\item{MT.Perc}{A numeric variable. The percentage of total count in any given cell attributed to counts of mitochondrial genes. Cells with MT.Perc greater than the specified number will be filtered out (Default is 100, so no MT filtering).}

\item{RP.Perc}{A numeric variable. The percentage of total count in any given cell attributed to counts of ribosomal genes. Cells with RP.Perc greater than the specified number will be filtered out (Default is 100, so no RP filtering).}

\item{verbose}{A logical variable. Specifies whether messages generated while running the function should be displayed (default is T)}
}
\value{
Generates a matrix.mtx.gz, a features.tsv, and a barcodes.tsv file in the same directory specified in the Path. The barcodes are modified to contain information about the folders. Eg: "Folder1Name_Barcode1", "Folder1Name_Barcode2",..,"Folder2Name_BarcodeX",...
}
\description{
This function combines the single-cell counts files (.mtx or .mtx.gz) obtained from different samples/studies/conditions etc
}
\examples{
\dontrun{
 CombineMTX(Path = "~/Documents/10X_PBMC_DataSets/") #Default
 CombineMTX(Path = "~/Documents/10X_PBMC_DataSets/",
 MinFeaturesPerCell = 100, MT.Perc = 50,RP.Perc = 90)
}
}
