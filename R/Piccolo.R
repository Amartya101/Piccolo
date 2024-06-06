
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
utils::globalVariables(c("."))
utils::globalVariables(c("UMAP1","UMAP2","UMAP_1", "UMAP_2", "Label","UMAP 1","UMAP 2","tSNE1","tSNE2","tSNE 1","tSNE 2","Index","Z Score"))
#' @import utils
#' @import stats
#' @import ggplot2
#' @import network
#' @import RSpectra
#' @import cluster
#' @importFrom ggplot2 aes
#' @importFrom methods is
#' @importFrom grDevices hcl
#' @import ggnetwork

#' @title  ConvertToMTX
#' @description  This function converts .csv or .txt counts files to .mtx format. Prepares the corresponding features.tsv and barcodes.tsv files as well.
#' @export
#' @param X A character variable. Specifies the name of the file that contains the raw counts data (should be in .csv, .txt, or .tsv format). Make sure it is in the genes along rows and cells along columns format.
#' @return Generates an .mtx.gz file (and features and barcodes files .tsv files) in the working directory that contains the input file.
#' @examples
#' \dontrun{
#' ConvertToMTX(X = "10X_PBMC3k.csv")
#' }

ConvertToMTX <- function(X){

  #Read in UMI counts csv or txt or tsv file
  UMI.Mat <- data.table::fread(X)
  colnames(UMI.Mat) <- c("Gene.ID",colnames(UMI.Mat)[-1])

  Gene.ID <- UMI.Mat$Gene.ID

  UMI.Mat <- as.matrix(UMI.Mat[,-1])

  rownames(UMI.Mat) <- Gene.ID

  #create sparse matrix
  sparse.Mat <- Matrix::Matrix(UMI.Mat,sparse = T)

  FileNameMTX <- paste0(substr(X,1,nchar(X) - 4),"_matrix.mtx.gz")

  #writeMMgz function by Kamil Slowikowski
  writeMMgz <- function(x,file){
    mtype <- "real"
    if (is(x,"dgCMatrix")){
      mtype <- "integer"
    }
    writeLines(c(sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))),gzfile(file))
    data.table::fwrite(x = Matrix::summary(x),file = file,append = T,sep = " ",row.names = F,col.names = F)
  }

  #write mtx file
  writeMMgz(x = sparse.Mat, file = FileNameMTX)

  #write feature names file

  FileNameFeatures <- paste0(substr(X,1,nchar(X)-4),"_features.tsv")
  write(x = rownames(UMI.Mat), file = FileNameFeatures)

  #write barcodes file
  FileNameBarcodes <- paste0(substr(X,1,nchar(X)-4),"_barcodes.tsv")
  write(x = colnames(UMI.Mat), file = FileNameBarcodes)

  closeAllConnections()

}

#' @title  CombineMTX
#' @description  This function combines the single-cell counts files (.mtx or .mtx.gz) obtained from different samples/studies/conditions etc
#' @export
#' @param Path A character variable. Specifies the full path of the directory that contains the folders with the .mtx or .mtx.gz files that need to be combined into one file.
#' @param MinFeaturesPerCell An integer variable. The minimum number of features with non-zero counts that any given cell is expected to contain. Cells with fewer than the specified number will be filtered out (Default is 10).
#' @param MT.Perc A numeric variable. The percentage of total count in any given cell attributed to counts of mitochondrial genes. Cells with MT.Perc greater than the specified number will be filtered out (Default is 100, so no MT filtering).
#' @param RP.Perc A numeric variable. The percentage of total count in any given cell attributed to counts of ribosomal genes. Cells with RP.Perc greater than the specified number will be filtered out (Default is 100, so no RP filtering).
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T)
#' @return Generates a matrix.mtx.gz, a features.tsv, and a barcodes.tsv file in the same directory specified in the Path. The barcodes are modified to contain information about the folders. Eg: "Folder1Name_Barcode1", "Folder1Name_Barcode2",..,"Folder2Name_BarcodeX",...
#' @examples
#' \dontrun{
#'  CombineMTX(Path = "~/Documents/10X_PBMC_DataSets/") #Default
#'  CombineMTX(Path = "~/Documents/10X_PBMC_DataSets/",
#'  MinFeaturesPerCell = 100, MT.Perc = 50,RP.Perc = 90)
#' }
CombineMTX <- function (Path, MinFeaturesPerCell = 10, MT.Perc = 100, RP.Perc = 100,verbose = T)
{
  setwd(Path)
  FolderNames <- list.dirs()[-1]
  mat1 <- Matrix::Matrix(0, nrow = 2, ncol = 2)
  colnames(mat1) <- c("A", "B")
  for (i in 1:length(FolderNames)) {
    Folder <- gsub("./", "", FolderNames[i])
    Files <- list.files(FolderNames[i])
    mtxcheck <- unique(grep(".mtx", x = Files, fixed = T),grep(".mtx.gz", x = Files, fixed = T))
    featuresbarcodescheck <- unique(grep(".tsv", x = Files,fixed = T), grep(".tsv.gz", x = Files, fixed = T))
    if (length(mtxcheck) == 1 & length(featuresbarcodescheck) == 2) {
      barcodes.file <- Files[grep("barcodes", Files, fixed = T)]
      barcodes.path <- paste0(FolderNames[i], "/", barcodes.file)
      features.file <- Files[grep("features", Files, fixed = T)]
      if (length(features.file) == 0) {
        features.file <- Files[grep("genes", Files, fixed = T)]
      }
      features.path <- paste0(FolderNames[i], "/", features.file)
      matrix.file <- Files[grep("matrix", Files, fixed = T)]
      matrix.path <- paste0(FolderNames[i], "/", matrix.file)
      message(paste0("Importing file ", i, "..."))
      mat <- Matrix::readMM(file = matrix.path)
      if (i == 1) {
        mat1 <- Matrix::Matrix(0, nrow = nrow(mat), ncol = 2)
        colnames(mat1) <- c("A", "B")
      }
      feature.names <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
      barcode.names <- read.delim(barcodes.path, header = FALSE, stringsAsFactors = FALSE)
      if (ncol(feature.names) > 1) {
        List.For.GS.Col <- vector(mode = "list", length = ncol(feature.names))
        for (j in 1:ncol(feature.names)) {
          List.For.GS.Col[[j]] <- which(substr(toupper(feature.names[,j]), 1, 3) == "RPL" | substr(toupper(feature.names[,
                                                                                                                        j]), 1, 3) == "RPS")
        }
        GS.Col <- which(unlist(lapply(List.For.GS.Col,length)) != 0)
        if (length(GS.Col) == 1) {
          Duplicated.Features <- which(duplicated(feature.names[,1]) == T)
          if (length(Duplicated.Features) != 0) {
            mat <- Matrix::t(mat)
            mat <- mat[, -Duplicated.Features]
            mat <- Matrix::t(mat)
            feature.names <- feature.names[-Duplicated.Features,]
            rownames(mat) <- feature.names[,1]
            if (i == 1) {
              mat1 <- Matrix::t(mat1)
              mat1 <- mat1[, -Duplicated.Features]
              mat1 <- Matrix::t(mat1)
            }
          }
          else {
            rownames(mat) <- feature.names[,1]
          }
        }
        else {
          warning("Features file does not contain gene symbols. MT genes and RP genes filtering will not be performed.")
          Duplicated.Features <- which(duplicated(feature.names[,1]) == T)
          if (length(Duplicated.Features) != 0) {
            mat <- Matrix::t(mat)
            mat <- mat[, -Duplicated.Features]
            mat <- Matrix::t(mat)
            feature.names <- feature.names[-Duplicated.Features,]
            rownames(mat) <- feature.names[,1]
            if (i == 1) {
              mat1 <- Matrix::t(mat1)
              mat1 <- mat1[,-Duplicated.Features]
              mat1 <- Matrix::t(mat1)
            }
          }
          else {
            rownames(mat) <- feature.names[,1]
          }
        }
      }
      else {
        Duplicated.Features <- which(duplicated(feature.names[,1]) == T)
        if (length(Duplicated.Features) != 0) {
          mat <- Matrix::t(mat)
          mat <- mat[, -Duplicated.Features]
          mat <- Matrix::t(mat)
          feature.names <- feature.names[-Duplicated.Features,]
          rownames(mat) <- feature.names[,1]
          if (i == 1) {
            mat1 <- Matrix::t(mat1)
            mat1 <- mat1[, -Duplicated.Features]
            mat1 <- Matrix::t(mat1)
            feature.names <- feature.names[-Duplicated.Features,]
            rownames(mat1) <- feature.names[,1]
          }
        }
        else {
          rownames(mat) <- feature.names[,1]
        }
      }
      if (nrow(mat1) != nrow(mat) & i != 1) {
        if (nrow(mat1) > nrow(mat)) {
          Matching.Features <- intersect(toupper(rownames(mat1)),toupper(rownames(mat)))
          Matching.Features.Mat1 <- which(toupper(rownames(mat1)) %in% Matching.Features)
          Matching.Features.Mat1 <- Matching.Features.Mat1[!is.na(Matching.Features.Mat1)]
          mat1 <- Matrix::t(mat1)
          mat1 <- mat1[, Matching.Features.Mat1]
          mat1 <- Matrix::t(mat1)
          Matching.Features.Mat <- which(toupper(rownames(mat)) %in% Matching.Features)
          Matching.Features.Mat <- Matching.Features.Mat[!is.na(Matching.Features.Mat)]
          mat <- Matrix::t(mat)
          mat <- mat[, Matching.Features.Mat]
          mat <- Matrix::t(mat)
        }
        else if (nrow(mat1) < nrow(mat)) {
          Matching.Features <- intersect(toupper(rownames(mat)), toupper(rownames(mat1)))
          Matching.Features.Mat <- which(toupper(rownames(mat)) %in% Matching.Features)
          Matching.Features.Mat <- Matching.Features.Mat[!is.na(Matching.Features.Mat)]
          mat <- Matrix::t(mat)
          mat <- mat[, Matching.Features.Mat]
          mat <- Matrix::t(mat)
          Matching.Features.Mat1 <- which(toupper(rownames(mat1)) %in% Matching.Features)
          Matching.Features.Mat1 <- Matching.Features.Mat1[!is.na(Matching.Features.Mat1)]
          mat1 <- Matrix::t(mat1)
          mat1 <- mat1[, Matching.Features.Mat1]
          mat1 <- Matrix::t(mat1)
        }
      }
      colnames(mat) <- paste0(Folder, "_", barcode.names[,1])
      rownames(mat1) <- rownames(mat)
      mat1 <- cbind(mat1, mat)
    }
  }
  feature.names <- feature.names[feature.names[,1] %in% rownames(mat1),]
  mat1 <- mat1[, -c(1, 2)]
  Col.Sums.Vec <- Matrix::colSums(mat1)
  FeatureCounts.Per.Cell <- Matrix::diff(mat1@p)
  if (verbose == T){
    message("Filtering...")
  }

  if (is.character(feature.names) != T) {
    MT.Features <- grep("MT-", toupper(feature.names[, GS.Col]),fixed = T)
  } else {
    MT.Features <- grep("MT-", toupper(feature.names),fixed = T)
  }

  if (length(MT.Features) != 0) {
    MT.mat <- mat1[MT.Features,]
    MT.Col.Sums <- Matrix::colSums(MT.mat)
    MT.In.Prop.Total.Sum <- MT.Col.Sums/Col.Sums.Vec
  } else {
    if (verbose == T){
      message("MT genes not detected in features.")
    }
    MT.In.Prop.Total.Sum <- c()
  }

  if (is.character(feature.names) != T) {
    RP.Features <- which(substr(toupper(feature.names[, GS.Col]), 1, 3) == "RPL" | substr(toupper(feature.names[, GS.Col]), 1, 3) == "RPS" | substr(toupper(feature.names[, GS.Col]),1, 3) %in% c("FAU", "UBA52"))
  } else {
    RP.Features <- which(substr(toupper(feature.names), 1, 3) == "RPL" | substr(toupper(feature.names), 1, 3) == "RPS" | substr(toupper(feature.names),1, 3) %in% c("FAU", "UBA52"))
  }

  if (length(RP.Features) != 0) {
    RP.mat <- mat1[RP.Features,]
    RP.Col.Sums <- Matrix::colSums(RP.mat)
    RP.In.Prop.Total.Sum <- RP.Col.Sums/Col.Sums.Vec
  }
  else {
    if (verbose == T){
      message("RP genes not detected in features.")
    }

    RP.In.Prop.Total.Sum <- c()
  }
  Cells.To.Remove <- unique(c(which(MT.In.Prop.Total.Sum >
                                      MT.Perc/100), which(RP.In.Prop.Total.Sum > RP.Perc/100),
                              which(FeatureCounts.Per.Cell < MinFeaturesPerCell)))
  if (length(Cells.To.Remove) != 0) {
    mat1 <- mat1[, -Cells.To.Remove]
  }
  writeMMgz <- function(x, file) {
    mtype <- "real"
    if (methods::is(x, "dgCMatrix")) {
      mtype <- "integer"
    }
    writeLines(c(sprintf("%%%%MatrixMarket matrix coordinate %s general",
                         mtype), sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))),
               gzfile(file))
    data.table::fwrite(x = Matrix::summary(x), file = file,
                       append = TRUE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  message(paste0("There are ", nrow(mat1), " features and ",
                 ncol(mat1)), " cells in the filtered matrix.")
  message("Writing the mtx and tsv files...")
  suppressWarnings(writeMMgz(x = mat1, file = "matrix.mtx.gz"))
  write.table(x = feature.names, file = "features.tsv", sep = "\t",
              row.names = F, col.names = F, quote = F)
  write(x = colnames(mat1), file = "barcodes.tsv")
  closeAllConnections()
  message("Successfully generated the output files.")
}


#' @title  CreatePiccoloList Function
#' @description  This function creates a list object containing the counts matrix, the features (genes) list, and the barcodes.
#' @export
#' @param MTX A character variable. Specifies the name of the .mtx or .mtx.gz file that contains the counts.
#' @param Genes A character variable. Specifies the name of the features (genes) file (.tsv format).
#' @param Barcodes A character variable. Specifies the name of the barcodes file (.tsv format).
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @return A list containing the counts matrix, the genes list, and the barcodes.
#' @examples
#' \dontrun{
#'  pbmc3k <- CreatePiccoloList(X = "10X_PBMC3k_matrix.mtx.gz",
#'  Gene = "10X_PBMC3k_features.tsv",
#'  Barcode = "10X_PBMC3k_barcodes.tsv")
#'  pbmc3k <- CreatePiccoloList(X = "10X_PBMC3k_matrix.mtx.gz",
#'  Gene = "10X_PBMC3k_features.tsv",
#'  Barcode = "10X_PBMC3k_barcodes.tsv"
#'  MinFeaturesPerCell = 100, MT.Perc = 20,RP.Perc = 90) #changing the filtering criteria
#' }
CreatePiccoloList <- function(MTX, Genes, Barcodes,verbose = T)
{
  if (verbose == T){
    message("Importing files...")
  }
  UMI.Mat <- Matrix::readMM(file = MTX)
  UMI.Mat <- methods::as(UMI.Mat, "CsparseMatrix")
  Gene.IDs <- read.delim(Genes, header = F, stringsAsFactors = F)
  if (ncol(Gene.IDs) == 1){
    colnames(Gene.IDs) <- c("V1")
    Gene.IDs <- as.vector(Gene.IDs$V1)
  }
  Barcodes <- read.delim(Barcodes, header = F, stringsAsFactors = F)
  Barcodes <- Barcodes$V1

  PiccoloList <- list(CountsOriginal = UMI.Mat,GenesOriginal = Gene.IDs,BarcodesOriginal = Barcodes,Counts = UMI.Mat,Genes = Gene.IDs,Barcodes = Barcodes)
  if (verbose == T){
    message("Successfully imported.")
  }
  return(PiccoloList)
}

#' @title  FilterCells Function
#' @description  This function is used for filtering cells based on criterias like the minimum number of unique features per cell, the contribution of the mitochondrial and the ribosomal genes to the total counts, as well as identification of outliers based on the total counts.
#' @export
#' @param PiccoloList A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function.
#' @param MinFeaturesPerCell An integer variable. Specifies the minimum number of genes with non-zero counts within each cell (default is 50).
#' @param MT.Perc A numeric variable. Specifies the maximum percentage of total counts within each cell contributed by mitochondrial genes. Cells with MT.Perc greater than specified threshold will be filtered out (default MT.Perc = 100). Can only be used when gene symbols are available in the features file.
#' @param RP.Perc A numeric variable. Specifies the maximum percentage of total counts within each cell contributed by ribosomal genes. Cells with RP.Perc greater than specified threshold will be filtered out (default RP.Perc = 100). Can only be used when gene symbols are available in the features file.
#' @param TotalCountsMADHigh A numeric variable. Specifies the median absolute deviation (MAD) of the total count above which the cells will be filtered out. Default is NULL, so no filtering based on total count.
#' @param TotalCountsMADLow A numeric variable. Specifies the median absolute deviation (MAD) of the total count above which the cells will be filtered out. Default is NULL, so no filtering based on total count.
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @return A list containing the filtered counts matrix, the genes list, and the barcodes.
#' @examples
#' \dontrun{
#'  pbmc3k <- FilterCells(PiccoloList = pbmc3k)
#'  pbmc3k <- FilterCells(PiccoloList = pbmc3k,
#'  MinFeaturesPerCell = 100, MT.Perc = 50,
#'  RP.Perc = 70,TotalCountsMADHigh = 3.5,
#'  TotalCountsMADLow = 3.5) #changing the filtering criteria
#' }
FilterCells <- function(PiccoloList,MinFeaturesPerCell = 50,MT.Perc = 100,RP.Perc = 100,TotalCountsMADHigh = NULL,TotalCountsMADLow = NULL,verbose = T){
  UMI.Mat <- PiccoloList$CountsOriginal
  Gene.IDs <- PiccoloList$GenesOriginal
  Barcodes <- PiccoloList$BarcodesOriginal
  if (is.vector(Gene.IDs)){
    List.For.GS.Col <- vector(mode = "list", length = length(Gene.IDs))
    for (j in 1:length(Gene.IDs)) {
      List.For.GS.Col[[j]] <- which(substr(toupper(Gene.IDs[j]),1,3) == "RPL" | substr(toupper(Gene.IDs[j]),1,3) == "RPS")
    }
    if (length(List.For.GS.Col) != 0) {
      GS.Col <- 1
    }
    Duplicated.Features <- which(duplicated(Gene.IDs) == T)
    if (length(Duplicated.Features) != 0) {
      UMI.Mat <- Matrix::t(UMI.Mat)
      UMI.Mat <- UMI.Mat[, -Duplicated.Features]
      UMI.Mat <- Matrix::t(UMI.Mat)
      Gene.IDs <- Gene.IDs[-Duplicated.Features]
      rownames(UMI.Mat) <- Gene.IDs
    }
    else {
      rownames(UMI.Mat) <- Gene.IDs
    }
  } else if (ncol(Gene.IDs) > 1) {
    List.For.GS.Col <- vector(mode = "list", length = ncol(Gene.IDs))
    for (j in 1:ncol(Gene.IDs)) {
      List.For.GS.Col[[j]] <- which(substr(toupper(Gene.IDs[,j]),1,3) == "RPL" | substr(toupper(Gene.IDs[,j]),1,3) == "RPS")
    }
    GS.Col <- which(unlist(lapply(List.For.GS.Col,length)) !=  0)
    if (length(GS.Col) == 1) {
      Duplicated.Features <- which(duplicated(Gene.IDs[,1]) == T)
      if (length(Duplicated.Features) != 0) {
        UMI.Mat <- Matrix::t(UMI.Mat)
        UMI.Mat <- UMI.Mat[, -Duplicated.Features]
        UMI.Mat <- Matrix::t(UMI.Mat)
        Gene.IDs <- Gene.IDs[-Duplicated.Features,]
        rownames(UMI.Mat) <- Gene.IDs[, 1]
      } else {
        rownames(UMI.Mat) <- Gene.IDs[, 1]
      }
    } else {
      warning("Features file does not contain gene symbols. MT genes and RP genes filtering will not be performed.")
      Gene.IDs <- Gene.IDs[,1]
      Duplicated.Features <- which(duplicated(Gene.IDs[,1]) == T)
      if (length(Duplicated.Features) != 0) {
        UMI.Mat <- Matrix::t(UMI.Mat)
        UMI.Mat <- UMI.Mat[, -Duplicated.Features]
        UMI.Mat <- Matrix::t(UMI.Mat)
        Gene.IDs <- Gene.IDs[-Duplicated.Features,]
        rownames(UMI.Mat) <- Gene.IDs[,1]
      } else {
        rownames(UMI.Mat) <- Gene.IDs[,1]
      }
    }
  } else {
    Gene.IDs <- Gene.IDs$V1
    List.For.GS.Col <- vector(mode = "list", length = length(Gene.IDs))
    for (j in 1:length(Gene.IDs)) {
      List.For.GS.Col[[j]] <- which(substr(toupper(Gene.IDs[j]),1,3) == "RPL" | substr(toupper(Gene.IDs[j]),1,3) == "RPS")
    }

    if (length(List.For.GS.Col) != 0) {
      GS.Col <- 1
    }

    Duplicated.Features <- which(duplicated(Gene.IDs) ==  T)
    if (length(Duplicated.Features) != 0) {
      UMI.Mat <- Matrix::t(UMI.Mat)
      UMI.Mat <- UMI.Mat[,-Duplicated.Features]
      UMI.Mat <- Matrix::t(UMI.Mat)
      Gene.IDs <- Gene.IDs[-Duplicated.Features]
      rownames(UMI.Mat) <- Gene.IDs
    } else {
      rownames(UMI.Mat) <- Gene.IDs
    }
  }

  Col.Sums.Vec <- Matrix::colSums(UMI.Mat)
  FeatureCounts.Per.Cell <- Matrix::diff(UMI.Mat@p)
  if (verbose == T){
    message("Filtering...")
  }
  if (is.null(ncol(Gene.IDs)) != T) {
    MT.Features <- grep("MT-", toupper(Gene.IDs[,GS.Col]),fixed = T)
  } else {
    MT.Features <- grep("MT-", toupper(Gene.IDs), fixed = T)
  }

  if (length(MT.Features) > 1) {
    MT.mat <- UMI.Mat[MT.Features, ]
    MT.Col.Sums <- Matrix::colSums(MT.mat)
    MT.In.Prop.Total.Sum <- MT.Col.Sums/Col.Sums.Vec
  } else {
    if (verbose == T){
      message("MT genes not detected in features.")
    }
    MT.In.Prop.Total.Sum <- c()
  }

  if (is.null(ncol(Gene.IDs)) != T) {
    RP.Features <- which(substr(toupper(Gene.IDs[,GS.Col]),1, 3) == "RPL" | substr(toupper(Gene.IDs[,GS.Col]),1,3) == "RPS" | substr(toupper(Gene.IDs[,GS.Col]),1,3) %in% c("FAU","UBA52"))
  } else {
    RP.Features <- which(substr(toupper(Gene.IDs),1,3) == "RPL" | substr(toupper(Gene.IDs), 1, 3) == "RPS" | substr(toupper(Gene.IDs),1,3) %in% c("FAU","UBA52"))
  }

  if (length(RP.Features) != 0) {
    RP.mat <- UMI.Mat[RP.Features,]
    RP.Col.Sums <- Matrix::colSums(RP.mat)
    RP.In.Prop.Total.Sum <- RP.Col.Sums/Col.Sums.Vec
  } else {
    if (verbose == T){
      message("RP genes not detected in features.")
    }
    RP.In.Prop.Total.Sum <- c()
  }

  Cells.To.Remove <- unique(c(which(MT.In.Prop.Total.Sum > MT.Perc/100), which(RP.In.Prop.Total.Sum > RP.Perc/100), which(FeatureCounts.Per.Cell < MinFeaturesPerCell)))
  if (length(Cells.To.Remove) != 0) {
    UMI.Mat <- UMI.Mat[,-Cells.To.Remove]
    Barcodes <- Barcodes[-Cells.To.Remove]
  }

  TotalCountsCells <- Matrix::colSums(UMI.Mat)
  if (is.null(TotalCountsMADHigh) != T){
    Cells.To.Remove.High <- which(TotalCountsCells > median(TotalCountsCells) + TotalCountsMADHigh*mad(TotalCountsCells))
  } else {
    Cells.To.Remove.High <- c()
  }

  if (is.null(TotalCountsMADLow) != T){
    Cells.To.Remove.Low <- which(TotalCountsCells < median(TotalCountsCells) - TotalCountsMADLow*mad(TotalCountsCells))
  } else {
    Cells.To.Remove.Low <- c()
  }

  Cells.To.Remove <- c(Cells.To.Remove.Low,Cells.To.Remove.High)

  if (length(Cells.To.Remove) != 0) {
    UMI.Mat <- UMI.Mat[,-Cells.To.Remove]
    Barcodes <- Barcodes[-Cells.To.Remove]
  }

  PiccoloList$Counts <-  UMI.Mat
  PiccoloList$Genes <- Gene.IDs
  PiccoloList$Barcodes = Barcodes

  if (verbose == T){
    message("Done.")
  }
  return(PiccoloList)
}

#' @title  PrepareCountsForSeurat function
#' @description  This function adds gene names and barcodes to the counts matrix in PiccoloList in order to use it as input for the CreateSeuratObject function in Seurat.
#' @export
#' @param PiccoloList A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function and optionally processed through the \link[Piccolo]{FilterCells}
#' @return An updated PiccoloList object containing the counts matrix with the gene names and the barcodes.
#' @examples
#' \dontrun{
#' pbmc3k <- PrepareCountsForSeurat(PiccoloList = pbmc3k)
#' }
PrepareCountsForSeurat <- function(PiccoloList){
  if(is.null(dim(PiccoloList$Genes)) != T & length(PiccoloList$Genes$V2) != 0){
    Genes <- PiccoloList$Genes$V2
    Genes <- make.unique(Genes,sep = "//")
  } else if (length(PiccoloList$Genes$V1) != 0){
    Genes <- PiccoloList$Genes$V1
    Genes <- make.unique(Genes,sep = "//")
  } else {
    Genes <- PiccoloList$Genes
    Genes <- make.unique(Genes,sep = "//")
  }
  rownames(PiccoloList$Counts) <- Genes
  colnames(PiccoloList$Counts) <- PiccoloList$Barcodes

  return(PiccoloList)
}

#' @title  AddMetadataForSeurat function
#' @description  This function can be used to add cell metadata information in the Seurat object.
#' @export
#' @param Obj A Seurat object. This is created using the CreateSeuratObject in Seurat.
#' @param FieldName A character input. Specifies the name of the field to be added to the metadata
#' @param Values A numeric or character variable. Specifies the labels for the cells. The labels should be in the order as the cells in the data set.
#' @return An updated Seurat object containing the metadata with the new field.
#' @examples
#' \dontrun{
#' Field <- "CellTypes"
#' CellLabels <- c(rep("B-cells",900),rep("T-cells",900),rep("Monocytes",900))
#' pbmc3k <- AddMetadataForSeurat(PiccoloList = pbmc3k,FieldName = Field,Values = CellLabels)
#' }
AddMetadataForSeurat <- function(Obj,FieldName,Values){
  MetaData.df <- Obj@meta.data
  New.MetaData.df <- data.frame(MetaData.df,X1 = Values)
  colnames(New.MetaData.df) <- c(colnames(MetaData.df),FieldName)
  Obj@meta.data <- New.MetaData.df
  return(Obj)
}

#' @title  SelectFeatures function
#' @description  This function performs feature selection by identifying highly variable genes and stable genes.
#' @export
#' @param PiccoloList A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function and processed through the \link[Piccolo]{FilterCells}
#' @param NoOfHVG A numeric (integer) variable. The number of variable features to be shortlisted. If left NULL (default), will shortlist based on the threshold of ReferenceLevel.
#' @param Batch An optional character vector. Specifies the batch labels for the cells. The order of the batch labels should match the order of the cells in the barcodes file.
#' @param MinPercNonZeroCells A numeric variable. Specifies the minimum percentage of cells that must have non-zero counts for each gene in the data set. Default is 0.5 (\%).
#' @param Reference A numeric variable (value should be greater than 0 but less than 1). Specifies the reference level against which features are identified as variable. Default is the 10th quantile for each bin (Reference = 0.1).
#' @param Out A logical variable. Specifies whether to return output files (.csv) with the HVGs and the stable genes (if set to T), or not (if set to F). Default is FALSE.
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @return An updated PiccoloList object containing data frames with the HVGs and the stable genes.
#' @examples
#' \dontrun{
#' pbmc3k <- SelectFeatures(PiccoloList = pbmc3k)
#' pbmc3k <- SelectFeatures(PiccoloList = pbmc3k,
#' NoOfHVG = 3000, Out = T)
#' }
SelectFeatures <- function(PiccoloList,NoOfHVG = NULL,Batch = NULL,MinPercNonZeroCells = 0.5,Reference = NULL,Out = F,verbose = T){

  if (is.null(Reference)){
    Reference <- 0.1
  } else if (Reference > 0.5){
    message("Reference should not be greater than 0.5. Resetting it to 0.1 (default)")
    Reference <- 0.1
  } else if (Reference <= 0){
    message("Reference has to be greater than 0. Resetting it to 0.1 (default)")
    Reference <- 0.1
  }

  if (is.null(Batch)){
    UMI.Mat <- Matrix::t(PiccoloList$Counts)

    Gene.IDs <- PiccoloList$Genes

    if(verbose == T){
      message("Filtering features...")
    }

    colVarsSPM <- function(X) {
      stopifnot( methods::is(X, "CsparseMatrix"))
      ans <- sapply(base::seq.int(X@Dim[2]),function(j) {
        if(X@p[j+1] == X@p[j]) { return(0) } # all entries are 0: var is 0
        mean <- base::sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]
        sum((X@x[(X@p[j]+1):X@p[j+1] ] - mean)^2) +
          mean^2 * (X@Dim[1] - (X@p[j+1] - X@p[j]))})/(X@Dim[1] - 1)
      names(ans) <- X@Dimnames[[2]]
      ans
    }

    colOverdispQPCoef <- function(X,alternative = "greater"){
      stopifnot( methods::is(X,"CsparseMatrix"))
      ans <- sapply( base::seq.int(X@Dim[2]),function(j){
        if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
        #mean <- exp(sum(log(X@x[(X@p[j]+1):X@p[j+1]]+1))/X@Dim[1]) - 1
        est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]

        aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                    X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean

        mean(aux) + 1})
    }

    Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
    Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)

    Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)

    No.of.Non.Zero.Per.Feature <- diff(UMI.Mat@p)

    Perc.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature/nrow(UMI.Mat) * 100

    Irrelevant.Features <- unique(c(Irrelevant.Features,which(Perc.Non.Zero.Per.Feature <= MinPercNonZeroCells)))

    if (length(Irrelevant.Features) != 0){
      UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
      if (is.null(ncol(Gene.IDs)) != T){
        Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
      } else {
        Gene.IDs <- Gene.IDs[-Irrelevant.Features]
      }
      Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
      Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
      No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
      Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
    }

    if (verbose == T){
      message("Estimating dispersion coefficients...")
    }

    Alpha.QP <- colOverdispQPCoef(UMI.Mat)

    Irrelevant.Features <- which(Alpha.QP <= 1)

    if (length(Irrelevant.Features) != 0){
      UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
      if (is.null(ncol(Gene.IDs)) != T){
        Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
      } else {
        Gene.IDs <- Gene.IDs[-Irrelevant.Features]
      }
      Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
      Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
      Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
      No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
      Alpha.QP <- Alpha.QP[-Irrelevant.Features]
    }

    Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature

    #Binning based approach

    if (verbose == T){
      message("Shortlisting variable features...")
    }

    Mean.Quantiles <- quantile(Mean.Arith.Per.Feature,probs = seq(0.001,1,0.001))
    Diff.AlphaQP.AlphaQPFit <- vector(mode = "numeric",length = length(Gene.IDs))
    Features.In.Bins <- vector(mode = "list",length = length(Mean.Quantiles))
    for (i in 1:length(Mean.Quantiles)){
      if (i == 1){
        Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
        Features.In.Bins[[1]] <- Features.In.Bin
      } else {
        Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
        if (length(intersect(Features.In.Bin,unlist(Features.In.Bins))) < length(Features.In.Bin)){
          Features.In.Bin <- Features.In.Bin[!Features.In.Bin %in% unlist(Features.In.Bins)]
          Features.In.Bins[[i]] <- Features.In.Bin
        } else {
          Features.In.Bins[[i]] <- Features.In.Bins[[(i-1)]]
          Features.In.Bin <- Features.In.Bins[[(i-1)]]
        }
      }
      Reference.AlphaQP.Bin <- quantile(Alpha.QP[Features.In.Bin],probs = c(Reference))
      Diff.AlphaQP.AlphaQPFit[Features.In.Bin] <- Alpha.QP[Features.In.Bin] - Reference.AlphaQP.Bin
    }

    PiccoloList$RelevantGenes <- Gene.IDs

    if(is.null(ncol(PiccoloList$RelevantGenes)) == T){
      RelevantGenesLength <- length(PiccoloList$RelevantGenes)
    } else {
      RelevantGenesLength <- length(PiccoloList$RelevantGenes[,1])
    }

    RelevantGenes.Ser.Nos <- rep(0,RelevantGenesLength)
    for (i in 1:length(RelevantGenes.Ser.Nos)){
      if (is.null(ncol(PiccoloList$RelevantGenes)) == T){
        RelevantGenes.Ser.Nos[i] <- which(PiccoloList$Genes == PiccoloList$RelevantGenes[i])
      } else {
        RelevantGenes.Ser.Nos[i] <- which(PiccoloList$Genes[,1] == PiccoloList$RelevantGenes[,1][i])
      }
    }

    PiccoloList$RelevantGenes.Ser.Nos <- RelevantGenes.Ser.Nos

    PiccoloList$DiffAlpha <- Diff.AlphaQP.AlphaQPFit

    #PiccoloList$VG.Ser.Nos <- Top.Features

    if (is.null(NoOfHVG)){
      Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
      Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]
      NoOfHVG <- length(Top.Features)
    } else if (is.numeric(NoOfHVG)){
      Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
      Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]
    } else {
      NoOfHVG <- RelevantGenesLength
      Top.Features <- order(Diff.AlphaQP.AlphaQPFit,decreasing = T)
    }

    if (NoOfHVG > length(Top.Features)){
      if (verbose == T){
        message("Number of HVGs shortlisted with given reference level is lesser than desired number of HVGs. Try lowering Referencelevel.")
      }
    } else {
      Top.Features <- Top.Features[1:NoOfHVG]
    }

    if (is.null(ncol(Gene.IDs)) != T){
      Top.Gene.IDs <- Gene.IDs[Top.Features,]
      Top.Genes <- data.frame(Top.Gene.IDs,Alpha.QP[Top.Features],Alpha.NB.Est[Top.Features],Diff.AlphaQP.AlphaQPFit[Top.Features])
      colnames(Top.Genes) <- c(colnames(Top.Gene.IDs),"AlphaQP","AlphaNB","DiffAlpha")
    } else {
      Top.Gene.IDs <- Gene.IDs[Top.Features]
      Top.Genes <- data.frame(Top.Gene.IDs,Alpha.QP[Top.Features],Alpha.NB.Est[Top.Features],Diff.AlphaQP.AlphaQPFit[Top.Features])
      colnames(Top.Genes) <- c("V1","AlphaQP","AlphaNB","DiffAlpha")
    }

    #Identify least variable features
    Default.Features <- which(Diff.AlphaQP.AlphaQPFit < 0)
    Bottom.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features])]

    if (is.null(ncol(Gene.IDs)) != T){
      Bottom.Gene.IDs <- Gene.IDs[Bottom.Features,]
      Bottom.Genes <- data.frame(Bottom.Gene.IDs,Alpha.QP[Bottom.Features],Alpha.NB.Est[Bottom.Features],Diff.AlphaQP.AlphaQPFit[Bottom.Features])
      colnames(Bottom.Genes) <- c(colnames(Bottom.Gene.IDs),"AlphaQP","AlphaNB","DiffAlpha")
    } else {
      Bottom.Gene.IDs <- Gene.IDs[Bottom.Features]
      Bottom.Genes <- data.frame(Bottom.Gene.IDs,Alpha.QP[Bottom.Features],Alpha.NB.Est[Bottom.Features],Diff.AlphaQP.AlphaQPFit[Bottom.Features])
      colnames(Bottom.Genes) <- c("V1","AlphaQP","AlphaNB","DiffAlpha")
    }

    PiccoloList$DispCoef <- data.frame(AlphaQP = Alpha.QP,AlphaNB = Alpha.NB.Est)

    PiccoloList$HVG <- Top.Genes

    Top.Features.Ser.Nos <- rep(0,length(PiccoloList$HVG$V1))
    for (i in 1:length(Top.Features.Ser.Nos))
    {
      if (is.null(ncol(PiccoloList$Genes)) == T){
        Top.Features.Ser.Nos[i] <- which(PiccoloList$Genes == PiccoloList$HVG$V1[i])
      } else {
        Top.Features.Ser.Nos[i] <- which(PiccoloList$Genes[,1] == PiccoloList$HVG$V1[i])
      }
    }

    PiccoloList$HVG.Ser.Nos <- Top.Features.Ser.Nos

    PiccoloList$StableGenes <- Bottom.Genes

    Bottom.Features.Ser.Nos <- rep(0,length(PiccoloList$StableGenes$V1))
    for (i in 1:length(Bottom.Features.Ser.Nos))
    {
      if (is.null(ncol(PiccoloList$Genes)) == T){
        Bottom.Features.Ser.Nos[i] <- which(PiccoloList$Genes == PiccoloList$StableGenes$V1[i])
      } else {
        Bottom.Features.Ser.Nos[i] <- which(PiccoloList$Genes[,1] == PiccoloList$StableGenes$V1[i])
      }
    }

    PiccoloList$Stable.Ser.Nos <- Bottom.Features.Ser.Nos

    if (verbose == T){
      message("Shortlisted highly variable genes (HVG) and stable genes.")
    }

    if (Out == T){
      write.csv(PiccoloList$HVG, file = paste0("Top",dim(PiccoloList$HVG)[1],"Features", ".csv"),row.names = F)
      message("Successfully prepared .csv file containing list of highly variable features.")
      write.csv(PiccoloList$StableGenes, file = paste0("Bottom",dim(PiccoloList$StableGenes)[1],"Features", ".csv"),row.names = F)
      message("Successfully prepared .csv file containing list of stable features.")
    }

  } else if (is.null(Batch) != T){

    UMI.Mat <- PiccoloList$Counts

    NoOfHVGOrig <- NoOfHVG
    stopifnot(length(Batch) == ncol(UMI.Mat))
    Batch <- as.factor(Batch)
    BatchLevels <- levels(Batch)
    Zero.Count.Features <- vector(mode = "list", length = length(BatchLevels))
    Top.Genes.List <- vector(mode = "list",length = length(BatchLevels))
    Top.Genes.Ser.Nos.List <- vector(mode = "list",length = length(BatchLevels))
    Bottom.Genes.List <- vector(mode = "list",length = length(BatchLevels))
    Bottom.Genes.Ser.Nos.List <- vector(mode = "list",length = length(BatchLevels))
    Diff.Alpha.List <- vector(mode = "list",length = length(BatchLevels))
    RelevantGenes.List <- vector(mode = "list",length = length(BatchLevels))
    RelevantGenes.Ser.Nos.List <- vector(mode = "list",length = length(BatchLevels))
    #VG.Ser.Nos.List <- vector(mode = "list",length = length(BatchLevels))
    Disp.Coef.List <- vector(mode = "list",length = length(BatchLevels))
    for (i in 1:length(BatchLevels)){

      BatchIndex <- which(Batch == BatchLevels[i])

      if(verbose == T){
        message(paste0("Batch ",i))
      }

      Temp.Mat <- UMI.Mat[,BatchIndex]

      Temp.Mat <- Matrix::t(Temp.Mat)

      Gene.IDs <- PiccoloList$Genes

      if (verbose == T){
        message("Filtering features...")
      }

      colVarsSPM <- function(X) {
        stopifnot( methods::is(X, "CsparseMatrix"))
        ans <- sapply(base::seq.int(X@Dim[2]),function(j) {
          if(X@p[j+1] == X@p[j]) { return(0) } # all entries are 0: var is 0
          mean <- base::sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]
          sum((X@x[(X@p[j]+1):X@p[j+1] ] - mean)^2) +
            mean^2 * (X@Dim[1] - (X@p[j+1] - X@p[j]))})/(X@Dim[1] - 1)
        names(ans) <- X@Dimnames[[2]]
        ans
      }

      colOverdispQPCoef <- function(X,alternative = "greater"){
        stopifnot( methods::is(X,"CsparseMatrix"))
        ans <- sapply( base::seq.int(X@Dim[2]),function(j){
          if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
          #mean <- exp(sum(log(X@x[(X@p[j]+1):X@p[j+1]]+1))/X@Dim[1]) - 1
          est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]

          aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                      X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean

          mean(aux) + 1})
      }


      Var.Arith.Per.Feature <- colVarsSPM(Temp.Mat)
      Mean.Arith.Per.Feature <- Matrix::colMeans(Temp.Mat)

      Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)

      No.of.Non.Zero.Per.Feature <- diff(Temp.Mat@p)

      Perc.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature/nrow(Temp.Mat) * 100

      Irrelevant.Features <- unique(c(Irrelevant.Features,which(Perc.Non.Zero.Per.Feature <= MinPercNonZeroCells)))

      if (length(Irrelevant.Features) != 0){
        Temp.Mat <- Temp.Mat[,-Irrelevant.Features]
        if (is.null(ncol(Gene.IDs)) != T){
          Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
        } else {
          Gene.IDs <- Gene.IDs[-Irrelevant.Features]
        }
        Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
        Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
        No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
        Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
      }

      if (verbose == T){
        message("Estimating dispersion coefficients...")
      }

      Alpha.QP <- colOverdispQPCoef(Temp.Mat)

      Irrelevant.Features <- which(Alpha.QP <= 1)

      if (length(Irrelevant.Features) != 0){
        Temp.Mat <- Temp.Mat[,-Irrelevant.Features]
        if (is.null(ncol(Gene.IDs)) != T){
          Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
        } else {
          Gene.IDs <- Gene.IDs[-Irrelevant.Features]
        }
        Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
        Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
        Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
        No.of.Non.Zero.Per.Feature<- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
        Alpha.QP <- Alpha.QP[-Irrelevant.Features]
      }

      Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature

      RelevantGenes.List[[i]] <- Gene.IDs

      if(is.null(ncol(RelevantGenes.List[[i]])) == T){
        RelevantGenesLength <- length(RelevantGenes.List[[i]])
      } else {
        Temp.df <- RelevantGenes.List[[i]]
        RelevantGenesLength <- length(Temp.df[,1])
      }

      RelevantGenes.Ser.Nos <- rep(0,RelevantGenesLength)
      for (k in 1:length(RelevantGenes.Ser.Nos))
      {
        if (is.null(ncol(PiccoloList$RelevantGenes)) == T){
          RelevantGenes.Ser.Nos[k] <- which(PiccoloList$Genes == RelevantGenes.List[[i]][k])
        } else {
          Temp.df <- RelevantGenes.List[[i]]
          RelevantGenes.Ser.Nos[k] <- which(PiccoloList$Genes[,1] == Temp.df[,1][k])
        }
      }

      RelevantGenes.Ser.Nos.List[[i]] <- RelevantGenes.Ser.Nos

      #Binning based approach

      if (verbose == T){
        message("Shortlisting variable features...")
      }

      Mean.Quantiles <- quantile(Mean.Arith.Per.Feature,probs = seq(0.001,1,0.001))
      Diff.AlphaQP.AlphaQPFit <- vector(mode = "numeric",length = length(Gene.IDs))
      Features.In.Bins <- vector(mode = "list",length = length(Mean.Quantiles))
      for (k in 1:length(Mean.Quantiles))
      {
        if (k == 1){
          Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[k])
          Features.In.Bins[[1]] <- Features.In.Bin
        } else {
          #Features.In.Bin <- which(Mean.Arith.Per.Feature > Mean.Quantiles[k-1] & Mean.Arith.Per.Feature <= Mean.Quantiles[k])
          Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[k])
          if (length(intersect(Features.In.Bin,unlist(Features.In.Bins))) < length(Features.In.Bin)){
            Features.In.Bin <- Features.In.Bin[!Features.In.Bin %in% unlist(Features.In.Bins)]
            Features.In.Bins[[k]] <- Features.In.Bin
          } else {
            Features.In.Bins[[k]] <- Features.In.Bins[[(k-1)]]
            Features.In.Bin <- Features.In.Bins[[(k-1)]]
          }
        }

        Reference.AlphaQP.Bin <- quantile(Alpha.QP[Features.In.Bin],probs = c(Reference))
        Diff.AlphaQP.AlphaQPFit[Features.In.Bin] <- Alpha.QP[Features.In.Bin] - Reference.AlphaQP.Bin
      }

      Diff.Alpha.List[[i]] <- Diff.AlphaQP.AlphaQPFit

      #VG.Ser.Nos.List[[i]] <- Top.Features

      if (is.null(NoOfHVGOrig)){
        Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
        Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]
        NoOfHVG <- length(Top.Features)
      } else if (is.numeric(NoOfHVG)){
        Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
        Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]
      } else {
        NoOfHVG <- RelevantGenesLength
        Top.Features <- order(Diff.AlphaQP.AlphaQPFit,decreasing = T)
      }

      # if (is.na(NoOfHVGOrig)){
      #   NoOfHVG <- RelevantGenesLength
      #   Top.Features <- order(Diff.AlphaQP.AlphaQPFit,decreasing = T)
      # } else if (is.null(NoOfHVGOrig)){
      #   Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
      #   Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]
      #   NoOfHVG <- length(Top.Features)
      #   Top.Features <- Top.Features[1:NoOfHVG]
      # } else {
      #   Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
      #   Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]
      # }

      if (NoOfHVG > length(Top.Features)){
        if (verbose == T){
          message("Number of HVGs shortlisted with given reference level is lesser than desired number of HVGs. Try lowering HVGlevel.")
        }
      } else {
        Top.Features <- Top.Features[1:NoOfHVG]
      }

      if (is.null(ncol(Gene.IDs)) != T){
        Top.Gene.IDs <- Gene.IDs[Top.Features,]
        Top.Genes <- data.frame(Top.Gene.IDs,Alpha.QP[Top.Features],Alpha.NB.Est[Top.Features],Diff.AlphaQP.AlphaQPFit[Top.Features])
        colnames(Top.Genes) <- c(colnames(Top.Gene.IDs),"AlphaQP","AlphaNB","DiffAlpha")
      } else {
        Top.Gene.IDs <- Gene.IDs[Top.Features]
        Top.Genes <- data.frame(Top.Gene.IDs,Alpha.QP[Top.Features],Alpha.NB.Est[Top.Features],Diff.AlphaQP.AlphaQPFit[Top.Features])
        colnames(Top.Genes) <- c("V1","AlphaQP","AlphaNB","DiffAlpha")
      }

      #Identify stable features
      Default.Features <- which(Diff.AlphaQP.AlphaQPFit < 0)
      Bottom.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features])]

      if (is.null(ncol(Gene.IDs)) != T){
        Bottom.Gene.IDs <- Gene.IDs[Bottom.Features,]
        Bottom.Genes <- data.frame(Bottom.Gene.IDs,Alpha.QP[Bottom.Features],Alpha.NB.Est[Bottom.Features],Diff.AlphaQP.AlphaQPFit[Bottom.Features])
        #colnames(Bottom.Genes) <- paste0("V",seq(1,ncol(Bottom.Genes),1))
        colnames(Bottom.Genes) <- c(colnames(Bottom.Gene.IDs),"AlphaQP","AlphaNB","DiffAlpha")
      } else {
        Bottom.Gene.IDs <- Gene.IDs[Bottom.Features]
        Bottom.Genes <- data.frame(Bottom.Gene.IDs,Alpha.QP[Bottom.Features],Alpha.NB.Est[Bottom.Features],Diff.AlphaQP.AlphaQPFit[Bottom.Features])
        colnames(Bottom.Genes) <- c("V1","AlphaQP","AlphaNB","DiffAlpha")
      }


      Disp.Coef.List[[i]] <- data.frame(AlphaQP = Alpha.QP,AlphaNB = Alpha.NB.Est)

      Top.Genes.List[[i]] <- Top.Genes

      Temp.df <- Top.Genes.List[[i]]
      Top.Features.Ser.Nos <- rep(0,length(Temp.df$V1))
      for (k in 1:length(Top.Features.Ser.Nos))
      {
        if (is.null(ncol(PiccoloList$Genes)) == T){
          Top.Features.Ser.Nos[k] <- which(PiccoloList$Genes == Temp.df$V1[k])
        } else {
          Top.Features.Ser.Nos[k] <- which(PiccoloList$Genes[,1] == Temp.df$V1[k])
        }
      }

      Top.Genes.Ser.Nos.List[[i]] <- Top.Features.Ser.Nos

      Bottom.Genes.List[[i]] <- Bottom.Genes

      Temp.df <- Bottom.Genes.List[[i]]
      Bottom.Features.Ser.Nos <- rep(0,length(Temp.df$V1))
      for (k in 1:length(Bottom.Features.Ser.Nos))
      {
        if (is.null(ncol(PiccoloList$Genes)) == T){
          Bottom.Features.Ser.Nos[k] <- which(PiccoloList$Genes == Temp.df$V1[k])
        } else {
          Bottom.Features.Ser.Nos[k] <- which(PiccoloList$Genes[,1] == Temp.df$V1[k])
        }
      }
      Bottom.Genes.Ser.Nos.List[[i]] <- Bottom.Features.Ser.Nos
    }

    PiccoloList$RelevantGenes <- RelevantGenes.List

    PiccoloList$RelevantGenes.Ser.Nos <- RelevantGenes.Ser.Nos.List

    PiccoloList$DiffAlpha <- Diff.Alpha.List

    PiccoloList$DispCoef <- Disp.Coef.List

    k <- 1
    Consensus.Top.Ser.Nos <- Top.Genes.Ser.Nos.List[[k]]
    while (k < length(BatchLevels)){
      k <- k + 1
      Consensus.Top.Ser.Nos <- intersect(Consensus.Top.Ser.Nos,Top.Genes.Ser.Nos.List[[k]])
    }

    PiccoloList1 <- Piccolo::SelectStableFeaturesAllCells(PiccoloList = PiccoloList,verbose = T)
    Stable.Genes.AllCells <- PiccoloList1$StableGenes.AllCells
    Stable.Genes.AllCells.SerNos <- PiccoloList1$Stable.Ser.Nos.AllCells

    Consensus.Bottom.Ser.Nos <- Stable.Genes.AllCells.SerNos

    Top.Genes <- vector(mode = "list",length = length(BatchLevels))
    Consensus.Top.Genes.Ser.Nos <- vector(mode = "list",length = length(BatchLevels))
    Bottom.Genes <- vector(mode = "list",length = length(BatchLevels))
    Consensus.Bottom.Genes.Ser.Nos <- vector(mode = "list",length = length(BatchLevels))
    for (k in 1:length(BatchLevels))
    {
      Ser.Nos.Temp <- match(Consensus.Top.Ser.Nos,Top.Genes.Ser.Nos.List[[k]])
      Temp.df <- Top.Genes.List[[k]]
      Top.Genes[[k]] <- Temp.df[Ser.Nos.Temp,]

      Temp.Top.Ser.Nos <- vector(mode = "numeric",length = length(Top.Genes[[k]]$V1))
      for (l in 1:length(Top.Genes[[k]]$V1))
      {
        if (is.null(ncol(PiccoloList$Genes)) == T){
          Temp.Top.Ser.Nos[l] <- which(PiccoloList$Genes == Top.Genes[[k]]$V1[l])
        } else {
          Temp.Top.Ser.Nos[l] <- which(PiccoloList$Genes[,1] == Top.Genes[[k]]$V1[l])
        }
      }
      Consensus.Top.Genes.Ser.Nos[[k]] <- Temp.Top.Ser.Nos

      Ser.Nos.Temp <- match(Consensus.Bottom.Ser.Nos,Bottom.Genes.Ser.Nos.List[[k]])
      Ser.Nos.Temp <- Ser.Nos.Temp[!Ser.Nos.Temp %in% NA]
      Temp.df <- Bottom.Genes.List[[k]]
      Bottom.Genes[[k]] <- Temp.df[Ser.Nos.Temp,]

      Temp.Bottom.Ser.Nos <- vector(mode = "numeric",length = length(Bottom.Genes[[k]]$V1))
      for (l in 1:length(Bottom.Genes[[k]]$V1))
      {
        if (is.null(ncol(PiccoloList$Genes)) == T){
          Temp.Bottom.Ser.Nos[l] <- which(PiccoloList$Genes == Bottom.Genes[[k]]$V1[l])
        } else {
          Temp.Bottom.Ser.Nos[l] <- which(PiccoloList$Genes[,1] == Bottom.Genes[[k]]$V1[l])
        }
      }
      Consensus.Bottom.Genes.Ser.Nos[[k]] <- Temp.Bottom.Ser.Nos
    }

    PiccoloList$HVG <- Top.Genes
    PiccoloList$HVG.Ser.Nos <- Consensus.Top.Genes.Ser.Nos
    PiccoloList$Stable.Genes <- Bottom.Genes
    PiccoloList$Stable.Ser.Nos <- Consensus.Bottom.Genes.Ser.Nos

    if (Out == T){
      for (i in 1:length(BatchLevels))
      {
        write.csv(PiccoloList$HVG[[i]], file = paste0("Top",dim(PiccoloList$HVG[[i]])[1],"Features","_Batch",i,".csv"),row.names = F)
        message("Successfully prepared .csv file containing list of highly variable features.")
        write.csv(PiccoloList$StableGenes[[i]], file = paste0("Bottom",dim(PiccoloList$StableGenes[[i]])[1],"Features","_Batch",i,".csv"),row.names = F)
        message("Successfully prepared .csv file containing list of stable features.")
      }
    }
  }
  return(PiccoloList)
}

#' @title  SelectFeaturesForSeurat function
#' @description  This function performs feature selection by identifying highly variable genes and stable genes using a Seurat object.
#' @export
#' @param Obj A Seurat object. This should be the list created using the CreateSeuratObject function in Seurat
#' @param NoOfHVG A numeric (integer) variable. The number of variable features to be shortlisted. If left NULL (default), will shortlist based on the threshold of ReferenceLevel.
#' @param MinPercNonZeroCells A numeric variable. Specifies the minimum percentage of cells that must have non-zero counts for each gene in the data set. Default is 0.5 (\%).
#' @param Reference A numeric variable (value should be greater than 0 but less than 1). Specifies the reference level against which features are identified as variable. Default is the 10th quantile for each bin (Reference = 0.1).
#' @param Out A logical variable. Specifies whether to return output files (.csv) with the HVGs and the stable genes (if set to T), or not (if set to F). Default is FALSE.
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @return An updated Seurat object containing the information of the shortlisted variable genes in the misc list entities of the assays present in the Seurat object.
#' @examples
#' \dontrun{
#' pbmc3kSeurat <- SelectFeaturesForSeurat(Obj = pbmc3kSeurat)
#' pbmc3kSeurat <- SelectFeaturesForSeurat(Obj = pbmc3kSeurat,
#' NoOfHVG = 3000, Out = T)
#' }
SelectFeaturesForSeurat <- function (Obj, NoOfHVG = NULL, MinPercNonZeroCells = 0.5, Reference = NULL, 
          Out = F, verbose = T) {
  if (is.null(Reference)) {
    Reference <- 0.1
  }
  else if (Reference > 0.5) {
    message("Reference should not be greater than 0.5. Resetting it to 0.1 (default)")
    Reference <- 0.1
  }
  else if (Reference <= 0) {
    message("Reference has to be greater than 0. Resetting it to 0.1 (default)")
    Reference <- 0.1
  }
  
  if (is(Obj@assays$RNA)[1] == "Assay5"){
    UMI.Mat <- Matrix::t(Obj@assays$RNA$counts)
    Gene.IDs <- rownames(Obj@assays$RNA$counts)
  } else {
    UMI.Mat <- Matrix::t(Obj@assays$RNA@counts)
    Gene.IDs <- rownames(Obj@assays$RNA@counts)
  }
    
  if (verbose == T) {
    message("Filtering features...")
  }
  colVarsSPM <- function(X) {
    stopifnot(methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
      if (X@p[j + 1] == X@p[j]) {
        return(0)
      }
      mean <- base::sum(X@x[(X@p[j] + 1):X@p[j + 1]])/X@Dim[1]
      sum((X@x[(X@p[j] + 1):X@p[j + 1]] - mean)^2) + mean^2 * 
        (X@Dim[1] - (X@p[j + 1] - X@p[j]))
    })/(X@Dim[1] - 1)
    names(ans) <- X@Dimnames[[2]]
    ans
  }
  colOverdispQPCoef <- function(X, alternative = "greater") {
    stopifnot(methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
      if (X@p[j + 1] == X@p[j]) {
        return(0)
      }
      est.mean <- sum(X@x[(X@p[j] + 1):X@p[j + 1]])/X@Dim[1]
      aux <- c(((X@x[(X@p[j] + 1):X@p[j + 1]] - est.mean)^2 - 
                  X@x[(X@p[j] + 1):X@p[j + 1]]), rep(est.mean^2, 
                                                     X@Dim[1] - length(X@x[(X@p[j] + 1):X@p[j + 1]])))/est.mean
      mean(aux) + 1
    })
  }
  Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
  Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)
  Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)
  No.of.Non.Zero.Per.Feature <- diff(UMI.Mat@p)
  Perc.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature/nrow(UMI.Mat) * 
    100
  Irrelevant.Features <- unique(c(Irrelevant.Features, which(Perc.Non.Zero.Per.Feature <= 
                                                               MinPercNonZeroCells)))
  if (length(Irrelevant.Features) != 0) {
    UMI.Mat <- UMI.Mat[, -Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T) {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features, ]
    }
    else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
    Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
  }
  if (verbose == T) {
    message("Estimating dispersion coefficients...")
  }
  Alpha.QP <- colOverdispQPCoef(UMI.Mat)
  Irrelevant.Features <- which(Alpha.QP <= 1)
  if (length(Irrelevant.Features) != 0) {
    UMI.Mat <- UMI.Mat[, -Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T) {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features, ]
    }
    else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
    No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
    Alpha.QP <- Alpha.QP[-Irrelevant.Features]
  }
  Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature
  if (verbose == T) {
    message("Shortlisting variable features...")
  }
  Mean.Quantiles <- quantile(Mean.Arith.Per.Feature, probs = seq(0.001,1,0.001))
  Diff.AlphaQP.AlphaQPFit <- vector(mode = "numeric", length = length(Gene.IDs))
  Features.In.Bins <- vector(mode = "list", length = length(Mean.Quantiles))
  for (i in 1:length(Mean.Quantiles)) {
    if (i == 1) {
      Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
      Features.In.Bins[[1]] <- Features.In.Bin
    }
    else {
      Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
      if (length(intersect(Features.In.Bin, unlist(Features.In.Bins))) < 
          length(Features.In.Bin)) {
        Features.In.Bin <- Features.In.Bin[!Features.In.Bin %in% unlist(Features.In.Bins)]
        Features.In.Bins[[i]] <- Features.In.Bin
      }
      else {
        Features.In.Bins[[i]] <- Features.In.Bins[[(i - 1)]]
        Features.In.Bin <- Features.In.Bins[[(i - 1)]]
      }
    }
    Reference.AlphaQP.Bin <- quantile(Alpha.QP[Features.In.Bin], probs = c(Reference))
    Diff.AlphaQP.AlphaQPFit[Features.In.Bin] <- Alpha.QP[Features.In.Bin] - 
      Reference.AlphaQP.Bin
  }
  if (is(Obj@assays$RNA)[1] == "Assay5"){
    Obj.List <- list(Genes = rownames(Obj@assays$RNA$counts))
  } else {
    Obj.List <- list(Genes = rownames(Obj@assays$RNA@counts))
  }

  Obj.List$RelevantGenes <- Gene.IDs
  if (is.null(ncol(Obj.List$RelevantGenes)) == T) {
    RelevantGenesLength <- length(Obj.List$RelevantGenes)
  }
  else {
    RelevantGenesLength <- length(Obj.List$RelevantGenes[,1])
  }
  RelevantGenes.Ser.Nos <- rep(0, RelevantGenesLength)
  for (i in 1:length(RelevantGenes.Ser.Nos)) {
    if (is.null(ncol(Obj.List$Genes)) == T) {
      RelevantGenes.Ser.Nos[i] <- which(Obj.List$Genes == Obj.List$RelevantGenes[i])
    }
    else {
      RelevantGenes.Ser.Nos[i] <- which(Obj.List$Genes[,1] == Obj.List$RelevantGenes[,1][i])
    }
  }
  Obj.List$RelevantGenes.Ser.Nos <- RelevantGenes.Ser.Nos
  Obj.List$DiffAlpha <- Diff.AlphaQP.AlphaQPFit
  if (is.null(NoOfHVG)) {
    Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
    Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features], 
                                           decreasing = T)]
    NoOfHVG <- length(Top.Features)
  }
  else if (is.numeric(NoOfHVG)) {
    Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
    Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features], 
                                           decreasing = T)]
  }else {
    NoOfHVG <- RelevantGenesLength
    Top.Features <- order(Diff.AlphaQP.AlphaQPFit, decreasing = T)
  }
  if (NoOfHVG > length(Top.Features)) {
    if (verbose == T) {
      message("Number of HVGs shortlisted with given reference level is lesser than desired number of HVGs. Try lowering Referencelevel.")
    }
  }else {
    Top.Features <- Top.Features[1:NoOfHVG]
  }
  if (is.null(ncol(Gene.IDs)) != T) {
    Top.Gene.IDs <- Gene.IDs[Top.Features, ]
    Top.Genes <- data.frame(Top.Gene.IDs, Alpha.QP[Top.Features], 
                            Alpha.NB.Est[Top.Features], Diff.AlphaQP.AlphaQPFit[Top.Features])
    colnames(Top.Genes) <- c(colnames(Top.Gene.IDs), "AlphaQP", 
                             "AlphaNB", "DiffAlpha")
  }else {
    Top.Gene.IDs <- Gene.IDs[Top.Features]
    Top.Genes <- data.frame(Top.Gene.IDs, Alpha.QP[Top.Features], 
                            Alpha.NB.Est[Top.Features], Diff.AlphaQP.AlphaQPFit[Top.Features])
    colnames(Top.Genes) <- c("V1", "AlphaQP", "AlphaNB", 
                             "DiffAlpha")
  }
  Default.Features <- which(Diff.AlphaQP.AlphaQPFit < 0)
  Bottom.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features])]
  if (is.null(ncol(Gene.IDs)) != T) {
    Bottom.Gene.IDs <- Gene.IDs[Bottom.Features, ]
    Bottom.Genes <- data.frame(Bottom.Gene.IDs, Alpha.QP[Bottom.Features], 
                               Alpha.NB.Est[Bottom.Features], Diff.AlphaQP.AlphaQPFit[Bottom.Features])
    colnames(Bottom.Genes) <- c(colnames(Bottom.Gene.IDs), 
                                "AlphaQP", "AlphaNB", "DiffAlpha")
  }else {
    Bottom.Gene.IDs <- Gene.IDs[Bottom.Features]
    Bottom.Genes <- data.frame(Bottom.Gene.IDs, Alpha.QP[Bottom.Features], 
                               Alpha.NB.Est[Bottom.Features], Diff.AlphaQP.AlphaQPFit[Bottom.Features])
    colnames(Bottom.Genes) <- c("V1", "AlphaQP", "AlphaNB", 
                                "DiffAlpha")
  }
  Obj.List$DispCoef <- data.frame(AlphaQP = Alpha.QP, AlphaNB = Alpha.NB.Est)
  Obj.List$HVG <- Top.Genes
  Top.Features.Ser.Nos <- rep(0, length(Obj.List$HVG$V1))
  for (i in 1:length(Top.Features.Ser.Nos)) {
    if (is.null(ncol(Obj.List$Genes)) == T) {
      Top.Features.Ser.Nos[i] <- which(Obj.List$Genes == Obj.List$HVG$V1[i])
    }
    else {
      Top.Features.Ser.Nos[i] <- which(Obj.List$Genes[,1] == Obj.List$HVG$V1[i])
    }
  }
  Obj.List$HVG.Ser.Nos <- Top.Features.Ser.Nos
  Obj.List$StableGenes <- Bottom.Genes
  Bottom.Features.Ser.Nos <- rep(0, length(Obj.List$StableGenes$V1))
  for (i in 1:length(Bottom.Features.Ser.Nos)) {
    if (is.null(ncol(Obj.List$Genes)) == T) {
      Bottom.Features.Ser.Nos[i] <- which(Obj.List$Genes == 
                                            Obj.List$StableGenes$V1[i])
    }else {
      Bottom.Features.Ser.Nos[i] <- which(Obj.List$Genes[,1] == Obj.List$StableGenes$V1[i])
    }
  }
  Obj.List$Stable.Ser.Nos <- Bottom.Features.Ser.Nos
  Obj@assays$RNA@misc$PiccoloInfo <- Obj.List
  Obj@assays$SCT@misc$PiccoloInfo <- Obj.List
  if (is(Obj@assays$RNA)[1] == "Assay5"){
    Obj@assays$SCT@var.features <- Obj.List$HVG$V1
  } else {
    Obj@assays$RNA@var.features <- Obj.List$HVG$V1
    Obj@assays$SCT@var.features <- Obj.List$HVG$V1
  }
  
  
  if (verbose == T) {
    message("Shortlisted highly variable genes (HVG) and stable genes.")
  }
  if (Out == T) {
    write.csv(Obj.List$HVG, file = paste0("Top", dim(Obj.List$HVG)[1], 
                                          "Features", ".csv"), row.names = F)
    message("Successfully prepared .csv file containing list of highly variable features.")
    write.csv(Obj.List$StableGenes, file = paste0("Bottom", 
                                                  dim(Obj.List$StableGenes)[1], "Features", ".csv"), 
              row.names = F)
    message("Successfully prepared .csv file containing list of stable features.")
  }
  return(Obj)
}


#' @title  SelectStableFeaturesAllCells function
#' @description  This function performs feature selection by identifying stable genes across all cells.
#' @export
#' @param PiccoloList A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function and processed through the \link[Piccolo]{FilterCells}.
#' @param Reference A numeric variable (value should be greater than 0 but less than 1). Specifies the reference level against which features are identified as variable. Default is the 10th quantile for each bin (Reference = 0.1).
#' @param MinPercNonZeroCells A numeric variable. Specifies the minimum percentage of cells that must have non-zero counts for each gene in the data set. Default is 0.5 (\%).
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @return An updated PiccoloList object containing the data frame with the stable genes.
#' @examples
#' \dontrun{
#' pbmc3k <- SelectStableFeaturesAllCells(PiccoloList = pbmc3k)
#' pbmc3k <- SelectStableFeaturesAllCells(PiccoloList = pbmc3k, verbose = T)
#' }
SelectStableFeaturesAllCells <- function (PiccoloList, Reference = NULL, MinPercNonZeroCells = 0.5, verbose = F) {
  if (is.null(Reference)) {
    Reference <- 0.1
  }
  else if (Reference > 0.5) {
    message("Reference should not be greater than 0.5. Resetting it to 0.1 (default)")
    Reference <- 0.1
  }
  else if (Reference <= 0) {
    message("Reference has to be greater than 0. Resetting it to 0.1 (default)")
    Reference <- 0.1
  }
  UMI.Mat <- Matrix::t(PiccoloList$Counts)
  Gene.IDs <- PiccoloList$Genes
  if (verbose == T) {
    message("Filtering features...")
  }
  colVarsSPM <- function(X) {
    stopifnot(methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
      if (X@p[j + 1] == X@p[j]) {
        return(0)
      }
      mean <- base::sum(X@x[(X@p[j] + 1):X@p[j + 1]])/X@Dim[1]
      sum((X@x[(X@p[j] + 1):X@p[j + 1]] - mean)^2) + mean^2 * 
        (X@Dim[1] - (X@p[j + 1] - X@p[j]))
    })/(X@Dim[1] - 1)
    names(ans) <- X@Dimnames[[2]]
    ans
  }
  colOverdispQPCoef <- function(X, alternative = "greater") {
    stopifnot(methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
      if (X@p[j + 1] == X@p[j]) {
        return(0)
      }
      est.mean <- sum(X@x[(X@p[j] + 1):X@p[j + 1]])/X@Dim[1]
      aux <- c(((X@x[(X@p[j] + 1):X@p[j + 1]] - est.mean)^2 - 
                  X@x[(X@p[j] + 1):X@p[j + 1]]), rep(est.mean^2, X@Dim[1] - length(X@x[(X@p[j] + 1):X@p[j + 1]])))/est.mean
      mean(aux) + 1
    })
  }
  Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
  Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)
  Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)
  No.of.Non.Zero.Per.Feature <- diff(UMI.Mat@p)
  Perc.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature/nrow(UMI.Mat) * 100
  Irrelevant.Features <- unique(c(Irrelevant.Features, which(Perc.Non.Zero.Per.Feature <= MinPercNonZeroCells)))
  if (length(Irrelevant.Features) != 0) {
    UMI.Mat <- UMI.Mat[, -Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T) {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features, ]
    }
    else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
    Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
  }
  if (verbose == T) {
    message("Estimating dispersion coefficients...")
  }
  Alpha.QP <- colOverdispQPCoef(UMI.Mat)
  Irrelevant.Features <- which(Alpha.QP <= 1)
  if (length(Irrelevant.Features) != 0) {
    UMI.Mat <- UMI.Mat[, -Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T) {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features, ]
    }
    else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    Perc.Non.Zero.Per.Feature <- Perc.Non.Zero.Per.Feature[-Irrelevant.Features]
    No.of.Non.Zero.Per.Feature <- No.of.Non.Zero.Per.Feature[-Irrelevant.Features]
    Alpha.QP <- Alpha.QP[-Irrelevant.Features]
  }
  Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature
  if (verbose == T) {
    message("Shortlisting variable features...")
  }
  Mean.Quantiles <- quantile(Mean.Arith.Per.Feature, probs = seq(0.001,1,0.001))
  Diff.AlphaQP.AlphaQPFit <- vector(mode = "numeric", length = length(Gene.IDs))
  Features.In.Bins <- vector(mode = "list", length = length(Mean.Quantiles))
  for (i in 1:length(Mean.Quantiles)) {
    if (i == 1) {
      Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
      Features.In.Bins[[1]] <- Features.In.Bin
    }
    else {
      Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
      if (length(intersect(Features.In.Bin, unlist(Features.In.Bins))) < length(Features.In.Bin)) {
        Features.In.Bin <- Features.In.Bin[!Features.In.Bin %in% unlist(Features.In.Bins)]
        Features.In.Bins[[i]] <- Features.In.Bin
      }
      else {
        Features.In.Bins[[i]] <- Features.In.Bins[[(i - 1)]]
        Features.In.Bin <- Features.In.Bins[[(i - 1)]]
      }
    }
    Reference.AlphaQP.Bin <- quantile(Alpha.QP[Features.In.Bin], probs = c(Reference))
    Diff.AlphaQP.AlphaQPFit[Features.In.Bin] <- Alpha.QP[Features.In.Bin] - Reference.AlphaQP.Bin
  }
  PiccoloList$RelevantGenes.AllCells <- Gene.IDs
  if (is.null(ncol(PiccoloList$RelevantGenes.AllCells)) == T) {
    RelevantGenesLength <- length(PiccoloList$RelevantGenes.AllCells)
  }
  else {
    RelevantGenesLength <- length(PiccoloList$RelevantGenes.AllCells[,1])
  }
  RelevantGenes.Ser.Nos <- rep(0, RelevantGenesLength)
  for (i in 1:length(RelevantGenes.Ser.Nos)) {
    if (is.null(ncol(PiccoloList$Genes)) == T) {
      RelevantGenes.Ser.Nos[i] <- which(PiccoloList$Genes == PiccoloList$RelevantGenes.AllCells[i])
    }
    else {
      RelevantGenes.Ser.Nos[i] <- which(PiccoloList$Genes[,1] == PiccoloList$RelevantGenes.AllCells[, 1][i])
    }
  }
  PiccoloList$RelevantGenes.Ser.Nos.AllCells <- RelevantGenes.Ser.Nos
  PiccoloList$DiffAlpha.AllCells <- Diff.AlphaQP.AlphaQPFit
  Default.Features <- which(Diff.AlphaQP.AlphaQPFit < 0)
  Bottom.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features])]
  if (is.null(ncol(Gene.IDs)) != T) {
    Bottom.Gene.IDs <- Gene.IDs[Bottom.Features, ]
    Bottom.Genes <- data.frame(Bottom.Gene.IDs, Alpha.QP[Bottom.Features], 
                               Alpha.NB.Est[Bottom.Features], Diff.AlphaQP.AlphaQPFit[Bottom.Features])
    colnames(Bottom.Genes) <- c(colnames(Bottom.Gene.IDs), "AlphaQP", "AlphaNB", "DiffAlpha")
  }
  else {
    Bottom.Gene.IDs <- Gene.IDs[Bottom.Features]
    Bottom.Genes <- data.frame(Bottom.Gene.IDs, Alpha.QP[Bottom.Features], 
                               Alpha.NB.Est[Bottom.Features], Diff.AlphaQP.AlphaQPFit[Bottom.Features])
    colnames(Bottom.Genes) <- c("V1", "AlphaQP", "AlphaNB", "DiffAlpha")
  }
  PiccoloList$DispCoef.AllCells <- data.frame(AlphaQP = Alpha.QP, 
                                              AlphaNB = Alpha.NB.Est)
  if (verbose == T) {
    message("Shortlisted stable genes.")
  }
  PiccoloList$StableGenes.AllCells <- Bottom.Genes
  Bottom.Features.Ser.Nos <- rep(0, length(PiccoloList$StableGenes$V1))
  for (i in 1:length(Bottom.Features.Ser.Nos)) {
    if (is.null(ncol(PiccoloList$Genes)) == T) {
      Bottom.Features.Ser.Nos[i] <- which(PiccoloList$Genes == PiccoloList$StableGenes$V1[i])
    }
    else {
      Bottom.Features.Ser.Nos[i] <- which(PiccoloList$Genes[,1] == PiccoloList$StableGenes$V1[i])
    }
  }
  PiccoloList$Stable.Ser.Nos.AllCells <- Bottom.Features.Ser.Nos
  return(PiccoloList)
}

#' @title  Normalize function
#' @description  This function performs normalization for the counts of the genes shortlisted using the \link[Piccolo]{SelectFeatures} function.
#' @export
#' @param PiccoloList A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function and processed through the \link[Piccolo]{FilterCells} and the \link[Piccolo]{SelectFeatures} functions.
#' @param Transform A character variable. Specifies the variance stabilizing transformation that will be applied to the counts. The default is the log transform (Transform = "log"). Other options include the Sqrt transform (Transform = "sqrt") and the Box-Cox power law transform (Transform = "bc").
#' @param SizeFactors A numeric variable. Can be used to specify size factors per cell obtained from another method. Should be specified in the same order as the cells listed in the barcodes file..
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @param Out A logical variable. Specifies whether to return output file (.csv) with the z-scores (if set to T), or not (if set to F). Default is FALSE.
#' @return An updated PiccoloList object containing the matrix with the residuals (z-scores) obtained from the counts matrix.
#' @examples
#' \dontrun{
#' pbmc3k <- Normalize(PiccoloList = pbmc3k)
#' pbmc3k <- Normalize(PiccoloList = pbmc3k, Transform = "bc")
#' }
Normalize <- function (PiccoloList, Transform = "log",SizeFactors = NULL, verbose = T, Out = F){

  if (length(PiccoloList$BatchLabels) != 0){
    Batch <- PiccoloList$BatchLabels
  } else {
    Batch <- NULL
  }

  UMI.Mat <- PiccoloList$Counts
  Gene.IDs <- PiccoloList$Genes
  Barcodes <- PiccoloList$Barcodes
  TransformType <- Transform

  if (is.null(Batch)) {

    UMI.Mat <- Matrix::t(UMI.Mat)

    SF.Per.Cell <- SizeFactors

    PiccoloList$SizeFactors <- SF.Per.Cell

    if (is.null(SizeFactors)){
      #Use all stable features and filtered genes
      Seq.Nos <- 1:dim(PiccoloList$Counts)[1]
      Temp.UMI.Mat <- UMI.Mat[,c(PiccoloList$Stable.Ser.Nos,Seq.Nos[!Seq.Nos %in% PiccoloList$FilteredGenes.Ser.Nos])]
      Zero.Count.Features <- which(Matrix::colSums(Temp.UMI.Mat) == 0)
      if (length(Zero.Count.Features) != 0){
        Temp.UMI.Mat <- Temp.UMI.Mat[,-Zero.Count.Features]
      }

      Total.UMI.Counts <- Matrix::rowSums(Temp.UMI.Mat)

      SF.Per.Cell <- Total.UMI.Counts/mean(Total.UMI.Counts)

      TempI <- 100
      while (length(which(SF.Per.Cell == 0)) != 0){
        Temp.UMI.Mat <- UMI.Mat[,c(PiccoloList$Stable.Ser.Nos[[i]],Seq.Nos[!Seq.Nos %in% PiccoloList$FilteredGenes.Ser.Nos],rev(PiccoloList$HVG.Ser.Nos)[1:TempI])]
        Zero.Count.Features <- which(Matrix::colSums(UMI.Mat) == 0)
        if (length(Zero.Count.Features) != 0){
          Temp.UMI.Mat <- Temp.UMI.Mat[,-Zero.Count.Features]
        }

        Total.UMI.Counts <- Matrix::rowSums(Temp.UMI.Mat)

        SF.Per.Cell <- Total.UMI.Counts/mean(Total.UMI.Counts)

        TempI <- TempI + 100
      }
      PiccoloList$SizeFactors <- SF.Per.Cell
    }

    UMI.Mat <- UMI.Mat[,PiccoloList$HVG.Ser.Nos]

    Temp.SF.Mat <- UMI.Mat/SF.Per.Cell

    if (verbose == T){
      message("Revising dispersion coefficients...")
    }

    colOverdispQPCoef <- function(X,alternative = "greater"){
      stopifnot( methods::is(X,"CsparseMatrix"))
      ans <- sapply( base::seq.int(X@Dim[2]),function(j){
        if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
        #mean <- exp(sum(log(X@x[(X@p[j]+1):X@p[j+1]]+1))/X@Dim[1]) - 1
        est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]

        aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                    X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean

        mean(aux) + 1})
    }

    AlphaRevised <- colOverdispQPCoef(X = Temp.SF.Mat)

    #Find which features exhibit underdispersion
    TempSerNos <- which(AlphaRevised <= 1)
    if (length(TempSerNos) != 0){
      UMI.Mat <- UMI.Mat[,-TempSerNos]
      PiccoloList$HVG <- PiccoloList$HVG[-TempSerNos,]
      PiccoloList$HVG.Ser.Nos <- PiccoloList$HVG.Ser.Nos[-TempSerNos]
      AlphaRevised <- AlphaRevised[-TempSerNos]
    }

    PiccoloList$RevisedAlpha <- AlphaRevised

    if (Transform == "bc" | Transform == "BC"){
      #Box-Cox transform function for sparse CsparseMatrix
      colBCLambda <- function(X,lower = -1, upper = 1, eps = 0.001){
        stopifnot( methods::is(X,"CsparseMatrix"))
        ans <- sapply( base::seq.int(X@Dim[2]),function(j){
          if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
          n <- X@Dim[1]
          yj_LL <- function(lambda){
            if (abs(lambda) < eps){
              x_t <- log(c(X@x[(X@p[j]+1):X@p[j+1]] + 1),rep(1,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))
            } else {
              x_t <- (c((X@x[(X@p[j]+1):X@p[j+1]] + 1),rep(1,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))^lambda - 1)/lambda
            }
            x_t_bar <- mean(x_t)
            x_t_var <- var(x_t) * (n-1)/n
            constant <- sum(sign(X@x[(X@p[j]+1):X@p[j+1]]) * log(abs(X@x[(X@p[j]+1):X@p[j+1]] ) + 1))
            -0.5 * n * log(x_t_var) + (lambda - 1) * constant
          }
          results <- suppressWarnings(optimize(yj_LL, lower = lower, upper = upper, maximum = TRUE, tol = 1e-04))
          results$maximum
        })
      }

      if (verbose == T){
        message("Estimating lambdas for power transform...")
      }
      Lambdas <- colBCLambda(X = UMI.Mat)
      PiccoloList$Lambda <- Lambdas

      Std.Mat <- ResidualSF(X = UMI.Mat,SF = SF.Per.Cell,Transform = TransformType,Lambda = Lambdas)

    } else {

      Std.Mat <- ResidualSF(X = UMI.Mat,SF = SF.Per.Cell,Transform = TransformType,Lambda = NULL)

    }

    if (Out == T) {
      FileName <- paste0(Transform, "Residuals.csv")
      data.table::fwrite(data.frame(Std.Mat), file = FileName,row.names = F,col.names = F, sep = ",")
    }
    PiccoloList$NormCounts <- Std.Mat
    return(PiccoloList)

  } else {
    stopifnot(length(Batch) == ncol(UMI.Mat))
    Batch <- as.factor(Batch)

    BatchLevels <- levels(Batch)
    SF.Per.Cell.List <- vector(mode = "list",length = length(BatchLevels))
    Batch.Cells.Ser.Nos <- vector(mode = "list",length = length(BatchLevels))
    names(SF.Per.Cell.List) <- as.character(BatchLevels)
    Norm.Counts <- vector(mode = "list",length = length(BatchLevels))
    AlphaRevisedList <- vector(mode = "list",length = length(BatchLevels))
    if (Transform == "bc" | Transform == "BC"){
      LambdaList <- vector(mode = "list",length = length(BatchLevels))
    }
    UnderDispersedFeatures <- vector(mode = "list",length = length(BatchLevels))
    for (i in 1:length(BatchLevels)) {
      BatchIndex <- which(Batch == BatchLevels[i])

      Batch.Cells.Ser.Nos[[i]] <- BatchIndex

      if (verbose == T){
        message(paste0("Batch ",i))
      }

      Temp.Mat <- UMI.Mat[,BatchIndex]
      Temp.Mat <- Matrix::t(Temp.Mat)

      Seq.Nos <- 1:dim(PiccoloList$Counts)[1]
      Temp.UMI.Mat <- Temp.Mat[,c(PiccoloList$Stable.Ser.Nos[[i]],Seq.Nos[!Seq.Nos %in% PiccoloList$FilteredGenes.Ser.Nos[[i]]])]
      Zero.Count.Features <- which(Matrix::colSums(Temp.UMI.Mat) == 0)
      if (length(Zero.Count.Features) != 0){
        Temp.UMI.Mat <- Temp.UMI.Mat[,-Zero.Count.Features]
      }

      Total.UMI.Counts <- Matrix::rowSums(Temp.UMI.Mat)

      SF.Per.Cell.Batch <- Total.UMI.Counts/mean(Total.UMI.Counts)

      TempI <- 100
      while (length(which(SF.Per.Cell.Batch == 0)) != 0){
        Temp.UMI.Mat <- Temp.Mat[,c(PiccoloList$Stable.Ser.Nos[[i]],Seq.Nos[!Seq.Nos %in% PiccoloList$FilteredGenes.Ser.Nos[[i]]],rev(PiccoloList$HVG.Ser.Nos[[i]])[1:TempI])]
        Zero.Count.Features <- which(Matrix::colSums(Temp.UMI.Mat) == 0)
        if (length(Zero.Count.Features) != 0){
          Temp.UMI.Mat <- Temp.UMI.Mat[,-Zero.Count.Features]
        }

        Total.UMI.Counts <- Matrix::rowSums(Temp.UMI.Mat)

        SF.Per.Cell.Batch <- Total.UMI.Counts/mean(Total.UMI.Counts)

        TempI <- TempI + 100
      }

      SF.Per.Cell.List[[i]] <- SF.Per.Cell.Batch

      Temp.Mat <- Temp.Mat[,PiccoloList$HVG.Ser.Nos[[i]]]

      Temp.SF.Mat <- Temp.Mat/SF.Per.Cell.Batch

      if (verbose == T){
        message("Revising dispersion coefficients...")
      }

      colOverdispQPCoef <- function(X,alternative = "greater"){
        stopifnot( methods::is(X,"CsparseMatrix"))
        ans <- sapply( base::seq.int(X@Dim[2]),function(j){
          if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
          #mean <- exp(sum(log(X@x[(X@p[j]+1):X@p[j+1]]+1))/X@Dim[1]) - 1
          est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]

          aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                      X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean

          mean(aux) + 1})
      }

      AlphaRevised <- colOverdispQPCoef(X = Temp.SF.Mat)

      UnderDispersedFeatures[[i]] <- which(AlphaRevised <= 1)

      AlphaRevisedList[[i]] <- AlphaRevised

      if (Transform == "bc" | Transform == "BC"){

        #Box-Cox transform function for sparse CsparseMatrix
        colBCLambda <- function(X,lower = -1, upper = 1, eps = 0.001){
          stopifnot( methods::is(X,"CsparseMatrix"))
          ans <- sapply( base::seq.int(X@Dim[2]),function(j){
            if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
            n <- X@Dim[1]
            yj_LL <- function(lambda){
              if (abs(lambda) < eps){
                x_t <- log(c(X@x[(X@p[j]+1):X@p[j+1]] + 1),rep(1,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))
              } else {
                x_t <- (c((X@x[(X@p[j]+1):X@p[j+1]] + 1),rep(1,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))^lambda - 1)/lambda
              }
              x_t_bar <- mean(x_t)
              x_t_var <- var(x_t) * (n-1)/n
              constant <- sum(sign(X@x[(X@p[j]+1):X@p[j+1]]) * log(abs(X@x[(X@p[j]+1):X@p[j+1]] ) + 1))
              -0.5 * n * log(x_t_var) + (lambda - 1) * constant
            }
            results <- suppressWarnings(optimize(yj_LL, lower = lower, upper = upper, maximum = TRUE, tol = 1e-04))
            results$maximum
          })
        }

        if (verbose == T){
          message("Estimating lambdas for power transform...")
        }

        Lambdas <- colBCLambda(X = Temp.Mat)

        LambdaList[[i]] <- Lambdas

        Std.Mat <- ResidualSF(X = Temp.Mat,SF = SF.Per.Cell.Batch,Transform = TransformType,Lambda = Lambdas)

      } else {

        Std.Mat <- ResidualSF(X = Temp.Mat,SF = SF.Per.Cell.Batch,Transform = TransformType,Lambda = NULL)

      }

      Norm.Counts[[i]] <- Std.Mat
    }

    UniqueUnderdispersedFeatures <- unique(unlist(UnderDispersedFeatures))
    if (length(UniqueUnderdispersedFeatures) != 0){
      for (i in 1:length(BatchLevels))
      {
        TempHVG <- PiccoloList$HVG[[i]]
        PiccoloList$HVG[[i]] <- TempHVG[-UniqueUnderdispersedFeatures,]

        TempHVGSerNos <- PiccoloList$HVG.Ser.Nos[[i]]
        PiccoloList$HVG.Ser.Nos[[i]] <- TempHVGSerNos[-UniqueUnderdispersedFeatures]

        TempAlphaRevised <- AlphaRevisedList[[i]]
        AlphaRevisedList[[i]] <- TempAlphaRevised[-UniqueUnderdispersedFeatures]

        TempNormCounts <- Norm.Counts[[i]]
        Norm.Counts[[i]] <- TempNormCounts[-UniqueUnderdispersedFeatures,]
      }
    }

    if (Transform == "bc" | Transform == "BC"){
      PiccoloList$Lambda <- LambdaList
    }
    PiccoloList$Batch.Cells.Ser.Nos <- Batch.Cells.Ser.Nos
    names(PiccoloList$Batch.Cells.Ser.Nos) <- as.character(BatchLevels)
    names(PiccoloList$Batch.Cells.Ser.Nos) <- as.character(BatchLevels)
    names(PiccoloList$HVG) <- as.character(BatchLevels)
    names(PiccoloList$HVG.Ser.Nos) <- as.character(BatchLevels)
    names(AlphaRevisedList) <- as.character(BatchLevels)
    names(Norm.Counts) <- as.character(BatchLevels)

    PiccoloList$SizeFactors <- SF.Per.Cell.List
    PiccoloList$RevisedAlpha <- AlphaRevisedList

    PiccoloList$NormCounts.Batch <- Norm.Counts
    #NormCounts <- do.call(cbind, PiccoloList$NormCounts.Batch)
    NormCounts <- matrix(0,ncol = length(PiccoloList$BatchLabels),nrow = nrow(PiccoloList$NormCounts.Batch[[1]]))
    for (i in 1:length(PiccoloList$Batch.Cells.Ser.Nos)){
        NormCounts[,as.numeric(PiccoloList$Batch.Cells.Ser.Nos[[i]])] <- PiccoloList$NormCounts.Batch[[i]]
    }
    PiccoloList$NormCounts <- NormCounts
    if (Out == T) {
      if (is.null(PiccoloList$BatchLabels) != T){
        FileName <- paste0("BatchesCombined_",Transform, "Residuals.csv")
        data.table::fwrite(data.frame(PiccoloList$NormCounts), file = FileName,row.names = F,col.names = F,sep = ",")
      }
    }
    return(PiccoloList)
  }
}

#' @title  NormalizeForSeurat function
#' @description  This function performs normalization for the counts of the genes shortlisted using the \link[Piccolo]{SelectFeaturesForSeurat} function.
#' @export
#' @param Obj A Seurat object. This should be the list created using the CreateSeuratObject function of Seurat.
#' @param Transform A character variable. Specifies the variance stabilizing transformation that will be applied to the counts. The default is the log transform (Transform = "log"). Other options include the Sqrt transform (Transform = "sqrt") and the Box-Cox power law transform (Transform = "bc").
#' @param SizeFactors A numeric variable. Can be used to specify size factors per cell obtained from another method. Should be specified in the same order as the cells listed in the barcodes file..
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @param Out A logical variable. Specifies whether to return output file (.csv) with the z-scores (if set to T), or not (if set to F). Default is FALSE.
#' @return An updated Seurat object containing the matrix with the residuals (z-scores) obtained from the counts matrix in the scale.data slot in the SCT assay.
#' @examples
#' \dontrun{
#' pbmc3kSeurat <- NormalizeForSeurat(Obj = pbmc3kSeurat)
#' pbmc3kSeurat <- NormalizeForSeurat(Obj = pbmc3kSeurat,
#' Transform = "bc")
#' }
NormalizeForSeurat <- function (Obj, Transform = "log", SizeFactors = NULL, verbose = T, 
          Out = F) {
  if (is(Obj@assays$RNA)[1] == "Assay5"){
    UMI.Mat <- Matrix::t(Obj@assays$RNA$counts)
    Gene.IDs <- rownames(Obj@assays$RNA$counts)
    Barcodes <- colnames(Obj@assays$RNA$counts)
  } else {
    UMI.Mat <- Matrix::t(Obj@assays$RNA@counts)
    Gene.IDs <- rownames(Obj@assays$RNA@counts)
    Barcodes <- colnames(Obj@assays$RNA@counts)
  }
  
  TransformType <- Transform
  PiccoloList <- Obj@assays$RNA@misc$PiccoloInfo
  SF.Per.Cell <- SizeFactors
  PiccoloList$SizeFactors <- SF.Per.Cell
  if (is.null(SizeFactors)) {
    Seq.Nos <- 1:dim(UMI.Mat)[2]
    Temp.UMI.Mat <- UMI.Mat[, c(PiccoloList$Stable.Ser.Nos, 
                                Seq.Nos[!Seq.Nos %in% PiccoloList$FilteredGenes.Ser.Nos])]
    Zero.Count.Features <- which(Matrix::colSums(Temp.UMI.Mat) == 
                                   0)
    if (length(Zero.Count.Features) != 0) {
      Temp.UMI.Mat <- Temp.UMI.Mat[, -Zero.Count.Features]
    }
    Total.UMI.Counts <- Matrix::rowSums(Temp.UMI.Mat)
    SF.Per.Cell <- Total.UMI.Counts/mean(Total.UMI.Counts)
    TempI <- 100
    while (length(which(SF.Per.Cell == 0)) != 0) {
      Temp.UMI.Mat <- UMI.Mat[, c(PiccoloList$Stable.Ser.Nos, 
                                  Seq.Nos[!Seq.Nos %in% PiccoloList$FilteredGenes.Ser.Nos], 
                                  rev(PiccoloList$HVG.Ser.Nos)[1:TempI])]
      Zero.Count.Features <- which(Matrix::colSums(UMI.Mat) == 
                                     0)
      if (length(Zero.Count.Features) != 0) {
        Temp.UMI.Mat <- Temp.UMI.Mat[, -Zero.Count.Features]
      }
      Total.UMI.Counts <- Matrix::rowSums(Temp.UMI.Mat)
      SF.Per.Cell <- Total.UMI.Counts/mean(Total.UMI.Counts)
      TempI <- TempI + 100
    }
    PiccoloList$SizeFactors <- SF.Per.Cell
  }
  UMI.Mat <- UMI.Mat[,PiccoloList$HVG.Ser.Nos]
  Temp.SF.Mat <- UMI.Mat/SF.Per.Cell
  if (verbose == T) {
    message("Revising dispersion coefficients...")
  }
  colOverdispQPCoef <- function(X, alternative = "greater") {
    stopifnot(methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
      if (X@p[j + 1] == X@p[j]) {
        return(0)
      }
      est.mean <- sum(X@x[(X@p[j] + 1):X@p[j + 1]])/X@Dim[1]
      aux <- c(((X@x[(X@p[j] + 1):X@p[j + 1]] - est.mean)^2 - 
                  X@x[(X@p[j] + 1):X@p[j + 1]]), rep(est.mean^2, 
                                                     X@Dim[1] - length(X@x[(X@p[j] + 1):X@p[j + 1]])))/est.mean
      mean(aux) + 1
    })
  }
  AlphaRevised <- colOverdispQPCoef(X = Temp.SF.Mat)
  TempSerNos <- which(AlphaRevised <= 1)
  if (length(TempSerNos) != 0) {
    UMI.Mat <- UMI.Mat[, -TempSerNos]
    PiccoloList$HVG <- PiccoloList$HVG[-TempSerNos, ]
    PiccoloList$HVG.Ser.Nos <- PiccoloList$HVG.Ser.Nos[-TempSerNos]
    AlphaRevised <- AlphaRevised[-TempSerNos]
  }
  PiccoloList$RevisedAlpha <- AlphaRevised
  if (Transform == "bc" | Transform == "BC") {
    colBCLambda <- function(X, lower = -1, upper = 1, eps = 0.001) {
      stopifnot(methods::is(X, "CsparseMatrix"))
      ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
        if (X@p[j + 1] == X@p[j]) {
          return(0)
        }
        n <- X@Dim[1]
        yj_LL <- function(lambda) {
          if (abs(lambda) < eps) {
            x_t <- log(c(X@x[(X@p[j] + 1):X@p[j + 1]] + 
                           1), rep(1, X@Dim[1] - length(X@x[(X@p[j] + 
                                                               1):X@p[j + 1]])))
          }
          else {
            x_t <- (c((X@x[(X@p[j] + 1):X@p[j + 1]] + 
                         1), rep(1, X@Dim[1] - length(X@x[(X@p[j] + 
                                                             1):X@p[j + 1]])))^lambda - 1)/lambda
          }
          x_t_bar <- mean(x_t)
          x_t_var <- var(x_t) * (n - 1)/n
          constant <- sum(sign(X@x[(X@p[j] + 1):X@p[j + 
                                                      1]]) * log(abs(X@x[(X@p[j] + 1):X@p[j + 1]]) + 
                                                                   1))
          -0.5 * n * log(x_t_var) + (lambda - 1) * constant
        }
        results <- suppressWarnings(optimize(yj_LL, lower = lower, 
                                             upper = upper, maximum = TRUE, tol = 1e-04))
        results$maximum
      })
    }
    if (verbose == T) {
      message("Estimating lambdas for power transform...")
    }
    Lambdas <- colBCLambda(X = UMI.Mat)
    PiccoloList$Lambda <- Lambdas
    Std.Mat <- ResidualSF(X = UMI.Mat, SF = SF.Per.Cell, 
                          Transform = TransformType, Lambda = Lambdas)
  }else {
    Std.Mat <- ResidualSF(X = UMI.Mat, SF = SF.Per.Cell, 
                          Transform = TransformType, Lambda = NULL)
  }
  rownames(Std.Mat) <- PiccoloList$HVG$V1
  colnames(Std.Mat) <- Barcodes
  if (Out == T) {
    FileName <- paste0(Transform, "Residuals.csv")
    data.table::fwrite(data.frame(Std.Mat), file = FileName, 
                       row.names = F, col.names = F, sep = ",")
  }
  Obj@assays$SCT@scale.data <- Std.Mat
  Obj@assays$SCT@misc$PiccoloInfo <- PiccoloList
  if (is(Obj@assays$RNA)[1] == "Assay5"){
    UMI.Mat <- Matrix::t(Obj@assays$RNA$counts)
  } else {
    UMI.Mat <- Matrix::t(Obj@assays$RNA@counts)
  }
  
  colVarsSPM <- function(X) {
    stopifnot(methods::is(X, "CsparseMatrix"))
    ans <- sapply(base::seq.int(X@Dim[2]), function(j) {
      if (X@p[j + 1] == X@p[j]) {
        return(0)
      }
      mean <- base::sum(X@x[(X@p[j] + 1):X@p[j + 1]])/X@Dim[1]
      sum((X@x[(X@p[j] + 1):X@p[j + 1]] - mean)^2) + mean^2 * 
        (X@Dim[1] - (X@p[j + 1] - X@p[j]))
    })/(X@Dim[1] - 1)
    names(ans) <- X@Dimnames[[2]]
    ans
  }
  GeneMeans <- Matrix::colMeans(UMI.Mat[, PiccoloList$HVG.Ser.Nos]/PiccoloList$SizeFactors)
  GeneVars <- colVarsSPM(UMI.Mat[, PiccoloList$HVG.Ser.Nos]/PiccoloList$SizeFactors)
  PiccoloList$GeneMeans <- GeneMeans
  PiccoloList$GeneVars <- GeneVars
  CorrectedCountsMat <- Std.Mat * (sqrt(PiccoloList$GeneVars)/(PiccoloList$GeneMeans + 
                                                                 1)) + log(PiccoloList$GeneMeans + 1)
  rownames(CorrectedCountsMat) <- rownames(Std.Mat)
  colnames(CorrectedCountsMat) <- colnames(Std.Mat)
  Obj@assays$SCT@data <- methods::as(CorrectedCountsMat, "CsparseMatrix")
  
  if (is(Obj@assays$RNA)[1] == "Assay5"){
    Obj@assays$SCT@counts <- Obj@assays$RNA$counts[PiccoloList$HVG.Ser.Nos,]
  } else {
    Obj@assays$SCT@counts <- Obj@assays$RNA@counts[PiccoloList$HVG.Ser.Nos,]
  }
  
  Obj@assays$SCT@misc$PiccoloInfo <- PiccoloList
  return(Obj)
}



ResidualSF <- function(X,Transform,SF,Lambda = NULL,verbose = T){

  if (Transform == "log" | Transform == "Log"){

    #For log transform
    if (verbose == T){
      message("Log transforming counts and computing residuals...")
    }

    Seq.Nos <- seq(1,ncol(X),500)

    Trans.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
    for (i in 1:length(Seq.Nos))
    {
      if (i < length(Seq.Nos)){
        Start.Point <- Seq.Nos[i]
        End.Point <- Seq.Nos[i+1] - 1
      } else {
        Start.Point <- Seq.Nos[i]
        End.Point <- ncol(X)
      }

      Temp.UMI.Mat <- as.matrix(X[,Start.Point:End.Point])
      Temp.mat <- matrix(0,nrow = ncol(Temp.UMI.Mat),ncol = nrow(Temp.UMI.Mat))
      StartEndVec <- Start.Point:End.Point

      for (j in 1:ncol(Temp.UMI.Mat))
      {
        RowMean <- mean(Temp.UMI.Mat[,j]/SF)
        RowVar <- var(Temp.UMI.Mat[,j]/SF)

        EstimatedMeans <- SF*RowMean
        Numerator <- log1p(Temp.UMI.Mat[,j]) - log1p(EstimatedMeans)

        Variance <- (SF^2)*RowVar/((EstimatedMeans + 1)^2)
        Denominator <- sqrt(Variance)

        Temp.mat[j,] <- Numerator/Denominator
        Temp.mat[j,] <- Temp.mat[j,] - mean(Temp.mat[j,])
      }

      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  } else if (Transform == "sqrt" | Transform == "Sqrt"){

    #For sqrt transform
    if (verbose == T){
      message("Sqrt transforming counts and computing residuals...")
    }

    Seq.Nos <- seq(1,ncol(X),500)

    Trans.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
    for (i in 1:length(Seq.Nos))
    {
      if (i < length(Seq.Nos)){
        Start.Point <- Seq.Nos[i]
        End.Point <- Seq.Nos[i+1] - 1
      } else {
        Start.Point <- Seq.Nos[i]
        End.Point <- ncol(X)
      }

      Temp.UMI.Mat <- as.matrix(X[,Start.Point:End.Point])
      Temp.mat <- matrix(0,nrow = ncol(Temp.UMI.Mat),ncol = nrow(Temp.UMI.Mat))
      StartEndVec <- Start.Point:End.Point
      for (j in 1:ncol(Temp.UMI.Mat))
      {
        RowMean <- mean(Temp.UMI.Mat[,j]/SF)
        RowVar <- var(Temp.UMI.Mat[,j]/SF)
        EstimatedMeans <- SF*RowMean

        SigmaSqSF <- (SF^2)*RowVar

        Numerator <- sqrt(Temp.UMI.Mat[,j]) - sqrt(EstimatedMeans)
        Denominator <-  sqrt(1/(4 * EstimatedMeans) * SigmaSqSF)

        Temp.mat[j,] <- Numerator/Denominator
        Temp.mat[j,] <- Temp.mat[j,] - mean(Temp.mat[j,])
      }
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  } else if (Transform == "linear" | Transform == "lin" | Transform == "Linear" | Transform == "Lin"){

    #For no VST transform
    if (verbose == T){
      message("Computing residuals for raw counts...")
    }

    Seq.Nos <- seq(1,ncol(X),500)

    #Trans.Mat <- methods::as(matrix(0,nrow = 2,ncol = nrow(X)),"CsparseMatrix")
    Trans.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
    for (i in 1:length(Seq.Nos))
    {
      if (i < length(Seq.Nos)){
        Start.Point <- Seq.Nos[i]
        End.Point <- Seq.Nos[i+1] - 1
      } else {
        Start.Point <- Seq.Nos[i]
        End.Point <- ncol(X)
      }

      Temp.UMI.Mat <- as.matrix(X[,Start.Point:End.Point])
      Temp.mat <- matrix(0,nrow = ncol(Temp.UMI.Mat),ncol = nrow(Temp.UMI.Mat))
      StartEndVec <- Start.Point:End.Point
      for (j in 1:ncol(Temp.UMI.Mat))
      {
        RowMean <- mean(Temp.UMI.Mat[,j]/SF)
        RowVar <- var(Temp.UMI.Mat[,j]/SF)
        EstimatedMeans <- SF*RowMean

        SigmaSqSF <- (SF^2)*RowVar

        Numerator <- Temp.UMI.Mat[,j] - EstimatedMeans
        Denominator <- sqrt(SigmaSqSF)

        Temp.mat[j,] <- Numerator/Denominator
        Temp.mat[j,] <- Temp.mat[j,] - mean(Temp.mat[j,])
      }
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  } else if (Transform == "AnalyticPearson"){

    #For analytic Pearson transform
    if (verbose == T){
      message("Calculating Analytic Pearson residuals...")
    }
    Seq.Nos <- seq(1,ncol(X),500)

    Trans.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
    for (i in 1:length(Seq.Nos))
    {
      if (i < length(Seq.Nos)){
        Start.Point <- Seq.Nos[i]
        End.Point <- Seq.Nos[i+1] - 1
      } else {
        Start.Point <- Seq.Nos[i]
        End.Point <- ncol(X)
      }

      Temp.UMI.Mat <- as.matrix(X[,Start.Point:End.Point])
      Temp.mat <- matrix(0,nrow = ncol(Temp.UMI.Mat),ncol = nrow(Temp.UMI.Mat))
      StartEndVec <- Start.Point:End.Point
      for (j in 1:ncol(Temp.UMI.Mat))
      {
        EstimatedMeans <- SF*mean(Temp.UMI.Mat[,j])
        Numerator <- Temp.UMI.Mat[,j] - EstimatedMeans
        Denominator <- sqrt(EstimatedMeans + 0.01*(EstimatedMeans)^2)

        Temp.mat[j,] <- Numerator/Denominator
        Temp.mat[j,] <- Temp.mat[j,] - mean(Temp.mat[j,])
      }
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  } else if (Transform == "logSF" | Transform == "LogSF"){

    #For log transform
    if (verbose == T){
      message("LogSF transforming counts...")
    }


    Seq.Nos <- seq(1,ncol(X),500)

    #Trans.Mat <- methods::as(matrix(0,nrow = 2,ncol = nrow(X)),"CsparseMatrix")
    Trans.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
    for (i in 1:length(Seq.Nos))
    {
      if (i < length(Seq.Nos)){
        Start.Point <- Seq.Nos[i]
        End.Point <- Seq.Nos[i+1] - 1
      } else {
        Start.Point <- Seq.Nos[i]
        End.Point <- ncol(X)
      }

      Temp.UMI.Mat <- as.matrix(X[,Start.Point:End.Point])
      Temp.mat <- matrix(0,nrow = ncol(Temp.UMI.Mat),ncol = nrow(Temp.UMI.Mat))
      StartEndVec <- Start.Point:End.Point
      for (j in 1:ncol(Temp.UMI.Mat))
      {
        Temp.mat[j,] <- log1p(Temp.UMI.Mat[,j]/SF)
      }
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  }  else if (Transform == "bc" | Transform == "BC"){

    Seq.Nos <- seq(1,ncol(X),500)

    Trans.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
    for (i in 1:length(Seq.Nos))
    {
      if (i < length(Seq.Nos)){
        Start.Point <- Seq.Nos[i]
        End.Point <- Seq.Nos[i+1] - 1
      } else {
        Start.Point <- Seq.Nos[i]
        End.Point <- ncol(X)
      }

      Temp.UMI.Mat <- as.matrix(X[,Start.Point:End.Point])
      Temp.mat <- matrix(0,nrow = ncol(Temp.UMI.Mat),ncol = nrow(Temp.UMI.Mat))

      for (j in 1:ncol(Temp.UMI.Mat))
      {
        UMI.Vec <- Temp.UMI.Mat[,j] + 1

        TempJ <- Start.Point + (j-1)

        EstimatedMeans <- SF*mean(Temp.UMI.Mat[,j]/SF)
        RowVar <- var(Temp.UMI.Mat[,j]/SF)

        Numerator <- ((UMI.Vec^Lambda[TempJ]) - 1)/Lambda[TempJ] - ((EstimatedMeans+1)^Lambda[TempJ] - 1)/Lambda[TempJ]

        Denominator <- sqrt((SF^2)*RowVar) * abs((EstimatedMeans+1)^(Lambda[TempJ] - 1))

        Temp.mat[j,] <- Numerator/Denominator
        Temp.mat[j,] <- Temp.mat[j,] - mean(Temp.mat[j,])
      }
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  }
  return(Trans.Mat)
}

#' @title  ComputePC function
#' @description  This function performs principal components analysis (PCA) using the residuals matrix obtained with the \link[Piccolo]{Normalize} function.
#' @export
#' @param PiccoloList A PiccoloList object containing NormCounts data frame in it.
#' @param NoOfPC An integer variable. Specifies the number of principal components (PCs) to shortlist. Default is 50.
#' @param Out A logical variable. Specifies whether to return output file (.csv) with the coordinates for each cell (if set to T), or not (if set to F). Default is FALSE.
#' @param Center A logical variable. Specifies whether to center the rows of the residuals matrix for PCA (default is F).
#' @param Scale A logical variable. Specifies whether to scale the rows of the residuals matrix for PCA (default is F).
#' @return An updated PiccoloList containing a list with the output of PCA.
#' @examples
#' \dontrun{
#' pbmc3k <- ComputePC(PiccoloList = pbmc3k)
#' pbmc3k <- ComputePC(PiccoloList = pbmc3k, NoOfPC = 20)
#' }
ComputePC <- function (PiccoloList, NoOfPC = 50, Out = F, Center = F, Scale = F){
  if (is.null(PiccoloList$BatchLabels) != T){
    Std.Mat <- matrix(0,ncol = length(PiccoloList$BatchLabels),nrow = nrow(PiccoloList$NormCounts.Batch[[1]]))
    for (i in 1:length(PiccoloList$Batch.Cells.Ser.Nos)){
      Std.Mat[,as.numeric(PiccoloList$Batch.Cells.Ser.Nos[[i]])] <- PiccoloList$NormCounts.Batch[[i]]
    }
  } else {
    Std.Mat <- PiccoloList$NormCounts
  }

  if (is.null(dim(PiccoloList$HVG))){
    if (is.list(PiccoloList$HVG)){
      Features <- PiccoloList$HVG[[1]]$V1
    } else {
      Features <- PiccoloList$HVG
    }
  } else {
    Features <- PiccoloList$HVG$V1
  }

  Barcodes <- PiccoloList$Barcodes
  rownames(Std.Mat) <- Features
  colnames(Std.Mat) <- Barcodes

  MatSVD <- RSpectra::svds(
    A    = t(Std.Mat),
    k    = NoOfPC,
    opts = list(center = Center, scale = Scale)
  )
  MatSVD$sdev <- MatSVD$d * sqrt(nrow(MatSVD$u) - 1)
  MatSVD$x <- MatSVD$u %*% diag(MatSVD$d)
  if (Out == T) {
    PC.df <- data.frame(MatSVD$x)
    FileName <- paste0("Top", NoOfPC, "PrinComp", ".csv")
    data.table::fwrite(PC.df, file = FileName, row.names = T,col.names = F,sep = ",")
  }
  PiccoloList$PCA <- MatSVD
  return(PiccoloList)
}

#' @title  ScreePlot function
#' @description  This function creates the Scree plot to show the percentage of total variation explained by each principal component (PC).
#' @export
#' @param PiccoloList A PiccoloList object after using the \link[Piccolo]{ComputePC} function.
#' @param MaxPC An integer variable. Specifies the maximum number of PCs to use for the Scree plot.
#' @param cex A numeric variable. Specifies the size of the elements in the plot.
#' @return A Scree plot showing the percentage of the overall variance explained by each PC.
#' @examples
#' \dontrun{
#' ScreePlot(PiccoloList = pbmc3k)
#' ScreePlot(PiccoloList = pbmc3k, MaxPC = 15)
#' }
ScreePlot <- function(PiccoloList,MaxPC = NULL,cex = 1.25){
  if (is.null(MaxPC)){
    if (dim(PiccoloList$PCA$x)[2] > 20){
      MaxPC <- 20
    } else {
      MaxPC <- dim(PiccoloList$PCA$x)[2]
    }
  }
  plot((PiccoloList$PCA$sdev^2/sum(PiccoloList$PCA$sdev^2))*100, xlim = c(0, MaxPC), type = "b", pch = 15, xlab = "Principal Components", ylab = "Variance Explained (%)",
       cex.lab = cex,
       cex.axis = cex,
       cex.main = cex,
       cex.sub = cex)
}

#' @title  UMAP coordinates
#' @description  This function generates the 2-dimensional coordinates for the cells using UMAP.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} and the \link[Piccolo]{ComputePC} functions.
#' @param NoOfPC An integer variable. Specifies the number of PCs to use for the UMAP. The default is NULL, which corresponds to the use of all the shortlisted principal components.
#' @param k An integer variable. Specifies the number of nearest neighbors. The default is 10.
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the UMAP coordinates (if set to T), or not (if set to F). Default is F.
#' @return An updated PiccoloList with a data frame containing the UMAP coordinates of the cells.
#' @examples
#' \dontrun{
#' pbmc3k <- UMAPcoords(PiccoloList = pbmc3k,
#' Out = T)
#' }
UMAPcoords <- function (PiccoloList, NoOfPC = NULL, k = 10, Out = F) {
  if (is.null(NoOfPC)) {
    NoOfPC <- dim(PiccoloList$PCA$x)[2]
  }
  PC.Mat <- as.matrix(PiccoloList$PCA$x[,1:NoOfPC])
  rownames(PC.Mat) <- PiccoloList$Barcodes
  x <- uwot::umap(PC.Mat,n_neighbors = k+1,ret_nn = T,pca = NULL,verbose = T)
  UMAP.df <- data.frame(CellID = rownames(x$embedding), UMAP1 = x$embedding[,1], UMAP2 = x$embedding[,2])
  PiccoloList$UMAP <- UMAP.df
  kNN <- list()
  kNN$dist <- x$nn$euclidean$dist[,-1]
  kNN$id <- x$nn$euclidean$idx[,-1]
  kNN$k <- k
  kNN$sort <- TRUE
  PiccoloList$kNN <- kNN
  
  if (Out == T) {
    FileName <- c("UMAPcoords.csv")
    data.table::fwrite(UMAP.df, file = FileName, row.names = F, sep = ",")
    FileName <- c("UMAP_NearestNeighborsID.csv")
    data.table::fwrite(PiccoloList$kNN$id, file = FileName, row.names = T, sep = ",")
    FileName <- c("UMAP_NearestNeighborsDist.csv")
    data.table::fwrite(PiccoloList$kNN$dist, file = FileName, row.names = T, sep = ",")
  }
  
  return(PiccoloList)
}

#' @title  tSNE coordinates
#' @description  This function generates the 2-dimensional coordinates for the cells using tSNE.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} and the \link[Piccolo]{ComputePC} functions.
#' @param NoOfPC An integer variable. Specifies the number of PCs to use for the UMAP. The default is NULL, which corresponds to the use of all the shortlisted principal components.
#' @param k An integer variable. Specifies the number of nearest neighbors. The default is 10.
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the UMAP coordinates (if set to T), or not (if set to F). Default is F.
#' @return A data frame containing the coordinates of the cells in the first 2 UMAP dimensions.
#' @examples
#' \dontrun{
#' pbmc3k <- tSNEcoords(PiccoloList = pbmc3k,
#' Out = T)
#' }
tSNEcoords <- function (PiccoloList,k = 10,NoOfPC = NULL, Out = F) {
  if (is.null(NoOfPC)) {
    NoOfPC <- dim(PiccoloList$PCA$x)[2]
  }
  PC.Mat <- as.matrix(PiccoloList$PCA$x[,1:NoOfPC])
  rownames(PC.Mat) <- PiccoloList$Barcodes
  x <- Rtsne::Rtsne(PC.Mat, pca = F, pca_scale = F)
  tSNE.df <- data.frame(CellID = rownames(PC.Mat), tSNE1 = x$Y[,1], tSNE2 = x$Y[,2])
  x$Y <- tSNE.df
  PiccoloList$kNN <- dbscan::kNN(x = PiccoloList$PCA$x, k = k, 
                                 query = NULL, sort = T, search = "kdtree", bucketSize = 10, 
                                 splitRule = "suggest", approx = 0)
  rownames(PiccoloList$kNN$id) <- PiccoloList$Barcodes
  rownames(PiccoloList$kNN$dist) <- PiccoloList$Barcodes
  if (Out == T) {
    FileName <- c("tSNEcoords.csv")
    data.table::fwrite(tSNE.df, file = FileName, row.names = F, sep = ",")
    FileName <- c("tSNE_NearestNeighborsID.csv")
    data.table::fwrite(PiccoloList$kNN$id, file = FileName, row.names = T, sep = ",")
    FileName <- c("tSNE_NearestNeighborsDist.csv")
    data.table::fwrite(PiccoloList$kNN$dist, file = FileName, row.names = T, sep = ",")
  }
  PiccoloList$tSNE <- x
  return(PiccoloList)
}

#' @title  Find k nearest neighbors
#' @description  This function finds the k nearest neighbors for every cell based on the coordinates in the PC space.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function and the \link[Piccolo]{ComputePC} function.
#' @param k An integer variable. Specifies the number of nearest neighbours we wish to shortlist for every cell. Default is 10.
#' @param query A data matrix with the points to query. If query is not specified, the NN for all the points in x is returned. If query is specified then x needs to be a data matrix.
#' @param sort Sort the neighbors by distance. Note that some search methods already sort the results. Sorting is expensive and sort = FALSE may be much faster for some search methods. kNN objects can be sorted using sort().
#' @param search  Nearest neighbor search strategy (one of "kdtree", "linear" or "dist").
#' @param bucketSize Maximum size of the kd-tree leafs.
#' @param splitRule rule to split the kd-tree. One of "STD", "MIDPT", "FAIR", "SL_MIDPT", "SL_FAIR" or "SUGGEST" (SL stands for sliding). "SUGGEST" uses ANNs best guess.
#' @param approx use approximate nearest neighbors. All NN up to a distance of a factor of 1 + approx eps may be used. Some actual NN may be omitted leading to spurious clusters and noise points. However, the algorithm will enjoy a significant speedup.
#' @return An updated PiccoloList containing nn object of class kNN (subclass of NN) containing a list with the following components: dist(a matrix of distances), id (a matrix with ids), k (number of nearest neighbors used).
#' @examples
#' \dontrun{
#' pbmc3k <- KNearestNeighbors(PiccoloList = pbmc3k)
#' pbmc3k <- KNearestNeighbors(PiccoloList = pbmc3k, k = 15,
#' Out = F)
#' }
KNearestNeighbors <- function(PiccoloList,k = 10,query = NULL,sort = TRUE,search = "kdtree",bucketSize = 10,splitRule = "suggest",approx = 0){
  PiccoloList$kNN <- dbscan::kNN(x = PiccoloList$PCA$x,k = k,query = query,sort = sort,search = search,bucketSize = bucketSize,splitRule = splitRule,approx = approx)
  return(PiccoloList)
}


#' @title  Perform Louvain Clustering
#' @description  This function performs Louvain clustering on the graph based on the k nearest neighbors identified using \link[Piccolo]{KNearestNeighbors} function.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{KNearestNeighbors} function.
#' @param Resolution An numeric variable. An optional parameter that allows the user to adjust the resolution parameter of the modularity function that the algorithm uses internally. Lower values typically yield fewer, larger clusters. The original definition of modularity is recovered when the resolution parameter is set to 1.
#' @return An updated PiccoloList containing a communities object which contains information about the cluster memberships of the cells.
#' @examples
#' \dontrun{
#' pbmc3k <- LouvainClustering(PiccoloList = pbmc3k)
#' pbmc3k <- LouvainClustering(PiccoloList = pbmc3k,
#' Resolution = 1.5)
#' }
LouvainClustering <- function(PiccoloList,Resolution = 1){
  g.Col2 <- apply(PiccoloList$kNN$id,1, as.list)
  g.Col1 <- rep(1:length(g.Col2),each = length(g.Col2[[1]]))
  graph <- igraph::simplify(igraph::graph_from_edgelist(cbind(g.Col1,unlist(g.Col2)),directed = F))
  # louvain clustering in igraph - multi-level modularity optimization algorithm for finding community structure
  louvain <- igraph::cluster_louvain(graph,resolution = Resolution)
  PiccoloList$Louvain <- louvain
  PiccoloList$ClusterLabels <- paste0("Cluster ",louvain$membership)
  return(PiccoloList)
}

#' @title  Perform Leiden Clustering
#' @description  This function performs Leiden clustering on the graph based on the k nearest neighbors identified using \link[Piccolo]{KNearestNeighbors} function.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{KNearestNeighbors} function.
#' @param Resolution An numeric variable. An optional parameter that allows the user to adjust the resolution parameter:. Higher resolutions lead to more smaller communities, while lower resolutions lead to fewer larger communities.
#' @param ObjectiveFunction A character variable. Specified whether to use the Constant Potts Model (CPM) or modularity. Must be either "CPM" or "modularity".
#' @param Weights A positive numeric vector. The weights of the edges. If it is NULL and the input graph has a ‘weight’ edge attribute, then that attribute will be used. If NULL and no such attribute is present, then the edges will have equal weights. Set this to NA if the graph was a ‘weight’ edge attribute, but you don't want to use it for community detection. A larger edge weight means a stronger connection for this function.
#' @param Beta A numeric variable. Specifies the parameter affecting the randomness in the Leiden algorithm. This affects only the refinement step of the algorithm.
#' @param InitialMembership If provided, the Leiden algorithm will try to improve this provided membership. If no argument is provided, the aglorithm simply starts from the singleton partition.
#' @param N_Iterations An integer variable. Specifies	the number of iterations to iterate the Leiden algorithm. Each iteration may improve the partition further.
#' @param VertexWeights The vertex weights used in the Leiden algorithm. If this is not provided, it will be automatically determined on the basis of the objective_function. Please see the details of this function how to interpret the vertex weights.
#' @return An updated PiccoloList containing a communities object which contains information about the cluster memberships of the cells.
#' @examples
#' \dontrun{
#' pbmc3k <- LeidenClustering(PiccoloList = pbmc3k)
#' pbmc3k <- LeidenClustering(PiccoloList = pbmc3k,
#' Resolution = 1.5)
#' }
LeidenClustering <- function(PiccoloList,Resolution = 1,ObjectiveFunction = "modularity",Weights = NULL,Beta = 0.01,InitialMembership = NULL,N_Iterations = 2, VertexWeights = NULL){
  g.Col2 <- apply(PiccoloList$kNN$id,1, as.list)
  g.Col1 <- rep(1:length(g.Col2),each = length(g.Col2[[1]]))
  graph <- igraph::simplify(igraph::graph_from_edgelist(cbind(g.Col1,unlist(g.Col2)),directed = F))

  leiden <- igraph::cluster_leiden(graph,resolution = Resolution,objective_function = ObjectiveFunction,weights = Weights,beta = Beta,initial_membership = InitialMembership,n_iterations = N_Iterations,vertex_weights = VertexWeights)
  PiccoloList$Leiden <- leiden
  PiccoloList$ClusterLabels <- paste0("Cluster ",leiden$membership)
  return(PiccoloList)
}

#' @title  UMAP for Louvain Clusters
#' @description  This function generates the UMAP plot with the cells labeled by the cluster labels identitied by Louvain clustering.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{LouvainClustering} function.
#' @param Levels A character variable. Specifies the order in which the cluster labels should be listed in the UMAP plot.
#' @param Alpha A numeric variable. Specifies the transparency of the dots in the UMAP plot. Smaller values lead to greater transparency. Default is 0.7.
#' @param Size A numeric variable. Specifies the size of the dots in the UMAP plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the UMAP plot. Default is 28.
#' @param Title A character variable. Specifies the title of the UMAP plot.
#' @param LegendPosition A character variable. Specifies the position in the plot where the legend should be placed. Default is "bottom".
#' @return A UMAP plot with the cells colored according to the clusters they belong to as identified using the Louvain algorithm.
#' @examples
#' \dontrun{
#' LouvainUMAP(PiccoloList = pbmc3k)
#' pbmc3k <- LeidenClustering(PiccoloList = pbmc3k,
#' Resolution = 1.5)
#' }
LouvainUMAP <- function (PiccoloList, Levels = NULL, Alpha = 0.7, Size = 1.4, 
                         BaseSize = 28, Title = "Piccolo", LegendPosition = "bottom") {
  UMAP.Coord.df <- PiccoloList$UMAP
  Labels <- as.character(PiccoloList$Louvain$membership)
  if (length(Labels) != length(UMAP.Coord.df$CellID)) {
    stop("The length of the Labels vector provided does not match the number of cells in the UMAP.")
  }
  plot_data <- data.frame(UMAP.Coord.df, Labels)
  colnames(plot_data) <- c("CellID", "UMAP 1", "UMAP 2", "Label")
  if (is.null(Levels) == F) {
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }
  GetClusterMedoidsForLabels <- function(PiccoloList, Clusters) {
    ClusterLabels <- as.factor(Clusters)
    LevelsLabels <- levels(ClusterLabels)
    SerNos <- 1:length(levels(ClusterLabels))
    ClusterLabelSerNos <- rep(0, length(ClusterLabels))
    for (i in 1:length(SerNos)) {
      TempI <- which(ClusterLabels == LevelsLabels[i])
      ClusterLabelSerNos[TempI] <- SerNos[i]
    }
    Medoids <- cluster::medoids(as.matrix(PiccoloList$UMAP[, -1]), clustering = ClusterLabelSerNos)
    names(Medoids) <- LevelsLabels
    return(Medoids)
  }
  Medoids <- GetClusterMedoidsForLabels(PiccoloList = PiccoloList, 
                                        Clusters = PiccoloList$Louvain$membership)
  Trajectory.df <- data.frame(UMAP.Coord.df[Medoids, ], names(Medoids))
  colnames(Trajectory.df) <- c("CellID", "UMAP 1", "UMAP 2", 
                               "Index")
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  if (length(table(Labels)) < 10) {
    cpal <- c("#4E79A7FF", "#F28E2BFF", "#B82E2EFF", "#59A14FFF", 
              "#EDC948FF", "#76B7B2FF", "#B07AA1FF", "#FF9DA7FF", 
              "#9C755FFF", "#BAB0ACFF")
  }
  else {
    cpal <- gg_color_hue(length(table(Labels)))
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(`UMAP 1`, `UMAP 2`)) + 
    ggplot2::geom_point(ggplot2::aes(color = Label), alpha = Alpha, 
                        size = Size) + ggrepel::geom_text_repel(data = Trajectory.df, 
                                                                ggplot2::aes(label = Index, fontface = "bold"), size = 7) + 
    ggplot2::scale_colour_manual(values = cpal) + ggplot2::ggtitle(Title) + 
    ggplot2::theme_bw(base_size = BaseSize) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                                                             panel.grid.minor = ggplot2::element_blank()) + ggplot2::theme(legend.position = LegendPosition)
  return(p)
}


#' @title  UMAP for Leiden Clusters
#' @description  This function generates the UMAP plot with the cells labeled by the cluster labels identitied by Leiden clustering.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{LeidenClustering} function.
#' @param Levels A character variable. Specifies the order in which the cluster labels should be listed in the UMAP plot.
#' @param Alpha A numeric variable. Specifies the transparency of the dots in the UMAP plot. Smaller values lead to greater transparency. Default is 0.7.
#' @param Size A numeric variable. Specifies the size of the dots in the UMAP plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the UMAP plot. Default is 28.
#' @param Title A character variable. Specifies the title of the UMAP plot.
#' @param LegendPosition A character variable. Specifies the position in the plot where the legend should be placed. Default is "bottom".
#' @return A UMAP plot with the cells colored according to the clusters they belong to as identified using the Leiden algorithm.
#' @examples
#' \dontrun{
#' LeidenUMAP(PiccoloList = pbmc3k)
#' }
LeidenUMAP <- function (PiccoloList,Levels = NULL, Alpha = 0.7, Size = 1.4,BaseSize = 28,Title = "Piccolo",LegendPosition = "bottom"){

  UMAP.Coord.df <- PiccoloList$UMAP
  Labels <- as.character(PiccoloList$Leiden$membership)

  if (length(Labels) != length(UMAP.Coord.df$CellID)) {
    stop("The length of the Labels vector provided does not match the number of cells in the UMAP.")
  }
  plot_data <- data.frame(UMAP.Coord.df, Labels)
  colnames(plot_data) <- c("CellID", "UMAP 1", "UMAP 2", "Label")
  if (is.null(Levels) == F) {
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }

  GetClusterMedoidsForLabels <- function(PiccoloList,Clusters){
    ClusterLabels <- as.factor(Clusters)
    LevelsLabels <- levels(ClusterLabels)
    SerNos <- 1:length(levels(ClusterLabels))
    ClusterLabelSerNos <- rep(0,length(ClusterLabels))
    for (i in 1:length(SerNos))
    {
      TempI <- which(ClusterLabels == LevelsLabels[i])
      ClusterLabelSerNos[TempI] <- SerNos[i]
    }

    Medoids <- cluster::medoids(as.matrix(PiccoloList$UMAP[,-1]),clustering = ClusterLabelSerNos)
    names(Medoids) <- LevelsLabels

    return(Medoids)
  }

  Medoids <- GetClusterMedoidsForLabels(PiccoloList = PiccoloList,Clusters = PiccoloList$Leiden$membership)


  Trajectory.df <- data.frame(UMAP.Coord.df[Medoids,],names(Medoids))
  colnames(Trajectory.df) <- c("CellID", "UMAP 1", "UMAP 2", "Index")

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (length(table(Labels)) < 10){
    cpal <- c("#4E79A7FF","#F28E2BFF","#B82E2EFF","#59A14FFF","#EDC948FF","#76B7B2FF","#B07AA1FF","#FF9DA7FF","#9C755FFF","#BAB0ACFF")
  } else {
    cpal <- gg_color_hue(length(table(Labels)))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(`UMAP 1`, `UMAP 2`)) +
    ggplot2::geom_point(ggplot2::aes(color = Label), alpha = Alpha,size = Size) +
    ggrepel::geom_text_repel(data = Trajectory.df,ggplot2::aes(label = Index,fontface = "bold"),size = 7) +
    ggplot2::scale_colour_manual(values=cpal) +
    ggplot2::ggtitle(Title) +
    ggplot2::theme_bw(base_size = BaseSize) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) +
    #ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank()) +
    ggplot2::theme(legend.position = LegendPosition)

  return(p)
}

#' @title  tSNE Plot for Louvain Clusters
#' @description  This function generates the tSNE plot with the cells labeled by the cluster labels identitied by Louvain clustering.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{LouvainClustering} function.
#' @param Levels A character variable. Specifies the order in which the cluster labels should be listed in the tSNE plot.
#' @param Alpha A numeric variable. Specifies the transparency of the dots in the tSNE plot. Smaller values lead to greater transparency. Default is 0.7.
#' @param Size A numeric variable. Specifies the size of the dots in the tSNE plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the tSNE plot. Default is 28.
#' @param Title A character variable. Specifies the title of the tSNE plot.
#' @param LegendPosition A character variable. Specifies the position in the plot where the legend should be placed. Default is "bottom".
#' @return A tSNE plot with the cells colored according to the clusters they belong to, as identified using the Louvain algorithm.
#' @examples
#' \dontrun{
#' LouvaintSNE(PiccoloList = pbmc3k)
#' }
LouvaintSNE <- function (PiccoloList,Levels = NULL, Alpha = 0.7, Size = 1.4,BaseSize = 28,Title = "Piccolo",LegendPosition = "bottom"){

  tSNE.Coord.df <- PiccoloList$tSNE$Y
  Labels <- as.character(PiccoloList$Louvain$membership)

  if (length(Labels) != length(tSNE.Coord.df$CellID)) {
    stop("The length of the Labels vector provided does not match the number of cells in the tSNE.")
  }
  plot_data <- data.frame(tSNE.Coord.df, Labels)
  colnames(plot_data) <- c("CellID", "tSNE 1", "tSNE 2", "Label")
  if (is.null(Levels) == F) {
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }

  GetClusterMedoidsForLabels <- function(PiccoloList,Clusters){
    ClusterLabels <- as.factor(Clusters)
    LevelsLabels <- levels(ClusterLabels)
    SerNos <- 1:length(levels(ClusterLabels))
    ClusterLabelSerNos <- rep(0,length(ClusterLabels))
    for (i in 1:length(SerNos))
    {
      TempI <- which(ClusterLabels == LevelsLabels[i])
      ClusterLabelSerNos[TempI] <- SerNos[i]
    }


    Medoids <- cluster::medoids(as.matrix(PiccoloList$tSNE$Y[,-1]),clustering = ClusterLabelSerNos)
    names(Medoids) <- LevelsLabels

    return(Medoids)
  }

  Medoids <- GetClusterMedoidsForLabels(PiccoloList = PiccoloList,Clusters = PiccoloList$Louvain$membership)


  #Trajectory data frame
  Trajectory.df <- data.frame(tSNE.Coord.df[Medoids,],names(Medoids))
  colnames(Trajectory.df) <- c("CellID", "tSNE 1", "tSNE 2", "Index")

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (length(table(Labels)) < 10){
    cpal <- c("#4E79A7FF","#F28E2BFF","#B82E2EFF","#59A14FFF","#EDC948FF","#76B7B2FF","#B07AA1FF","#FF9DA7FF","#9C755FFF","#BAB0ACFF")
  } else {
    cpal <- gg_color_hue(length(table(Labels)))
  }


   p <- ggplot2::ggplot(plot_data, ggplot2::aes(`tSNE 1`, `tSNE 2`)) +
    ggplot2::geom_point(ggplot2::aes(color = Label), alpha = Alpha,size = Size) +
    ggrepel::geom_text_repel(data = Trajectory.df,ggplot2::aes(label = Index,fontface = "bold"),size = 7) +
    ggplot2::scale_colour_manual(values=cpal) +
    ggplot2::ggtitle(Title) +
    ggplot2::theme_bw(base_size = BaseSize) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = LegendPosition)

  return(p)
}

#' @title  tSNE Plot for Leiden Clusters
#' @description  This function generates the tSNE plot with the cells labeled by the cluster labels identitied by Leiden clustering.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{LeidenClustering} function.
#' @param Levels A character variable. Specifies the order in which the cluster labels should be listed in the tSNE plot.
#' @param Alpha A numeric variable. Specifies the transparency of the dots in the tSNE plot. Smaller values lead to greater transparency. Default is 0.7.
#' @param Size A numeric variable. Specifies the size of the dots in the tSNE plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the tSNE plot. Default is 28.
#' @param Title A character variable. Specifies the title of the tSNE plot.
#' @param LegendPosition A character variable. Specifies the position in the plot where the legend should be placed. Default is "bottom".
#' @return A tSNE plot with the cells colored according to the clusters they belong to, as identified using the Leiden algorithm.
#' @examples
#' \dontrun{
#' LeidentSNE(PiccoloList = pbmc3k)
#' }
LeidentSNE <- function (PiccoloList,Levels = NULL, Alpha = 0.7, Size = 1.4,BaseSize = 28,Title = "Piccolo",LegendPosition = "bottom"){

  tSNE.Coord.df <- PiccoloList$tSNE$Y
  Labels <- as.character(PiccoloList$Leiden$membership)

  if (length(Labels) != length(tSNE.Coord.df$CellID)) {
    stop("The length of the Labels vector provided does not match the number of cells in the tSNE.")
  }
  plot_data <- data.frame(tSNE.Coord.df, Labels)
  colnames(plot_data) <- c("CellID", "tSNE 1", "tSNE 2", "Label")
  if (is.null(Levels) == F) {
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }

  GetClusterMedoidsForLabels <- function(PiccoloList,Clusters){
    ClusterLabels <- as.factor(Clusters)
    LevelsLabels <- levels(ClusterLabels)
    SerNos <- 1:length(levels(ClusterLabels))
    ClusterLabelSerNos <- rep(0,length(ClusterLabels))
    for (i in 1:length(SerNos))
    {
      TempI <- which(ClusterLabels == LevelsLabels[i])
      ClusterLabelSerNos[TempI] <- SerNos[i]
    }


    Medoids <- cluster::medoids(as.matrix(PiccoloList$tSNE$Y[,-1]),clustering = ClusterLabelSerNos)
    names(Medoids) <- LevelsLabels

    return(Medoids)
  }

  Medoids <- GetClusterMedoidsForLabels(PiccoloList = PiccoloList,Clusters = PiccoloList$Leiden$membership)


  #Trajectory data frame
  Trajectory.df <- data.frame(tSNE.Coord.df[Medoids,],names(Medoids))
  colnames(Trajectory.df) <- c("CellID", "tSNE 1", "tSNE 2", "Index")

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (length(table(Labels)) < 10){
    cpal <- c("#4E79A7FF","#F28E2BFF","#B82E2EFF","#59A14FFF","#EDC948FF","#76B7B2FF","#B07AA1FF","#FF9DA7FF","#9C755FFF","#BAB0ACFF")
  } else {
    cpal <- gg_color_hue(length(table(Labels)))
  }


  p <- ggplot2::ggplot(plot_data, ggplot2::aes(`tSNE 1`, `tSNE 2`)) +
    ggplot2::geom_point(ggplot2::aes(color = Label), alpha = Alpha,size = Size) +
    ggrepel::geom_text_repel(data = Trajectory.df,ggplot2::aes(label = Index,fontface = "bold"),size = 7) +
    ggplot2::scale_colour_manual(values=cpal) +
    ggplot2::ggtitle(Title) +
    ggplot2::theme_bw(base_size = BaseSize) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = LegendPosition)

  return(p)
}

#' @title  UMAP Plot With Cell Labels
#' @description  This function generates the UMAP plot with the cells labeled by the cluster labels identitied by Leiden clustering.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{UMAPcoords} function.
#' @param Labels A character variable. Specifies the cell-type/group/batch labels for the cells in the same order that they are listed in the barcodes file.
#' @param Levels A character variable. Specifies the order in which the cluster labels should be listed in the UMAP plot.
#' @param Alpha A numeric variable. Specifies the transparency of the dots in the UMAP plot. Smaller values lead to greater transparency. Default is 0.7.
#' @param Size A numeric variable. Specifies the size of the dots in the UMAP plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the UMAP plot. Default is 28.
#' @param Title A character variable. Specifies the title of the UMAP plot.
#' @param LegendPosition A character variable. Specifies the position in the plot where the legend should be placed. Default is "right".
#' @return A UMAP plot with the cells colored according to the cell labels provided by the user.
#' @examples
#' \dontrun{
#' CellLabels <- c("b-cells","b-cells",..,"cd14 monocytes",..,"NK cells",..)
#' LabelUMAP(PiccoloList = pbmc3k,
#' Labels = CellLabels,
#' Levels = c("b-cells","cd14 monocytes","dendritic","NK cells","naive cytotoxic"),
#' Title = "PBMC3k")
#' }
LabelUMAP <- function (PiccoloList, Labels, Levels = NULL, Alpha = 0.7, Size = 1.4,BaseSize = 28,Title = "Piccolo",LegendPosition = "right"){
  UMAP.Coord.df <- PiccoloList$UMAP
  if (length(Labels) != length(UMAP.Coord.df$CellID)) {
    stop("The length of the Labels vector provided does not match the number of cells in the UMAP.")
  }
  plot_data <- data.frame(UMAP.Coord.df, Labels)
  colnames(plot_data) <- c("CellID", "UMAP 1", "UMAP 2", "Label")
  if (is.null(Levels) == F) {
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (length(table(Labels)) < 10){
    cpal <- c("#4E79A7FF","#F28E2BFF","#B82E2EFF","#59A14FFF","#EDC948FF","#76B7B2FF","#B07AA1FF","#FF9DA7FF","#9C755FFF","#BAB0ACFF")
  } else {
    cpal <- gg_color_hue(length(table(Labels)))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(`UMAP 1`, `UMAP 2`)) +
    ggplot2::geom_point(ggplot2::aes(color = Label), alpha = Alpha,size = Size) +
    ggplot2::scale_colour_manual(values=cpal) +
    ggplot2::ggtitle(Title) +
    ggplot2::theme_bw(base_size = BaseSize) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) +
    #ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank()) +
    ggplot2::theme(legend.position = LegendPosition)

  return(p)
}


#' @title  tSNE Plot With Cell Labels
#' @description  This function generates the tSNE plot with the cells labeled by the cluster labels identitied by Leiden clustering.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{tSNEcoords} function.
#' @param Labels A character variable. Specifies the cell-type/group/batch labels for the cells in the same order that they are listed in the barcodes file.
#' @param Levels A character variable. Specifies the order in which the cluster labels should be listed in the tSNE plot.
#' @param Alpha A numeric variable. Specifies the transparency of the dots in the tSNE plot. Smaller values lead to greater transparency. Default is 0.7.
#' @param Size A numeric variable. Specifies the size of the dots in the tSNE plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the tSNE plot. Default is 28.
#' @param Title A character variable. Specifies the title of the tSNE plot.
#' @param LegendPosition A character variable. Specifies the position in the plot where the legend should be placed. Default is "right".
#' @return A tSNE plot with the cells colored according to the cell labels provided by the user.
#' @examples
#' \dontrun{
#' LabeltSNE(PiccoloList = pbmc3k)
#' LabeltSNE(PiccoloList = pbmc3k, Title = "PBMC3k")
#' }
LabeltSNE <- function (PiccoloList, Labels, Levels = NULL, Alpha = 0.7, Size = 1.4,BaseSize = 24,Title = "Piccolo",LegendPosition = "right"){
  tSNE.Coord.df <- PiccoloList$tSNE$Y
  if (length(Labels) != length(tSNE.Coord.df$CellID)) {
    stop("The length of the Labels vector provided does not match the number of cells in the tSNE.")
  }
  plot_data <- data.frame(tSNE.Coord.df, Labels)
  colnames(plot_data) <- c("CellID", "tSNE 1", "tSNE 2", "Label")
  if (is.null(Levels) == F) {
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (length(table(Labels)) < 10){
    cpal <- c("#4E79A7FF","#F28E2BFF","#B82E2EFF","#59A14FFF","#EDC948FF","#76B7B2FF","#B07AA1FF","#FF9DA7FF","#9C755FFF","#BAB0ACFF")
  } else {
    cpal <- gg_color_hue(length(table(Labels)))
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(`tSNE 1`, `tSNE 2`)) +
    ggplot2::geom_point(ggplot2::aes(color = Label), alpha = Alpha,size = Size) +
    ggplot2::scale_colour_manual(values=cpal) +
    ggplot2::ggtitle(Title) +
    ggplot2::theme_bw(base_size = BaseSize) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank()) +
    #ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank()) +
    ggplot2::theme(legend.position = LegendPosition)

  return(p)
}


#' @title  Calculate z-scores for gene set
#' @description  This function combines and calculates z-scores for gene sets (even works for single genes).
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function.
#' @param Genes A character variable. Specifies gene IDs or names. In case of features.tsv or genes.tsv file that contain multiple columns, the input IDs should be present in IDs listed in the first column of the file.
#' @return An updated PiccoloList with a list item containing the combined z-scores for the input gene set.
#' @examples
#' \dontrun{
#' GeneSet <- c("ENSG00000123374","ENSG00000135446",
#' "ENSG00000124762","ENSG00000111276","ENSG00000118007")
#' pbmc3k <- GeneSetZScores(PiccoloList = pbmc3k,Genes = GeneSet)
#' GeneSet <- c("CDK2","CDK4","CDKN1A","CDKN1B","STAG1")
#' pbmc3k <- GeneSetZScores(PiccoloList = pbmc3k,Genes = GeneSet)
#' }
GeneSetZScores <- function(PiccoloList,Genes){

  if (is.null(dim(PiccoloList$HVG))){
    GenesSerNos <- which(PiccoloList$HVG[[1]]$V1 %in% Genes)
  } else {
    GenesSerNos <- which(PiccoloList$HVG$V1 %in% Genes)
    if (length(GenesSerNos) == 0){
      stop("None of the genes provided were shortlisted as variable genes.")
    }
  }
  ZScoresMat <- matrix(0,nrow = length(GenesSerNos),ncol = ncol(PiccoloList$NormCounts))
  for (i in 1:nrow(ZScoresMat))
  {
    TempVec <- PiccoloList$NormCounts[GenesSerNos[i],]
    ZScoresMat[i,] <- (TempVec - mean(TempVec))/sd(TempVec)
  }
  StoufferZScores <- colSums(ZScoresMat)/sqrt(nrow(ZScoresMat))
  GenesZScores <- (StoufferZScores - mean(StoufferZScores))/sd(StoufferZScores)
  PiccoloList$GeneSetZScore <- GenesZScores
  return(PiccoloList)
}

#' @title  UMAP plot with z-scores for gene sets
#' @description  This function generates the UMAP plot showing the z-scores for individual genes or gene sets.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{GeneSetZScores} function.
#' @param Name A character variable. Specifies the name of the gene or gene set that you would like to display as the title.
#' @param xLabel A logical variable. Specifies whether the x-axis label should be displayed (T) or not (F). Default is T.
#' @param Size A numeric variable. Specifies the size of the dots in the UMAP plot. Default is 1.4.
#' @param yLabel A logical variable. Specifies whether the y-axis label should be displayed (T) or not (F). Default is T.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the UMAP plot. Default is 28.
#' @param UpperLowerSDThreshold A numeric variable. Specifies the number of standard deviations (SD) above or below which the z-scores should be capped to the value at the specified SD. Default is 2.5.
#' @param col_pal A character variable or vector. Color palette for z-scores. The default is Shuksan, another option is "viridis". Users can themselves also specify a vector containing a gradient of colors.
#' @return A UMAP plot with the cells colored according to the z-scores based on input gene sets specified in \link[Piccolo]{GeneSetZScores}.
#' @examples
#' \dontrun{
#' UMAPZScores(PiccoloList = pbmc3k,Name = "Cell Cycle")
#' UMAPZScores(PiccoloList = pbmc3k,Name = "Cell Cycle",
#' UpperLowerSDThreshold = 3.5)
#' }

UMAPZScores <- function (PiccoloList, Name, xLabel = T, yLabel = T, Size = 1.4, 
                         BaseSize = 28, UpperLowerSDThreshold = 2.5,col_pal = NULL) {
  GeneSetZScores <- PiccoloList$GeneSetZScore
  Outliers.High <- which(GeneSetZScores > mean(GeneSetZScores) + UpperLowerSDThreshold * sd(GeneSetZScores))
  if (length(Outliers.High) != 0) {
    GeneSetZScores[Outliers.High] <- mean(GeneSetZScores) + UpperLowerSDThreshold * sd(GeneSetZScores)
  }
  Outliers.Low <- which(GeneSetZScores < mean(GeneSetZScores) - UpperLowerSDThreshold * sd(GeneSetZScores))
  if (length(Outliers.Low) != 0) {
    GeneSetZScores[Outliers.Low] <- mean(GeneSetZScores) - 
      UpperLowerSDThreshold * sd(GeneSetZScores)
  }
  df <- data.frame(PiccoloList$UMAP[,2:3], GeneSetZScores)
  colnames(df) <- c("UMAP1", "UMAP2", "Z Score")
  if (xLabel == T) {
    xlabeltext <- "UMAP 1"
  }
  else {
    xlabeltext <- ""
  }
  if (yLabel == T) {
    ylabeltext <- "UMAP 2"
  }
  else {
    ylabeltext <- ""
  }
  
  if (is.null(col_pal)){
    ggplot2::ggplot(data = df, ggplot2::aes(x = UMAP1, y = UMAP2)) + 
      ggplot2::geom_point(ggplot2::aes(UMAP1, UMAP2, color = `Z Score`), 
                          size = Size) + 
      ggplot2::scale_color_gradientn(colors = rev(c("grey28", "#74677EFF", "#AC8EABFF", "#D7B1C5FF", "#EBBDC8FF", "#F2CEC7FF","#F8E3D1FF","#FEFBE9FF"))) + #"#33271EFF"
      
      ggplot2::theme_bw(base_size = BaseSize, base_line_size = 0.4) + 
      ggplot2::ggtitle(Name) + ggplot2::xlab(xlabeltext) + 
      ggplot2::ylab(ylabeltext) + 
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                     axis.ticks.y = ggplot2::element_blank())
  } else if (col_pal == "viridis"){
      ggplot2::ggplot(data = df, ggplot2::aes(x = UMAP1, y = UMAP2)) + 
        ggplot2::geom_point(ggplot2::aes(UMAP1, UMAP2, color = `Z Score`), 
                            size = Size) + 
        viridis::scale_color_viridis() + 
        ggplot2::theme_bw(base_size = BaseSize, base_line_size = 0.4) + 
        ggplot2::ggtitle(Name) + ggplot2::xlab(xlabeltext) + 
        ggplot2::ylab(ylabeltext) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                       axis.ticks.y = ggplot2::element_blank())
      
  } else {
      ggplot2::ggplot(data = df, ggplot2::aes(x = UMAP1, y = UMAP2)) + 
        ggplot2::geom_point(ggplot2::aes(UMAP1, UMAP2, color = `Z Score`), 
                            size = Size) + 
        ggplot2::scale_color_gradientn(colors = col_pal) + 
        
        ggplot2::theme_bw(base_size = BaseSize, base_line_size = 0.4) + 
        ggplot2::ggtitle(Name) + ggplot2::xlab(xlabeltext) + 
        ggplot2::ylab(ylabeltext) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                       axis.ticks.y = ggplot2::element_blank())
  }
    
}

#' @title  tSNE plot with z-scores for gene sets
#' @description  This function generates the tSNE plot showing the z-scores for individual genes or gene sets.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{GeneSetZScores} function.
#' @param Name A character variable. Specifies the name of the gene or gene set that you would like to display as the title.
#' @param xLabel A logical variable. Specifies whether the x-axis label should be displayed (T) or not (F). Default is T.
#' @param yLabel A logical variable. Specifies whether the y-axis label should be displayed (T) or not (F). Default is T.
#' @param Size A numeric variable. Specifies the size of the dots in the tSNE plot. Default is 1.4.
#' @param BaseSize A numeric variable. Specifies the base size of the text elements in the tSNE plot. Default is 28.
#' @param UpperLowerSDThreshold A numeric variable. Specifies the number of standard deviations (SD) above or below which the z-scores should be capped to the value at the specified SD.
#' @param col_pal A character variable or vector. Color palette for z-scores. The default is Shuksan, another option is "viridis". Users can themselves also specify a vector containing a gradient of colors.
#' @return A tSNE plot with the cells colored according to the z-scores based on input gene sets specified in \link[Piccolo]{GeneSetZScores}.
#' @examples
#' \dontrun{
#' tSNEZScores(PiccoloList = pbmc3k,Name = "Cell Cycle")
#' tSNEZScores(PiccoloList = pbmc3k,Name = "Cell Cycle",
#' UpperLowerSDThreshold = 3.5)
#' }
tSNEZScores <- function (PiccoloList, Name, xLabel = T, yLabel = T, Size = 1.4, 
                         BaseSize = 28, UpperLowerSDThreshold = 2.5) {
  GeneSetZScores <- PiccoloList$GeneSetZScore
  Outliers.High <- which(GeneSetZScores > mean(GeneSetZScores) + 
                           UpperLowerSDThreshold * sd(GeneSetZScores))
  if (length(Outliers.High) != 0) {
    GeneSetZScores[Outliers.High] <- mean(GeneSetZScores) + 
      UpperLowerSDThreshold * sd(GeneSetZScores)
  }
  Outliers.Low <- which(GeneSetZScores < mean(GeneSetZScores) - 
                          UpperLowerSDThreshold * sd(GeneSetZScores))
  if (length(Outliers.Low) != 0) {
    GeneSetZScores[Outliers.Low] <- mean(GeneSetZScores) - 
      UpperLowerSDThreshold * sd(GeneSetZScores)
  }
  df <- data.frame(PiccoloList$tSNE$Y[, 2:3], GeneSetZScores)
  colnames(df) <- c("tSNE1", "tSNE2", "Z Score")
  if (xLabel == T) {
    xlabeltext <- "tSNE 1"
  }
  else {
    xlabeltext <- ""
  }
  if (yLabel == T) {
    ylabeltext <- "tSNE 2"
  }
  else {
    ylabeltext <- ""
  }
  if (is.null(col_pal)){
    ggplot2::ggplot(data = df, ggplot2::aes(x = UMAP1, y = UMAP2)) + 
      ggplot2::geom_point(ggplot2::aes(UMAP1, UMAP2, color = `Z Score`), size = Size) + 
      ggplot2::scale_color_gradientn(colors = rev(c("grey28", "#74677EFF", "#AC8EABFF", "#D7B1C5FF", "#EBBDC8FF", "#F2CEC7FF","#F8E3D1FF","#FEFBE9FF"))) +
      
      ggplot2::theme_bw(base_size = BaseSize, base_line_size = 0.4) + 
      ggplot2::ggtitle(Name) + ggplot2::xlab(xlabeltext) + 
      ggplot2::ylab(ylabeltext) + 
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                     axis.ticks.y = ggplot2::element_blank())
    } else if (col_pal == "viridis"){
      ggplot2::ggplot(data = df, ggplot2::aes(x = UMAP1, y = UMAP2)) + 
        ggplot2::geom_point(ggplot2::aes(UMAP1, UMAP2, color = `Z Score`), 
                            size = Size) + 
        viridis::scale_color_viridis() + 
        ggplot2::theme_bw(base_size = BaseSize, base_line_size = 0.4) + 
        ggplot2::ggtitle(Name) + ggplot2::xlab(xlabeltext) + 
        ggplot2::ylab(ylabeltext) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                       axis.ticks.y = ggplot2::element_blank())
      
    } else {
      ggplot2::ggplot(data = df, ggplot2::aes(x = UMAP1, y = UMAP2)) + 
        ggplot2::geom_point(ggplot2::aes(UMAP1, UMAP2, color = `Z Score`), 
                            size = Size) + 
        ggplot2::scale_color_gradientn(colors = col_pal) + 
        
        ggplot2::theme_bw(base_size = BaseSize, base_line_size = 0.4) + 
        ggplot2::ggtitle(Name) + ggplot2::xlab(xlabeltext) + 
        ggplot2::ylab(ylabeltext) + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                       axis.ticks.y = ggplot2::element_blank())
    }
}


#' @title  Identify differentially expressed genes between 2 groups of cells
#' @description  This function performs differential expression analysis between 2 groups of cells provided by the user.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function.
#' @param Group1 A numeric (integers) vector. Specifies the serial numbers of cells in group 1 (serial numbers based on the order of barcodes).
#' @param Group2 A numeric (integers) vector. Specifies the serial numbers of cells in group 2 (serial numbers based on the order of barcodes).
#' @param Transform A character variable. Should be the same as the one used during normalization. Else, the default of "log" will be used.
#' @param Method A character variable. Specifies the method used for testing differences between the expression levels in the 2 groups. Currently, two tests are available: the Student's t-test ("t.test) and the Wilcoxon rank-sum test ("wilcoxon").
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the differential expression result (if set to T). Default is F.
#' @return An updated PiccoloList containing a data frame with the results of the differential expression analysis.
#' @examples
#' \dontrun{
#' Group1.vec <- 1:200
#' Group2.vec <- 301:500
#' pbmc3k <- PerformDiffExp(PiccoloList = pbmc3k,
#' Group1 = Group1.vec,
#' Group2 = Group2.vec,
#' Out = T)
#' }
PerformDiffExp <- function(PiccoloList,Group1,Group2,Transform = "log",Method = "t.test",Out = F){
  Relevant.Counts.Mat <- Matrix::t(PiccoloList$Counts)
  if (is.null(dim(PiccoloList$Genes))) {
    if (is.list(PiccoloList$HVG)) {
      Features <- PiccoloList$HVG[[1]]
      FeaturesSerNos <- PiccoloList$HVG.Ser.Nos[[1]]
    } else {
      Features <- PiccoloList$HVG
      FeaturesSerNos <- PiccoloList$HVG.Ser.Nos
    }
  } else {
    Features <- data.frame(V1 = PiccoloList$HVG$V1,V2 = PiccoloList$HVG$V2)
    FeaturesSerNos <- PiccoloList$HVG.Ser.Nos
  }

  Relevant.Counts.Mat <- Relevant.Counts.Mat[,FeaturesSerNos]
  if (is.list(PiccoloList$SizeFactors)){
    Scaled.Counts.Mat <- Relevant.Counts.Mat/unlist(PiccoloList$SizeFactors)
  } else {
    Scaled.Counts.Mat <- Relevant.Counts.Mat/PiccoloList$SizeFactors
  }  
  
  Relevant.Counts.Mat.Group1 <- Matrix::t(Scaled.Counts.Mat)
  Relevant.Counts.Mat.Group1 <- Relevant.Counts.Mat.Group1[,Group1]
  Relevant.Counts.Mat.Group1 <- Matrix::t(Relevant.Counts.Mat.Group1)
  PercNonZero.Group1 <- diff(Relevant.Counts.Mat.Group1@p)/length(Group1) * 100
  Relevant.Counts.Mat.Group2 <- Matrix::t(Scaled.Counts.Mat)
  Relevant.Counts.Mat.Group2 <- Relevant.Counts.Mat.Group2[,Group2]
  Relevant.Counts.Mat.Group2 <- Matrix::t(Relevant.Counts.Mat.Group2)
  PercNonZero.Group2 <- diff(Relevant.Counts.Mat.Group2@p)/length(Group2) * 100
  GeneMeans.Group1 <- Matrix::colMeans(Relevant.Counts.Mat.Group1)
  GeneMeans.Group2 <- Matrix::colMeans(Relevant.Counts.Mat.Group2)
  Log2FC.Means <- log2((GeneMeans.Group1 + 1/dim(Relevant.Counts.Mat)[1])/(GeneMeans.Group2 + 1/dim(Relevant.Counts.Mat)[1]))
  Barcodes <- PiccoloList$Barcodes
  Seq.Nos <- seq(1, ncol(Relevant.Counts.Mat), 500)
  if (Method == "t.test"){
    Res.df <- data.frame(obs.x = c(1,2), obs.y = c(1,2), obs.tot = c(1,2), mean.x = c(1,2), mean.y = c(1,2), mean.diff = c(1,2), var.x = c(1,2), var.y = c(1,2), var.pooled = c(1,2), stderr = c(1,2), df = c(1,2), statistic = c(1,2), pvalue = c(1,2), conf.low = c(1,2), conf.high = c(1,2))
  } else if (Method == "wilcoxon"){
    Res.df <- data.frame(obs.x = c(1,2), obs.y = c(1,2), obs.tot = c(1,2), mean.x = c(1,2), mean.y = c(1,2), mean.diff = c(1,2), statistic = c(1,2), pvalue = c(1,2))
  }

  if (Method == "t.test"){
      Temp.df <- matrixTests::row_t_equalvar(PiccoloList$NormCounts[,Group1],PiccoloList$NormCounts[,Group2])
      Temp.df <- Temp.df[,1:15]
  } else if (Method == "wilcoxon") {
    Temp.Res.df <- matrixTests::row_wilcoxon_twosample(PiccoloList$NormCounts[,Group1],PiccoloList$NormCounts[,Group2])
    Mean.x <- rowMeans(PiccoloList$NormCounts[,Group1])
    Mean.y <- rowMeans(PiccoloList$NormCounts[,Group2])
    Mean.Diff <- Mean.x - Mean.y
    Temp.Res.df <- Temp.Res.df[,1:5]
    Temp.df <- data.frame(obs.x = Temp.Res.df$obs.x, obs.y = Temp.Res.df$obs.y, obs.tot = Temp.Res.df$obs.tot, mean.x = Mean.x, mean.y = Mean.y, mean.diff = Mean.Diff, statistic = Temp.Res.df$statistic, pvalue = Temp.Res.df$pvalue)
  }

  Res.df <- rbind(Res.df,Temp.df)
  Res.df <- Res.df[-c(1,2),]

  p.val.vec <- Res.df$pvalue
  p.adj.vec <- p.adjust(p.val.vec, method = "BH")
  Genes <- Features
  SortOrder <- order(p.adj.vec)
  Res.df <- Res.df[SortOrder,]
  log2FC.vec <- Log2FC.Means[SortOrder]
  p.adj.vec <- p.adj.vec[SortOrder]
  GeneMeans.Group1 <- GeneMeans.Group1[SortOrder]
  GeneMeans.Group2 <- GeneMeans.Group2[SortOrder]
  PercNonZero.Group1 <- PercNonZero.Group1[SortOrder]
  PercNonZero.Group2 <- PercNonZero.Group2[SortOrder]
  if (is.null(dim(Genes))){
    Genes <- Genes[SortOrder]
    DE.Res.df <- data.frame(Genes, Res.df, p.adj.vec, log2FC.vec, GeneMeans.Group1, GeneMeans.Group2, PercNonZero.Group1, PercNonZero.Group2)
    colnames(DE.Res.df) <- c("Gene.ID", colnames(Res.df), "p.adj(BH)", "Log2FC.Means", "Mean.Group1", "Mean.Group2", "PercNonZero.Group1", "PercNonZero.Group2")
  }else {
    Genes <- Genes[SortOrder,]
    DE.Res.df <- data.frame(Genes, Res.df, p.adj.vec, log2FC.vec, GeneMeans.Group1, GeneMeans.Group2, PercNonZero.Group1, PercNonZero.Group2)
    colnames(DE.Res.df) <- c(colnames(Genes), colnames(Res.df), "p.adj(BH)", "Log2FC.Means", "Mean.Group1", "Mean.Group2", "PercNonZero.Group1", "PercNonZero.Group2")
  }

  if (Out == T) {
    FileName <- paste0("DEResults", ".csv")
    data.table::fwrite(DE.Res.df, file = FileName, row.names = F, col.names = T, sep = ",")
  }
  PiccoloList$DE.Results <- DE.Res.df
  return(PiccoloList)
}



#' @title  Differential expression analysis between cells in each cluster vs the rest.
#' @description  This function performs differential expression analysis between cells belonging to each cluster (as identified by Leiden or UMAP clustering) and the rest of the cells.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{LeidenClustering} or \link[Piccolo]{LouvainClustering} functions.
#' @param Transform A character variable. Should be the same as the one used during normalization. Else, the default of "log" will be used.
#' @param Method A character variable. Specifies the method used for testing differences between the expression levels in the 2 groups. Currently, two tests are available: the default Student's t-test ("t.test"), and the Wilcoxon rank-sum test ("wilcoxon").
#' @param Out A logical variable. Specifies whether to return an output files (.csv) with the differential expression results (if set to T) or not (set to F). Default is F.
#' @param verbose A logical variable. Specifies whether messages generated while running the function should be displayed (default is T).
#' @return An updated PiccoloList containing a data frame with the results of the differential expression analysis.
#' @examples
#' \dontrun{
#' pbmc3k <- PerformDiffExpClusterwise(PiccoloList = pbmc3k,
#' Out = T)
#' pbmc3k <- PerformDiffExpClusterwise(PiccoloList = pbmc3k,
#' Method = "wilcoxon", Out = T)
#' }
PerformDiffExpClusterwise <- function (PiccoloList,Transform = "log",Method = "t.test",Out = F,verbose = T){
  ClusterLevels <- table(PiccoloList$ClusterLabels)
  Cluster.DE <- vector(mode = "list", length = length(ClusterLevels))
  for (i in 1:length(ClusterLevels)){
    if (verbose == T) {
      message(paste0("Cluster ", i, " vs Others..."))
    }

    CellSerNos <- 1:length(PiccoloList$Barcodes)
    Group1.Vec <- which(PiccoloList$ClusterLabels == paste0("Cluster ",i))
    Group2.Vec <- CellSerNos[!CellSerNos %in% Group1.Vec]
    PiccoloList1 <- PerformDiffExp(PiccoloList = PiccoloList,Group1 = Group1.Vec,Group2 = Group2.Vec,Transform = Transform,Method = Method,Out = F)
    if (verbose == T) {
      message("Done.")
    }
    Cluster.DE[[i]] <- PiccoloList1$DE.Results
    if (Out == T) {
      FileName <- paste0("Cluster", i, "vsRest_DEResults", ".csv")
      data.table::fwrite(PiccoloList1$DE.Results, file = FileName, row.names = F, col.names = T, sep = ",")
    }
  }
  names(Cluster.DE) <- paste0("Cluster ", 1:length(ClusterLevels))
  PiccoloList$DE.Results <- Cluster.DE
  return(PiccoloList)
}

#' @title  Get marker genes for each cluster
#' @description  This function shortlists marker gene sets for each cluster.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after performing clustering (Leiden or Louvain).
#' @param MeanDiffSD A numeric variable. Specifies the minimum number of standard deviations that the difference of mean residuals between cluster cells and non-cluster cells should exhibit in order to be identified as a marker.
#' @param pValcutoff A numeric variable. Specifies the minimum raw (uncorrected) p-value that a gene must have to be identified as differentially expressed between the cluster cells and the non-cluster cells. Default is NULL.
#' @param FDRcutoff A numeric variable. Specifies the minimum FDR value that a gene must exhibit for it to be considered differentially expressed between the cluster cells and the non-cluster cells. Default is 0.005.
#' @param MaxSharedClusters An integer variable. The maximum number of clusters that a given gene may be identified as differentially expressed in. Default is 2.
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the lists of marker genes.
#' @return An updated PiccoloList containing a list with the differential expression result tables containing the marker genes for each cluster.
#' @examples
#' \dontrun{
#' pbmc3k <- GetClusterMarkers(PiccoloList = pbmc3k,
#' Out = T)
#' pbmc3k <- GetClusterMarkers(PiccoloList = pbmc3k,
#' MaxSharedClusters = 1,Out = T)
#' }
GetClusterMarkers <- function(PiccoloList,MeanDiffSD = 2.5,pValcutoff = NULL,FDRcutoff = 0.005,MaxSharedClusters = 2,Out = F){

  ClusterLevels <- table(PiccoloList$ClusterLabels)
  Mean.Diff <- vector(mode = "list",length = length(PiccoloList$DE.Results))
  for (i in 1:length(Mean.Diff)){
    Mean.Diff[[i]] <- PiccoloList$DE.Results[[i]]$mean.diff
  }
  SD.Overall <- sd(unlist(Mean.Diff))
  ClusterDEGenes.Up <- vector(mode = "list",length = length(ClusterLevels))
  ClusterGenes.Up <- vector(mode = "list",length = length(ClusterLevels))
  for (i in 1:length(ClusterLevels))
  {
    ClusterI <- paste0("Cluster ",i)

    DE.Results <- PiccoloList$DE.Results[[which(names(PiccoloList$DE.Results) == ClusterI)]]
    colnames(DE.Results) <- c("Gene.ID",colnames(DE.Results)[-1])

    if (is.null(pValcutoff) != T){
      RelevantUpGenes <- which(DE.Results$mean.diff > mean(DE.Results$mean.diff) + MeanDiffSD*SD.Overall & DE.Results$pvalue < pValcutoff)
      Up.Genes <- DE.Results[RelevantUpGenes,]
    } else {
      RelevantUpGenes <- which(DE.Results$mean.diff >= mean(DE.Results$mean.diff) + MeanDiffSD*SD.Overall & DE.Results$`p.adj(BH)` < FDRcutoff)
      Up.Genes <- DE.Results[RelevantUpGenes,]
    }
    ClusterDEGenes.Up[[i]] <- Up.Genes
    ClusterGenes.Up[[i]] <- Up.Genes$Gene.ID
  }

  Table.Cluster.Markers <- table(unlist(ClusterGenes.Up))

  Temp.Ser.Genes.To.Be.Removed <- which(Table.Cluster.Markers > MaxSharedClusters)

  if (length(Temp.Ser.Genes.To.Be.Removed) != 0){
    Genes.To.Be.Removed <- names(Table.Cluster.Markers)[Temp.Ser.Genes.To.Be.Removed]
  }

  for (k in 1:length(ClusterDEGenes.Up)){
    Temp.df <- ClusterDEGenes.Up[[k]]
    if (length(Genes.To.Be.Removed) != 0){
      Rows.To.Remove <- which(Temp.df$Gene.ID %in% Genes.To.Be.Removed)
      Temp.df <- Temp.df[-Rows.To.Remove,]
      if (nrow(Temp.df) == 0){
        warning("No marker genes shortlisted for ","Cluster ",k," with the given thresholds.")
      }
      ClusterDEGenes.Up[[k]] <- Temp.df
    }
  }

  names(ClusterDEGenes.Up) <- paste0("Cluster ",1:length(ClusterLevels))

  PiccoloList$ClusterMarkers <- ClusterDEGenes.Up

  if (Out == T){
    if (is.null(pValcutoff)){
      DirectoryName <- paste0("/ClusterMarkers","_MeanDiffSD",MeanDiffSD,"_FDR",FDRcutoff,"_","SharedClusters",MaxSharedClusters)
    } else {
      DirectoryName <- paste0("/ClusterMarkers","_MeanDiffSD",MeanDiffSD,"_pval",pValcutoff,"_","SharedClusters",MaxSharedClusters)
    }

    OrigDir <- getwd()

    dir.create(paste0(getwd(),DirectoryName))
    setwd(paste0(getwd(),DirectoryName))
    for (i in 1:length(PiccoloList$ClusterMarkers)){
      if (is.null(pValcutoff)){
        FileName <- paste0("Cluster",i,".csv")
      } else {
        FileName <- paste0("Cluster",i,".csv")
      }
      write.csv(PiccoloList$ClusterMarkers[[i]]$Gene.ID,file = FileName,row.names = F)
    }
    setwd(OrigDir)
  }

  return(PiccoloList)
}




