
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
utils::globalVariables(c("UMAP_1", "UMAP_2", "Label"))
#' @import utils
#' @import stats
#' @import ggplot2
#' @import network
#' @importFrom ggplot2 aes
#' @import ggnetwork

#' @title  ConvertToMTX
#' @description  This function converts .csv or .txt counts files to .mtx format. Prepares the corresponding features.tsv and barcodes.tsv files as well.
#' @export
#' @param X A character variable. Specifies the name of the file that contains the raw counts data (should be in .csv, .txt, or .tsv format). Make sure it is in the genes along rows and samples along columns format.
#' @return Generates an .mtx file (and features and barcodes files .tsv files) in the working directory that contains the input file.
#' @examples
#' \dontrun{
#' ConvertToMTX(X = "10X_PBMC3k.csv")
#' }

ConvertToMTX <- function(X){

  #Read in the UMI counts file
  UMI.Mat <- data.table::fread(X)
  colnames(UMI.Mat) <- c("Gene.ID",colnames(UMI.Mat)[-1])

  Gene.ID <- UMI.Mat$Gene.ID

  UMI.Mat <- as.matrix(UMI.Mat[,-1])

  rownames(UMI.Mat) <- Gene.ID

  #Make sparse matrix
  sparse.Mat <- Matrix::Matrix(UMI.Mat, sparse = T)

  FileNameMTX <- paste0(substr(X,1,nchar(X)-4),"_matrix.mtx.gz")

  #Nice writeMMgz function by Kamil Slowikowski
  writeMMgz <- function(x, file) {
    mtype <- "real"
    if (methods::is(x, "dgCMatrix")) {
      mtype <- "integer"
    }
    writeLines(
      c(
        sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
        sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
      ),
      gzfile(file)
    )
    data.table::fwrite(
      x = Matrix::summary(x),
      file = file,
      append = TRUE,
      sep = " ",
      row.names = FALSE,
      col.names = FALSE
    )
  }

  #write .mtx file
  suppressWarnings(writeMMgz(x = sparse.Mat, file=FileNameMTX))

  FileNameFeatures <- paste0(substr(X,1,nchar(X)-4),"_features.tsv")

  #write gene names file
  write(x = rownames(UMI.Mat), file = FileNameFeatures)


  FileNameBarcodes <- paste0(substr(X,1,nchar(X)-4),"_barcodes.tsv")

  #write barcodes file
  write(x = colnames(UMI.Mat), file = FileNameBarcodes)

}

#' @title  CombineMTX
#' @description  This function combines the single-cell counts files (.mtx or .mtx.gz) obtained from different samples/studies/conditions etc
#' @export
#' @param Path A character variable. Specifies the full path of the directory that contains the folders with the .mtx or .mtx.gz files that need to be combined into one file.
#' @param MinFeaturesPerCell A numeric variable. The minimum number of features with non-zero counts that any given cell is expected to contain. Cells with fewer than the specified number will be filtered out (Default is 200).
#' @param MT.Perc A numeric variable. The percentage of total count in any given cell attributed to counts of mitochondrial genes. Cells with MT.Perc greater than the specified number will be filtered out (Default is 10).
#' @param RP.Perc A numeric variable. The percentage of total count in any given cell attributed to counts of ribosomal genes. Cells with RP.Perc greater than the specified number will be filtered out (Default is 70).
#' @return Generates a matrix.mtx, a features.tsv, and a barcodes.tsv file in the same directory specified in the Path. The barcodes are modified to contain information about the folders. Eg: "Folder1Name_Barcode1", "Folder1Name_Barcode2",..,"Folder2Name_BarcodeX",...
#' @examples
#' \dontrun{
#'  CombineMTX(Path = "~/Documents/10X_PBMC_DataSets/") #Default
#'  CombineMTX(Path = "~/Documents/10X_PBMC_DataSets/",
#'  MinFeaturesPerCell = 100, MT.Perc = 20,RP.Perc = 90)
#' }
CombineMTX <- function(Path,MinFeaturesPerCell = 200, MT.Perc = 10,RP.Perc = 70){

  setwd(Path)

  FolderNames <- list.dirs()[-1]

  mat1 <- Matrix::Matrix(0,nrow = 2,ncol = 2)
  colnames(mat1) <- c("A","B")
  for (i in 1:length(FolderNames))
  {
    Folder <- gsub("./","",FolderNames[i])

    Files <- list.files(FolderNames[i])

    mtxcheck <- unique(grep(".mtx",x = Files,fixed = T),grep(".mtx.gz",x = Files,fixed = T))

    featuresbarcodescheck <- unique(grep(".tsv",x = Files,fixed = T),grep(".tsv.gz",x = Files,fixed = T))

    if (length(mtxcheck) == 1 & length(featuresbarcodescheck) == 2){
      barcodes.file <- Files[grep("barcodes",Files,fixed = T)]
      barcodes.path <- paste0(FolderNames[i],"/",barcodes.file)
      features.file <- Files[grep("features",Files,fixed = T)]
      if (length(features.file) == 0){
        features.file <- Files[grep("genes",Files,fixed = T)]
      }
      features.path <- paste0(FolderNames[i],"/",features.file)
      matrix.file <- Files[grep("matrix",Files,fixed = T)]
      matrix.path <- paste0(FolderNames[i],"/",matrix.file)

      message(paste0("Importing file ",i,"..."))

      mat <- Matrix::readMM(file = matrix.path)
      if (i == 1){
        mat1 <- Matrix::Matrix(0,nrow = nrow(mat),ncol = 2)
        colnames(mat1) <- c("A","B")
      }

      feature.names <- read.delim(features.path,header = FALSE,stringsAsFactors = FALSE)

      barcode.names <- read.delim(barcodes.path, header = FALSE,stringsAsFactors = FALSE)

      if (ncol(feature.names) > 1){
        List.For.GS.Col <- vector(mode = "list",length = ncol(feature.names))
        for (j in 1:ncol(feature.names))
        {
          List.For.GS.Col[[j]] <- which(substr(toupper(feature.names[,j]),1,3) == "RPL" | substr(toupper(feature.names[,j]),1,3) == "RPS")
        }

        GS.Col <- which(unlist(lapply(List.For.GS.Col,length)) != 0)

        if (length(GS.Col) == 1){
          Duplicated.Features <- which(duplicated(feature.names[,1]) == T)
          if (length(Duplicated.Features) != 0){
            mat <- Matrix::t(mat)
            mat <- mat[,-Duplicated.Features]
            mat <- Matrix::t(mat)
            feature.names <- feature.names[-Duplicated.Features,]
            rownames(mat) <- feature.names[,1]
            if (i == 1){
              mat1 <- Matrix::t(mat1)
              mat1 <- mat1[,-Duplicated.Features]
              mat1 <- Matrix::t(mat1)
            }
          } else {
            rownames(mat) <- feature.names[,1]
          }
        } else {
          warning("Features file does not contain gene symbols. MT genes and RP genes filtering will not be performed.")
          feature.names <- feature.names[,1]
          Duplicated.Features <- which(duplicated(feature.names[,1]) == T)
          if (length(Duplicated.Features) != 0){
            mat <- Matrix::t(mat)
            mat <- mat[,-Duplicated.Features]
            mat <- Matrix::t(mat)
            feature.names <- feature.names[-Duplicated.Features,]
            rownames(mat) <- feature.names[,1]
            if (i == 1){
              mat1 <- Matrix::t(mat1)
              mat1 <- mat1[,-Duplicated.Features]
              mat1 <- Matrix::t(mat1)
            }
          } else {
            rownames(mat) <- feature.names[,1]
          }
        }

      } else {
        feature.names <- feature.names[,1]
        Duplicated.Features <- which(duplicated(feature.names[,1]) == T)
        if (length(Duplicated.Features) != 0){
          mat <- Matrix::t(mat)
          mat <- mat[,-Duplicated.Features]
          mat <- Matrix::t(mat)
          feature.names <- feature.names[-Duplicated.Features,]
          rownames(mat) <- feature.names[,1]
          if (i == 1){
            mat1 <- Matrix::t(mat1)
            mat1 <- mat1[,-Duplicated.Features]
            mat1 <- Matrix::t(mat1)
            feature.names <- feature.names[-Duplicated.Features,]
            rownames(mat1) <- feature.names[,1]
          }
        } else {
          rownames(mat) <- feature.names[,1]
        }
      }

      if (nrow(mat1) != nrow(mat) & i != 1){
        if (nrow(mat1) > nrow(mat)){
          Matching.Features <- intersect(toupper(rownames(mat1)),toupper(rownames(mat)))
          Matching.Features.Mat1 <- which(toupper(rownames(mat1)) %in% Matching.Features)
          Matching.Features.Mat1 <- Matching.Features.Mat1[!is.na(Matching.Features.Mat1)]
          mat1 <- Matrix::t(mat1)
          mat1 <- mat1[,Matching.Features.Mat1]
          mat1 <- Matrix::t(mat1)
          Matching.Features.Mat <- which(toupper(rownames(mat)) %in% Matching.Features)
          Matching.Features.Mat <- Matching.Features.Mat[!is.na(Matching.Features.Mat)]
          mat <- Matrix::t(mat)
          mat <- mat[,Matching.Features.Mat]
          mat <- Matrix::t(mat)
        } else if (nrow(mat1) < nrow(mat)) {
          Matching.Features <- intersect(toupper(rownames(mat)),toupper(rownames(mat1)))
          Matching.Features.Mat <- which(toupper(rownames(mat)) %in% Matching.Features)
          Matching.Features.Mat <- Matching.Features.Mat[!is.na(Matching.Features.Mat)]
          mat <- Matrix::t(mat)
          mat <- mat[,Matching.Features.Mat]
          mat <- Matrix::t(mat)
          Matching.Features.Mat1 <- which(toupper(rownames(mat1)) %in% Matching.Features)
          Matching.Features.Mat1 <- Matching.Features.Mat1[!is.na(Matching.Features.Mat1)]
          mat1 <- Matrix::t(mat1)
          mat1 <- mat1[,Matching.Features.Mat1]
          mat1 <- Matrix::t(mat1)
        }

      }

      colnames(mat) <- paste0(Folder,"_",barcode.names[,1])

      rownames(mat1) <- rownames(mat)

      mat1 <- cbind(mat1,mat)
    }
  }

  feature.names <- feature.names[feature.names[,1] %in% rownames(mat1),]

  mat1 <- mat1[,-c(1,2)]

  Col.Sums.Vec <- Matrix::colSums(mat1)

  FeatureCounts.Per.Cell <- Matrix::diff(mat1@p)

  message("Filtering...")

  MT.Features <- grep("MT-",toupper(feature.names[,GS.Col]),fixed = T)
  if (length(MT.Features) != 0){
    MT.mat <- mat1[MT.Features,]
    MT.Col.Sums <- Matrix::colSums(MT.mat)
    MT.In.Prop.Total.Sum <- MT.Col.Sums/Col.Sums.Vec
  } else {
    warning("MT genes not detected in features.")
    MT.In.Prop.Total.Sum <- c()
  }

  RP.Features <- which(substr(toupper(feature.names[,GS.Col]),1,3) == "RPL" | substr(toupper(feature.names[,GS.Col]),1,3) == "RPS" |  substr(toupper(feature.names[,GS.Col]),1,3) %in% c("FAU","UBA52"))
  if (length(RP.Features) != 0){
    RP.mat <- mat1[RP.Features,]
    RP.Col.Sums <- Matrix::colSums(RP.mat)
    RP.In.Prop.Total.Sum <- RP.Col.Sums/Col.Sums.Vec
  } else {
    warning("RP genes not detected in features.")
    RP.In.Prop.Total.Sum <- c()
  }

  Cells.To.Remove <- unique(c(which(MT.In.Prop.Total.Sum > MT.Perc/100),which(RP.In.Prop.Total.Sum > RP.Perc/100),which(FeatureCounts.Per.Cell < MinFeaturesPerCell)))

  if (length(Cells.To.Remove) != 0){
    mat1 <- mat1[,-Cells.To.Remove]
  }

  #NonZeroProp.Per.Feature <- Matrix::diff(Matrix::t(mat1)@p)/ncol(mat1)

  #Nice writeMMgz function by Kamil Slowikowski
  writeMMgz <- function(x, file) {
    mtype <- "real"
    if (methods::is(x, "dgCMatrix")) {
      mtype <- "integer"
    }
    writeLines(
      c(
        sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
        sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
      ),
      gzfile(file)
    )
    data.table::fwrite(
      x = Matrix::summary(x),
      file = file,
      append = TRUE,
      sep = " ",
      row.names = FALSE,
      col.names = FALSE
    )
  }

  message(paste0("There are ",nrow(mat1)," features and ",ncol(mat1))," cells in the filtered matrix.")

  message("Writing the mtx and tsv files...")

  #write .mtx file
  suppressWarnings(writeMMgz(x = mat1, file="matrix.mtx.gz"))

  #write gene names file
  write.table(x = feature.names, file = "features.tsv",sep = "\t",row.names = F,col.names = F,quote = F)

  #write barcodes file
  write(x = colnames(mat1), file = "barcodes.tsv")

  message("Successfully generated the output files.")
}

#' @title  CreatePiccoloList Function
#' @description  This function creates a list object containing the counts matrix, the features (genes) list, and the barcodes
#' @export
#' @param X A character variable. Specifies the name of the .mtx or .mtx.gz file that contains the counts.
#' @param Gene A character variable. Specifies the name of the features (genes) file (.tsv format)
#' @param Barcode A character variable. Specifies the name of the barcodes file (.tsv format)
#' @param MinFeaturesPerCell A numeric variable. The minimum number of features with non-zero counts that any given cell is expected to contain. Cells with fewer than the specified number will be filtered out (Default is 200).
#' @param MT.Perc A numeric variable. The percentage of total count in any given cell attributed to counts of mitochondrial genes. Cells with MT.Perc greater than the specified number will be filtered out (Default is 10).
#' @param RP.Perc A numeric variable. The percentage of total count in any given cell attributed to counts of ribosomal genes. Cells with RP.Perc greater than the specified number will be filtered out (Default is 70).
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
CreatePiccoloList <-  function (X, Gene, Barcode, MinFeaturesPerCell = 200, MT.Perc = 10, RP.Perc = 70)
{
  message("Importing files...")
  UMI.Mat <- Matrix::readMM(file = X)
  UMI.Mat <- methods::as(UMI.Mat, "dgCMatrix")
  Gene.IDs <- read.delim(Gene, header = F, stringsAsFactors = F)
  Barcodes <- read.delim(Barcode, header = F, stringsAsFactors = F)
  Barcodes <- Barcodes$V1
  if (ncol(Gene.IDs) > 1) {
    List.For.GS.Col <- vector(mode = "list", length = ncol(Gene.IDs))
    for (j in 1:ncol(Gene.IDs)) {
      List.For.GS.Col[[j]] <- which(substr(toupper(Gene.IDs[,
                                                            j]), 1, 3) == "RPL" | substr(toupper(Gene.IDs[,
                                                                                                          j]), 1, 3) == "RPS")
    }
    GS.Col <- which(unlist(lapply(List.For.GS.Col, length)) !=
                      0)
    if (length(GS.Col) == 1) {
      Duplicated.Features <- which(duplicated(Gene.IDs[,
                                                       1]) == T)
      if (length(Duplicated.Features) != 0) {
        UMI.Mat <- Matrix::t(UMI.Mat)
        UMI.Mat <- UMI.Mat[, -Duplicated.Features]
        UMI.Mat <- Matrix::t(UMI.Mat)
        Gene.IDs <- Gene.IDs[-Duplicated.Features, ]
        rownames(UMI.Mat) <- Gene.IDs[, 1]
      }
      else {
        rownames(UMI.Mat) <- Gene.IDs[, 1]
      }
    }
    else {
      warning("Features file does not contain gene symbols. MT genes and RP genes filtering will not be performed.")
      Gene.IDs <- Gene.IDs[, 1]
      Duplicated.Features <- which(duplicated(Gene.IDs[,
                                                       1]) == T)
      if (length(Duplicated.Features) != 0) {
        UMI.Mat <- Matrix::t(UMI.Mat)
        UMI.Mat <- UMI.Mat[, -Duplicated.Features]
        UMI.Mat <- Matrix::t(UMI.Mat)
        Gene.IDs <- Gene.IDs[-Duplicated.Features, ]
        rownames(UMI.Mat) <- Gene.IDs[, 1]
      }
      else {
        rownames(UMI.Mat) <- Gene.IDs[, 1]
      }
    }
  }
  else {
    Gene.IDs <- Gene.IDs$V1
    List.For.GS.Col <- vector(mode = "list", length = length(Gene.IDs))
    for (j in 1:length(Gene.IDs)) {
      List.For.GS.Col[[j]] <- which(substr(toupper(Gene.IDs[j]), 1, 3) == "RPL" | substr(toupper(Gene.IDs[j]), 1, 3) == "RPS")
    }

    if (length(List.For.GS.Col) != 0){
      GS.Col <- 1
    }


    Duplicated.Features <- which(duplicated(Gene.IDs) ==
                                   T)
    if (length(Duplicated.Features) != 0) {
      UMI.Mat <- Matrix::t(UMI.Mat)
      UMI.Mat <- UMI.Mat[, -Duplicated.Features]
      UMI.Mat <- Matrix::t(UMI.Mat)
      Gene.IDs <- Gene.IDs[-Duplicated.Features, ]
      rownames(UMI.Mat) <- Gene.IDs
    }
    else {
      rownames(UMI.Mat) <- Gene.IDs
    }
  }
  Col.Sums.Vec <- Matrix::colSums(UMI.Mat)
  FeatureCounts.Per.Cell <- Matrix::diff(UMI.Mat@p)
  message("Filtering...")
  if (is.null(ncol(Gene.IDs)) != T) {
    MT.Features <- grep("MT-", toupper(Gene.IDs[, GS.Col]), fixed = T)
  } else {
    MT.Features <- grep("MT-", toupper(Gene.IDs), fixed = T)
  }

  if (length(MT.Features) != 0) {
    MT.mat <- UMI.Mat[MT.Features, ]
    MT.Col.Sums <- Matrix::colSums(MT.mat)
    MT.In.Prop.Total.Sum <- MT.Col.Sums/Col.Sums.Vec
  }
  else {
    warning("MT genes not detected in features.")
    MT.In.Prop.Total.Sum <- c()
  }

  if (is.null(ncol(Gene.IDs)) != T) {
    RP.Features <- which(substr(toupper(Gene.IDs[, GS.Col]),
                                1, 3) == "RPL" | substr(toupper(Gene.IDs[, GS.Col]),
                                                        1, 3) == "RPS" | substr(toupper(Gene.IDs[, GS.Col]),
                                                                                1, 3) %in% c("FAU", "UBA52"))
  } else {
    RP.Features <- which(substr(toupper(Gene.IDs),
                                1, 3) == "RPL" | substr(toupper(Gene.IDs),
                                                        1, 3) == "RPS" | substr(toupper(Gene.IDs),
                                                                                1, 3) %in% c("FAU", "UBA52"))
  }

  if (length(RP.Features) != 0) {
    RP.mat <- UMI.Mat[RP.Features, ]
    RP.Col.Sums <- Matrix::colSums(RP.mat)
    RP.In.Prop.Total.Sum <- RP.Col.Sums/Col.Sums.Vec
  }
  else {
    warning("RP genes not detected in features.")
    RP.In.Prop.Total.Sum <- c()
  }
  Cells.To.Remove <- unique(c(which(MT.In.Prop.Total.Sum >
                                      MT.Perc/100), which(RP.In.Prop.Total.Sum > RP.Perc/100),
                              which(FeatureCounts.Per.Cell < MinFeaturesPerCell)))

  if (length(Cells.To.Remove) != 0) {
    UMI.Mat <- UMI.Mat[, -Cells.To.Remove]
    Barcodes <- Barcodes[-Cells.To.Remove]
  }
  PiccoloList <- list(Counts = UMI.Mat, Genes = Gene.IDs, Barcodes = Barcodes)
  return(PiccoloList)
}

#' @title  Normalization and feature selection
#' @description  This function performs feature selection and prepares standardized counts matrix for the shortlisted features
#' @export
#' @param PiccoloList A list object. This should be the list created using the \link[Piccolo]{CreatePiccoloList} function
#' @param VarFeatures A numeric (integer) variable. The number of variable features to be shortlisted. If left NULL, will shortlist based on the threshold of ReferenceLevel
#' @param Transform A character variable. Specifies the non-linear transformation that will be applied to the counts. Currently, we offer the log transform (default) and Yeo-Johnson transform. These can be specified as "log" and "yj", respectively. The default is "log".
#' @param Batch An optional character vector. Specifies the batch labels for the cells. The order of the batch labels should match the order of the cells in the counts (or barcodes) file.
#' @param ReferenceLevel A numeric variable (value should be greater than 0 but less than 1). Specifies the reference level against which features are identified as variable. Default is the median (ReferenceLevel = 0.5).
#' @param MinPercNonZero A numeric variable. Specifies the minimum percentage of cells that must have non-zero counts for each gene in the data set. Default is 1 (\%).
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the standardized values (if set to T), or not (if set to F). Default is F
#' @return A list object containing the normalized counts and the variable features.
#' @examples
#' \dontrun{
#' pbmc3k <- Normalize(PiccoloList = pbmc3k)
#' pbmc3k <- Normalize(PiccoloList = pbmc3k,
#' ReferenceLevel = 0.3,MinPercNonZero = 0.5,Out = T)
#' }
Normalize <- function(PiccoloList,VarFeatures = NULL,Transform = "log", Batch=NULL,ReferenceLevel = NULL,MinPercNonZero = 1,Out = F){

  message("Importing files...")

  UMI.Mat <- PiccoloList$Counts

  Gene.IDs <- PiccoloList$Genes

  Barcodes <- PiccoloList$Barcodes

  TransformType <- Transform

  if(is.null(Batch)){

    Total.UMI.Counts <- Matrix::colSums(UMI.Mat)

    if(TransformType == "log"){
      message("Calculating size factors per cell...")

      LogTrans.TotalCounts <- log(Total.UMI.Counts)

      SF.Per.Cell <- LogTrans.TotalCounts/mean(LogTrans.TotalCounts)

    } else if(TransformType == "yj"){
      yeojohnson <- function (x, eps = 0.001)
      {
        stopifnot(is.numeric(x))
        lambda <- est_yj_lambda(x, eps = eps)
        x.t <- x
        na_idx <- is.na(x)
        x.t[!na_idx] <- yj_trans(x[!na_idx], lambda, eps)
        x.t
      }

      est_yj_lambda <- function (x, lower = -5, upper = 5, eps = 0.001)
      {
        n <- length(x)
        ccID <- !is.na(x)
        x <- x[ccID]
        yj_LL <- function(lambda) {
          x_t <- yj_trans(x, lambda, eps)
          x_t_bar <- mean(x_t)
          x_t_var <- var(x_t) * (n - 1)/n
          constant <- sum(sign(x) * log(abs(x) + 1))
          -0.5 * n * log(x_t_var) + (lambda - 1) * constant
        }
        results <- optimize(yj_LL, lower = lower, upper = upper,
                            maximum = TRUE, tol = 1e-04)
        results$maximum
      }

      yj_trans <- function (x, lambda, eps = 0.001)
      {
        pos_idx <- x >= 0
        neg_idx <- x < 0
        if (any(pos_idx)) {
          if (abs(lambda) < eps) {
            x[pos_idx] <- log(x[pos_idx] + 1)
          }
          else {
            x[pos_idx] <- ((x[pos_idx] + 1)^lambda - 1)/lambda
          }
        }
        if (any(neg_idx)) {
          if (abs(lambda - 2) < eps) {
            x[neg_idx] <- -log(-x[neg_idx] + 1)
          }
          else {
            x[neg_idx] <- -((-x[neg_idx] + 1)^(2 - lambda) -
                              1)/(2 - lambda)
          }
        }
        x
      }

      YJ.Trans.TotalCounts <- yeojohnson(Total.UMI.Counts)

      SF.Per.Cell <- YJ.Trans.TotalCounts/mean(YJ.Trans.TotalCounts)
    }

    Top.Features <- FeatureSelect(X = UMI.Mat,Gene = Gene.IDs,Reference = ReferenceLevel,Min.Perc.Non.Zero.Cells = MinPercNonZero)

    if (is.null(VarFeatures)){
      if (is.null(ncol(Top.Features)) != T){
        VarFeatures <- length(Top.Features[,1])
      } else {
        VarFeatures <- length(Top.Features)
      }
    }

    if (is.null(ncol(Top.Features)) != T){
      if (dim(Top.Features)[1] >= VarFeatures){
        Top.Features <- Top.Features[1:VarFeatures,]
      }
      Top.Features.Ser.Nos <- vector(mode = "numeric",length = dim(Top.Features)[1])
      for (i in 1:length(Top.Features.Ser.Nos)){
        Top.Features.Ser.Nos[i] <- which(Gene.IDs[,1] == Top.Features[i,1])
      }
    } else {
      if (length(Top.Features) >= VarFeatures){
        Top.Features <- Top.Features[1:VarFeatures]
      }
      Top.Features.Ser.Nos <- vector(mode = "numeric",length = length(Top.Features))
      for (i in 1:length(Top.Features.Ser.Nos)){
        Top.Features.Ser.Nos[i] <- which(Gene.IDs == Top.Features[i])
      }
    }

    write.csv(Top.Features,file = paste0("Top",length(Top.Features.Ser.Nos),"Features",".csv"),row.names = F)
    message("Successfully prepared .csv file containing list of highly variable features.")

    UMI.Mat <- Matrix::t(UMI.Mat)

    UMI.Mat <- UMI.Mat[,Top.Features.Ser.Nos]

    Std.Mat <- Standardize(X = UMI.Mat,Transform = TransformType,SF = SF.Per.Cell)

    if (Out == T){
      FileName <- paste0(Transform,"TransformedStandardizedCounts.csv")
      data.table::fwrite(data.frame(Std.Mat),file = FileName,row.names = F,col.names = F,sep = ",")
    }
    PiccoloList$NormCounts <- Std.Mat
    PiccoloList$VariableFeatures <- Top.Features
    return(PiccoloList)
  } else { #For batches
    stopifnot(length(Batch) == ncol(UMI.Mat))

    Batch <- as.factor(Batch)

    Zero.Count.Features <- vector(mode = "list",length = length(levels(Batch)))
    for(i in levels(Batch)){
      BatchIndex <- which(Batch == i)
      Temp.Mat <- UMI.Mat[,BatchIndex]
      Zero.Count.Features[[i]] <- which(Matrix::rowSums(Temp.Mat) == 0)
    }

    Zero.Count.Features <- unique(unlist(Zero.Count.Features))

    Zero.Count.Features <- Gene.IDs[Zero.Count.Features,1]

    Top.Features <- FeatureSelect(X = UMI.Mat,Gene = Gene.IDs,Reference  = ReferenceLevel,Min.Perc.Non.Zero.Cells = MinPercNonZero)

    if (is.null(VarFeatures)){
      if (is.null(ncol(Top.Features)) != T){
        VarFeatures <- length(Top.Features[,1])
      } else {
        VarFeatures <- length(Top.Features)
      }
    }

    if (is.null(ncol(Top.Features)) != T){
      Top.Features <- Top.Features[!Top.Features[,1] %in% Zero.Count.Features,]
    } else {
      Top.Features <- Top.Features[!Top.Features %in% Zero.Count.Features]
    }

    if (is.null(ncol(Top.Features)) != T){
      if (dim(Top.Features)[1] > VarFeatures){
        Top.Features <- Top.Features[1:VarFeatures,]
      }
      Top.Features.Ser.Nos <- vector(mode = "numeric",length = dim(Top.Features)[1])
      for (i in 1:length(Top.Features.Ser.Nos)){
        Top.Features.Ser.Nos[i] <- which(Gene.IDs[,1] == Top.Features[i,1])
      }
    } else {
      if (length(Top.Features) >= VarFeatures){
        Top.Features <- Top.Features[1:VarFeatures]
      }
      Top.Features.Ser.Nos <- vector(mode = "numeric",length = length(Top.Features))
      for (i in 1:length(Top.Features.Ser.Nos)){
        Top.Features.Ser.Nos[i] <- which(Gene.IDs == Top.Features[i])
      }
    }

    write.csv(Top.Features,file = paste0("Top",length(Top.Features.Ser.Nos),"Features",".csv"),row.names = F)
    message("Successfully prepared .csv file containing list of highly variable features.")

    Std.Mat <- matrix(0, nrow = length(Top.Features.Ser.Nos), ncol = ncol(UMI.Mat))

    for(i in levels(Batch)){
      BatchIndex <- which(Batch == i)
      Temp.Mat <- UMI.Mat[,BatchIndex]

      Total.UMI.Counts <- Matrix::colSums(Temp.Mat)

      if(TransformType == "log"){
        message("Calculating size factors per cell...")

        LogTrans.TotalCounts <- log(Total.UMI.Counts)

        SF.Per.Cell <- LogTrans.TotalCounts/mean(LogTrans.TotalCounts)

      } else if(TransformType == "yj"){
        yeojohnson <- function (x, eps = 0.001)
        {
          stopifnot(is.numeric(x))
          lambda <- est_yj_lambda(x, eps = eps)
          x.t <- x
          na_idx <- is.na(x)
          x.t[!na_idx] <- yj_trans(x[!na_idx], lambda, eps)
          x.t
        }

        est_yj_lambda <- function (x, lower = -5, upper = 5, eps = 0.001)
        {
          n <- length(x)
          ccID <- !is.na(x)
          x <- x[ccID]
          yj_LL <- function(lambda) {
            x_t <- yj_trans(x, lambda, eps)
            x_t_bar <- mean(x_t)
            x_t_var <- var(x_t) * (n - 1)/n
            constant <- sum(sign(x) * log(abs(x) + 1))
            -0.5 * n * log(x_t_var) + (lambda - 1) * constant
          }
          results <- optimize(yj_LL, lower = lower, upper = upper,
                              maximum = TRUE, tol = 1e-04)
          results$maximum
        }

        yj_trans <- function (x, lambda, eps = 0.001)
        {
          pos_idx <- x >= 0
          neg_idx <- x < 0
          if (any(pos_idx)) {
            if (abs(lambda) < eps) {
              x[pos_idx] <- log(x[pos_idx] + 1)
            }
            else {
              x[pos_idx] <- ((x[pos_idx] + 1)^lambda - 1)/lambda
            }
          }
          if (any(neg_idx)) {
            if (abs(lambda - 2) < eps) {
              x[neg_idx] <- -log(-x[neg_idx] + 1)
            }
            else {
              x[neg_idx] <- -((-x[neg_idx] + 1)^(2 - lambda) -
                                1)/(2 - lambda)
            }
          }
          x
        }

        YJ.Trans.TotalCounts <- yeojohnson(Total.UMI.Counts)

        SF.Per.Cell <- YJ.Trans.TotalCounts/mean(YJ.Trans.TotalCounts)
      }

      Temp.Mat <- Matrix::t(Temp.Mat)

      Temp.Mat <- Temp.Mat[,Top.Features.Ser.Nos]

      Std.Mat[,BatchIndex] <- Standardize(X = Temp.Mat,SF = SF.Per.Cell,Transform = TransformType)
    }

    if (Out == T){
      FileName <- paste0(Transform,"TransformedStandardizedCounts_BatchCorrected.csv")
      data.table::fwrite(data.frame(Std.Mat),file = FileName,row.names = F,col.names = T,sep = ",")

    }
    PiccoloList$NormCounts <- Std.Mat
    PiccoloList$VariableFeatures <- Top.Features
    return(PiccoloList)
  }
}


#Core functions
colVarsSPM <- function(X) {
  stopifnot( methods::is(X, "dgCMatrix"))
  ans <- sapply( base::seq.int(X@Dim[2]),function(j) {
    if(X@p[j+1] == X@p[j]) { return(0) } # all entries are 0: var is 0
    mean <- base::sum(X@x[ (X@p[j]+1):X@p[j+1]])/X@Dim[1]
    sum((X@x[(X@p[j]+1):X@p[j+1] ] - mean)^2) +
      mean^2 * (X@Dim[1] - (X@p[j+1] - X@p[j]))})/(X@Dim[1] - 1)
  names(ans) <- X@Dimnames[[2]]
  ans
}

colOverdispQPCoef <- function(X,alternative = "greater"){
  stopifnot( methods::is(X,"dgCMatrix"))
  ans <- sapply( base::seq.int(X@Dim[2]),function(j){
    if(X@p[j+1] == X@p[j]){return(0)} # all entries are 0: var is 0
    #mean <- exp(sum(log(X@x[(X@p[j]+1):X@p[j+1]]+1))/X@Dim[1]) - 1
    est.mean <- sum(X@x[(X@p[j]+1):X@p[j+1]])/X@Dim[1]

    aux <- c(((X@x[(X@p[j]+1):X@p[j+1]] - est.mean)^2 -
                X@x[(X@p[j]+1):X@p[j+1]]),rep(est.mean^2,X@Dim[1] - length(X@x[(X@p[j]+1):X@p[j+1]])))/est.mean

    mean(aux) + 1})
}

#Function to shortlist highly variable features
FeatureSelect <- function(X,Gene,Reference,Min.Perc.Non.Zero.Cells){

  if (is.null(Reference)){
    Reference <- 0.5
  } else if (Reference > 0.5){
    message("Reference should not be greater than 0.5. Resetting it to 0.5 (default)")
    Reference <- 0.5
  }

  UMI.Mat <- Matrix::t(X)

  Gene.IDs <- Gene

  Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
  Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)

  message("Filtering features...")
  Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)

  No.of.Non.Zero.Per.Row <- diff(UMI.Mat@p)

  Perc.Non.Zero.Per.Row <- No.of.Non.Zero.Per.Row/nrow(UMI.Mat) * 100

  Irrelevant.Features <- unique(c(Irrelevant.Features,which(Perc.Non.Zero.Per.Row <= Min.Perc.Non.Zero.Cells)))

  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (is.null(ncol(Gene.IDs)) != T){
      Gene.IDs <- Gene.IDs[-Irrelevant.Features,]
    } else {
      Gene.IDs <- Gene.IDs[-Irrelevant.Features]
    }
    Mean.Arith.Per.Feature <- Mean.Arith.Per.Feature[-Irrelevant.Features]
    Var.Arith.Per.Feature <- Var.Arith.Per.Feature[-Irrelevant.Features]
    No.of.Non.Zero.Per.Row <- No.of.Non.Zero.Per.Row[-Irrelevant.Features]
    Perc.Non.Zero.Per.Row <- Perc.Non.Zero.Per.Row[-Irrelevant.Features]
  }

  message("Estimating overdispersion coefficients...")
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
    Perc.Non.Zero.Per.Row <- Perc.Non.Zero.Per.Row[-Irrelevant.Features]
    No.of.Non.Zero.Per.Row <- No.of.Non.Zero.Per.Row[-Irrelevant.Features]
    Alpha.QP <- Alpha.QP[-Irrelevant.Features]
  }

  Alpha.NB.Est <- (Alpha.QP - 1)/Mean.Arith.Per.Feature

  #Binning based approach

  message("Shortlisting variable features...")

  Mean.Quantiles <- quantile(Mean.Arith.Per.Feature,probs = seq(0.001,1,0.001))
  Diff.AlphaQP.AlphaQPFit <- vector(mode = "numeric",length = length(Gene.IDs))
  for (i in 1:length(Mean.Quantiles))
  {
    if (i == 1){
      Features.In.Bin <- which(Mean.Arith.Per.Feature <= Mean.Quantiles[i])
    } else {
      Features.In.Bin <- which(Mean.Arith.Per.Feature > Mean.Quantiles[i-1] & Mean.Arith.Per.Feature <= Mean.Quantiles[i])
    }

    Median.AlphaQP.Bin <- quantile(Alpha.QP[Features.In.Bin],probs = c(Reference))
    Diff.AlphaQP.AlphaQPFit[Features.In.Bin] <- Alpha.QP[Features.In.Bin] - Median.AlphaQP.Bin
  }

  #Identify top features
  Default.Features <- which(Diff.AlphaQP.AlphaQPFit > 0)
  Top.Features <- Default.Features[order(Diff.AlphaQP.AlphaQPFit[Default.Features],decreasing = T)]

  if (is.null(ncol(Gene.IDs)) != T){
    Gene.IDs <- Gene.IDs[Top.Features,]
  } else {
    Gene.IDs <- Gene.IDs[Top.Features]
  }

  return(Gene.IDs)
}

#Function to prepare standardized counts
Standardize <- function(X,Transform,SF){

  if (Transform == "log"){

    #For log transform
    message("Log transforming counts...")

    Seq.Nos <- seq(1,ncol(X),500)

    Trans.Mat <- methods::as(matrix(0,nrow = 2,ncol = nrow(X)),"dgCMatrix")
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
      #Res.List <- vector(mode = "list",length = ncol(Temp.UMI.Mat))

      for (j in 1:ncol(Temp.UMI.Mat))
      {
        Temp.mat[j,] <- log1p(Temp.UMI.Mat[,j])
      }
      Temp.mat <- methods::as(Temp.mat,"dgCMatrix")
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  } else if (Transform == "yj"){

    #Create Yeo-Johnson transform function for sparse dgCMatrix
    colYJLambda <- function(X,lower = -5, upper = 5, eps = 0.001){
      stopifnot( methods::is(X,"dgCMatrix"))
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

    message("Estimating lambdas for power transform...")

    Lambda <- colYJLambda(X)

    Seq.Nos <- seq(1,ncol(X),500)

    Trans.Mat <- methods::as(matrix(0,nrow = 2,ncol = nrow(X)),"dgCMatrix")
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
      #Res.List <- vector(mode = "list",length = ncol(Temp.UMI.Mat))

      for (j in 1:ncol(Temp.UMI.Mat))
      {
        UMI.Vec <- Temp.UMI.Mat[,j] + 1

        TempJ <- Start.Point + (j-1)

        Temp.mat[j,] <- ((UMI.Vec^Lambda[TempJ]) - 1)/Lambda[TempJ]

      }
      Temp.mat <- methods::as(Temp.mat,"dgCMatrix")
      Trans.Mat <- rbind(Trans.Mat,Temp.mat)
      if (i == 1){
        Trans.Mat <- Trans.Mat[-c(1,2),]
      }
    }
  }

  #No cell-size scaling
  #SF.Per.Cell <- rep(1,nrow(UMI.Mat))

  Trans.Mat <- Matrix::t(Trans.Mat)

  #Function to standardize
  StandardizeCounts <- function(x,sf){
    EstimatedMeans <- mean(x)*sf
    #Alpha <- (sum((x - EstimatedMeans)^2)/(length(x)-1))/EstimatedMeans
    #(x-EstimatedMeans)/sqrt(EstimatedMeans*Alpha)#QP
    (x-EstimatedMeans)/sd(x-EstimatedMeans)
  }

  message("Preparing standardized counts matrix...")

  #Sequence of numbers separated by 500
  Seq.Nos <- seq(1,ncol(X),500)

  Std.Mat <- matrix(0,nrow = 2,ncol = nrow(X))
  for (i in 1:length(Seq.Nos))
  {
    if (i < length(Seq.Nos)){
      Start.Point <- Seq.Nos[i]
      End.Point <- Seq.Nos[i+1] - 1
    } else {
      Start.Point <- Seq.Nos[i]
      End.Point <- ncol(Trans.Mat)
    }

    Temp.mat <- as.matrix(Trans.Mat[,Start.Point:End.Point])
    Temp.Std.Mat <- matrix(0,nrow = ncol(Temp.mat),ncol = nrow(Temp.mat))
    for (j in 1:ncol(Temp.mat))
    {
      UMI.Vec <- Temp.mat[,j]

      TempJ <- Start.Point + (j-1)

      Temp.Std.Mat[j,] <- StandardizeCounts(x = UMI.Vec,sf = SF)
    }
    Std.Mat <- rbind(Std.Mat,Temp.Std.Mat)
    if (i == 1){
      Std.Mat <- Std.Mat[-c(1,2),]
    }
  }
  return(Std.Mat)
}

#' @title  Max-Min Normalization Function
#' @description  This function will normalize the standardized values to normalized values in the range 0-1
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T). Default is F.
#' @return A numeric matrix containing the normalized values.
#' @examples
#' \dontrun{
#' NormMat <- MaxMinNormMat(PiccoloList = pbmc3k,
#' Out = F)
#' }
MaxMinNormMat <- function(PiccoloList,Out = F){
  Std.Mat <- PiccoloList$NormCounts

  Norm.Mat <- t(apply(Std.Mat,1,function(x) (x- min(x))/(max(x) - min(x))))
  if (Out == T){
    Norm.df <- data.frame(Norm.Mat)
    FileName <- paste0("MaxMinNormMat.csv")
    data.table::fwrite(Norm.df,file = FileName,row.names = F,col.names = F,sep = ",")
  }
  return(Norm.Mat)
}

#' @title  Identify differentially expressed between 2 groups of cells
#' @description  This function performs differential expression analysis (using the Mann-Whitney test) between 2 groups of cells provided by the user
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function
#' @param Group1 A numeric (integers) vector. Specifies the serial numbers of cells in group 1 (serial numbers based on the order of barcodes)
#' @param Group2 A numeric (integers) vector. Specifies the serial numbers of cells in group 2 (serial numbers based on the order of barcodes)
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the differential expression result (if set to T). Default is F.
#' @return A data frame containing the gene IDs, the log2 fold change (FC) of normalized values between group 1 and group 2 (positive FC indicates higher expression in group 1), the p-values from the Mann-Whitney test, and the adjusted p-values (p-adj) after Benjamini-Hochberg correction
#' @examples
#' \dontrun{
#' Group1.vec <- 1:200
#' Group2.vec <- 300:500
#' DE.Genes.df <- DEfeatures(PiccoloList = pbmc3k,
#' Group1 = Group1.vec,
#' Group2 = Group2.vec,
#' Out = T)
#' }
DEfeatures <- function(PiccoloList,Group1,Group2,Out = F){
  Std.Mat <- PiccoloList$NormCounts

  Features <- PiccoloList$VariableFeatures

  Barcodes <- PiccoloList$Barcodes

  Norm.Mat <- t(apply(Std.Mat,1,function(x) (x- min(x))/(max(x) - min(x))))

  log2FC.vec <- rep(1,nrow(Norm.Mat))
  p.val.vec <- rep(1,nrow(Norm.Mat))
  Base.Mean.vec <- rep(1,nrow(Norm.Mat))
  for (i in 1:length(p.val.vec))
  {
    log2FC.vec[i] <- log2(mean(Norm.Mat[i,Group1])/mean(Norm.Mat[i,Group2]))
    Base.Mean.vec[i] <- mean(Norm.Mat[i,])
    p.val.vec[i] <- wilcox.test(Norm.Mat[i,Group1],Norm.Mat[i,Group2])$p.val
  }

  p.val.vec <- p.val.vec[order(log2FC.vec,decreasing = T)]
  p.adj.vec <- p.adjust(p.val.vec,method = "BH")
  if (length(dim(Features)) > 1){
    Genes <- Features[order(log2FC.vec,decreasing = T),]
  } else {
    Genes <- Features[,order(log2FC.vec,decreasing = T)]
  }

  Base.Mean.vec <- Base.Mean.vec[order(log2FC.vec,decreasing = T)]
  log2FC.vec <- log2FC.vec[order(log2FC.vec,decreasing = T)]

  DE.Res.df <- data.frame(Genes,Base.Mean.vec,log2FC.vec,p.val.vec,p.adj.vec)

  if (Out == T){
    FileName <- paste0("DEGenes",".csv")
    data.table::fwrite(DE.Res.df,file = FileName,row.names = F,col.names = F,sep = ",")
  }
  return(DE.Res.df)
}

#' @title  Compute Principal Components
#' @description  This function will calculate the principal components for the matrix with the standardized values obtained from the \link[Piccolo]{Normalize} function.
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function
#' @param NoOfPC A numeric (integer) variable. No of principal components to be retained in the output. Default is 50.
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T), or not (if set to F). Default is T
#' @return A numeric matrix containing the standardized values.
#' @examples
#' \dontrun{
#' pbmc3k <- ComputePC(PiccoloList = pbmc3k)
#' pbmc3k <- ComputePC(PiccoloList = pbmc3k,NoOfPC = 20,Out = T)
#' }
ComputePC <- function(PiccoloList,NoOfPC =  50,Out = F){
  Std.Mat <- PiccoloList$NormCounts

  Features <- PiccoloList$VariableFeatures

  if (length(dim(Features)) > 1){
    Features <- Features[,1]
  }

  Barcodes <- PiccoloList$Barcodes

  rownames(Std.Mat) <- Features
  colnames(Std.Mat) <- Barcodes

  res.pca <- stats::prcomp(t(Std.Mat))

  if (Out == T){
    PC.df <- data.frame(res.pca$x[,1:NoOfPC])
    FileName <- paste0("Top",NoOfPC,"PrinComp",".csv")
    data.table::fwrite(PC.df,file = FileName,row.names = T,col.names = F,sep = ",")
  }
  PiccoloList$PCmat <- res.pca$x[,1:NoOfPC]
  return(PiccoloList)
}

#' @title  UMAP coordinates
#' @description  This function will run UMAP for the matrix with the principal components
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} and the \link[Piccolo]{ComputePC} functions
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T), or not (if set to F). Default is T
#' @return A data frame containing the coordinates of the cells in the first 2 UMAP dimensions.
#' @examples
#' \dontrun{
#' pbmc3k <- UMAPcoords(PiccoloList = pbmc3k,
#' Out = T)
#' }
UMAPcoords <- function(PiccoloList,Out = F){

  PC.Mat <- as.matrix(PiccoloList$PCmat)

  x <- umap::umap(PC.Mat)

  UMAP.df <- data.frame(CellID = rownames(x$layout),UMAP1  = x$layout[,1], UMAP2 = x$layout[,2])

  if (Out == T){
    FileName <- c("UMAPcoords.csv")
    data.table::fwrite(UMAP.df,file = FileName,row.names = F,sep = ",")
  }
  PiccoloList$UMAPcoord <- UMAP.df
  return(PiccoloList)
}

#' @title  Label cells on UMAP
#' @description  This function can be used to label (color) cells on the UMAP plot
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize}, the \link[Piccolo]{ComputePC}, and the \link[Piccolo]{UMAPcoords} functions
#' @param Labels A character vector. Should contain the character labels for cells in the same order as the cells in the counts matrix
#' @param Levels An optional character vector. Should specify the unique labels in the order in which the labels of the cells should be presented
#' @param Alpha A numeric variable (strictly greater than 0). Specified the transparency of the dots on the UMAP plot. Default Alpha = 0.7
#' @param Size A numeric variable (strictly greater than 0). Specified the size of the dots on the UMAP plot. Default Size = 0.9
#' @return A ggplot2 object
#' @examples
#' \dontrun{
#' p <- LabelUMAP(PiccoloList = pbmc3k,
#' Labels = c("b-cells","b-cells",..,"cd14 monocytes",..,"NK cells",..),
#' Levels = c("b-cells","cd14 monocytes","dendritic","NK cells","naive cytotoxic"))
#' }

LabelUMAP <- function(PiccoloList,Labels,Levels = NULL,Alpha = 0.7,Size = 0.9){

  UMAP.Coord.df <- PiccoloList$UMAPcoord

  if (length(Labels) != length(UMAP.Coord.df$CellID)){
    stop("The length of the Labels vector provided does not match the number of cells in the UMAP.")
  }

  plot_data <- data.frame(UMAP.Coord.df,Labels)
  colnames(plot_data) <- c("CellID","UMAP_1","UMAP_2","Label")

  if (is.null(Levels) == F){
    #Specify the factor levels in the order you want
    plot_data$Label <- factor(plot_data$Label, levels = Levels)
  }

  p <- ggplot2::ggplot(plot_data,ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(ggplot2::aes(color=Label),alpha = Alpha,size = Size) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(p)
}

#' @title  GenePairsSingleCell function
#' @description  This function identifies the gene-pairs with significant overlaps between their percentile sets
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function
#' @param PercSetSize A numeric variable. Specifies the size of the top percentile set (in percent). Default value is 5.
#' @param JcdInd A numeric variable. Specifies the minimum extent of the overlap between percentile sets of any given gene-pair.
#' @param Stop A numeric variable. Specifies the stop point for the sliding window. For use as TuBA, leave this NULL (default).
#' @param Out A logical variable. If set to T (default is F), will prepare .csv files containing the information about the gene-pairs, cells in percentile sets etc.
#' @return A list object containing the gene-pairs data frame, the genes percentile sets data frame, and a data frame with cell IDs.
#' @examples
#' \dontrun{
#' pbmc3k <- GenePairsSingleCell(PiccoloList = pbmc3k,PercSetSize = 5,JcdInd = 0.3)
#' }

#Function that generates p-values
GenePairsSingleCell <- function(PiccoloList,PercSetSize = NULL,JcdInd,Stop = NULL,Out = F)
{
  Norm.Mat <- PiccoloList$NormCounts

  Sample.IDs <- PiccoloList$Barcodes

  Gene.Names <- PiccoloList$VariableFeatures
  if (is.null(ncol(Gene.Names)) != T){
    Gene.Names <- Gene.Names[,1]
  }

  if (is.null(PercSetSize)){
    PercSetSize <- 5
  }

  if(is.null(Stop)){
    Stop <- PercSetSize
  } else if (Stop > 30){
    message("Stop point can not be larger than 30. Setting it to 30 (maximum permitted).")
    Stop <- 30
  }

  CutOffPerc <- PercSetSize/100

  WindowShift <- PercSetSize/2

  Stop.Point <- 100 - Stop

  Stop.Point <- round(Stop.Point/100 * ncol(Norm.Mat))

  Window.Indices <- c(ncol(Norm.Mat))
  k <- 2
  while (Window.Indices[k-1] >= Stop.Point){
    Window.Indices <- c(Window.Indices,Window.Indices[k-1] - ceiling(CutOffPerc*ncol(Norm.Mat)/2))
    k <- k + 1
  }

  List.With.Window.Ranges <- vector(mode = "list")
  k <- 1
  Start.Point <- 1
  while(Start.Point+2 <= length(Window.Indices)){
    List.With.Window.Ranges[[k]] <- c(Start.Point,Start.Point+2)
    Start.Point <- Start.Point + 1
    k <- k + 1
  }

  message("Finding gene-pairs...")

  NonInformativeSamples <- vector(mode = "list",length(List.With.Window.Ranges))
  Scores.Mat <- matrix("0",nrow = nrow(Norm.Mat),ncol = nrow(Norm.Mat))
  SerNos.Reverse <- length(List.With.Window.Ranges):1
  for (i in 1:length(List.With.Window.Ranges))
  {
    #Binary matrix with 1 for samples that are in the quantile set for any given gene
    List.Samples.In.Window <- vector(mode = "list",length = nrow(Norm.Mat))
    Row.Index.List <- vector(mode = "list",length = nrow(Norm.Mat))
    Values.List <- vector(mode = "list",length = nrow(Norm.Mat))
    for (j in 1:nrow(Norm.Mat))
    {
      Samples.In.Order <-  order(Norm.Mat[j,])

      Samples.In.Window <- Samples.In.Order[Window.Indices[List.With.Window.Ranges[[i]][1]]:(Window.Indices[List.With.Window.Ranges[[i]][2]]+1)]

      if (length(Samples.In.Window) != 0){
        List.Samples.In.Window[[j]] <- Samples.In.Window
        Row.Index.List[[j]] <- rep(j,length(Samples.In.Window))
        Values.List[[j]] <- rep(1,length(Samples.In.Window))
      }
    }

    Binary.Mat.For.Genes.Window <- Matrix::sparseMatrix(i = unlist(Row.Index.List),j = unlist(List.Samples.In.Window),x = unlist(Values.List))

    #Find the frequencies for the samples
    Sample.Frequencies <- Matrix::colSums(Binary.Mat.For.Genes.Window)

    #Identify samples that only show up in one percentile set
    ExclusionSamples <- which(Sample.Frequencies == 1)
    if (length(ExclusionSamples) != 0){
      Binary.Mat.For.Genes.Window <- Binary.Mat.For.Genes.Window[,-ExclusionSamples]
      NonInformativeSamples[[i]] <- c(NonInformativeSamples[[i]],ExclusionSamples)
    }

    #Find the threshold for maximum frequency based on percentile set size
    MaxThreshold <- max(Sample.Frequencies)
    n <- choose(n = nrow(Binary.Mat.For.Genes.Window),k = 2)
    a <- choose(n = MaxThreshold,k = 2)
    p <- a/n
    while (p > CutOffPerc) {
      MaxThreshold <- MaxThreshold -1
      a <- choose(n = MaxThreshold,k = 2)
      p <- a/n
    }

    ExclusionSamples <- which(Sample.Frequencies > MaxThreshold)
    if (length(ExclusionSamples) != 0){
      Binary.Mat.For.Genes.Window <- Binary.Mat.For.Genes.Window[,-ExclusionSamples]
      NonInformativeSamples[[i]] <- c(NonInformativeSamples[[i]],ExclusionSamples)
    }

    SampleOverlaps.Matrix <- Matrix::tcrossprod(Binary.Mat.For.Genes.Window,Binary.Mat.For.Genes.Window)

    SampleOverlaps.Matrix <- as.matrix(SampleOverlaps.Matrix)

    SampleUnion <- 2*ceiling(CutOffPerc*ncol(Norm.Mat))
    SampleUnion.Matrix <- SampleUnion - SampleOverlaps.Matrix

    rm(Binary.Mat.For.Genes.Window)

    Jaccard.Ind.Mat <- SampleOverlaps.Matrix/SampleUnion.Matrix
    rm(SampleUnion.Matrix)
    rm(SampleOverlaps.Matrix)

    #Find Jaccard distance between gene-pairs
    Relevant.Gene.Pairs <- which(Jaccard.Ind.Mat >= JcdInd,arr.ind = T)
    rm(Jaccard.Ind.Mat)

    if (length(Relevant.Gene.Pairs) != 0){
      Relevant.Gene.Pairs <- Relevant.Gene.Pairs[which(Relevant.Gene.Pairs[,1] < Relevant.Gene.Pairs[,2]),]
      Scores.Mat[Relevant.Gene.Pairs] <- paste0(as.character(length(List.With.Window.Ranges) - SerNos.Reverse[i]),"/",Scores.Mat[Relevant.Gene.Pairs])
    }
  }

  Scores.Summary <- table(Scores.Mat[Scores.Mat != "0"])

  Names.Scores <- names(Scores.Summary)

  Names.Scores <- substr(Names.Scores,1,nchar(Names.Scores)-1)

  List.Scores <- vector(mode = "list",length = length(Names.Scores))
  for (i in 1:length(Names.Scores))
  {
    List.Scores[[i]] <- unlist(strsplit(Names.Scores[i],"/",fixed = T))
  }

  if (length(List.Scores) == 1){
    Relevant.Indices <- paste0(Names.Scores,"0")
  } else {
    Threshold.Vec <- rep(0,length(List.With.Window.Ranges))
    for (i in 1:length(List.With.Window.Ranges))
    {
      CharI <- as.character(i-1)
      TempVec <- c()
      for (j in 1:length(List.Scores))
      {
        TempNo <- which(List.Scores[[j]] %in% CharI)
        if (length(TempNo) != 0){
          TempVec <- c(TempVec,j)
        }
      }
      Threshold.Vec[i] <- sum(Scores.Summary[TempVec])
    }

    Relevant.Ind <- paste0("/",which(Threshold.Vec < 2*Threshold.Vec[1])-1,"/")

    Matches.List <- vector(mode = "list",length(Relevant.Ind))
    for (i in 1:length(Matches.List))
    {
      Matches.List[[i]] <- grep(Relevant.Ind[i],x = Names.Scores)
    }

    Matches.Vec <- unique(unlist(Matches.List))

    Relevant.Indices <- paste0(Names.Scores[Matches.Vec],"0")

  }

  if (length(Relevant.Indices) < length(names(Scores.Summary))){
    Relevant.Gene.Pairs <- arrayInd(which(Scores.Mat %in% Relevant.Indices),dim(Scores.Mat))
    Relevant.Scores <- Scores.Mat[Relevant.Gene.Pairs]
    Relevant.Scores <- substr(Relevant.Scores,1,nchar(Relevant.Scores)-2)
  } else {
    Relevant.Gene.Pairs <- arrayInd(which(Scores.Mat != "0"),dim(Scores.Mat))
    Relevant.Scores <- Scores.Mat[Relevant.Gene.Pairs]
    Relevant.Scores <- substr(Relevant.Scores,1,nchar(Relevant.Scores)-2)
  }

  if (length(Relevant.Gene.Pairs) == 0){
    stop("No gene-pairs found for given choice of JcdInd.")
  }

  Gene.Pairs.Col1 <- Relevant.Gene.Pairs[,1]

  Gene.Pairs.Col2 <- Relevant.Gene.Pairs[,2]

  Relevant.Feature.Ser.Nos <- unique(c(Gene.Pairs.Col1,Gene.Pairs.Col2))

  Relevant.Feature.Names <- Gene.Names[Relevant.Feature.Ser.Nos]

  Relevant.Scores.Mat <- Scores.Mat[Relevant.Feature.Ser.Nos,Relevant.Feature.Ser.Nos]

  rm(Scores.Mat)

  Gene.Pairs.df <- data.frame(Gene.Pairs.Col1,Gene.Pairs.Col2,Relevant.Scores,rep(CutOffPerc,length(Relevant.Scores)),rep(JcdInd,length(Relevant.Scores)))
  colnames(Gene.Pairs.df) <- c("Gene.1","Gene.2","Window.Indices","PercSet","JcdInd")
  Gene.Pairs.df$Window.Indices <- as.character(Gene.Pairs.df$Window.Indices)

  message("Preparing gene-pairs file...")

  if (Out == T){
    FileName <- paste0("GenePairsAndWindowIndices","_H",CutOffPerc,"_JcdInd",JcdInd,".csv")
    write.table(Gene.Pairs.df,file = FileName, row.names = F,col.names = T,sep = ",")
    message(paste0("Successfully generated .csv file containing edgelist with ",length(Relevant.Scores)," gene-pairs. Maximum window index is ",Gene.Pairs.df$Window.Indices[which(nchar(Gene.Pairs.df$Window.Indices) == max(nchar(Gene.Pairs.df$Window.Indices)))[1]]))
  }

  #Make the reduced scores matrix symmetrical
  Relevant.Scores.Mat[lower.tri(Relevant.Scores.Mat)] = t(Relevant.Scores.Mat)[lower.tri(Relevant.Scores.Mat)]

  Max.Score.Indices <- which(nchar(Relevant.Scores.Mat) == max(nchar(Relevant.Scores.Mat)),arr.ind = T)[1,]

  Max.Relevant.Score.Exp <- as.numeric(unlist(strsplit((Relevant.Scores.Mat[Max.Score.Indices[1],Max.Score.Indices[2]]),"/",fixed = T))[1])

  message("Preparing samples info...")

  Matrix.With.Collated.Samples <- matrix(0,nrow = nrow(Norm.Mat),ncol = (Max.Relevant.Score.Exp+1))
  for (i in 1:nrow(Relevant.Scores.Mat))
  {
    TempI <- Relevant.Feature.Ser.Nos[i]
    Samples.In.Order <-  order(Norm.Mat[TempI,])
    Temp.Max.Index <- which(nchar(Relevant.Scores.Mat[i,]) == max(nchar(Relevant.Scores.Mat[i,])))[1]
    if (Temp.Max.Index != 0){
      Max.Index <- as.numeric(unlist(strsplit(Relevant.Scores.Mat[i,Temp.Max.Index],"/",fixed = T))[1]) + 1
      for (j in 1:Max.Index)
      {
        Samples.In.Window <- Samples.In.Order[Window.Indices[List.With.Window.Ranges[[j]][1]]:(Window.Indices[List.With.Window.Ranges[[j]][2]]+1)]

        NonInformativeSamplesJ <- NonInformativeSamples[[j]]
        if (length(NonInformativeSamplesJ) != 0){
          Samples.In.Window <- Samples.In.Window[!Samples.In.Window %in% NonInformativeSamplesJ]
        }

        if (length(Samples.In.Window) != 0){
          Matrix.With.Collated.Samples[TempI,j] <- paste(Samples.In.Window,collapse = "/")
        }
      }
    } else {
      Matrix.With.Collated.Samples[TempI,] <- 0
    }
  }

  Col.Names <- paste0("W",as.character(1:(Max.Relevant.Score.Exp+1)))

  colnames(Matrix.With.Collated.Samples) <- Col.Names

  Sample.IDs.df <- data.frame(Sample.ID = Sample.IDs)

  if (Out == T){
    FileName <- paste0("CellIDs","_H",CutOffPerc,"_JcdInd",JcdInd,".csv")
    write.csv(Sample.IDs.df,file = FileName,row.names = F)
  }

  #Prepare Genes-Samples window indices matrix
  Genes.Samples.df <- data.frame(Gene.Names,Matrix.With.Collated.Samples)
  colnames(Genes.Samples.df) <- c("Gene.ID",colnames(Matrix.With.Collated.Samples))

  if (Out == T){
    FileName <- paste0("GenesSamplesMatrix","_H",CutOffPerc,"_JcdInd",JcdInd,".csv")
    write.table(Genes.Samples.df,file = FileName,row.names = F,col.names = T,sep = ",")
    message("Successfully generated .csv file containing the genes-sample window indices matrix.")
  }

  PiccoloList$GenePairs <- Gene.Pairs.df
  PiccoloList$CellIDsForBiclustering <- Sample.IDs.df
  PiccoloList$GenesPercentileSets <- Genes.Samples.df

  return(PiccoloList)

}

#' @title  GenePairsSingleCell function
#' @description  This function identifies the gene-pairs with significant overlaps between their percentile sets
#' @export
#' @param PiccoloList A list object. Piccolo list object obtained after applying the \link[Piccolo]{Normalize} function
#' @param Window A numeric (integer) variable. Specifies which window to perform biclustering on. For TuBA, it is set to 1 (default).
#' @param MinGenes A numeric variable. Specifies the minimum number of genes that a bicluster should contain. Cannot be less than 3 (default).
#' @param MinSamples A numeric variable. Specifies the minimum number of samples that a bicluster must contain. If left NULL (default), no minimum limit is imposed.
#' @param SampleEnrichment A numeric variable. Specifies the enrichment level that samples must exhibit in order to belong to a bicluster.
#' @param Out A logical variable. If set to T (default is F), will prepare a .csv file containing the bicluster results.
#' @return A list object containing the bicluster results data frame.
#' @examples
#' \dontrun{
#' pbmc3k <- Bicluster(PiccoloList = pbmc3k,SampleEnrichment = 0.05)
#' }

Bicluster <- function(PiccoloList,Window = 1,MinGenes = 3,MinSamples = NULL,SampleEnrichment = NULL,Out = F){
  if(MinGenes < 3){
    MinGenes <- 3
    message("MinGenes cannot be less than 3. Minimum of 3 set as default.")
  }

  if(is.null(MinSamples)){
    MinSamples <- 2
  }

  #Import data frame that contains the Node pairs found by the significant node pairs function
  Node.Pairs.df <- PiccoloList$GenePairs
  colnames(Node.Pairs.df) <- c("Node.1","Node.2","WndInd","CutOffPerc","JcdInd")

  Node.Pairs.df$WndInd <- as.character(Node.Pairs.df$WndInd)

  MaxInd <- which(nchar(Node.Pairs.df$WndInd) == max(nchar(Node.Pairs.df$WndInd)))[1]

  MaxWindowSize <- as.numeric(unlist(strsplit(Node.Pairs.df$WndInd[MaxInd],"/",fixed = T))[1]) + 1

  if (Window > MaxWindowSize){
    stop(paste0("Window size choice exceeds maximum value of ",MaxWindowSize, " in the input file."))
  }

  List.WndInd <- vector(mode = "list",length = length(Node.Pairs.df$WndInd))
  for (i in 1:length(List.WndInd))
  {
    List.WndInd[[i]] <- as.numeric(unlist(strsplit(Node.Pairs.df$WndInd[i],"/",fixed = T)))
  }

  JaccardInd <- round(min(Node.Pairs.df$JcdInd),digits = 2)

  message("Importing files..")

  #Import data frame that contains the window indices matrix
  Nodes.Samples.Window.df <- PiccoloList$GenesPercentileSets
  colnames(Nodes.Samples.Window.df) <- c("Node.ID",colnames(Nodes.Samples.Window.df[,-1]))

  #Names of Nodes
  Node.Names <- as.character(Nodes.Samples.Window.df$Node.ID)

  #Binary matrix of Nodes and percentile set samples
  Matrix.For.Samples.Windows <- as.matrix(Nodes.Samples.Window.df[,-1])

  rm(Nodes.Samples.Window.df)

  No.of.Samples.In.Windows <- as.numeric(unlist(lapply(sapply(as.character(Matrix.For.Samples.Windows[,Window]),function(x) strsplit(x,split = "/",fixed = T)),length)))

  N.PercentileSet <- max(No.of.Samples.In.Windows)

  Eligible.Column.For.Window <- which(colnames(Matrix.For.Samples.Windows) == paste0("W",Window))

  List.Samples.Per.Node <- sapply(as.character(Matrix.For.Samples.Windows[,Window]),function(x) strsplit(x,split = "/",fixed = T))

  #IDs of samples (or conditions)
  Sample.IDs <- PiccoloList$CellIDsForBiclustering$Sample.ID

  Eligible.NodePairs <- which(sapply(List.WndInd,function(x) x[order(x)][Window]) == Window-1)

  #Column 1 of node pairs data frame with eligible window indices
  Qualified.Node1 <- Node.Pairs.df$Node.1[Eligible.NodePairs]

  #Column 2 of node pairs data frame with eligible window indices
  Qualified.Node2 <- Node.Pairs.df$Node.2[Eligible.NodePairs]

  if(length(Qualified.Node1) == 0){
    stop("Please check input file. No node pairs found.")
  } else {
    message(paste0("There are a total of ",length(Qualified.Node1)," edges in the graph."))
  }

  #Summarize graph info in a table
  Node.Summary.Info <- table(c(Qualified.Node1,Qualified.Node2))

  #All Nodes in graph
  All.Nodes.In.Graph <- as.numeric(names(Node.Summary.Info))

  rm(Node.Pairs.df)

  List.Samples.Per.Node <- List.Samples.Per.Node[All.Nodes.In.Graph]

  message("Preparing graph...")

  #Vector to convert node serial numbers to the corresponding serial number in the reduced adjacency matrix
  Node.Annotation.Conversion.Vec <- vector(mode = "numeric", length = length(Node.Names))
  Node.Ser.Nos <- 1:length(Node.Names)
  Matching.Nodes.Indices <- match(All.Nodes.In.Graph,Node.Ser.Nos)
  Node.Annotation.Conversion.Vec[Matching.Nodes.Indices] <- 1:length(All.Nodes.In.Graph)

  #Prepare adjacency matrix for the graph
  Adjacency.Mat.Nodes <- matrix(0,nrow = length(All.Nodes.In.Graph),ncol = length(All.Nodes.In.Graph))
  Adjacency.Mat.Nodes[cbind(Node.Annotation.Conversion.Vec[Qualified.Node1],Node.Annotation.Conversion.Vec[Qualified.Node2])] <- 1

  #Make the adjacency matrix symmetric
  Adjacency.Mat.Nodes <- Adjacency.Mat.Nodes + t(Adjacency.Mat.Nodes)
  rownames(Adjacency.Mat.Nodes) <- All.Nodes.In.Graph
  colnames(Adjacency.Mat.Nodes) <- All.Nodes.In.Graph

  message("Finding biclusters...")

  if (length(Qualified.Node1) > 200000)
    message("This may take several minutes due to the large size of the graph.")

  Total.No.of.Edges.In.Unpruned.Graph <- length(Qualified.Node1)

  Nodes.Per.Sample <- Nodes.Per.Sample <- vector(mode = "list",length = length(Sample.IDs))
  for (i in 1:length(All.Nodes.In.Graph))
  {
    Node.Samples <- as.numeric(List.Samples.Per.Node[[i]])
    for (j in 1:length(Node.Samples))
    {
      Nodes.Per.Sample[[Node.Samples[j]]] <- c(Nodes.Per.Sample[[Node.Samples[j]]],All.Nodes.In.Graph[i])
    }
  }

  #Find sample background counts along edges in graph
  Samples.Background.Frequencies <- vector(mode = "numeric", length = length(Sample.IDs))
  for (i in 1:length(Samples.Background.Frequencies))
  {
    Temp.Nodes.Per.Sample <- Nodes.Per.Sample[[i]]
    if (length(Temp.Nodes.Per.Sample) >= 2){
      Sub.Adj.Mat <- Adjacency.Mat.Nodes[Node.Annotation.Conversion.Vec[Temp.Nodes.Per.Sample],Node.Annotation.Conversion.Vec[Temp.Nodes.Per.Sample]]
      Samples.Background.Frequencies[i] <- sum(Sub.Adj.Mat)/2
    } else{
      Samples.Background.Frequencies[i] <- 0
    }
  }

  #Find the number of triangular cliques associated with each node-pair in graph
  No.of.Nodes.Associated <- vector(mode = "numeric",length = length(Qualified.Node1))
  for (i in 1:length(Qualified.Node1))
  {
    TempI1 <- Qualified.Node1[i]
    TempI2 <- Qualified.Node2[i]

    Col.Sums.Vec <- colSums(Adjacency.Mat.Nodes[c(Node.Annotation.Conversion.Vec[TempI1],Node.Annotation.Conversion.Vec[TempI2]),])
    No.of.Nodes.Associated[i] <- length(which(Col.Sums.Vec == 2))
  }

  #Non-triangular gene-pairs
  Node.Pairs.For.Filtering <- which(No.of.Nodes.Associated == 0)

  if (length(Node.Pairs.For.Filtering) != 0){
    No.of.Nodes.Associated <- No.of.Nodes.Associated[-Node.Pairs.For.Filtering]

    Qualified.Node1 <- Qualified.Node1[-Node.Pairs.For.Filtering]
    Qualified.Node2 <- Qualified.Node2[-Node.Pairs.For.Filtering]
  }

  if(max(No.of.Nodes.Associated) == 0){
    stop("No clique of size 3 found in input graph.")
  }

  Decreasing.No.of.Nodes <- order(No.of.Nodes.Associated,decreasing = T)

  #Sort column 1 and column2 based on decreasing order of nodes associated
  Qualified.Node1 <- Qualified.Node1[Decreasing.No.of.Nodes]
  Qualified.Node2 <- Qualified.Node2[Decreasing.No.of.Nodes]

  #Find dense subgraphs
  i <- 1
  Temp.Vec1 <- Qualified.Node1
  Temp.Vec2 <- Qualified.Node2
  Nodes.In.Subgraph <- vector(mode = "list")
  Temp.Nodes.Vec <- c()
  while (length(Temp.Vec1) != 0){
    Sub.Adj.Mat <- Adjacency.Mat.Nodes[c(Node.Annotation.Conversion.Vec[Temp.Vec1[1]],Node.Annotation.Conversion.Vec[Temp.Vec2[1]]),]
    Nodes.In.Subgraph[[i]] <- c(Temp.Vec1[1],Temp.Vec2[1],All.Nodes.In.Graph[colSums(Sub.Adj.Mat) == 2])
    Nodes.In.Subgraph[[i]] <- Nodes.In.Subgraph[[i]][!Nodes.In.Subgraph[[i]] %in% Temp.Nodes.Vec]
    Temp.Nodes.Vec <- c(Temp.Nodes.Vec,Nodes.In.Subgraph[[i]])

    #Remove those edges that contain the nodes in the dense subgraph
    Edges.To.Be.Removed <- which(Temp.Vec1 %in% Nodes.In.Subgraph[[i]] | Temp.Vec2 %in% Nodes.In.Subgraph[[i]])

    Temp.Vec1 <- Temp.Vec1[-Edges.To.Be.Removed]
    Temp.Vec2 <- Temp.Vec2[-Edges.To.Be.Removed]

    i <- i + 1
  }

  if (length(Nodes.In.Subgraph) != 0){
    message("Dense subgraphs identified.")
  } else {
    stop("No dense subgraph with at least 3 nodes identified.")
  }

  #Reintroduce dense subgraphs back in original graph and add nodes that share edges with at least 2 nodes in the dense subgraphs
  Nodes.In.Bicluster <- vector(mode = "list",length = length(Nodes.In.Subgraph))
  for (i in 1:length(Nodes.In.Subgraph))
  {
    Col.Ser.Nos.For.Bicluster <- Node.Annotation.Conversion.Vec[Nodes.In.Subgraph[[i]]]
    Sub.Adj.Mat <- Adjacency.Mat.Nodes[,Col.Ser.Nos.For.Bicluster]
    Nodes.In.Bicluster[[i]] <- unique(c(Nodes.In.Subgraph[[i]],as.numeric(rownames(Sub.Adj.Mat)[which(rowSums(Sub.Adj.Mat) >= 2)])))
  }

  #Find biclusters that have nodes that are subsets of other biclusters
  Nested.Biclusters <- vector(mode = "numeric")
  for (i in 1:length(Nodes.In.Bicluster))
  {
    TempI <- length(Nodes.In.Bicluster) - (i-1)
    Temp.Bicluster.I <- Nodes.In.Bicluster[[TempI]]
    j <- 1
    Temp.Intersection <- 1
    while(Temp.Intersection != 0 & TempI != j){
      Temp.Bicluster.J <- Nodes.In.Bicluster[[j]]
      Temp.Intersection <- length(Temp.Bicluster.I) - length(intersect(Temp.Bicluster.I,Temp.Bicluster.J))
      if (Temp.Intersection == 0 & TempI != j & length(Temp.Bicluster.J) >= length(Temp.Bicluster.I)){
        Nested.Biclusters <- c(Nested.Biclusters,TempI)
      }
      j <- j + 1
    }
  }

  #Remove biclusters that are nested within larger biclusters
  if (length(Nested.Biclusters) != 0){
    Nodes.In.Bicluster <- Nodes.In.Bicluster[-Nested.Biclusters]
  }

  if (length(Nodes.In.Bicluster) > 1){
    Bicluster.Size.Order <- order(unlist(lapply(Nodes.In.Bicluster,length)),decreasing = T)
  } else if (length(Nodes.In.Bicluster) < 2 & length(Nodes.In.Bicluster[[1]]) != 0){
    Bicluster.Size.Order <- 1
  } else {
    stop("No bicluster found with given choice of parameters.")
  }

  Nodes.In.Bicluster <- Nodes.In.Bicluster[Bicluster.Size.Order]

  #Find samples preferentially associated with biclusters
  n <- Total.No.of.Edges.In.Unpruned.Graph
  Samples.In.Bicluster <- vector(mode = "list",length = length(Nodes.In.Bicluster))
  for (i in 1:length(Nodes.In.Bicluster))
  {
    Temp.Nodes.In.Bicluster <- Nodes.In.Bicluster[[i]]

    Sub.Adj.Mat <- Adjacency.Mat.Nodes[Node.Annotation.Conversion.Vec[Temp.Nodes.In.Bicluster],Node.Annotation.Conversion.Vec[Temp.Nodes.In.Bicluster]]
    No.of.Edges.In.Bicluster <- sum(Sub.Adj.Mat)/2

    Bicluster.Samples.Frequencies <- vector(mode = "numeric",length = length(Sample.IDs))
    for (j in 1:length(Bicluster.Samples.Frequencies))
    {
      Temp.Nodes.Per.Sample <- Nodes.Per.Sample[[j]]
      if (length(Temp.Nodes.Per.Sample) != 0){
        Nodes.Per.Sample.In.Bicluster <- intersect(Temp.Nodes.Per.Sample,Temp.Nodes.In.Bicluster)
      } else {
        Nodes.Per.Sample.In.Bicluster <- NULL
      }
      if (length(Nodes.Per.Sample.In.Bicluster) > 0){
        Sub.Adj.Mat <- Adjacency.Mat.Nodes[Node.Annotation.Conversion.Vec[Nodes.Per.Sample.In.Bicluster],Node.Annotation.Conversion.Vec[Nodes.Per.Sample.In.Bicluster]]
        Bicluster.Samples.Frequencies[j] <- sum(Sub.Adj.Mat)/2
      } else {
        Bicluster.Samples.Frequencies[j] <- 0
      }
    }

    Valid.Sample.Ser.Nos <- which(Bicluster.Samples.Frequencies != 0)
    Bicluster.Samples.Frequencies <- Bicluster.Samples.Frequencies[Valid.Sample.Ser.Nos]
    Order.Bicluster.Samples.Frequencies <- order(Bicluster.Samples.Frequencies,decreasing = T)
    Bicluster.Samples.Frequencies <- Bicluster.Samples.Frequencies[Order.Bicluster.Samples.Frequencies]
    Valid.Sample.Ser.Nos <- Valid.Sample.Ser.Nos[Order.Bicluster.Samples.Frequencies]

    if (is.null(SampleEnrichment)){
      Sample.Ratios.Per.Bicluster <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      for (k in 1:length(Valid.Sample.Ser.Nos))
      {
        t <- Bicluster.Samples.Frequencies[k]
        Sample.Count.In.Graph <- Samples.Background.Frequencies[Valid.Sample.Ser.Nos[k]]

        Sample.Ratios.Per.Bicluster[k] <- (t/Sample.Count.In.Graph)/(No.of.Edges.In.Bicluster/n)
      }

      Temp.ser.nos <- which(Sample.Ratios.Per.Bicluster >= 1)
      Samples.In.Bicluster[[i]] <- Valid.Sample.Ser.Nos[Temp.ser.nos]
    } else {
      p.value.sample <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      Sample.Ratios.Per.Bicluster <- vector(mode = "numeric",length = length(Valid.Sample.Ser.Nos))
      for (k in 1:length(Valid.Sample.Ser.Nos))
      {
        Temp.Sample.Frequency <- Bicluster.Samples.Frequencies[k]
        t <- Temp.Sample.Frequency

        Sample.Count.In.Graph <- Samples.Background.Frequencies[Valid.Sample.Ser.Nos[k]]
        if (No.of.Edges.In.Bicluster <= Sample.Count.In.Graph){
          b <- No.of.Edges.In.Bicluster
          a <- Sample.Count.In.Graph
        } else {
          a <- No.of.Edges.In.Bicluster
          b <- Sample.Count.In.Graph
        }
        p.value.sample[k] <- sum(dhyper(t:b,a,n-a,b))
        Sample.Ratios.Per.Bicluster[k] <- (t/Sample.Count.In.Graph)/(No.of.Edges.In.Bicluster/n)
      }
      Temp.ser.nos <- which(p.value.sample <= SampleEnrichment & Sample.Ratios.Per.Bicluster >= 1)
      Samples.In.Bicluster[[i]] <- Valid.Sample.Ser.Nos[Temp.ser.nos]
    }
  }

  #Filter out biclusters that have fewer genes than the specified threshold (MinGenes)
  if (length(which(unlist(lapply(Nodes.In.Bicluster,length)) >= MinGenes)) == 0){
    stop("No biclusters with ",MinGenes," genes found. Please choose different parameters.")
  } else {
    SatisfactoryMinSizeGenes.Biclusters <- which(unlist(lapply(Nodes.In.Bicluster,length)) >= MinGenes)
    Nodes.In.Bicluster <- Nodes.In.Bicluster[SatisfactoryMinSizeGenes.Biclusters]
    Samples.In.Bicluster <- Samples.In.Bicluster[SatisfactoryMinSizeGenes.Biclusters]
  }

  #Find biclusters that contain at least the minimum of samples specified (MinSamples)
  Biclusters.With.Some.Samples <- which(unlist(lapply(Samples.In.Bicluster,length)) >= MinSamples)

  if (length(Biclusters.With.Some.Samples) == 0 & is.null(SampleEnrichment) == T){
    stop("No biclusters were found with the given choice of SampleEnrichment. Try with SampleEnrichment = NULL or SampleEnrichment = 1")
  } else if (length(Biclusters.With.Some.Samples) == 0 & (SampleEnrichment == 1 | is.null(SampleEnrichment) == F)){
    stop("No bicluster found with at least 3 genes")
  } else {
    message("Samples enriched in biclusters identified.")
  }

  #Filter nodes in biclusters based on the samples found enriched in each bicluster - Approach 1 (Remove nodes based on the overlap of their percentile sets with samples found enriched in the bicluster)
  Nodes.In.Bicluster <- Nodes.In.Bicluster[Biclusters.With.Some.Samples]
  Samples.In.Bicluster <- Samples.In.Bicluster[Biclusters.With.Some.Samples]

  #Filter nodes based on their association with samples in subgraphs
  Nodes.In.Final.Biclusters <- vector(mode = "list",length = length(Nodes.In.Bicluster))
  for (i in 1:length(Nodes.In.Bicluster))
  {
    Temp.Nodes.In.Bicluster <- Nodes.In.Bicluster[[i]]
    JInd <- vector(mode = "numeric",length = length(Temp.Nodes.In.Bicluster))
    for (j in 1:length(Temp.Nodes.In.Bicluster))
    {
      TempJ <- which(All.Nodes.In.Graph == Temp.Nodes.In.Bicluster[j])
      Samples.In.WindowJ <- as.numeric(List.Samples.Per.Node[[TempJ]])
      Intersecting.Samples <- intersect(Samples.In.Bicluster[[i]],Samples.In.WindowJ)
      JInd[j] <- length(Intersecting.Samples)/(2*N.PercentileSet -  length(Intersecting.Samples))
    }
    Temp.Filter.Index <- which(JInd < JaccardInd)
    if (length(Temp.Filter.Index != 0)){
      Nodes.In.Final.Biclusters[[i]] <- Nodes.In.Bicluster[[i]][-Temp.Filter.Index]
    } else {
      Nodes.In.Final.Biclusters[[i]] <- Nodes.In.Bicluster[[i]]
    }
  }

  #Find biclusters that have nodes that are subsets of other biclusters
  Nested.Biclusters <- vector(mode = "numeric")
  for (i in 1:length(Nodes.In.Final.Biclusters))
  {
    TempI <- length(Nodes.In.Final.Biclusters) - (i-1)
    Temp.Bicluster.I <- Nodes.In.Final.Biclusters[[TempI]]
    j <- 1
    Temp.Intersection <- 1
    while(Temp.Intersection != 0 & TempI != j){
      Temp.Bicluster.J <- Nodes.In.Final.Biclusters[[j]]
      Temp.Intersection <- length(Temp.Bicluster.I) - length(intersect(Temp.Bicluster.I,Temp.Bicluster.J))
      if (Temp.Intersection == 0 & TempI != j & length(Temp.Bicluster.J) >= length(Temp.Bicluster.I)){
        Nested.Biclusters <- c(Nested.Biclusters,TempI)
      }
      j <- j + 1
    }
  }

  #Remove biclusters that are nested within larger biclusters
  if (length(Nested.Biclusters) != 0){
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[-Nested.Biclusters]
    Samples.In.Bicluster <- Samples.In.Bicluster[-Nested.Biclusters]
  }

  if (length(Nodes.In.Final.Biclusters) > 1){
    Bicluster.Size.Order <- order(unlist(lapply(Nodes.In.Final.Biclusters,length)),decreasing = T)
  } else if (length(Nodes.In.Final.Biclusters) < 2 & length(Nodes.In.Final.Biclusters[[1]]) != 0){
    Bicluster.Size.Order <- 1
  } else {
    stop("No bicluster found with given choice of parameters.")
  }

  Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[Bicluster.Size.Order]
  Samples.In.Bicluster <- Samples.In.Bicluster[Bicluster.Size.Order]

  #Filter out biclusters that have fewer genes than the specified threshold (MinGenes)
  if (length(which(unlist(lapply(Nodes.In.Final.Biclusters,length)) >= MinGenes)) == 0){
    stop("No biclusters with ",MinGenes," genes found. Please choose different parameters.")
  } else {
    SatisfactoryMinSizeGenes.Biclusters <- which(unlist(lapply(Nodes.In.Final.Biclusters,length)) >= MinGenes)
    Nodes.In.Final.Biclusters <- Nodes.In.Final.Biclusters[SatisfactoryMinSizeGenes.Biclusters]
    Samples.In.Bicluster <- Samples.In.Bicluster[SatisfactoryMinSizeGenes.Biclusters]
  }

  message("Nodes in final biclusters identified.")

  #Make bicluster-samples matrix
  Bicluster.Samples.List <- vector(mode = "list",length = length(Samples.In.Bicluster))
  for (i in 1:length(Samples.In.Bicluster))
  {
    Bicluster.Samples.List[[i]] <- paste(Samples.In.Bicluster[[i]],collapse = "/")
  }

  if (length(lapply(Bicluster.Samples.List,length)) != 0){
    MinBiclusterSamples <- min(unlist(lapply(Samples.In.Bicluster,length)))
    No.of.Samples <- length(Sample.IDs)
  } else {
    MinBiclusterSamples <- 0
  }

  Temp.No.of.Nodes.In.Bicluster <- unlist(lapply(Nodes.In.Final.Biclusters,length))
  if (length(Temp.No.of.Nodes.In.Bicluster) != 0){
    MinBiclusterGenes <- min(Temp.No.of.Nodes.In.Bicluster)
  } else {
    MinBiclusterGenes <- 0
  }

  OriginalMinSamples <- MinSamples
  if (MinBiclusterSamples > MinSamples)
    MinSamples <- MinBiclusterSamples

  if (OriginalMinSamples != 2){
    message(paste0("Found ",length(Nodes.In.Final.Biclusters)," biclusters with at least ",MinGenes," genes and ",OriginalMinSamples," samples."))
  } else {
    message(paste0("Found ",length(Nodes.In.Final.Biclusters)," biclusters with at least ",MinGenes," genes and ",MinSamples," samples."))
  }

  #Prepare the output files
  if (length(Nodes.In.Final.Biclusters) != 0){

    Nodes.Biclusters.Info.df <- data.frame(Node.Names[unlist(Nodes.In.Final.Biclusters)],rep(1:length(Nodes.In.Final.Biclusters),unlist(lapply(Nodes.In.Final.Biclusters,length))))
    colnames(Nodes.Biclusters.Info.df) <- c("Gene.ID","Bicluster.No")

    No.of.Samples.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Samples.Per.Gene.In.Bicluster <- vector(mode = "list",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    No.of.Samples.Per.Gene.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Proportion.of.Samples.In.Bicluster <- vector(mode = "numeric",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Genes.Bicluster.Samples.List <- vector(mode = "list",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    Samples.In.Bicluster.Vec <- vector(mode = "character",length = length(Nodes.Biclusters.Info.df$Gene.ID))
    for (i in 1:length(Nodes.Biclusters.Info.df$Gene.ID))
    {
      Gene.Ser.No <- which(Node.Names == Nodes.Biclusters.Info.df$Gene.ID[i])
      Bicluster.No <- Nodes.Biclusters.Info.df$Bicluster.No[i]

      No.of.Samples.In.Bicluster[i] <- length(Samples.In.Bicluster[[Bicluster.No]])
      Samples.Per.Gene.In.Bicluster[[i]] <- intersect(List.Samples.Per.Node[[which(All.Nodes.In.Graph == Gene.Ser.No)]],Samples.In.Bicluster[[Bicluster.No]])
      Samples.In.Bicluster.Vec[i] <- Bicluster.Samples.List[[Bicluster.No]]
      No.of.Samples.Per.Gene.In.Bicluster[i] <- length(Samples.Per.Gene.In.Bicluster[[i]])
      Proportion.of.Samples.In.Bicluster[i] <- No.of.Samples.Per.Gene.In.Bicluster[i]/No.of.Samples.In.Bicluster[i]
      Genes.Bicluster.Samples.List[[i]] <- paste(Samples.Per.Gene.In.Bicluster[[i]],collapse = "/")
    }

    Nodes.Biclusters.Info.df$Total.Samples.In.Bicluster <- No.of.Samples.In.Bicluster
    Nodes.Biclusters.Info.df$Samples.In.Bicluster <- Samples.In.Bicluster.Vec
    Nodes.Biclusters.Info.df$Total.Samples.Per.Gene<- No.of.Samples.Per.Gene.In.Bicluster
    Nodes.Biclusters.Info.df$Samples.Per.Gene.In.Bicluster <- unlist(Genes.Bicluster.Samples.List)
    Nodes.Biclusters.Info.df$Proportion.of.Samples <- Proportion.of.Samples.In.Bicluster

    if (Out == T){

      File.Name <- paste0("MinGenes",MinBiclusterGenes,"_MinSamples",MinBiclusterSamples,"_GenesInBiclusters","_Window",Window,".csv")
      write.table(Nodes.Biclusters.Info.df,file = File.Name,row.names = F,col.names = T,sep = ",")
    }
    PiccoloList$BiclusterOutput <- Nodes.Biclusters.Info.df
  }
  return(PiccoloList)
}


