
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

#' @title  Standardization Function
#' @description  This function performs feature selection and prepares standardized counts matrix for the shortlisted features
#' @export
#' @param X A character variable. Specifies the name of the .mtx or .mtx.gz file that contains the counts.
#' @param Gene A character variable. Specifies the name of the features (genes) file (.tsv format)
#' @param Barcode A character variable. Specifies the name of the barcodes file (.tsv format)
#' @param VarFeatures A numeric (integer) variable. The number of variable features to be shortlisted.
#' @param Transform A character variable. Specifies the non-linear transformation that will be applied to the counts. Currently, we offer the log transform (default) and Yeo-Johnson transform. These can be specified as "log" and "yj", respectively. The default is "log".
#' @param Batch An optional character vector. Specifies the batch labels for the cells. The order of the batch labels should match the order of the cells in the counts (or barcodes) file.
#' @param ReferenceLevel A numeric variable (value should be greater than 0 but less than 1). Specifies the reference level against which features are identified as variable. Default is the median (ReferenceLevel = 0.5).
#' @param MinPercNonZero A numeric variable. Specifies the minimum percentage of cells that must have non-zero counts for each gene in the data set. Default is 1 (\%).
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the standardized values (if set to T), or not (if set to F). Default is F
#' @return A numeric matrix containing the standardized values.
#' @examples
#' \dontrun{
#'  StandardizeMat(X = "10X_PBMC3k_matrix.mtx.gz",
#'  Gene = "10X_PBMC3k_features.tsv",
#'  BinaryMatrix = "10X_PBMC3k_barcodes.tsv",Out = T)
#'  StandardizeMat(X = "10X_PBMC3k_matrix.mtx.gz",
#'  Gene = "10X_PBMC3k_features.tsv",
#'  BinaryMatrix = "10X_PBMC3k_barcodes.tsv",
#'  ReferenceLevel = 0.3,MinPercNonZero = 0.5,Out = T)
#' }
StandardizeMat <- function(X,Gene,Barcode,VarFeatures = NULL,Transform = "log", Batch=NULL,ReferenceLevel = NULL,MinPercNonZero = 1,Out = F){

  message("Importing files...")

  UMI.Mat <- Matrix::readMM(file = X)

  UMI.Mat <- methods::as(UMI.Mat,"dgCMatrix")

  Gene.IDs <- read.delim(Gene,header = F,stringsAsFactors = F)

  Barcodes <- read.delim(Barcode,header = F,stringsAsFactors = F)
  Barcodes <- Barcodes$V1

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

    Top.Features <- FeatureSelect(X = UMI.Mat,Reference = ReferenceLevel,GeneID = Gene.IDs,Min.Perc.Non.Zero.Cells = MinPercNonZero)

    if (is.null(VarFeatures)){
      VarFeatures <- length(Top.Features)
    }

    if(length(dim(Top.Features)) > 1){
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
        Top.Features.Ser.Nos[i] <- which(Gene.IDs[,1] == Top.Features[i])
      }
    }

    write.csv(Top.Features,file = paste0("TopFeatures",unlist(strsplit(X,".",fixed = T))[1],".csv"),row.names = F)
    message("Successfully prepared .csv file containing list of highly variable features.")

    UMI.Mat <- Matrix::t(UMI.Mat)

    UMI.Mat <- UMI.Mat[,Top.Features.Ser.Nos]

    Std.Mat <- Standardize(X = UMI.Mat,Transform = TransformType,SF = SF.Per.Cell)

    if (Out == T){
      FileName <- paste0(Transform,"TransformedStandardizedCounts.csv")
      data.table::fwrite(data.frame(Std.Mat),file = FileName,row.names = F,col.names = F,sep = ",")
    }
    return(Std.Mat)
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

    Top.Features <- FeatureSelect(X = UMI.Mat,GeneID = Gene.IDs,Min.Perc.Non.Zero.Cells = 0.5)

    if (is.null(VarFeatures)){
      VarFeatures <- length(Top.Features)
    }

    if (length(dim(Top.Features)) > 1){
      Top.Features <- Top.Features[!Top.Features[,1] %in% Zero.Count.Features,]
    } else {
      Top.Features <- Top.Features[!Top.Features %in% Zero.Count.Features]
    }

    if(length(dim(Top.Features)) > 1){
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
        Top.Features.Ser.Nos[i] <- which(Gene.IDs[,1] == Top.Features[i])
      }
    }

    write.csv(Top.Features,file = paste0("TopFeatures",unlist(strsplit(X,".",fixed = T))[1],".csv"),row.names = F)
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
    return(Std.Mat)
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
FeatureSelect <- function(X,GeneID,Reference = NULL,Min.Perc.Non.Zero.Cells = 1){

  if (is.null(Reference)){
    Reference <- 0.5
  } else if (Reference > 0.5){
    message("Reference should not be greater than 0.5. Resetting it to 0.5 (default)")
    Reference <- 0.5
  }

  UMI.Mat <- Matrix::t(X)

  Gene.IDs <- GeneID

  Var.Arith.Per.Feature <- colVarsSPM(UMI.Mat)
  Mean.Arith.Per.Feature <- Matrix::colMeans(UMI.Mat)

  message("Filtering features...")
  Irrelevant.Features <- which(Var.Arith.Per.Feature <= Mean.Arith.Per.Feature)

  No.of.Non.Zero.Per.Row <- diff(UMI.Mat@p)

  Perc.Non.Zero.Per.Row <- No.of.Non.Zero.Per.Row/nrow(UMI.Mat) * 100

  Irrelevant.Features <- unique(c(Irrelevant.Features,which(Perc.Non.Zero.Per.Row <= Min.Perc.Non.Zero.Cells)))

  if (length(Irrelevant.Features) != 0){
    UMI.Mat <- UMI.Mat[,-Irrelevant.Features]
    if (length(dim(Gene.IDs)) > 1){
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
    if (length(dim(Gene.IDs)) > 1){
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

  if (length(dim(Gene.IDs)) > 1){
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
#' @param X A character variable. Specifies the name of the .csv file containing the standardized values obtained from the \link[Piccolo]{StandardizeMat} function
#' @param TopFeatures A character variable. Specifies the name of the .csv file containing the list of genes that were shortlisted as highly variable by the \link[Piccolo]{StandardizeMat} function
#' @param Barcode A character variable. Specifies the name of the barcodes file (.tsv format)
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T), or not (if set to F). Default is F
#' @return A numeric matrix containing the standardized values.
#' @examples
#' \dontrun{
#' MaxMinNormMat(X = "10X_PBMC3k_logTransformedStandardizedCounts.csv",
#' TopFeatures = "TopFeatures10X_PBMC3k_matrix.csv",
#' Barcode = "10X_PBMC3k_barcodes.tsv",Out = T)
#' }
MaxMinNormMat <- function(X,TopFeatures,Barcode,Out = F){
  Std.Mat <- as.matrix(data.table::fread(X))

  Features <- read.csv(TopFeatures)
  Features <- Features$x

  Barcodes <- read.delim(Barcode,header = F,stringsAsFactors = F)
  Barcodes <- Barcodes$V1

  Barcodes <- gsub(".","-",Barcodes,fixed = T)

  rownames(Std.Mat) <- Features
  colnames(Std.Mat) <- Barcodes

  Norm.Mat <- t(apply(Std.Mat,1,function(x) (x- min(x))/(max(x) - min(x))))
  if (Out == T){
    Norm.df <- data.frame(Norm.Mat)
    FileName <- paste0(substr(X,1,nchar(X)-4),"_MaxMinNormMat.csv")
    data.table::fwrite(Norm.df,file = FileName,row.names = F,col.names = F,sep = ",")
  }
  return(Norm.Mat)
}

#' @title  Identify differentially expressed between 2 groups of cells
#' @description  This function performs differential expression analysis (using the Mann-Whitney test) between 2 groups of cells provided by the user
#' @export
#' @param X A character variable. Specifies the name of the .csv file containing the standardized values obtained from the \link[Piccolo]{StandardizeMat} function
#' @param TopFeatures A character variable. Specifies the name of the .csv file containing the list of genes that were shortlisted as highly variable by the \link[Piccolo]{StandardizeMat} function
#' @param Barcode A character variable. Specifies the name of the barcodes file (.tsv format)
#' @param Group1 A character vector. Specifies the barcodes of cells in group 1
#' @param Group2 A character vector. Specifies the barcodes of cells in group 2
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the differential expression result (if set to T), or not (if set to F). Default is F
#' @return A data frame containing the gene IDs, the log2 fold change (FC) of normalized values between group 1 and group 2 (positive FC indicates higher expression in group 1), the p-values from the Mann-Whitney test, and the adjusted p-values (p-adj) after Benjamini-Hochberg correction
#' @examples
#' \dontrun{
#' DEfeatures(X = "10X_PBMC3k_logTransformedStandardizedCounts.csv",
#' TopFeatures = "TopFeatures10X_PBMC3k_matrix.csv",
#' Barcode = "10X_PBMC3k_barcodes.tsv",
#' Group1 = c("Barcode1","Barcode23","Barcode47",..),
#' Group2 = c("Barcode3,"Barcode7, "Barcode11",..)
#' Out = T)
#' }

DEfeatures <- function(X,TopFeatures,Barcode,Group1,Group2,Out = F){
  Std.Mat <- as.matrix(data.table::fread(X))

  Features <- read.csv(TopFeatures)
  Features <- Features$x

  Barcodes <- read.delim(Barcode,header = F,stringsAsFactors = F)
  Barcodes <- Barcodes$V1

  Barcodes <- gsub(".","-",Barcodes,fixed = T)

  Group1 <- which(Barcodes %in% Group1)
  Group2 <- which(Barcodes %in% Group2)

  rownames(Std.Mat) <- Features
  colnames(Std.Mat) <- Barcodes

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
  Genes <- Features[order(log2FC.vec,decreasing = T)]
  Base.Mean.vec <- Base.Mean.vec[order(log2FC.vec,decreasing = T)]
  log2FC.vec <- log2FC.vec[order(log2FC.vec,decreasing = T)]

  DE.Res.df <- data.frame(Genes,Base.Mean.vec,log2FC.vec,p.val.vec,p.adj.vec)
  colnames(DE.Res.df) <- c("Gene","BaseMean","log2FC","p.val","p.adj")

  if (Out == T){
    FileName <- paste0(substr(X,1,nchar(X)-4),"DEGenes",".csv")
    data.table::fwrite(DE.Res.df,file = FileName,row.names = F,col.names = F,sep = ",")
  }

  return(DE.Res.df)
}

#' @title  Compute Principal Components
#' @description  This function will calculate the principal components for the matrix with the standardized values obtained from the \link[Piccolo]{StandardizeMat} function.
#' @export
#' @param X A character variable. Specifies the name of the .csv file containing the standardized values obtained from the \link[Piccolo]{StandardizeMat} function
#' @param TopFeatures A character variable. Specifies the name of the .csv file containing the list of genes that were shortlisted as highly variable by the \link[Piccolo]{StandardizeMat} function
#' @param Barcode A character variable. Specifies the name of the barcodes file (.tsv format)
#' @param NoOfPC A numeric (integer) variable. No of principal components to be retained in the output. Default is 50.
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T), or not (if set to F). Default is T
#' @return A numeric matrix containing the standardized values.
#' @examples
#' \dontrun{
#' ComputePC(X = "10X_PBMC3k_logTransformedStandardizedCounts.csv",
#' TopFeatures = "TopFeatures10X_PBMC3k_matrix.csv",
#' Barcode = "10X_PBMC3k_barcodes.tsv",Out = T)
#' }
ComputePC <- function(X,TopFeatures,Barcode,NoOfPC =  50,Out = T){
  Std.Mat <- as.matrix(data.table::fread(X))

  Features <- read.csv(TopFeatures)
  Features <- Features$x

  Barcodes <- read.delim(Barcode,header = F,stringsAsFactors = F)
  Barcodes <- Barcodes$V1

  Barcodes <- gsub(".","-",Barcodes,fixed = T)

  rownames(Std.Mat) <- Features
  colnames(Std.Mat) <- Barcodes

  res.pca <- stats::prcomp(t(Std.Mat))

  if (Out == T){
    PC.df <- data.frame(res.pca$x[,1:NoOfPC])
    FileName <- paste0(substr(X,1,nchar(X)-4),"_Top",NoOfPC,"PrinComp",".csv")
    data.table::fwrite(PC.df,file = FileName,row.names = T,col.names = F,sep = ",")
  }
  return(res.pca$x[,1:NoOfPC])
}

#' @title  UMAP coordinates
#' @description  This function will run UMAP for the matrix with the principal components
#' @export
#' @param X A character variable. Specifies the name of the .csv file containing the principal components obtained using the \link[Piccolo]{ComputePC} function
#' @param Out A logical variable. Specifies whether to return an output file (.csv) with the normalized values (if set to T), or not (if set to F). Default is T
#' @return A data frame containing the coordinates of the cells for the first 3 UMAP dimensions.
#' @examples
#' \dontrun{
#' UMAPCoords(X = "10X_PBMC3k_logTransformedStandardizedCounts_Top50PrinComp.csv",
#' Out = T)
#' }
UMAPCoords <- function(X,Out = T){

  PC.df <- data.table::fread(X)

  PC.Mat <- as.matrix(PC.df[,-1])
  rownames(PC.Mat.TopPC) <- PC.df$V1

  rm(PC.df)

  x <- umap::umap(PC.Mat)

  UMAP.df <- data.frame(CellID = rownames(x$layout),UMAP1  = x$layout[,1], UMAP2 = x$layout[,2], UMAP3 = x$layout[,3])

  if (Out == T){
    FileName <- paste0(substr(X,1,nchar(X)-4),"_UMAPCoords.csv")
    data.table::fwrite(UMAP.df,file = FileName,row.names = F,sep = ",")
  }
  return(UMAP.df)
}

#' @title  Label cells on UMAP
#' @description  This function can be used to label (color) cells on the UMAP plot
#' @export
#' @param X A character variable. Specifies the name of the .csv file containing the UMAP coordinates obtained from the \link[Piccolo]{UMAPCoords} function
#' @param Labels A character vector. Should contain the character labels for cells in the same order as the cells in the counts matrix
#' @param Levels An optional character vector. Should specify the unique labels in the order in which the labels of the cells should be presented
#' @param Alpha A numeric variable (strictly greater than 0). Specified the transparency of the dots on the UMAP plot. Default Alpha = 0.7
#' @param Size A numeric variable (strictly greater than 0). Specified the size of the dots on the UMAP plot. Default Size = 0.9
#' @return A ggplot2 object
#' @examples
#' \dontrun{
#' LabelUMAP(X = "10X_PBMC3k_logTransformedStandardizedCounts_Top50PrinComp_UMAPCoords.csv",
#' Labels = c("b-cells","b-cells",..,"cd14 monocytes",..,"NK cells",..),
#' Levels = c("b-cells","cd14 monocytes","dendritic","NK cells","naive cytotoxic"))
#' }

LabelUMAP <- function(X,Labels,Levels = NULL,Alpha = 0.7,Size = 0.9){

  UMAP.Coord.df <- data.table::fread(X)

  UMAP.Coord.df <- UMAP.Coord.df[,-4]

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



