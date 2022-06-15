# Piccolo
An R package for processing and analyzing single-cell transcriptomic data (manuscript under preparation)

## Getting Started

You will need an R version 3.4.0 or more recent in order to use the package.

### Installation

Currently, the development version of **Piccolo** can be installed using *install_github* from the **devtools** package.

Install **devtools** using:
```
install.packages("devtools",dependencies = T)
```
Once **devtools** is installed, run the following to install **Piccolo**:
```
devtools::install_github("Amartya101/Piccolo")
```
## Instructions for Use

The central object in the Piccolo workflow is the PiccoloList. Given the single-cell gene expression data in the .mtx or .mtx.gz format (and the corresponding features and barcodes .tsv files), we first need to create a PiccoloList that contains the matrix and the information about the features and barcodes. This can be done with the *CreatePiccoloList* function.

### CreatePiccoloList
The *CreatePiccoloList* function requires the following inputs - *X*, *Gene*, *Barcode*, *MinFeaturesPerCell*, *MT.Perc*, and *RP.Perc*. Assuming all the files are located in the current working directory, *X* should be the full name of the .mtx or .mtx.gz file that contains the gene counts. *Gene* should be the full name of the .tsv file that contains the information of the genes/features. *Barcode* should be the full name of the .tsv file that contains the information of the barcodes. These first 3 inputs are necessary. The following 3 arguments are optional: *MinFeaturesPerCell* (Minimum number of unique features that a cell must contain), *MT.Perc* (Maximum percentage of total counts due to mitochondrial genes that any given cell can exhibit), *RP.Perc* (Maximum percentage of total counts due to ribosomal genes that any given cell can exhibit). If not explicitly specified, they will employ the default thresholds of *MinFeaturesPerCell = 200*, *MT.Perc = 10*, and *RP.Perc = 70*. 

Here is an example of a valid function call

```
pbmc3k <- CreatePiccoloList(X = "10X_PBMC3k_matrix.mtx", Gene = "10X_PBMC3k_features.tsv", Barcode = "10X_PBMC3k_barcodes.tsv", RP.Perc = 80)
```
This will create a list that contains the counts matrix, the features data frame, and the barcodes.

### Normalize

After the PiccoloList object has been created, we need to perform feature selection and normalization on the counts. This is done by the *Normalize* function.

The *Normalize* function has 7 arguments - *PiccoloList*, *VarFeatures*, *Transform*, *Batch*, *ReferenceLevel*, *MinPercNonZero*, and *Out*. *PiccoloList* should be the list created using the *CreatePiccoloList* function. *VarFeatures* specifies the number of highly variable features you wish to shortlist. If left NULL, the number of variable features shortlisted will depend on how many features exceed the median overdispersion coefficient for every expression level. *Transform* specifies the non-linear transformation to be applied to the counts (default is log). Currently, apart from the log transformation we offer the option of applying the Yeo-Johnson power transform (*Transform = "yj"*). *Batch* requires the specification of a character vector that provides the batch labels for the cells in the same order as the cells in the count data/barcodes file. If left NULL, no batch correction will be performed. *ReferenceLevel* offers the option to alter the reference threshold for the shortlisting of variable features (default is 0.5, corresponding to median as the reference). *MinPercNonZero* specifies the minimum percentage of cells that must exhibit non-zero counts for any given feature. *Out* is a logical variable to specify whether an output .csv file containing the normalized counts shold be generated (default is F) 

Examples of valid function calls are provided below:
```
pbmc3k <- Normalize(PiccoloList = pbmc3k)
pbmc3k <- Normalize(PiccoloList = pbmc3k, ReferenceLevel = 0.3,MinPercNonZero = 0.5,Out = T)
```

### Compute PCs and get UMAP coordinates

After the normalized counts have been obtained using the *Normalize* function, we can compute the principal components and generate the UMAP coordinates. The *ComputePC* function takes in the PiccoloList prepared by the *Normalize* function and by default gives the top 50 principal components. After the principal components are obtained and stored in the list, we can apply the *UMAPcoords* function that takes in the PiccoloList and calculates the UMAP coordinates. 

Examples of valid function calls:
```
pbmc3k <- ComputePC(PiccoloList = pbmc3k)
pbmc3k <- ComputePC(PiccoloList = pbmc3k,NoOfPC = 20,Out = T) # for Top 20 PCs, and will generate an output .csv file containing the PCs

pbmc3k <- UMAPcoords(PiccoloList = pbmc3k)
```

### Label cells on UMAP
The *LabelUMAP* function can be used to make the UMAP plots and color the cells based on labels provided by the user. *Labels* should contain the character labels for all the cells in the same order as the cells in the counts matrix.

Example of a valid function call is provided below:
```
p <- LabelUMAP(PiccoloList = pbmc3k,
Labels = c("b-cells","b-cells",..,"cd14 monocytes",..,"NK cells",..),
Levels = c("b-cells","cd14 monocytes","dendritic","NK cells","naive cytotoxic"))

p
```
### Perform differential expression analysis between 2 groups of cells
The *DEfeatures* function can be used to normalize the standardized values such that they lie in the 0 to 1 range. Subsequently, the Mann-Whitney U test is employed to determine whether the expression levels of any given features are differentially expressed between 2 groups of cells specified by the user. *Group1* and *Group2* should contain the serial numbers of cells that belong to the respective groups. The output is a data frame that contains the log2FC of the normalized standardized values of the *Group1* cells relative to *Group2* cells, the base mean expression level (in normalized values) of each feature, and the p-values obtained from the Mann-Whitney U test.

Example of a valid function call is provided below:
```
DE.Res.df <- DEfeatures(PiccoloList = pbmc3k,
Group1 = c(2,3,5,7,11,..),
Group2 = c(4,6,8,9,10,..))
```

## Authors

* **Amartya Singh** - [Amartya101](https://github.com/Amartya101/)

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
* Hossein Khiabanian - [Khiabanian Lab](https://khiabanian-lab.org)

