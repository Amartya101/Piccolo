# Piccolo
An R package for processing and analyzing single-cell counts data (Preprint can be found at [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.03.02.530891v1) and the peer-reviewed manuscript (published in BMC Bioinformatics) can be found [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05872-w))

## Getting Started

We recommend an R version 4.2.0 or more recent in order to use the package.

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
The *CreatePiccoloList* function requires the following inputs - *MTX*, *Genes*, *Barcodes*. Assuming all the files are located in the current working directory, *MTX* should be the full name of the .mtx or .mtx.gz file that contains the gene counts. *Genes* should be the full name of the .tsv file that contains the information of the genes/features. *Barcodes* should be the full name of the .tsv file that contains the information of the barcodes. 

Here is an example of a valid function call (the 10X PBMC3k data set used as an example here can be found in the [PBMC3k_Data](https://github.com/Amartya101/PBMC3k_Data) repository):

```
pbmc3k <- CreatePiccoloList(MTX = "10X_PBMC3k_matrix.mtx.gz", Genes = "10X_PBMC3k_features.tsv", Barcodes = "10X_PBMC3k_barcodes.tsv")
```
This will create a list that contains the counts matrix, the features data frame, and the barcodes.

### FilterCells 
The *FilterCells* function is used to perform basic cell filtering based on criterias such as the minimum number of genes with non-zero counts within each cell (specified through the argument *MinFeaturesPerCell*), the maximum percentage of total counts contributed by mitochondrial genes within each cell (specified by the argument *MT.Perc*), the the maximum percentage of total counts contributed by ribosomal genes within each cell (specified by the argument *RP.Perc*), and the maximum or minimum total counts that any cell can have based on how many median absolute deviation away from the median the total count of any given cell is.

Examples of valid function calls are given below:
```
pbmc3k <- FilterCells(PiccoloList = pbmc3k)

pbmc3k <- FilterCells(PiccoloList = pbmc3k,
 MinFeaturesPerCell = 100, MT.Perc = 50,
 RP.Perc = 70,TotalCountsMADHigh = 3.5,
 TotalCountsMADLow = 3.5) #changing the filtering criteria
```

### SelectFeatures
The *SelectFeatures* function performs feature selection by identifying highly variable genes and stable genes. It takes in the PiccoloList obtained after using the *FilterCells* function as input and first performs basic gene filtering based on the criteria that at least some percentage of cells have non-zero counts for each gene (this is specified through the argument *MinPercNonZeroCells*). It then identifies highly variable genes and stable genes. You can specify the number of highly variable genes that you wish to shortlist by using the argument *NoOfHVG*. Output files containing the list of HVGs and the list of stable genes can be generated by specifying the argument *Out* to be T. 

Examples of valid function calls are given below:
```
pbmc3k <- SelectFeatures(PiccoloList = pbmc3k)

pbmc3k <- SelectFeatures(PiccoloList = pbmc3k,
 NoOfHVG = 3000, Out = T)
```

### Normalize

The *Normalize* function computes the residuals for the counts of the HVGs. You can specify the variance stabilization transformation that you wish to apply (default is "log") by specifying the *Transform* argument. Apart from log, options include "sqrt" (for the the square-root transform) and "bc" (for the Box-Cox power law transform). You can even generate an output file containing the residuals by specifying *Out* as T.


Examples of valid function calls are given below:
```
pbmc3k <- Normalize(PiccoloList = pbmc3k)

pbmc3k <- Normalize(PiccoloList = pbmc3k, Transform = "bc")
```

### Compute PCs and get UMAP coordinates

After the normalized counts have been obtained using the *Normalize* function, we can compute the principal components and generate the UMAP coordinates. The *ComputePC* function takes in the PiccoloList prepared by the *Normalize* function and by default gives the top 50 principal components. After the principal components are obtained and stored in the PiccoloList, we can apply the *UMAPcoords* function that takes in the PiccoloList and calculates the UMAP coordinates. 

Examples of valid function calls:
```
pbmc3k <- ComputePC(PiccoloList = pbmc3k)

pbmc3k <- ComputePC(PiccoloList = pbmc3k,NoOfPC = 20,Out = T) # for Top 20 PCs, and will generate an output .csv file containing the PCs

pbmc3k <- UMAPcoords(PiccoloList = pbmc3k, Out = T)
```
### Clustering
After the principal components have been computed, we can use the *KNearestNeighbors* function to identify the *k* nearest-neighbors of every cell based on the PC coordinates (default *k* is 10). After obtaining the nearest-neighbors we can use graph-based partitioning algorithm such as Leiden to identify clusters of cells. This can be accomplished by using the *LeidenClustering* function. 

Examples of valid function calls:

```
pbmc3k <- KNearestNeighbors(PiccoloList = pbmc3k)

pbmc3k <- LeidenClustering(PiccoloList = pbmc3k)

pbmc3k <- LeidenClustering(PiccoloList = pbmc3k,
 Resolution = 1.5)
```

### Label cells on UMAP
The *LabelUMAP* function can be used to make the UMAP plots and color the cells based on labels provided by the user. *Labels* should contain the character labels for all the cells in the same order as the cells in the counts matrix.

Examples of valid function calls are provided below:

```
CellLabels <- c("b-cells","b-cells",..,"cd14 monocytes",..,"NK cells",..)
p <- LabelUMAP(PiccoloList = pbmc3k,
 Labels = CellLabels,
 Levels = c("b-cells","cd14 monocytes","dendritic","NK cells","naive cytotoxic"),
 Title = "PBMC3k")

p

p <- LabelUMAP(PiccoloList = pbmc3k, Labels = PiccoloList$ClusterLabels,Title = "PBMC3k")

p
```

### Perform differential expression analysis between 2 groups of cells
The *PerformDiffExp* function employs either the Student's t-test (by specifying *Method = "t.test"*) or the Wilcoxon rank-sum test (by specifying *Method = "wilcoxon"*) to determine whether any of the given features are differentially expressed between 2 groups of cells specified by the user. *Group1* and *Group2* should contain the serial numbers of cells that belong to the respective groups.

Example of a valid function call is provided below:

```
Group1.vec <- 1:200
Group2.vec <- 301:500
pbmc3k <- PerformDiffExp(PiccoloList = pbmc3k,
 Group1 = Group1.vec,
 Group2 = Group2.vec,
 Out = T)
```

## Authors

* **Amartya Singh** - [Amartya101](https://github.com/Amartya101/)

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
* [Khiabanian Lab](https://khiabanian-lab.org)

