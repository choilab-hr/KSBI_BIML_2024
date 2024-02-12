# 0. Installation of ArchR

ArchR is designed to be run on Unix-based operating systems such as macOS and linux. **ArchR is NOT supported on Windows or other operating systems.**  

ArchR installation currently requires devtools and BiocManager for installation of GitHub and Bioconductor packages.  

First, install devtools (for installing GitHub packages) if it isn’t already installed:

```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```
Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```
Then, install ArchR:
```R
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
```
Install all of the ArchR dependencies that arent installed by default:
```R
library(ArchR)
ArchR::installExtraPackages()
```
Set a working directory variable for the session:
```R
biml_dir <- 'your/directory'
setwd(biml_dir)
```

# 1. Load ArchR and Download tutorial Data
First, load ArchR and set a random seed.
```R
library(ArchR)
set.seed(2024)
```
Set the default number of threads for parallelized operations in ArchR functions.
```R
((ncore <- parallel::detectCores()))
addArchRThreads(threads = ncore-2)
```
Get tutorial data for the session.
```R
inputFiles <- getTutorialData("Hematopoiesis")
inputFiles
```
<img width="612" alt="image" src="https://github.com/choilab-hr/KSBI_BIML_2024/assets/159281429/fcf10bf2-6d7f-4b51-afb4-2efe71d5d5c4">

Lastly, add a reference genome annotation for ArchR.
```R
addArchRGenome("hg19")
```
# 2. Creating Arrow files and Quality Control
Now we will create our Arrow files which will take 10-15 minutes. For each sample, this step will:  

1. Read accessible fragments from the provided input files.
2. Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
3. Filter cells based on quality control parameters.
4. Create a genome-wide TileMatrix using 500-bp bins.
5. Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().

```R
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles
```

# 3. Inferring scATAC-seq Doublets with ArchR
```R
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
```

# 4. Creating an ArchRProject and Doublet removal
## 4.1. Creating an ArchRProject
```R
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "2024_BIML_archR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
getAvailableMatrices(proj)
```
## 4.2. Filter doublets
```R
proj <- filterDoublets(ArchRProj = proj)
paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")
```

# 5. Save and Load an ArchRProject
Save an ArchRProject:  
```R
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)
```
Load an ArchRProject
```R
proj <- loadArchRProject(paste0(biml_dir, '/Save-Proj'))
```

# 6. Dimensionality reduction
```R
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
# saveArchRProject(ArchRProj = proj, outputDirectory = "01_Dimensionality_reduction", load = FALSE) # Optional
```
If a dimensionality reduction was not successfully done, don't worry.  
1. Download the following link https://www.dropbox.com/scl/fo/uac5vx8qboppxvddacxp2/h?rlkey=5w4hmidncjkp237mohjy402y7&dl=0
2. Execute the following script:
```R
proj <- loadArchRProject(paste0(biml_dir, '/01_Dimensionality_reduction'))
```














