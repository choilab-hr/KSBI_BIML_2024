# Installation of ArchR

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

# Load ArchR and Download tutorial Data
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




















