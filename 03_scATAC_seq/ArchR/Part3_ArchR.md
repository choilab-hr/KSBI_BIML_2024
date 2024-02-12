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

