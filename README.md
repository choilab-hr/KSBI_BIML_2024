# KSBI-BIML 2024, Single-cell Multiomics
Welcome to the **Single-cell Multiomics** session.  
This is a GitHub repository for the session.  

- **장소**: 서울대학교 자연과학대학 28동 102호
- **시간**: 2024년 2월 27일 (화), 오전 9시 30분 ~ 오후 4시 50분
- **연자**: 최정민 교수
- **조교**: 김지현, 유광민, 이다준, 이문영, 이호진, 최승지, 천하림, 홍주현
- **강의 계획서**: [Syllabus](https://journal-home.s3.ap-northeast-2.amazonaws.com/site/biml2024/intro/off-12.pdf, "syllabus link")
- **강의 자료**: [강의자료](https://www.dropbox.com/scl/fi/0q5rerxe1w1ig356dadbn/20240213_BIML2024_Combined_-_v2.pdf?rlkey=dqtrvv830jv37ibl627lgcl0z&dl=0, "강의자료")

**This session will cover:**  
1. Spatial Transcriptomics analysis
2. Multi-omics analysis with Spatial Transcriptomics
3. scATAC-seq data analysis and GRN construction

# Prerequisites
원활한 진행을 위해, **실습일 (2024년 2월 27일 화요일)** 전 R 및 R studio를 설치하시고,  
아래의 과정에 따라 실습에 사용할 R package 및 필요한 자료를 다운로드해오시길 바랍니다.
  
## 0️⃣ devtools, BiocManager
First, install devtools (for installing GitHub packages) if it isn’t already installed:
```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```
Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```
## 1️⃣ Part1. Spatial Transcriptomics analysis
### R package
#### ❗️Seurat
```R
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
```
#### ❗️BayesSpace
```R
BiocManager::install("BayesSpace")
```
#### ❗️hdf5r
```R
install.packages('hdf5r')
```
#### ❗️SingleCellExperiment
```R
BiocManager::install("SingleCellExperiment")
```
### Materials
추후 업데이트 예정
## 2️⃣ Part2. Multi-omics analysis with Spatial Transcriptomics
### R package
#### ❗️scDblFinder
```R
BiocManager::install("scDblFinder")
```
#### ❗️spacexr
```R
options(timeout = 600000000) ### set this to avoid timeout error
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```
#### ❗️CellChat
```R
options(timeout = 600000000) ### set this to avoid timeout error
  remotes::install_github("sqjin/CellChat")
```
### Materials
추후 업데이트 예정
## 3️⃣ Part3. scATAC-seq data analysis and GRN construction
### R package
#### ❗️ArchR
> ```R
> devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
> ```
> Install all of the ArchR dependencies that aren't installed by default:
> ```R
> library(ArchR)
> ArchR::installExtraPackages()
> ```
#### ❗️figR
```R
devtools::install_github("buenrostrolab/FigR")
```
#### ❗️patchwork
```R
devtools::install_github("thomasp85/patchwork")
```
#### ❗️BSgenome.Hsapiens.UCSC.hg19
```R
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
```
#### ❗️ggrastr
```R
devtools::install_github('VPetukhov/ggrastr')
```
### Materials
추후 업데이트 예정
