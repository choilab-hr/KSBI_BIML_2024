# KSBI-BIML 2024, Single-cell Multiomics
Welcome to the **Single-cell Multiomics** session.  
This is a GitHub repository for the session.  

- **장소**: 서울대학교 자연과학대학 28동 102호
- **시간**: 2024년 2월 27일 (화), 오전 9시 30분 ~ 오후 4시 50분
- **연자**: 최정민 교수
- **조교**: 김지현, 유광민, 이다준, 이문영, 이호진, 최승지, 천하림, 홍주현
- **강의 계획서**: [Syllabus](https://journal-home.s3.ap-northeast-2.amazonaws.com/site/biml2024/intro/off-12.pdf, "syllabus link")
- **강의 자료**: [강의자료](https://www.dropbox.com/scl/fi/mt5qcbcyhj0ot7logmrxz/20240219_BIML2024_Combined_v3_.pdf?rlkey=8b4bte0s8suyj83q4yh7715ob&dl=0, "강의자료")
- **강의 시간표**: ![timetable_hr](https://github.com/choilab-hr/KSBI_BIML_2024/assets/159281429/9790db00-1090-443d-8598-9f1590993911)


**This session will cover:**  
1. Spatial Transcriptomics analysis
2. Multi-omics analysis with Spatial Transcriptomics
3. scATAC-seq data analysis and GRN construction

# Prerequisites
~~원활한 진행을 위해, **실습일 (2024년 2월 27일 화요일)** 전 R 및 R studio를 설치하시고,~~  
~~아래의 과정에 따라 실습에 사용할 R package 및 필요한 자료를 다운로드해오시길 바랍니다.~~  
❗️ KOBIC 서버를 활용한 실습을 진행할 예정이기때문에 아래의 R package와 Materials는 서버내에 미리 구축해둘 예정입니다. ❗️  
❗️ 설치 및 다운로드는 선택사항입니다. ❗️
  
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
https://www.dropbox.com/scl/fo/r9epz6dl6dga1e7gs6anm/h?rlkey=5n0pzcrvw0elt2fp4m80m17ph&dl=0
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
#### ❗️Cottrazm
Please refer to '2. Installation and requirement' in the vignette.  
https://github.com/Yelab2020/Cottrazm/blob/main/doc/my-vignette.pdf

### 10x Software
#### ❗️10x genomics xenium explorer
Please agree to the End User License Agreement and download 10x genomics xenium explorer.  
https://www.10xgenomics.com/support/software/xenium-explorer/latest

### Materials
https://www.dropbox.com/scl/fo/0xntiwvduozhbxqvgp65z/h?rlkey=3w5xzwrd8d313a7pk8hnklym2&dl=0

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
1. ArchR: https://www.dropbox.com/scl/fi/6ao0s2jkgtuf9r6e70ir8/2024_BIML_markerGS.rds?rlkey=i1nft4ugjkwpldnmx3zxixzt9&dl=0
2. FigR: https://www.dropbox.com/scl/fi/q63f4wr4jlwtva7z72i9g/FigR_stim.zip?rlkey=gibefa8gdtj4zto78rnvmuym1&dl=0
