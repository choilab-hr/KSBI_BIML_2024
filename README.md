# KSBI-BIML 2024
Welcome to the **Single-cell Multiomics** session.  
This is a GitHub repository for the session.  
[Syllabus](https://journal-home.s3.ap-northeast-2.amazonaws.com/site/biml2024/intro/off-12.pdf, "syllabus link")

- 장소: 서울대학교 자연과학대학 28동 102호
- 시간: 2024년 2월 27일 (화), 오전 9시 30분 ~ 오후 4시 50분

This session will cover:  
1. Spatial Transcriptomics analysis
2. Multi-omics analysis with Spatial Transcriptomics
3. scATAC-seq data analysis and GRN construction

# Prerequisites
실습일 (2024년 2월 27일 화요일)전 아래의 과정을 따라서 R package와 필요한 자료를 다운로드 받아오시길 바랍니다.  
  
## 0. devtools, BiocManager
First, install devtools (for installing GitHub packages) if it isn’t already installed:

```R
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```
Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```
