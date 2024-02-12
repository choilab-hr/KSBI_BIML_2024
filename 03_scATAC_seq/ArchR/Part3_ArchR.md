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
***!!! OPTIONAL !!!***  
If a dimensionality reduction was not successfully done, don't worry.  
1. Download a directory to your working directory with the following link https://www.dropbox.com/scl/fo/uac5vx8qboppxvddacxp2/h?rlkey=5w4hmidncjkp237mohjy402y7&dl=0
2. Execute the following script:
```R
proj <- loadArchRProject(paste0(biml_dir, '/01_Dimensionality_reduction'))
```

# 7. Clustering with ArchR
```R
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
# saveArchRProject(ArchRProj = proj, outputDirectory = "02_Clustering", load = FALSE)
```
***!!! OPTIONAL !!!***  
If a dimensionality reduction was not successfully done, don't worry.  
1. Download a directory to your working directory with the following link https://www.dropbox.com/scl/fo/6urh6s8p9pkm86vla45k4/h?rlkey=2rlrgehcwatb28wyi2s01waed&dl=0
2. Execute the following script:
```R
proj <- loadArchRProject(paste0(biml_dir, '/02_Clustering'))
```
# 8. Visualization in a Two-dimensional space
```R
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

# UMAP colored by the Sample
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

# UMAP colored by the Clusters
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")

# Save a plot
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
# saveArchRProject(ArchRProj = proj, outputDirectory = "03_UMAP", load = FALSE)
```
<img width="286" alt="image" src="https://github.com/choilab-hr/KSBI_BIML_2024/assets/159281429/1bff0418-b076-445f-940b-fac321c8c161">
  
***!!! OPTIONAL !!!***  
If a dimensionality reduction was not successfully done, don't worry.  
1. Download a directory to your working directory with the following link https://www.dropbox.com/scl/fo/oyfjfit64sj5z9gu490fi/h?rlkey=gakkfdnjro1adyb7ukjsn0orh&dl=0
2. Execute the following script:
```R
proj <- loadArchRProject(paste0(biml_dir, '/03_UMAP'))
```

# 9. Identifying Marker genes for each cluster
```R
# Do not run below
# markersGS <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = "GeneScoreMatrix", 
#   groupBy = "Clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# saveRDS(markersGS, file = paste0(biml_dir, '/2024_BIML_markerGS.rds')) # Optional
```
Download the pre-calculated markerGS:  
https://www.dropbox.com/scl/fi/6ao0s2jkgtuf9r6e70ir8/2024_BIML_markerGS.rds?rlkey=i1nft4ugjkwpldnmx3zxixzt9&dl=0  
```R
markersGS <- readRDS(file = paste0(biml_dir, '/2024_BIML_markerGS.rds'))
```
Get the marker genes for the Cluster 6
```R
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
markerList$C6 %>% as.data.frame() %>% dplyr::arrange(desc(Log2FC)) %>% dplyr::pull(name) %>% head(10)
```
## 9.1. Visualization with Heatmap
```R
markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8", 
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

heatmapGS <- markerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = markerGenes,transpose = TRUE)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)
```
<img width="506" alt="image" src="https://github.com/choilab-hr/KSBI_BIML_2024/assets/159281429/63b8aba7-5602-47d9-b705-7fe09ecf99a9">

## 9.2. Visulization on an embedding
```R
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)
p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP", imputeWeights = getImputeWeights(proj))
p$CD14
```
### 9.2.1. Plot all genes
To plot all genes we can use cowplot to arrange the various marker genes into a single plot.  
```R
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)
```
## 9.3. Imputation with MAGIC
```R
proj <- addImputeWeights(proj)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = 'CD14', 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p
```

# 10. Visualizing Genome Browser Tracks
```R
p <- plotBrowserTrack(
  ArchRProj = proj, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$CD14)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)

```

# 11. Integration with scRNA-seq data - DO NOT RUN
```R
if(!file.exists("scRNA-Hematopoiesis-Granja-2019.rds")){
  download.file(
    url = "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds",
    destfile = "scRNA-Hematopoiesis-Granja-2019.rds"
  )
}

seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")
seRNA

colnames(colData(seRNA))

table(colData(seRNA)$BioClassification)

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "BioClassification",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

```











