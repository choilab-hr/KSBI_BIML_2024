# 0. Prerequisites -------------------------------------------------------
# First, install devtools (for installing GitHub packages) if it isn’t already installed:
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

#Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Then, install ArchR:
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

#Lastly, install all of the ArchR dependencies that arent installed by default:
library(ArchR)
ArchR::installExtraPackages()

# 1. Load library and download tutorial data ------------------------------
biml_dir <- '~/BIML/Part3/ArchR/'
setwd(biml_dir)

## 1.1. Load library ####
library(ArchR)
set.seed(2024)

## 1.2. Set default number of threads ####
((ncore <- parallel::detectCores())) # 10
addArchRThreads(threads = ncore-2) 

## 1.3. Get tutorial data ####
inputFiles <- getTutorialData("Hematopoiesis")
inputFiles

addArchRGenome("hg19")

# 2. Creating Arrow files -------------------------------------------------
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles


# 3. Inferring scATAC-seq Doublets with ArchR -----------------------------------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

# 4. Create an ArchRProject -----------------------------------------------
## 4.1. Create an ArchRProject ####
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "2024_BIML_archR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

getAvailableMatrices(proj)

## 4.2. Filter doublets ####
proj <- filterDoublets(ArchRProj = proj)
paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")


# 5. Save and Load an ArchRProject ----------------------------------------
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-Proj", load = FALSE)
proj <- loadArchRProject(paste0(biml_dir, '/Save-Proj'))

# 6. Dimensionality reduction ---------------------------------------------
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
# saveArchRProject(ArchRProj = proj, outputDirectory = "01_Dimensionality_reduction", load = FALSE)
# proj <- loadArchRProject(paste0(biml_dir, '/01_Dimensionality_reduction'))

# 7. Clustering with ArchR ------------------------------------------------
proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)
# saveArchRProject(ArchRProj = proj, outputDirectory = "02_Clustering", load = FALSE)
# proj <- loadArchRProject(paste0(biml_dir, '/02_Clustering'))


# 8. Visualization in a 2D ------------------------------------------------
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

## 8.1. UMAP colored by the Sample ####
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

## 8.2. UMAP colored by the Clusters ####
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")

plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## 8.3. TSNE ####
proj <- addTSNE(ArchRProj = proj, reducedDims = "IterativeLSI", name = "TSNE")

# UMAP colored by the Sample
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNE")

# UMAP colored by the Clusters
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")

# Save a plot
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# saveArchRProject(ArchRProj = proj, outputDirectory = "03_UMAP", load = FALSE)
# proj <- loadArchRProject(paste0(biml_dir, '/03_UMAP'))

# 9. Identifying Marker genes for each cluster ----------------------------
# markersGS <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = "GeneScoreMatrix", 
#   groupBy = "Clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# saveRDS(markersGS, file = paste0(biml_dir, '/2024_BIML_markerGS.rds')) # Optional
# markersGS <- readRDS(file = paste0(biml_dir, '/2024_BIML_markerGS.rds'))

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6
markerList$C6 %>% as.data.frame() %>% dplyr::arrange(desc(Log2FC)) %>% dplyr::pull(name) %>% head(10)

## 9.1.  Visualization with Heatmap ####
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

## 9.2. Visulization on an embedding ####
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)

p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP", imputeWeights = getImputeWeights(proj))

p$CD14

### 9.2.1. Plot all genes ####
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

## 9.3. Imputation with MAGIC ####
proj <- addImputeWeights(proj)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = 'CD14', 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p


# 10. Visualizing Genome Browser Tracks -----------------------------------
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

# 11. Integration with scRNA-seq data - Do not run -------------------------------------
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





