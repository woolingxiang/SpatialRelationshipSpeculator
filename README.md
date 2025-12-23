### Inroduction
The SpatialRelationshipSpeculator (SRS) algorithm quantitatively assesses the spatial interplay between two objects of interest (e.g., a histological/feature domain or cell population) within a spatial transcriptomic dataset. Briefly, for a primary object of interest (Object a), the algorithm evaluates its spatial relationship with a target object (Object b) by computationally simulating the displacement of Object b in eight cardinal directions across varying step sizes. At each step, the algorithm calculates the Jaccard index to quantify the overlap of spatial spots between the two objects. The change in the Jaccard index (ΔJaccard) relative to the original state for every direction and step is computed. The resulting profile of ΔJaccard values is used to classify the spatial relationship, such as co-localizing, approaching, or distant (see Methods).The algorithm characterizes this relationship by identifying the maximum and minimum ΔJaccard values across all displacements. These extreme values define a vector for each step that encapsulates the directional tendency of spatial overlap change. A vector oriented towards quadrant 3 (Q3), where the maximum ΔJaccard decreases, indicates that Object b is initially near the core of Object a. Conversely, a vector oriented towards quadrant 4 (Q4) or even quadrant 1 (Q1), where the maximum or minimum ΔJaccard increases, suggests Object b originates at the periphery of or outside Object a. The magnitude (length) of this vector reflects the degree of spatial spot sharing, indicating the strength of the spatial association. 

<img width="1788" height="772" alt="image" src="https://github.com/user-attachments/assets/1ac6db65-0faf-455b-b489-bc36dc55180a" />

### Installation
#### option 1
```{r}
install.packages('devtools') # version ≥ 2.4.3
devtools::install_github('woolingxiang/SpatialRelationshipSpeculator')
```
#### option 2
```{r}
install.packages('renv') # version ≥ 1.0.5
renv::init('.')
renv::install('woolingxiang/SpatialRelationshipSpeculator')
renv::snapshot(type='all')
```
#### dependency
Seurat ≤ 4.4.0

---

### Quick Start
```{r}

library(SpatialRelationshipSpeculator)
library(Seurat)

dat1 = readRDS('./genomicX10/spatial_seuratObj_sample1.RDS')
dat2 = readRDS('./genomicX10/spatial_seuratObj_sample2.RDS')
dat = list(sample1=dat1,sample2=dat2)

# run for a single sample
rs1 = spatial_vector(dat1,target='HIF1A',features=c('VEGFA','PTPRC','CD68'))

# run for a list of samples
rs2 = spatial_vectorX(dat,target='HIF1A',features=c('VEGFA','PTPRC','CD68'))

# explore gene sets
cellmarkers = list(malignant = c('PTPRZ1','SOX2','EGFR'),  
                   BMDM = c('FPR3','ITGA4','TGFBI','KYNU','S100A11','IFITM2'), 
                   microglia = c('SLC1A3','CX3CR1','P2RY12','SIGLEC8','NAV3'),
                   tnkcell = c('PTPRC','CD3D','CD3E','CD8A','CD4','KLRD1','NCAM1'))
for(i in 1:length(dat)) dat[[i]] = AddModuleScore(object = dat[[i]], features = cellmarkers , name = names(cellmarkers), seed=666, nbin = 10)
features = paste0(cellmarkers,1:length(cellmarkers))
rs3 = spatial_vectorX(dat,target=features[1],features=features[-1])
rs4 = spatial_vectorX(dat,target=features[1],features=c('VEGFA','PTPRC','CD68'))


# visualization
par(mfrowc(2,4))
spatial_vecPlot(rs2) # spatial vector plot
spatial_magPlot(rs2) # spatial vector magnitude & projected score plot
sa1 = spatial_adjust(dat1,feature='HIF1A',plot = T) # spatial map 
sa2 = spatial_adjust(dat1,feature='VEGFA',plot = T) # spatial map
sc = spatial_cordstat(sa1,sa2,plot=T,plot_bin_col=c('blue','orange'),operator_steps = rep(0,4)) # spatial map for two specified features
spatial_vecPlot(rs3) # spatial vector plot
spatial_vecPlot(rs4) # spatial vector plot

```
