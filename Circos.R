#========================================================================================
#
#            Learning how to use 'circlize' package to draw Circos plot
#
#========================================================================================

#========================================================================================
#                             Start from Point plot
#========================================================================================


## First, let's hava a look at the testing datasets provided by circlize package
library(circlize)
bed <- generateRandomBed(nc = 2)
head(bed); dim(bed)
### In the 'circlize' input file, the first three columns are always the same in genomics data analysis
### The first three columns are chromosome, start, end.

## Points ##
### When we draw the plot, first we need to is to set parameter:
circos.par(gap.degree = 4)  ## 设置各个板块之间的间隔
circos.par(start.degree = 90)  ## 设置各个板块之间的起始角度

### Next, we need to initialize the plot
circos.initializeWithIdeogram(plotType = NULL) # This step is like ggplot()

### Draw the points in the first track
circos.genomicTrackPlotRegion(
  bed, track.height = 0.1, ylim = c(-1,1), numeric.column = 4,
  panel.fun = function(region, value,...){
    circos.genomicPoints(
      region, value, numeric.column = 1, pch = 16, cex = 0.3)
  }
)
circos.clear() ## always remember to clear so that you can draw another plot



### Draw lines
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicTrack(
  bed, track.height = 0.1, ylim = c(-1,1), numeric.column = 4,
  panel.fun = function(region, value,...){
    circos.genomicLines(
      region, value, numeric.column = 1)
  }
)
circos.clear()


### Draw Rectangle
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicTrack(
  bed, track.height = 0.1, ylim = c(-1,1), numeric.column = 4,
  panel.fun = function(region, value,...){
    circos.genomicRect(
      region, value, numeric.column = 1)
  }
)
circos.clear()



############### You can use your own data and have a try ##############
setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/circos')
beta <- read.csv('beta_lincRNA.csv', header = T, row.names = 1)
dim(beta); head(beta); 
colnames(beta)[4] <- 'value1'
## Draw points
circos.initializeWithIdeogram(species = 'danRer10',chromosome.index = paste0("chr", 1:25))
circos.genomicTrackPlotRegion(
  beta, track.height = 0.1, ylim = c(-0.5,0.5), numeric.column = 4,
  panel.fun = function(region, value,...){
    circos.genomicPoints(
      region, value, numeric.column = 1, pch = 16, cex = 0.3)
  }
)
circos.clear()
View(beta)


#========================================================================================
#
#            Learning how to use 'RCircos' package to draw Circos plot
#
#========================================================================================
## Using demo datasets
### Load packages and inner testing datasets

library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)   #导入内建人类染色体数据
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=10, tracks.outside=0 )
#chr.exclude <- NULL;           设置不显示的染色体，如 c(1,3)          
#cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 设置染色体数据
#tracks.inside <- 10;    设置内部track 个数
#tracks.outside <- 0;    设置外部track 个数 

#导入内建人类染色体数据 (注意此处更换的人类参考基因组)
a <- function(...){
data(UCSC.HG19.Human.CytoBandIdeogram);
head(UCSC.HG19.Human.CytoBandIdeogram);
## 这里换了个参考基因组版本，请注意
chr.exclude <- NULL;   #设置不显示的染色体，如 c(1,3)          
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 
#设置染色体数据
tracks.inside <- 10;   #设置内部track 个数
tracks.outside <- 0;   #设置外部track 个数
#导入上面四个基本参数
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside); 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()}
a()
