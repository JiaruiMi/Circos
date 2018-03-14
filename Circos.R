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
bed <- generateRandomBed(nc = 1)
head(bed)
### In the 'circlize' input file, the first three columns are always the same in genomics data analysis
### The first three columns are chromosome, start, end.

## Points ##
### When we draw the plot, first we need to is to set parameter:
circos.par(gap.degree = 4)  ## 设置各个板块之间的间隔
circos.par(start.degree = 90)  ## 设置各个板块之间的起始角度

### Next, we need to initialize the plot
circos.initializeWithIdeogram() # This step is like ggplot()

### draw the points in the first track
circos.genomicTrackPlotRegion(
  bed, track.height = 0.1, ylim = c(-1,1),
  panel.fun = function(region, value,...){
    circos.genomicPoints(
      region, value, pch = 16, cex = 0.1)
  }
)
circos.clear() ## always remember to clear so that you can draw another plot

## Second, Let's try a list of datasets
bedlist = list(bed01 = generateRandomBed(nc = 2, nr = 200),
               bed02 = generateRandomBed(nc = 2, nr = 200))
bedlist

circos.initializeWithIdeogram(plotType = NULL)
circos.genomicTrackPlotRegion(
  bedlist, track.height = 0.1, ylim = c(-1,1),
  numeric.column = c(4,5),
  panel.fun = function(region, value, ...){
    i = getI(...)
    circos.genomicPoints(
      region, value, pch = 16, cex = 0.1, col = i)
  }
)











