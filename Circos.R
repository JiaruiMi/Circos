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
### circos 本来是一个perl程序，专门用来genomic features在全基因组的分布，连接的。 所以输入数据必然要
### 包含染色体，起始终止坐标的（即3列内容）。 然后可以是一列value用来调色或者高度信息，可以加入基因
### 名字（比如标签）或数值（散点图，热图，柱状图，线图）。如果有两个genomic features的坐标信息，就是
### 连接图啦（共六列）。


## Using demo datasets
### Load packages and inner testing datasets

### 最基本的Circos图
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)   #导入内建人类染色体数据
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL,tracks.inside=10, tracks.outside=0 )
#chr.exclude <- NULL;           设置不显示的染色体，如 c(1,3)          
#cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 设置染色体数据
#tracks.inside <- 10;    设置内部track 个数
#tracks.outside <- 0;    设置外部track 个数 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot() #绘制染色体表意文字，染色体的默认方法亮点和染色体名称。


### 复杂化
####导入内建人类染色体数据 (注意此处更换的人类参考基因组)
a <- function(...){
data(UCSC.HG19.Human.CytoBandIdeogram);
head(UCSC.HG19.Human.CytoBandIdeogram);
#### 这里换了个参考基因组版本，请注意
chr.exclude <- NULL;   #设置不显示的染色体，如 c(1,3)          
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 
#### 设置染色体数据
tracks.inside <- 10;   #设置内部track 个数
tracks.outside <- 0;   #设置外部track 个数
#### 导入上面四个基本参数
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside); 
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()}
a()

### 在基因组上面根据坐标标记基因
b <- function(...){
  a()
  #Gene Labels and connectors on RCircos Plot
  #RCircos.Gene.Connector.Plot   绘制染色体表意文字和基因名称之间的连接
  #RCircos.Gene.Name.Plot 在数据轨道上绘制基因名称
  #从图上我们看出150-300个gene用于绘制是比较理想的，不会太稀疏也不会太密
  data(RCircos.Gene.Label.Data);
  name.col <- 4;
  side <- "in";
  track.num <- 1;
  RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side);
  track.num <- 2;
  RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col,track.num, side);
}
b()


### 用热图来标记基因组每个区域的数据，一般是拷贝数变异等数据，用RCircos.Heatmap.Plot函数绘制一个数据轨道
### 的热图
data(RCircos.Heatmap.Data)
head(RCircos.Heatmap.Data)
c <- function(...){
  b()
  data.col <- 6;
  track.num <- 5;  ## 给字多留一点空间，所以从track.num = 1直接跳到track.num = 5
  side <- "in";
  RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col, track.num, side);
  
}
c()



### 加上散点图，用RCircos.Scatter.Plot函数来添加一个数据轨迹的扫描图
data(RCircos.Scatter.Data)
head(RCircos.Scatter.Data)
d <- function(...){
  c()
  RCircos.Scatter.Data$chromosome=paste0('chr',RCircos.Scatter.Data$chromosome)
  head(RCircos.Scatter.Data);
  data.col <- 5;
  track.num <- 6;
  side <- "in";
  by.fold <- 1;
  RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col,track.num, side, by.fold);
}
d()



### 加上直线，用RCircos.Line.Plot函数绘制线为一个数据轨道
data(RCircos.Line.Data)
head(RCircos.Line.Data)
e <- function(...){
  d()
  RCircos.Line.Data$chromosome=paste0('chr',RCircos.Line.Data$chromosome)
  head(RCircos.Line.Data);
  data.col <- 5;
  track.num <- 7;
  side <- "in";
  RCircos.Line.Plot(RCircos.Line.Data, data.col, track.num, side);
  
}
e()



### 加上直方图,这里用RCircos.Histogram.Plot 函数绘制一个数据轨迹的直方图
data(RCircos.Histogram.Data)
head(RCircos.Histogram.Data)
f <- function(...){
  e()
  data.col <- 4;
  track.num <- 8;
  side <- "in";
  RCircos.Histogram.Plot(RCircos.Histogram.Data, data.col, track.num, side);
  
}
f()



### 最后把有连接关系的基因连线,用RCircos.Link.Plot函数绘制两个或多个基因组位置之间的链接线，请自行查看
### RCircos.Link.Data数据是什么，如何映射到这个circos图上的
data(RCircos.Link.Data)
head(RCircos.Link.Data)
data(RCircos.Ribbon.Data)
head(RCircos.Ribbon.Data)  ## ribbon与link的区别在于，ribbon适合于比较宽的基因组区域的展示
h <- function(...){
  g()
  track.num <- 11;
  RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE); ## 注意TRUE的含义，查看帮助文档得知，同一染色体内部和不同染色体之间的关系不同
  data(RCircos.Ribbon.Data);
  RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, by.chromosome=FALSE, twist=FALSE)
}
h()

