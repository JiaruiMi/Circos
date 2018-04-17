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


### When we draw the plot, first we need to is to set parameter (initialize circos layout):
### You cannot see anything after running following code because no track has been added yet.
### Do remember when you explicitly set circos.par(), you need to call circos.clear() to finish the plotting.
####circos.par(gap.degree = 4)  ## 设置各个板块之间的间隔
circos.par(start.degree = 90)  ## 设置各个板块之间的起始角度, 默认值是起始3点钟
circos.par("track.height" = 0.2) ## 设置板块的高度

### Next, we need to initialize the plot
circos.initializeWithIdeogram(plotType = NULL) 

## Rules for making circular plot
## It follows the sequence of initialize layout -> create track -> add graphics -> create track -> add graphics 
## -> ... -> clear. Graphics can be added at any time as long as the tracks are created.

## Create plotting regions for the new track and add graphics. The new track is created just inside the 
## previously created one. Only after the creation of the track can you add other graphics on it. There are 
## three ways to add graphics in cells.
## Use panel.fun argument in circos.track() to add graphics immediately after the creation of a certain 
## cell. panel.fun needs two arguments x and y which are x values and y values that are in the current cell. 
## This subset operation is applied automatically. This is the most recommended way. Section 2.7 gives detailed explanation of using panel.fun argument.

## Points ##
### Draw the points in the first track
### After the circular layout is initialized, graphics can bd added to the plot in a track-by-track manner.
### Before drawing anything, we need to know that all tracks should be first created by circos.trackPlotRegion()
### In this case, "circos.genomicTrackPlotRegion", then the low-level functions can be added afterwards.
### # This step is like ggplot(), or plot() in base Graphics. Then you can use functions like points() or 
### lines() to add graphics.
### We should know that if only the points are needed to put in thecells, it would be more convenient to use
### circos.trackPoints(), because it is quite straightforward to understand that this function needs a 
### categorical variable (df$factors), values on x direction and y direction (df$x and df$y)
circos.genomicTrackPlotRegion(
  bed, ylim = c(-1,1), numeric.column = 4,
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


### Draw Rectangle(draw heatmaps)
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
alpha <- read.table('alpha_linc.txt', header = T, sep = '\t', stringsAsFactors = F)
beta <- read.table('beta_linc.txt', header = T, sep = '\t', stringsAsFactors = F)
delta <- read.table('delta_linc.txt', header = T, sep = '\t', stringsAsFactors = F)
acinar <- read.table('acinal_linc.txt', header = T, sep = '\t', stringsAsFactors = F)
### For genomics data, we expect all input data follow a data.frame structure, with the variables as
### "chr", "start", "end", and "value", we consider this form as pseudo-bed file format
dim(beta); head(beta); 
alpha$exprs <- log2(alpha$exprs+1) / log2(max(alpha$exprs)) * 10000
beta$exprs <- log2(beta$exprs+1) / log2(max(beta$exprs)) * 10000
delta$exprs <- log2(delta$exprs+1) / log2(max(delta$exprs)) * 10000
acinar$exprs <- log2(acinar$exprs+1) / log2(max(acinar$exprs)) * 10000
## Initialize cytoband
### By default, circos.initializeWithIdeogram() initializes the plot with cytoband data of human genome hg19. 
### Users can also initialize with other species by specifying species argument and it will automatically 
### download cytoband files for corresponding species. When you are dealing rare species and there is no 
### cytoband data available yet, circos.initializeWithIdeogram() will try to continue to download the 
### “chromInfo” file form UCSC, which also contains lengths of chromosomes, but of course, there is no 
### ideogram track on the plot.
### When there is no cytoband data for the specified species, and when chromInfo data is used instead, there
### may be many many extra short contigs. chromosome.index can also be useful to remove unnecessary contigs.

### 设置默认参数
library(gridBase)
plot.new()
circle_size = unit(1, "snpc") # snpc unit gives you a square region

### For circos
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)

circos.par(start.degree = 90)  
circos.par("track.height" = 0.12)

### 图像初始化
circos.initializeWithIdeogram(species = 'danRer10',chromosome.index = paste0("chr", 1:25))

### The labels can only be located either inside or outside of cytoband, but always outside of other tracks
# bed <- read.csv('beta_lincRNA_DGE.csv', header = T)
# head(bed); dim(head)
# circos.genomicLabels(bed, labels.column = 4, side = "inside")

### Draw lineplots
circos.genomicTrack(
  alpha, ylim = c(min(beta$exprs),max(beta$exprs)), numeric.column = 4, 
  panel.fun = function(region, value,...){
    circos.genomicLines(
      region, value, col = 'blue')
  })

circos.genomicTrack(
  beta, ylim = c(min(beta$exprs),max(beta$exprs)), numeric.column = 4, 
  panel.fun = function(region, value,...){
    circos.genomicLines(
      region, value, col = 'orange')
  })

circos.genomicTrack(
  delta, ylim = c(min(beta$exprs),max(beta$exprs)), numeric.column = 4, 
  panel.fun = function(region, value,...){
    circos.genomicLines(
      region, value, col = 'skyblue')
  })


circos.genomicTrack(
  acinar, ylim = c(min(beta$exprs),max(beta$exprs)), numeric.column = 4, 
  panel.fun = function(region, value,...){
    circos.genomicLines(
      region, value, col = 'brown')
  })

### Legends cannot be automatically generated by circlize, by using "ComplexHeatmap" package, it is just
### a few more lines to really implement it.In "ComplexHeatmap" package, there is a legend() function which
### customizes legends with various styles. In the following codes, lgd_lines is "grob" object (graphical
### object) defined by grid package, which you can think as boxes which containall graphical elements for 
### legends and they can be added to the plot by grid.draw()
### circlize is implemented in the base graphic system while ComplexHeatmap is implemented by grid graphic
### system. However, these two systems can be mixed somehow. We can directly add grid graphics to the base
### graphics. (Actually they are two independent layers but drawn on a same graphic device  )
library(ComplexHeatmap)
lgd_lines = Legend(at = c("alpha", "beta","delta", "acinar"), type = "lines", 
                   legend_gp = gpar(col = c("blue", "orange", "skyblue", "brown"), lwd = 4), title_position = "topleft", 
                   title = "Cell types", nrow = 1)
pushViewport(viewport(x = 0.5, y = unit(1, "npc") - circle_size, width = grobWidth(lgd_lines), 
                      height = grobHeight(lgd_lines), just = c("center","bottom")))
grid.draw(lgd_lines)

upViewport()
circos.clear()



## Draw points
circos.genomicTrackPlotRegion(
  beta, ylim = c(min(beta$exprs),max(beta$exprs)), numeric.column = 4,
  panel.fun = function(region, value,...){
    circos.genomicPoints(
      region, value, pch = 16, cex = 0.3)
  }
)

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

summary(RCircos.Heatmap.Data$A498)

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


################ use your own datasets to draw Rcircos plot ###############
#### 1, 整理表格最重要，三张表（染色体长度文件，txt文件和link文件），染色体长度文件是chr，start，end，band(和chr一样)，stain
#### txt文件前四列是chromosome，start，end，geneid（name），后面可以跟一列或者多列，用于画heatmap，scatterplot，lineplot，histogram
#### link文件是chr1，start1，end1，chr2，start2，end2这样一共六列
#### 2, chromosome需要加chr标签
#### 3, 删除所有没有必要的染色体和基因信息(chrKN...)
#### 4, 基因标签要优中选优，50-100足够了

setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/circos')
library(RCircos)

#### To draw Circos plot, we need to prepare three files, the file with chromosome information (if without
#### a complete cytoband available in UCSC, you can edit it a little bit), bed file and link file
#### The chromosome information is stored in danRer10_ref.csv(with chromosome, chromstart, chromend, band
#### and stain, in total five columns);
#### It is better to prepare file with different information into different files

cytoBandIdelgram = read.csv('danRer10_ref.csv', header = T)
chr.exclude <- NULL
cyto.info <- cytoBandIdelgram
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
RCircos.List.Plot.Parameters() # show inner parameters, e.g. text.size, heatmap.color

### Modifying RCircos core components
rcircos.param <- RCircos.Get.Plot.Parameters()
rcircos.param$text.size <- 0.5
RCircos.Reset.Plot.Parameters(rcircos.param)  ## 确认重新设置，这个时候text.size = 0.5
RCircos.List.Plot.Parameters()

# pdf(file = 'RCircosOut.pdf', height = 800, weight = 800)
RCircos.Set.Plot.Area()  ##设置画布区域，每次重新画图，哪怕是修改都有重新铺画板
RCircos.Chromosome.Ideogram.Plot()  ## 画染色体的最外圈

### Gene labels, the number should be limited, no more than 100
#### There are some transcripts located in unknown gene areas, and we need to get rid of them in advance
#### We use dplyr to pick rows randomly
#### The chromosome information should have 'chr' label on it
####

library(dplyr)
RCircos.Gene.Label.Data <- read.csv('beta_lincRNA.csv', header = T)
RCircos.Gene.Label.Data <- sample_n(tbl = RCircos.Gene.Label.Data, size = 100, replace = F)
table(RCircos.Gene.Label.Data$Chr)
str(RCircos.Gene.Label.Data)
name.col <- 4
side <- 'in'
track.num <- 1
RCircos.Gene.Connector.Plot(genomic.data = RCircos.Gene.Label.Data, # 画连接器
                            track.num = track.num, side = side)
track.num <- 2
RCircos.Gene.Name.Plot(gene.data = RCircos.Gene.Label.Data, 
                       name.col = name.col, track.num = track.num, side = side)

#### 注意，对于每一条染色体，都有一个最大基因标注数目的，我们可以使用RCircos.Get.Gene.Name.Plot.Parameters()
#### 来查看，我们可以看到大部分基因只能标注2个，少数可以标注4个，只有4号染色体可以标注4个，所以如果一开始选择
#### 标注的基因比较多的话，软件会给出提示警告，告诉作者超出数目的基因不会被标注出来; 所以建议作者在每次输入
#### 数据的时候，需要对表格做出调整，只选择那些最重要的基因进行标注
RCircos.Get.Gene.Name.Plot.Parameters()
list <- RCircos.Get.Gene.Name.Plot.Parameters()
sum(list$maxLabels) ## 最多标注58个基因


### Heatmap, 所有的表格前四列非常明确，chromosome，start，end，geneid
RCircos.Heatmap.Data <- read.csv('beta_lincRNA_expr.csv', header = T)
data.col <- 5
track.num <- 5 # 空一些距离，以免与label重叠
side <- 'in'
summary(RCircos.Heatmap.Data$beta_mean)
RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col = data.col, track.num = track.num, side = side, is.sorted = F)


### Scatterplot
RCircos.Scatter.Data <- read.csv('beta_lincRNA_expr.csv', header = T)
data.col <- 5
track.num <- 6
side <- 'in'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Scatter.Plot(scatter.data = RCircos.Scatter.Data, data.col = data.col,
                     track.num = track.num, by.fold = by.fold, is.sorted = F)


### Line
RCircos.Line.Data <- read.csv('beta_lincRNA_expr.csv', header = T)
data.col <- 5
track.num <- 7
side <- 'in'
RCircos.Line.Plot(line.data = RCircos.Line.Data, data.col = data.col,
                     track.num = track.num, is.sorted = F)


### Histogram
RCircos.Histogram.Data <- read.csv('beta_lincRNA_expr.csv', header = T)
data.col <- 5
track.num <- 8
side <- 'in'
RCircos.Histogram.Plot(hist.data = RCircos.Histogram.Data, data.col = data.col, track.num = track.num, side = side)


### link
RCircos.Link.Data <- read.csv('link.csv', header = T)
track.num <- 10
RCircos.Link.Plot(RCircos.Link.Data, track.num = track.num, by.chromosome = F, is.sorted = F) ## by.chromosome = T, 染色体之间是蓝色，染色体内部是红色，建议false


### Ribbon (一般需要10MB以上才有必要用ribbon)
RCircos.Ribbon.Data <- read.csv('ribbon.csv', header = T)
RCircos.Ribbon.Plot(ribbon.data = RCircos.Ribbon.Data, track.num = 10, is.sorted = F)



#======================================================================================
#
# Let's include data from different cell types, the order is acinal, alpha, beta, delta
#                                       Dotplot
#
#======================================================================================

setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/circos')
library(RCircos)

#### To draw Circos plot, we need to prepare three files, the file with chromosome information (if without
#### a complete cytoband available in UCSC, you can edit it a little bit), bed file and link file
#### The chromosome information is stored in danRer10_ref.csv(with chromosome, chromstart, chromend, band
#### and stain, in total five columns);
#### It is better to prepare file with different information into different files

cytoBandIdelgram = read.csv('danRer10_ref.csv', header = T)
chr.exclude <- NULL
cyto.info <- cytoBandIdelgram
tracks.inside <- 3
tracks.outside <- 4
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
RCircos.List.Plot.Parameters() # show inner parameters, e.g. text.size, heatmap.color

### Modifying RCircos core components
rcircos.param <- RCircos.Get.Plot.Parameters()
rcircos.param$text.size <- 0.5
RCircos.Reset.Plot.Parameters(rcircos.param)  ## 确认重新设置，这个时候text.size = 0.5
RCircos.List.Plot.Parameters()

# pdf(file = 'RCircosOut.pdf', height = 800, weight = 800)
RCircos.Set.Plot.Area()  ##设置画布区域，每次重新画图，哪怕是修改都有重新铺画板
RCircos.Chromosome.Ideogram.Plot()  ## 画染色体的最外圈

### Gene labels, the number should be limited, no more than 100
#### There are some transcripts located in unknown gene areas, and we need to get rid of them in advance
#### We use dplyr to pick rows randomly
#### The chromosome information should have 'chr' label on it
####

library(dplyr)
RCircos.Gene.Label.Data <- read.csv('beta_lincRNA_DGE.csv', header = T)
table(RCircos.Gene.Label.Data$Chr)
str(RCircos.Gene.Label.Data)
name.col <- 4
side <- 'in'
track.num <- 1
RCircos.Gene.Connector.Plot(genomic.data = RCircos.Gene.Label.Data, # 画连接器
                            track.num = track.num, side = side)
track.num <- 2
RCircos.Gene.Name.Plot(gene.data = RCircos.Gene.Label.Data, 
                       name.col = name.col, track.num = track.num, side = side)

#### 注意，对于每一条染色体，都有一个最大基因标注数目的，我们可以使用RCircos.Get.Gene.Name.Plot.Parameters()
#### 来查看，我们可以看到大部分基因只能标注2个，少数可以标注4个，只有4号染色体可以标注4个，所以如果一开始选择
#### 标注的基因比较多的话，软件会给出提示警告，告诉作者超出数目的基因不会被标注出来; 所以建议作者在每次输入
#### 数据的时候，需要对表格做出调整，只选择那些最重要的基因进行标注
RCircos.Get.Gene.Name.Plot.Parameters()
list <- RCircos.Get.Gene.Name.Plot.Parameters()
sum(list$maxLabels) ## 最多标注58个基因



### beta-cell
RCircos.Scatter.Data <- read.csv('beta_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data$exprs <- log2(RCircos.Line.Data$exprs+1) / log2(max(RCircos.Line.Data$exprs)) * 10000
data.col <- 4
track.num <- 1
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Scatter.Plot(scatter.data = RCircos.Scatter.Data, data.col = data.col,
                     track.num = track.num, by.fold = by.fold, is.sorted = F, side = side)


### alpha
RCircos.scatter.Data1 <- read.table('alpha_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data1$exprs <- log2(RCircos.Line.Data1$exprs+1) / log2(max(RCircos.Line.Data1$exprs)) * 10000
data.col <- 4
track.num <- 2
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Scatter.Plot(scatter.data = RCircos.scatter.Data1, data.col = data.col,
                     track.num = track.num, by.fold = by.fold, is.sorted = F, side = side)


### delta
RCircos.scatter.Data2 <- read.table('delta_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data2$exprs <- log2(RCircos.Line.Data2$exprs+1) / log2(max(RCircos.Line.Data2$exprs)) * 10000
data.col <- 4
track.num <- 3
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Scatter.Plot(scatter.data = RCircos.scatter.Data2, data.col = data.col,
                     track.num = track.num, by.fold = by.fold, is.sorted = F, side = side)

### acinal
RCircos.scatter.Data3 <- read.table('acinal_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data3$exprs <- log2(RCircos.Line.Data3$exprs+1) / log2(max(RCircos.Line.Data3$exprs)) * 10000
data.col <- 4
track.num <- 4
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Scatter.Plot(scatter.data = RCircos.scatter.Data3, data.col = data.col,
                     track.num = track.num, by.fold = by.fold, is.sorted = F, side = side)







#======================================================================================
#
# Let's include data from different cell types, the order is acinal, alpha, beta, delta
#                                       Lineplot
#
#======================================================================================

setwd('/Users/mijiarui/R_bioinformatics_project/Master_thesis_project/circos')
library(RCircos)

#### To draw Circos plot, we need to prepare three files, the file with chromosome information (if without
#### a complete cytoband available in UCSC, you can edit it a little bit), bed file and link file
#### The chromosome information is stored in danRer10_ref.csv(with chromosome, chromstart, chromend, band
#### and stain, in total five columns);
#### It is better to prepare file with different information into different files

cytoBandIdelgram = read.csv('danRer10_ref.csv', header = T)
chr.exclude <- NULL
cyto.info <- cytoBandIdelgram
tracks.inside <- 1
tracks.outside <- 5
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)
RCircos.List.Plot.Parameters() # show inner parameters, e.g. text.size, heatmap.color

### Modifying RCircos core components
rcircos.param <- RCircos.Get.Plot.Parameters()
rcircos.param$text.size <- 0.6
RCircos.Reset.Plot.Parameters(rcircos.param)  ## 确认重新设置，这个时候text.size = 0.5
RCircos.List.Plot.Parameters()

# pdf(file = 'RCircosOut.pdf', height = 800, weight = 800)
RCircos.Set.Plot.Area()  ##设置画布区域，每次重新画图，哪怕是修改都有重新铺画板
RCircos.Chromosome.Ideogram.Plot()  ## 画染色体的最外圈

### Gene labels, the number should be limited, no more than 100
#### There are some transcripts located in unknown gene areas, and we need to get rid of them in advance
#### We use dplyr to pick rows randomly
#### The chromosome information should have 'chr' label on it
####

library(dplyr)
RCircos.Gene.Label.Data <- read.csv('beta_lincRNA_DGE.csv', header = T)
RCircos.Gene.Label.Data
table(RCircos.Gene.Label.Data$Chr)
str(RCircos.Gene.Label.Data)
name.col <- 4
side <- 'in'
track.num <- 1
RCircos.Gene.Connector.Plot(genomic.data = RCircos.Gene.Label.Data, # 画连接器
                            track.num = track.num, side = side)
side <- 'in'
track.num <- 2
RCircos.Gene.Name.Plot(gene.data = RCircos.Gene.Label.Data, 
                       name.col = name.col, track.num = track.num, side = side)

#### 注意，对于每一条染色体，都有一个最大基因标注数目的，我们可以使用RCircos.Get.Gene.Name.Plot.Parameters()
#### 来查看，我们可以看到大部分基因只能标注2个，少数可以标注4个，只有4号染色体可以标注4个，所以如果一开始选择
#### 标注的基因比较多的话，软件会给出提示警告，告诉作者超出数目的基因不会被标注出来; 所以建议作者在每次输入
#### 数据的时候，需要对表格做出调整，只选择那些最重要的基因进行标注
RCircos.Get.Gene.Name.Plot.Parameters()
list <- RCircos.Get.Gene.Name.Plot.Parameters()
sum(list$maxLabels) ## 最多标注58个基因



### beta-cell
RCircos.Line.Data <- read.table('beta_linc.txt', header = T, sep = '\t', quote = "")
max(RCircos.Line.Data$exprs)
RCircos.Line.Data$exprs <- log2(RCircos.Line.Data$exprs+1) / log2(max(RCircos.Line.Data$exprs)) * 10000
data.col <- 4
track.num <- 1
side <- 'out'
by.fold <- 1 ## by.fold是针对scatterplot的一个属性，在lineplot中没有
RCircos.Line.Plot(line.data = RCircos.Line.Data, data.col = data.col,
                     track.num = track.num, is.sorted = F, side = side)
max(RCircos.Line.Data3$V4 )

### alpha
RCircos.Line.Data1 <- read.table('alpha_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data1$exprs <- log2(RCircos.Line.Data1$exprs+1)/ log2(max(RCircos.Line.Data1$exprs)) * 10000
data.col <- 4
track.num <- 2
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Line.Plot(line.data = RCircos.Line.Data1, data.col = data.col,
                     track.num = track.num, is.sorted = F, side = side)


### delta
RCircos.Line.Data2 <- read.table('delta_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data2$exprs <- log2(RCircos.Line.Data2$exprs+1) / log2(max(RCircos.Line.Data2$exprs))* 10000
data.col <- 4
track.num <- 3
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Line.Plot(line.data = RCircos.Line.Data2, data.col = data.col,
                     track.num = track.num,  is.sorted = F, side = side)

### acinal
RCircos.Line.Data3 <- read.table('acinal_linc.txt', header = T, sep = '\t', quote = "")
RCircos.Line.Data3$exprs <- log2(RCircos.Line.Data3$exprs+1) / log2(max(RCircos.Line.Data3$exprs))* 10000
data.col <- 4
track.num <- 4
side <- 'out'
by.fold <- 1 ## 设定1为颜色区分, by.fold > 1 为红色，< -1为蓝色，-1~1为黑色，可以用来标出高表达；默认值是0
RCircos.Line.Plot(line.data = RCircos.Line.Data3, data.col = data.col,
                     track.num = track.num,  is.sorted = F, side = side)
summary(RCircos.Line.Data)
summary(RCircos.Line.Data1)
summary(RCircos.Line.Data2)
summary(RCircos.Line.Data3)
