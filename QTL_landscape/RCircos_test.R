library(RCircos)

data(UCSC.HG19.Human.CytoBandIdeogram)
chr.exclude <- NULL
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)




out.file <- "RCircosDemoHumanGenome.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data(RCircos.Gene.Label.Data)
name.col <- 4
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num, side)

data(RCircos.Heatmap.Data);
data.col <- 6;
track.num <- 5;
side <- "in";
RCircos.Heatmap.Plot(RCircos.Heatmap.Data, data.col,track.num, side);

data(RCircos.Scatter.Data);
data.col <- 5;
track.num <- 6;
side <- "in";
by.fold <- 1;
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold);


data(RCircos.Line.Data);
data.col <- 5;
track.num <- 7;
side <- "in";
RCircos.Line.Plot(RCircos.Line.Data, data.col,track.num, side);


data(RCircos.Histogram.Data);
data.col <- 4;
track.num <- 8;
side <- "in";
RCircos.Histogram.Plot(RCircos.Histogram.Data,data.col, track.num, side)

data(RCircos.Tile.Data);
track.num <- 9;
side <- "in";
RCircos.Tile.Plot(RCircos.Tile.Data, track.num, side);

data(RCircos.Link.Data);
track.num <- 11;
RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
data(RCircos.Ribbon.Data);
RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data,track.num=11, by.chromosome=FALSE, twist=FALSE);
dev.off()


library(RCircos);
data(UCSC.HG19.Human.CytoBandIdeogram);
data(RCircos.Line.Data);
RCircos.Line.Data$chromosome = paste0('chr',RCircos.Line.Data$chromosome )
RCircos.Set.Core.Components(
cyto.info=UCSC.HG19.Human.CytoBandIdeogram,
chr.exclude=c("chrX", "chrY"),
tracks.inside=10, tracks.outside=5)
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot()
RCircos.Vertical.Line.Plot(RCircos.Line.Data[!RCircos.Line.Data$chromosome %in% c('chrX','chrY'),], track.num=1, side="in")
RCircos.Vertical.Line.Plot(RCircos.Line.Data[!RCircos.Line.Data$chromosome %in% c('chrX','chrY'),], side="in",inside.pos=1.5, outside.pos=1.75)


line.data = RCircos.Line.Data[!RCircos.Line.Data$chromosome %in% c('chrX','chrY'),]
side = 'in'
inside.pos = NULL
outside.pos = NULL
track.num = 2
line.width = 1
genomic.columns = 3
RCircos.Vertical.Line.Plot <- function(line.data=NULL, track.num=NULL, 
        side=c("in", "out"), line.width=1, inside.pos=NULL, outside.pos=NULL,
        genomic.columns=3, is.sorted=TRUE)
{
    if(is.null(line.data)) 
        stop("Genomic data missing in RCircos.Vertical.Line.Plot.\n");
    if(is.null(genomic.columns) || genomic.columns < 2) 
        stop("Missing number of columns for genomic position.\n");

    boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                outside.pos, FALSE);
    outerPos <- boundary[1];
    innerPos  <- boundary[2];

    if(length(line.width) == 1)
        line.width <- rep(line.width, nrow(line.data));

    line.data <- RCircos.Get.Single.Point.Positions(line.data, 
                        genomic.columns);
    point.index <- as.numeric(line.data$Location);
    
    RCircos.Pos <- RCircos.Get.Plot.Positions();
    RCircos.Par <- RCircos.Get.Plot.Parameters();
    line.colors <- RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color);
    line.colors <- RCircos.Get.Plot.Colors(line.data, 'red');

    for(a.point in seq_len((nrow(line.data))))
    {
        lines(  c(RCircos.Pos[point.index[a.point] ,1]*innerPos,
                    RCircos.Pos[point.index[a.point], 1]*outerPos),
                c(  RCircos.Pos[point.index[a.point], 2]*innerPos,
                    RCircos.Pos[point.index[a.point], 2]*outerPos),
                col=line.colors[a.point], lwd=line.width[a.point]);
    }
}

