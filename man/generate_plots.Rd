\name{generate_plots}
\alias{generate_plots}
\alias{plot_genome}
\alias{plot_chromosomes}
\alias{plot_zones}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate plots of genome, chromosomes and zones.
}
\description{
  Generate the plot from the processed GenomicOZone dataset object, including genome plots, chromosome plots and zone plots.
}
\usage{
plot_genome(GOZ.ds, plot.file,
            alpha = 0.05, min.effect.size = 0.8,
            plot.width = NULL, plot.height = NULL)

plot_chromosomes(GOZ.ds, plot.file,
            alpha = 0.05, min.effect.size = 0.8,
            plot.width = NULL, plot.height = NULL)

plot_zones(GOZ.ds, plot.file,
            alpha = 0.05, min.effect.size = 0.8,
            log.exp = TRUE, plot.all.zones = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOZ.ds}{
    a GenomicOZong dataset object after running \code{GenomicOZone} function.
  }

  \item{plot.file}{
    a output file name. The file type is "pdf".
  }

  \item{alpha}{
    a cutoff for selecting adjuted \code{p}-values.
  }

  \item{min.effect.size}{
    the minimum effect size required for an outstanding zone. The effect size for ANOVA ranging from 0 to 1 is calculated by R package \pkg{sjstats} \insertCite{sjstats2019}{GenomicOZone}. Default to 0.8.
  }

  \item{plot.width}{
    a numerical number to specify the width of page in the plot. Using \code{NULL} will automatically determine the page width. Default is \code{NULL}.
  }

  \item{plot.height}{
    a numerical number to specify the height of page in the plot. Using \code{NULL} will automatically determine the page height. Default is \code{NULL}.
  }

  \item{log.exp}{
    a logical indicating whether to use log-scaled activity in the plot.
  }

  \item{plot.all.zones}{
    a logical indicating whether to plot all zones into the file. If \code{FALSE}, only outstanding genomic zones will be plotted.
  }
}
\details{
  The three functions plot visualizations of the genome, chromosomes and zones. The R packages \code{ggplot2} \insertCite{Wickham2016ggplot2}{GenomicOZone} and \code{ggbio} \insertCite{yin2012ggbio}{GenomicOZone} are used to generate the plots.

  The function \code{plot_genome} plots the genome-wide overviews with marked significant differential zones.

  The function \code{plot_chromosomes} plots the chromosome-wide heatmap of normalized and linearized activity between sorted zones and samples, visualizing the zones with significant ones marked.

  The function \code{plot_zones} plots the line chart and box-plot of the activity of the genes within each significant zone, visualizing gene activity changes over sample conditions.
}
\value{
  The function takes an input of a object, which has been created by \code{\link{GOZDataSet}} and and processed by \code{\link{GenomicOZone}}. The functions accesse the object and generate visualizations. See \code{\link{GOZDataSet}} for how to create the input object. See \code{\link{GenomicOZone}} for how to process the input object and perform the analysis.
}
\examples{
  # Create an example of GOZ.ds
  data <- matrix(c(1,5,2,6,5,1,6,2), ncol = 2, byrow = TRUE)
  rownames(data) <- paste("Gene", 1:4, sep='')
  colnames(data) <- paste("Sample", c(1:2), sep='')

  colData <- data.frame(Sample_name = paste("Sample", c(1:2), sep=''),
                        Condition = c("Cancer", "Normal"))

  design <- ~ Condition

  rowData.GRanges <- GRanges(seqnames = Rle(rep("chr1", 4)),
                             ranges = IRanges(start = c(1,2,3,4), end = c(5,6,7,8)))
  names(rowData.GRanges) <- paste("Gene", 1:4, sep='')

  ks <- c(2)
  names(ks) <- "chr1"

  GOZ.ds <- GOZDataSet(data, colData, design,
                       rowData.GRanges,
                       ks = ks)
  ####

  # Run the zoing process
  GOZ.ds <- GenomicOZone(GOZ.ds)
  ####

  # Generate plots
  plot_genome(GOZ.ds, plot.file = "Test_genome.pdf",
              plot.width = 15, plot.height = 4)

  plot_chromosomes(GOZ.ds, plot.file = "Test_chromosome.pdf",
                   plot.width = 20, plot.height = 4)

  plot_zones(GOZ.ds, plot.file = "Test_zone.pdf",
             plot.all.zones = FALSE)
}
\references{
  \insertAllCited{}
}
