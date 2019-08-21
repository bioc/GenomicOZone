\name{extract_outputs}
\alias{extract_outputs}
\alias{extract_genes}
\alias{extract_zones}
\alias{extract_outstanding_zones}
\alias{extract_zone_expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract annotation of genes, zones, outstanding zones, and expression of zones
}
\description{
  Extract informtation from the output of \code{GenomicOZone} function, including a gene annotation object, a zone annotation object, an outstanding zone annotation object, or a zone activity matrix, respectively.
}
\usage{
  extract_genes(GOZ.ds)

  extract_zones(GOZ.ds)

  extract_outstanding_zones(
                GOZ.ds,
                alpha = 0.05,
                min.effect.size = 0.8)

  extract_zone_expression(GOZ.ds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOZ.ds}{
    a object returned from the \code{GenomicOZone} function.
  }

  \item{alpha}{
    a cutoff for adjusted \code{p}-values. Default to 0.05.
  }

  \item{min.effect.size}{
    the minimum effect size required for an outstanding zone. The effect size for ANOVA ranging from 0 to 1 is calculated by R package \pkg{sjstats} \insertCite{sjstats2019}{GenomicOZone}. Default to 0.8.
  }
}
\details{
  These functions take the input of an object created by \code{\link{GOZDataSet}} and processed by \code{\link{GenomicOZone}}. The functions access the object and fetch the results.
  The function \code{extract_zone_expression} offers the zone activity matrix. The activity of a zone is the total activity of genes within the zone for each sample. The activity of genes without annotation is included as last rows in the zone activity matrix.
}
\value{
  The first three functions return an object of \code{GRanges} class \insertCite{lawrence2013granges}{GenomicOZone} for all genes, all zones and outstanding genomic zones only. The gene \code{GRanges} object includes genome annotation and the zones where the genes belong. The zone \code{GRanges} object includes zone positions and \code{p}-values of differential zone analysis. Outstanding genomic zones are a subset of all zones that satisfy required \code{p}-value and effect size.
}
\seealso{
  See \code{\link{GOZDataSet}} for how to create the input object before outstanding genomic zone analysis. The object must contain information obtained from outstanding zone analysis function \code{\link{GenomicOZone}}.
}
\examples{
  # Create an object of GOZ.ds
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

  # Run outstanding zone analysis
  GOZ.ds <- GenomicOZone(GOZ.ds)
  ####

  # Extract output in various formats
  Gene.GRanges <- extract_genes(GOZ.ds)
  head(Gene.GRanges)

  Zone.GRanges <- extract_zones(GOZ.ds)
  head(Zone.GRanges)

  OZone.GRanges <- extract_outstanding_zones(
                                GOZ.ds,
                                alpha = 0.05,
                                min.effect.size = 0.8)
  head(OZone.GRanges)

  Zone.exp.mat <- extract_zone_expression(GOZ.ds)
  head(Zone.exp.mat)
}
\references{
  \insertAllCited{}
}
