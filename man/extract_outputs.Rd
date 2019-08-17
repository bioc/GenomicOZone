\name{extract_outputs}
\alias{extract_outputs}
\alias{extract_genes}
\alias{extract_zones}
\alias{extract_outstanding_zones}
\alias{extract_zone_expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract the gene/zone GRanges object
}
\description{
  Extract the GRanges object from the GenomicOZone function output. The four functions are respectively used to extract the gene GRanges object, zone GRanges object, outstanding genomic zone GRanges object and the zone expression matrix.
}
\usage{
  extract_genes(GOZ.ds)

  extract_zones(GOZ.ds)

  extract_outstanding_zones(
                GOZ.ds,
                p.value.cutoff = 0.05,
                effect.size.rate = 0.05)

  extract_zone_expression(GOZ.ds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GOZ.ds}{
    a GenomicOZong dataset list after running \code{GenomicOZone} function.
  }

  \item{p.value.cutoff}{
    a cutoff for selecting adjuted \code{p}-values.
  }

  \item{effect.size.rate}{
    a top rate for selecting the cutoff of effect size.
  }
}
\details{
  The first three functions respectively return the GRanges \insertCite{lawrence2013granges}{GenomicOZone} object of all genes, all zones and outstanding genomic zones only. The gene GRanges object includes the gene annotations and the zones where the genes belong to. The zone GRanges object includes the zone positions and the p-values of differential expression analysis. The outstanding genomic zone GRanges is a subset of the zone GRanges, containing selected zones only by the cutoffs.

  The function \code{extract_zone_expression} offers the zone expression matrix. The zone expression uses the total abundance of the genes within a zone for each sample. The expression of genes without annotation is appended at the bottom of the zone expression matrix.
}
\value{
  The function takes an input of a list, which has been created by \code{\link{GOZDataSet}} and processed by \code{\link{GenomicOZone}}. The functions accesse the list and fetche the results. See \code{\link{GOZDataSet}} for how to create the input list. See \code{\link{GenomicOZone}} for how to process the input list and perform the analysis.
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

  # Extract outputs
  Gene.GRanges <- extract_genes(GOZ.ds)
  head(Gene.GRanges)

  Zone.GRanges <- extract_zones(GOZ.ds)
  head(Zone.GRanges)

  OZone.GRanges <- extract_outstanding_zones(GOZ.ds)
  head(OZone.GRanges)

  Zone.exp.mat <- extract_zone_expression(GOZ.ds)
  head(Zone.exp.mat)
}
\references{
  \insertAllCited{}
}
