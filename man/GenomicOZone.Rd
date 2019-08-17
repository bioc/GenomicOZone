\name{GenomicOZone}
\alias{GenomicOZone}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Perform genomic outstanding zoning process
}
\description{
  The function performs the outstanding genomic zoning along the chromosomes, and evaluates the differentially expressed zones with consistent expression pattern of the within genes among the samples.
}
\usage{
  GenomicOZone(GOZ.ds)
}
\arguments{
  \item{GOZ.ds}{
    an input list created by function \code{GOZDataSet}.
  }
}
\details{
  This is the main function of the package. The function integrates gene annotation, expression matrix preprocess, chromosome zoning, and zone fifferential analysis.

  The gene annotation can be specified either by the user, or obtained from the R package \code{biomaRt} \insertCite{smedley2015biomart}{GenomicOZone} to access ensembl annotation databases \insertCite{zerbino2017ensembl}{GenomicOZone}.

  The function performs the weighted univariate clustering method Ckmeans \insertCite{wang2011ckmeans}{GenomicOZone}. If \code{ks} is specified, fixed number of zones at each chromosome will be received. If \code{ks} is \code{NULL}, optimal number of clusters at each chromosome will be determined by dynamic programming.

  The function also conducts zone differential expression analysis by using one-way ANOVA \insertCite{chambers1992statistical}{GenomicOZone} on gene ranks. Given \code{p}-value cutoff and top rate of effect size, outstanding genomic zones will be selected.
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
  GenomicOZone(GOZ.ds)
  ####
}

\references{
  \insertAllCited{}
}
