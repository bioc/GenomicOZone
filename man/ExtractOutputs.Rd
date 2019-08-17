\name{Extract outputs}
\alias{Extract outputs}
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
\references{
  \insertAllCited{}
}
