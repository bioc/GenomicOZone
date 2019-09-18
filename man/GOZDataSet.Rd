\name{GOZDataSet}
\alias{GOZDataSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create an object for outstanding genomic zone analysis
}
\description{
  The function prepares an object for outstanding genomic zone analysis. It integrates data, annotation, and analysis parameters into the object and performs additional check on data integrity.
}
\usage{
  GOZDataSet(data, colData, design,
             clustering.method = "1C",
             rowData.GRanges = NULL,
             ks = NULL,
             genome = NULL,
             ensembl.mirror = "www",
             gene.ID.type = NULL,
             ncores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    a numerical \code{matrix} of gene activity data. Rows represent genes. Columns represent samples. The activity data can be gene expression, methylation beta value, copy number variation segment mean, or other gene-based omic data. Row and column names of the matrix must be specified.
  }

  \item{colData}{
    a \code{dataframe} of sample information. The first column must be sample names, corresponding to columns in the \code{data} matrix. Each additional column must contain at least two experimental conditions necessary for differential zone analysis.
  }

  \item{design}{
    a one-sided \code{formula} with only right-hand side variables. Only one variable is supported in this version. The formula describes on which variable in \code{colData} to apply differential zone analysis.
  }

  \item{clustering.method}{
    a chracter string. An option to choose either using \code{"1C"} to accumulate all channels of weight into one channel, or using \code{"MC"} to allow multi channel of weight in the clustering. Default is "1D".
  }

  \item{rowData.GRanges}{
    an optional genome annotation of \code{GRanges} class. Rows of \code{rowData.GRanges} correspond to rows of the \code{data} matrix. Row names of \code{rowData.GRanges} must be consistent with row names of \code{data}. Their orders are not necessarily the same. Only annotated genes in \code{data} will be used in genomic zone analysis. One of \code{rowData.GRanges} and \code{genome} must be specified.
  }

  \item{ks}{
    an optional numerical vector to specify the number of zones to divide each chromosome into. The names of the \code{ks} vector must be chromosome names. It is only used with user-specified \code{rowData.GRanges}. The seqlevels of \code{rowData.GRanges} must have a corresponding name in \code{ks}. If not specified, an optimal k value (1~400) will be determined for each chromosome. 400 is equivalent to cluster the longest human chromosome into zones averaging wider than 1 million base pairs. Default is \code{NULL}.
  }

 \item{genome}{
    an optional value of \code{character} type to select a genome from \pkg{biomaRt}. Available genomes can be found in the "version" column of the available ensembl datasets from \pkg{biomaRt} database by calling \code{listDatasets(useMart("ensembl"))}. One of \code{rowData.GRanges} and \code{genome} must be specified.
  }

  \item{ensembl.mirror}{
    an optional Ensembl mirror server to connect to. It is used only when \code{genome} is not NULL. The options are \code{"www"}, \code{"uswest"}, \code{"useast"} and \code{"asia"}. Default is \code{"www"}.
  }

  \item{gene.ID.type}{
    an optional value of \code{character} type to specify a gene ID type. Options are \code{"hgnc_symbol"}, \code{"mgi_symbol"}, \code{"ensembl_gene_id"} and \code{"ensembl_transcript_id"}. Only these four types are allowed in this version. This parameter only works with user-specified \code{genome}. If unspecified, all four types would be evaluated to choose the best one.
  }

  \item{ncores}{
    an optional integer to specify the number of cores to use parallely in outstanding genomic zone analysis. Default is 1.
  }
}
\details{
    The function collects all the input information, checks requirement completeness and integrates the inputs into a list, in preparation for function \code{GenomicOZone} to perform outstanding zone analysis.

    A genome annotation parameter of \code{GRanges} class \insertCite{lawrence2013granges}{GenomicOZone} or a genome version must be assigned by the user. The annotation is used to sort genes by their genomic coordinates. The \code{genome} parameter is for function \code{GenomicOZone} to obtain genome annotation from the R package \pkg{biomaRt} \insertCite{smedley2015biomart}{GenomicOZone} to access Ensembl annotation databases \insertCite{zerbino2017ensembl}{GenomicOZone}. Using \code{rowData.GRanges} is recommended over using \code{genome}.
}
\value{
  A list object with all relevant information for oustanding genomic zone analysis. It will be expanded by further analysis.
}
\examples{
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
                       rowData.GRanges = rowData.GRanges,
                       ks = ks)
}

\references{
  \insertAllCited{}
}
