\name{GOZDataSet}
\alias{GOZDataSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a list of all input information
}
\description{
  The function collects all the inputs together and check completeness.
}
\usage{
  GOZDataSet(data, colData, design,
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
    a numerical \code{matrix} of input data. The row variable represents the genes. The column variable represents the samples. The data type can be gene expression, methylation beta value, copy number variation segment mean, or others. The row and column names of the matrix must be specified.
  }

  \item{colData}{
    a \code{dataframe} of sample information. The first column must be the column of sample names, corresponding to columns in the \code{data} matrix. Each of other columns is one variable containing at least two sample conditions for differential analysis.
  }

  \item{design}{
    a one-sided \code{formula} with only right-hand side variables. Only one variable is supported in this version. The formula describes on which variable in \code{colData} the zoning process and differential analysis depend.
  }

  \item{rowData.GRanges}{
    an optional gene annotation in \code{GRanges} format. Rows of \code{rowData.GRanges} are corresponding to row of \code{data}. The names of \code{rowData.GRanges} must be consistent with rownames in \code{data}. Their orders are not necessarily being the same. Subset of \code{data} with only annotated genes will be used in the zoning process. One of the two parameters \code{rowData.GRanges} and \code{genome} must be specified. Using \code{rowData.GRanges} is more recommended than using \code{genome}.
  }

  \item{ks}{
    an optional numerical vector labeled with the chromosome names. It is only used along with user-specified \code{rowData.GRanges}. The seqlevels of the \code{rowData.GRanges} must have a corresponding labeled value in \code{ks}. If not specified, the zoning process with automatically calculate an optimal k value (1~400) for each chromosome. 400 is enough for the algorithm to cluster the longest human chromosome into zones averagely wider than 1M basepair. Default is \code{NULL}.
  }

 \item{genome}{
    an optional \code{character} type selection of the genome from \code{biomaRt}. The available genomes can be found in the "version" column of the available ensembl datasets from \code{biomaRt} database by calling \code{listDatasets(useMart("ensembl"))}.
  }

  \item{ensembl.mirror}{
    an optional Ensembl mirror server to connect to. It is used only when \code{genome} is specified. The options are \code{"www"}, \code{"uswest"}, \code{"useast"} and \code{"asia"}. Default is \code{"www"}.
  }

  \item{gene.ID.type}{
    an optional \code{character} type to specify the kind of the gene names. Options are \code{"hgnc_symbol"}, \code{"mgi_symbol"}, \code{"ensembl_gene_id"} and \code{"ensembl_transcript_id"}. Only the four kinds are allowed in this version. This parameter only works along with used-specified \code{"genome"}. If not specified, all of the four types would be evaluated to choose a correct one, which would cause computational inefficiency.
  }

  \item{ncores}{
    an optional integer to specify the number of cores to use parallely in the zoning process. Default is 1.
  }
}
\details{
    The function is used to collect all the input information, check requirement completeness and integrate the inputs into a list for function \code{GenomicOZone} to perform all process.
    A gene annotation file in \code{GRanges} format \insertCite{lawrence2013granges}{GenomicOZone} or a genome version must be assigned by the user, in order to annotate the genes and sort the genes by coordinates. The \code{genome} information is for function \code{GenomicOZone} to obtain annotations from the R package \code{biomaRt} \insertCite{smedley2015biomart}{GenomicOZone} to access ensembl annotation databases \insertCite{zerbino2017ensembl}{GenomicOZone}. Using \code{rowData.GRanges} is more recommended than using \code{genome}.

    The argument \code{ks} specifies the number of zones to divide the chromosomes into. The function \code{GenomicOZone} will perform the weighted univariate clustering method Ckmeans \insertCite{wang2011ckmeans}{GenomicOZone} with the \code{ks}. If \code{ks} is \code{NULL}, optimal k values will be determined by dynamic programming.
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
                       rowData.GRanges,
                       ks = ks)
}

\references{
  \insertAllCited{}
}
