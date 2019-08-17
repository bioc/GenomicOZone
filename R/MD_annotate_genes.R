# Package: GenomicOZone
# MD_annotated_genes.R
#   - Annotate the genes with given genome infomation
# Created:
#   Hua Zhong
#   Dec 25, 2018


#' @import biomaRt
find.ensembl.dataset <- function(Dataset = NULL, Genome = NULL, mirror = "www"){
  Genomes.available <- listDatasets(useEnsembl("ensembl", mirror = "uswest"))

  if(is.null(Genome) && is.null(Dataset)){
    stop("ERROR: at least one parameter is specified!")
  }

  Ensembl.dataset <- ""

  if(!is.null(Dataset) && is.null(Genome)){
    Dataset.index <- grepl(pattern = paste(Dataset, "[_gene_ensembl]*", sep=""),
                           x = Genomes.available$dataset, ignore.case = TRUE)
    if(any(Dataset.index)){
      Ensembl.dataset <- Genomes.available$dataset[Dataset.index][1]
    }else{
      stop("Please find a correct enssembl dataset/genome from biomaRt.
           Use command: listDatasets( useMart(\"ensembl\"))")
    }
  }else if(is.null(Dataset) && !is.null(Genome)){
    if(tolower(Genome) == tolower("GRCh38") || tolower(Genome) == tolower("hg38")){
      # Human hg38
      Ensembl.dataset <- "hsapiens_gene_ensembl"
    }else if(tolower(Genome) == tolower("GRCm38") || tolower(Genome) == tolower("mm10")){
      # Mouse mm10
      Ensembl.dataset <- "mmusculus_gene_ensembl"
    }else{
      Genome.index <- grepl(pattern = paste(Genome, "[\\..*]*", sep=""),
                            x = Genomes.available$version,
                            ignore.case = TRUE, fixed = FALSE)
      if(any(Genome.index)){
        Ensembl.dataset <- Genomes.available$dataset[Genome.index][1]
      }else{
        stop("Please find a correct enssembl dataset/genome from biomaRt.
             Use command: listDatasets( useMart(\"ensembl\"))")
      }
    }
  }else{
    Dataset.index <- grepl(pattern = Dataset, x = Genomes.available$dataset, ignore.case = TRUE)
    Genome.index <- grepl(pattern = paste(Genome, "\\..*", sep=""), x = Genomes.available$version,
                          ignore.case = TRUE, fixed = FALSE)
    Index <- Dataset.index & Genome.index
    if(any(Index)){
      Ensembl.dataset <- Genomes.available$dataset[Index][1]
    }else{
      stop("Please find a correct enssembl dataset/genome from biomaRt.
           Use command: listDatasets( useMart(\"ensembl\"))")
    }
  }

  ensembl.db <- useEnsembl("ensembl", mirror = mirror,
                           dataset = Ensembl.dataset,
                           host = "www.ensembl.org")
  return(ensembl.db)
}


#' @import biomaRt
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
Search.genes.from.ensembl <- function(gene.names, Dataset = NULL, Genome = NULL, ID.type = NULL, mirror = "www"){
  #message("Start matching gene names from ensembl.")
  ensembl.db <- find.ensembl.dataset(Genome = Genome, mirror = mirror)

  gene.names.type.colname <- NULL
  attri.list <- NULL
  Match.1 <- NULL

  if(is.null(ID.type)) ID.type <- c("Ensembl_gene_ID", "Ensembl_transcript_ID", "Gene_name")

  for (gene.names.type in ID.type) {
    if(tolower(gene.names.type) == "gene_name"){
      gene.names.type.colname.tmp <- "hgnc_symbol"
      attri.list.tmp <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                          "start_position", "end_position", "strand")
    }else if(tolower(gene.names.type) == "ensembl_gene_id"){
      gene.names.type.colname.tmp <- "ensembl_gene_id"
      attri.list.tmp <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                          "start_position", "end_position", "strand")
    }else if(tolower(gene.names.type) == "ensembl_transcript_id"){
      gene.names.type.colname.tmp <- "ensembl_transcript_id"
      attri.list.tmp <- c("ensembl_transcript_id", "hgnc_symbol", "chromosome_name",
                          "transcript_start", "transcript_end", "strand")
    }else if(tolower(gene.names.type) == "mgi_symbol"){
      gene.names.type.colname.tmp <- "mgi_symbol"
      attri.list.tmp <- c("ensembl_gene_id", "mgi_symbol", "chromosome_name",
                          "start_position", "end_position", "strand")
    }
    Match.tmp <- getBM(attributes = attri.list.tmp,
                       values = as.character(gene.names), uniqueRows =  TRUE,
                       filters = gene.names.type.colname.tmp, mart = ensembl.db)
    if(nrow(Match.tmp) >= length(gene.names)*0.8) {
      Match.1 <- Match.tmp
      gene.names.type.colname <- gene.names.type.colname.tmp
      attri.list <- attri.list.tmp
      break
    }
  }

  if(nrow(Match.1) == 0) stop("The gene names cannot be annotated!")
  if(nrow(Match.1) < length(gene.names)*0.8) stop("Only less than 80% of the genes can be annotated! STOP!")

  colnames(Match.1[,c("chromosome_name", "start_position",  "end_position", "strand")]) <- c("chromosome_name", "start", "end", "strand")

  Match.1.genes <- do.call("rbind", lapply(gene.names, function(x){
    Match.tmp <- Match.1[,gene.names.type.colname] == x
    if(sum(Match.tmp) == 0){
      return(NULL)
    }else if(sum(Match.tmp) > 1){
      Match.tmp.pick <- which(Match.tmp)[!grepl("CHR\\_.*\\_PATCH",
                                                Match.1$chromosome_name[Match.tmp],
                                                ignore.case = TRUE, fixed = FALSE)]
      if(length(Match.tmp.pick) == 0){
        return(NULL)
      }else{
        Match.tmp[-Match.tmp.pick[1]] <- FALSE
      }
    }

    return(Match.1[Match.tmp,])
    #return(unique(c(unlist(Match.1[Match.1$hgnc_symbol == x,-1], use.names = F))))
  }))

  Match.genes <- Match.1.genes
  rownames(Match.genes) <- as.character(Match.genes[,gene.names.type.colname])

  Match.genes.GRanges <- GRanges(seqnames = Rle(paste("chr", Match.genes$chromosome_name, sep="")),
                                 ranges = IRanges(start = Match.genes$start, end = Match.genes$end),
                                 strand = Rle(Match.genes$strand))
  elementMetadata(Match.genes.GRanges) <- Match.genes[,c(1:2)]
  names(Match.genes.GRanges) <- Match.genes[,gene.names.type.colname]

  #message("Finish matching gene names from ensembl.")

  Match.genes.GRanges <- Add.chr.info(Match.genes.GRanges, Genome)
  return(Match.genes.GRanges)
}


#' @importFrom GenomeInfoDb seqlevels seqlengths seqlevels<- seqlengths<-
Add.chr.info <- function(X.GRanges, Genome){
  chr.size <- NULL

  if(!is.null(Genome)){
    if(tolower(Genome) == tolower("GRCh38") ||
       tolower(Genome) == tolower("hg38")){
      # Human hg38
      chr.all <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22", "chrX", "chrY")
      chr.size <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                    170805979, 159345973, 145138636, 138394717, 133797422,
                    135086622, 133275309, 114364328, 107043718, 101991189,
                    90338345, 83257441, 80373285, 58617616, 64444167,
                    46709983, 50818468, 156040895, 57227415)
      names(chr.size) <- chr.all
    }else if(tolower(Genome) == tolower("GRCh38") ||
             tolower(Genome) == tolower("hg19")){
      # Human hg19
      chr.all <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22", "chrX", "chrY")
      chr.size <- c(249250621, 243199373, 198022430, 191154276, 180915260,
                    171115067, 159138663, 146364022, 141213431, 135534747,
                    135006516, 133851895, 115169878, 107349540, 102531392,
                    90354753, 81195210, 78077248, 59128983, 63025520,
                    48129895, 51304566, 155270560, 59373566)
      names(chr.size) <- chr.all
    }else if(tolower(Genome) == tolower("GRCm38") ||
             tolower(Genome) == tolower("mm10")){
      # mm10
      chr.all <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrM")
      chr.size <- c(195471971, 182113224, 160039680, 156508116, 151834684,
                    149736546, 145441459, 129401213, 124595110, 130694993,
                    122082543, 120129022, 120421639, 124902244, 104043685,
                    98207768, 94987271, 90702639, 61431566,
                    171031299, 91744698, 16299)
      names(chr.size) <- chr.all
    }else if(tolower(Genome) == tolower("mm9")){
      # mm9
      chr.all <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrM")
      chr.size <- c(197195432, 181748087, 159599783, 155630120, 152537259,
                    149517037, 152524553, 131738871, 124076172, 129993255,
                    121843856, 121257530, 120284312, 125194864, 103494974,
                    98319150, 95272651, 90772031, 61342430,
                    166650296, 15902555, 16299)
      names(chr.size) <- chr.all
    }else{
      warning("Not supported genome for plotting!")
    }
  }else{
    chr.all <- as.character(unique(seqnames(X.GRanges)))
    chr.size <- sapply(chr.all, function(chr){
      X.GRanges.chr <- X.GRanges[seqnames(X.GRanges) == chr]
      return(max(end(X.GRanges.chr)))
    })
    names(chr.size) <- chr.all
    warning("Chromosome length are estimated by the end of the last gene!")
  }

  X.GRanges <- X.GRanges[as.character(seqnames(X.GRanges)) %in% names(chr.size)]
  seqlevels(X.GRanges) <- names(chr.size)
  seqlengths(X.GRanges) <- chr.size
  return(X.GRanges)
}


#' @importFrom GenomeInfoDb seqlevels seqlengths
MD.Annotate.data <- function(GOZ.ds, ID.type = NULL){
  mirror <- GOZ.ds$input.data$mirror
  data.GRanges <- NULL

  if(is.null(GOZ.ds$input.data$rowData.GRanges)){
    data.GRanges <- Search.genes.from.ensembl(gene.names = rownames(GOZ.ds$input.data$data),
                                              Genome = GOZ.ds$input.data$genome,
                                              ID.type = ID.type,
                                              mirror = mirror)
  }else{
    data.GRanges <- GOZ.ds$input.data$rowData.GRanges
    data.GRanges <- data.GRanges[names(data.GRanges) %in% rownames(GOZ.ds$input.data$data)]
  }
  if(is.null(seqlengths(data.GRanges)) || any(is.na(seqlengths(data.GRanges)))){
    data.GRanges <- Add.chr.info(X.GRanges = data.GRanges,
                                 Genome = GOZ.ds$input.data$genome)
  }

  data.GRanges <- sort(data.GRanges, by = ~ seqnames + start + end)
  return(data.GRanges)
}
