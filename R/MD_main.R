# Package: GenomicOZone
# MD_main.R
#   - Wraping functions
# Created:
#   Hua Zhong
#   Dec 25, 2018

#######################################
#### Imports for the whole package ####
#' @import GenomicRanges
#' @importFrom Rdpack reprompt
#######################################


#' @importFrom plyr is.formula
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @export
GOZDataSet <- function(data, colData, design,
                       rowData.GRanges = NULL,
                       ks = NULL,
                       genome = NULL,
                       ensembl.mirror = "www",
                       gene.ID.type = NULL,
                       ncores = 1
                       ){

  test <- "ANOVA"
  method <- "1D"

  # Parameter example:
  #
  # data: normalized gene expression matrix
  #              [Sample_1]  [Sample_2]
  #   Gene_1  |  1           2
  #   Gene_2  |  2           4
  #   ...       ...         ...
  #
  # colData:
  #   Sample_name     [Condition]     [tissue]
  #   Sample_1        Cancer          Lung
  #   Sample_2        Normal          Liver
  #   ...             ...             ...
  #
  # design: formula
  #   ~ Condition
  #
  # ks: integer vecter with labels being chromosome names
  #   chr1  chr2  chr3  ... chr10
  #   10    8     4     ... 5
  #
  # rowData.GRanges: gene annotation, rownames must be consistent with rownames in 'data'
  #
  # test:
  #   "ANOVA" ("Chisq", "FunChisq" and "G-test" will be supported)
  #
  # method:
  #   "1D" or "MD"
  #
  # ensembl.mirror:
  #   "uswest", "useast", "www", "asia"
  #
  # gene.ID.type:
  #   "hgnc_symbol", "mgi_symbol", "ensembl_gene_id", "ensembl_transcript_id"
  #
  #

  if(ncol(data) != nrow(colData))
    stop("Sample numbers are not consistent in \'data\' and \'colData\'")

  if(is.null(design) || is.null(colData))
    stop("ERROR: \'design\' and \'colData\' are required!")

  if(any(table(colnames(colData)[-1]) > 1))
    stop("ERROR: colData have duplicate colnames!")

  if(!is.formula(design))
    stop("ERROR: design must be a formula!")

  if(length(all.vars(design)) > 1)
    stop("ERROR: only one variable is supported in \"design\" for now!")

  if(!all.vars(design)[1] %in% colnames(colData)[-1])
    stop("ERROR: \"design\" variables are not consistent with sample names in colData!!")

  if(is.null(rowData.GRanges) && is.null(genome))
    stop("ERROR: one of \"rowData.GRanges\" and \"genome\" must be specified!")

  if(!is.null(rowData.GRanges)){
    if(!any(names(rowData.GRanges) %in% rownames(data)))
      stop("ERROR: \"rowData.GRanges\" and \"data\" have no matched rownames (gene name)!")

    if(!is.null(ks)){
      if(!(all(names(ks) %in% seqlevels(rowData.GRanges)) &&
           all(seqlevels(rowData.GRanges) %in% names(ks))))
        stop("ERROR: the names of ks don't match the seqlevels of rowData.GRanges")

      if(!is.numeric(ks))
        stop("ERROR: ks is not a numerical vecter!")
    }
  }else{
    if(is.null(ks))
      stop("ERROR: ks can only be specified when rowData.GRanges is provided.
            If using ensembl annotations, this parameter is not available yet.")
  }

  if(!ensembl.mirror %in% c("uswest", "useast", "www", "asia"))
    stop("ERROR: The mirror server of ensembl must be one of \"uswest\", \"useast\", \"www\" and \"asia\"!")

  if(is.null(rowData.GRanges) && is.null(gene.ID.type))
    warning("\"gene.ID.type\" is not specified, exhaustively checking from the candidate types will cost more computational time!")

  if(method != "1D" && method != "MD") stop("ERROR: \"method\" must be \"1D\" or \"MD\".")


  if(ncores < 1) stop("ERROR: \"ncores\" must be an integer no less than 1!")


  if(!is.matrix(data)) data <- data.matrix(data)

  rownames(colData) <- colData[,1]

  return(list(input.data = list(data = data, colData = colData,
                                design = design,
                                rowData.GRanges = rowData.GRanges,
                                p.value.test = test,
                                genome = genome,
                                method = method,
                                ks.method = if(is.null(ks)) "optimal" else "ks",
                                ks = ks,
                                mirror = ensembl.mirror,
                                ncores = ncores),
              runtime.var = list(design.treatment.var = tail(all.vars(design), 1),
                                 design.treatment.vals = levels(colData[,tail(all.vars(design), 1)]))
              ))
}


#' @export
GenomicOZone <- function(GOZ.ds){
  GOZ.ds$runtime.var$data.GRanges <- MD.Annotate.data(GOZ.ds)
  GOZ.ds$runtime.var$Data.matrix <- MD.Reduce.dimension(GOZ.ds)
  GOZ.ds$runtime.var$Weight.matrix <- MD.Prepare.weight.matrix(GOZ.ds, log.value = FALSE)

  GOZ.ds$runtime.var$data.GRanges <- MD.Chr.zoning.Granges(GOZ.ds)
  GOZ.ds$runtime.var$zone.GRanges <- MD.Create.zone.GRanges(GOZ.ds)

  Rank.res <- MD.rank.statistic(GOZ.ds)
  GOZ.ds$runtime.var$zone.stat <- Rank.res$Zone.stat
  GOZ.ds$runtime.var$zone.GRanges <- Rank.res$Zone.GRanges
  return(GOZ.ds)
}


#' @export
extract_genes <- function(GOZ.ds){return(GOZ.ds$runtime.var$data.GRanges)}


#' @export
extract_zones <- function(GOZ.ds){return(GOZ.ds$runtime.var$zone.GRanges)}


#' @export
extract_outstanding_zones <- function(GOZ.ds, alpha = 0.05, min.effect.size = 0.8){
  Zone.GRanges <- GOZ.ds$runtime.var$zone.GRanges

  effect.size.threshold <- min.effect.size
        #sort(Zone.GRanges$effect.size, na.last = TRUE, decreasing = TRUE)[floor(length(Zone.GRanges) * effect.size.rate)]
  Zone.GRanges <- Zone.GRanges[Zone.GRanges$p.value.adj <= alpha &
                               Zone.GRanges$effect.size >= effect.size.threshold]
  return(Zone.GRanges)
}


#' @export
extract_zone_expression <- function(GOZ.ds){
  data.mat <- GOZ.ds$input.data$data
  Gene.GRanges <- GOZ.ds$runtime.var$data.GRanges
  zones.all <- names(GOZ.ds$runtime.var$zone.GRanges)
  zone.exp <- do.call("rbind", lapply(zones.all, function(zone){
    return(colSums(data.mat[names(Gene.GRanges[Gene.GRanges$zone == zone]),, drop = FALSE]))
  }))
  rownames(zone.exp) <- zones.all
  colnames(zone.exp) <- colnames(data.mat)

  data.mat.left <- data.mat[! rownames(data.mat) %in% names(Gene.GRanges),, drop = FALSE]
  if(nrow(data.mat.left) != 0) zone.exp <- rbind(zone.exp, data.mat.left)
  return(zone.exp)
}


#' @export
plot_genome <- function(GOZ.ds, plot.file,
                        alpha = 0.05, min.effect.size = 0.8,
                        plot.width = NULL, plot.height = NULL){

  MD.Genome.plots(GOZ.ds, plot.file = plot.file,
                  alpha = alpha,
                  min.effect.size = min.effect.size,
                  plot.width = plot.width,
                  plot.height = plot.height)
}


#' @export
plot_chromosomes <- function(GOZ.ds, plot.file,
                            alpha = 0.05, min.effect.size = 0.8,
                            plot.width = NULL, plot.height = NULL){

  MD.Chromosome.heatmap(GOZ.ds, plot.file = plot.file,
                        alpha = alpha,
                        min.effect.size = min.effect.size,
                        plot.width = plot.width,
                        plot.height = plot.height)
}


#' @export
plot_zones <- function(GOZ.ds, plot.file,
                      alpha = 0.05, min.effect.size = 0.8,
                      log.exp = TRUE, plot.all.zones = FALSE){

  MD.Zone.Gene.path.plots(GOZ.ds, plot.file = plot.file,
                          alpha = alpha,
                          min.effect.size = min.effect.size,
                          log.exp = log.exp,
                          plot.all.zones = plot.all.zones)
}

