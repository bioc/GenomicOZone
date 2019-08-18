# Package: GenomicOZone
# MD_perform_zoning.R
#   - Perform the zoning process
# Created:
#   Hua Zhong
#   Dec 25, 2018


MD.Reduce.dimension <- function(GOZ.ds){
  # Reduce sample dimensions
  X <- GOZ.ds$input.data$data[names(GOZ.ds$runtime.var$data.GRanges),] #GOZ.ds$runtime.var$data.scaled
  group.vars <- all.vars(GOZ.ds$input.data$design)

  X.new <- NULL
  if(GOZ.ds$input.data$method == "1D"){
    X.new <- matrix(rowSums(X, na.rm = TRUE), ncol = 1)
    rownames(X.new) <- rownames(X)
    colnames(X.new) <- "Weight"
  }else if(GOZ.ds$input.data$method == "MD"){
    # MD
    colData <- GOZ.ds$input.data$colData
    group.vars.combination <- expand.grid(lapply(group.vars, function(x){unique(colData[,x])}))
    colnames(group.vars.combination) <- group.vars

    X.new <- lapply(c(1:nrow(group.vars.combination)), function(x){
      x <- as.character(unlist(group.vars.combination[x,]))
      pick <- rep(TRUE, nrow(colData))

      for (i in c(1:length(x))) {
        pick <- pick & as.character(unlist(colData[,colnames(group.vars.combination)[i]])) == x[i]
      }

      pick <- as.character(unlist(colData[1, pick]))
      if(length(pick) == 1){
        return(as.numeric(X[,pick]))
      }else{
        #rowMedians(X[,pick], na.rm = TRUE)
        X.pick <- X[,pick, drop = FALSE]
        X.pick <- as.numeric(apply(X.pick, 1, median, na.rm = TRUE))
        return(X.pick)
      }
    })
    X.new <- do.call("cbind", X.new)
    X.new[is.nan(X.new)] <- NA
    rownames(X.new) <- rownames(X)
    colnames(X.new) <- apply(group.vars.combination, 1, paste, collapse='.')
  }
  return(X.new)
}


MD.Prepare.weight.matrix <- function(GOZ.ds, log.value = FALSE){
  weight.mat <- GOZ.ds$runtime.var$Data.matrix
  weight.mat.new <- apply(weight.mat, 2, function(x){
    x <- as.numeric(x)
    x.na <- is.na(x)
    if(any(x[!x.na] < 0)) x[!x.na] <- log2(2^x[!x.na] + 1)
    if(log.value) x[!x.na] <- log(x[!x.na] + 1)

    x[!x.na] <- x[!x.na] / sum(x[!x.na]) * length(x)
    x[x.na] <- 0
    return(x)
  })
  rownames(weight.mat.new) <- rownames(weight.mat)
  colnames(weight.mat.new) <- colnames(weight.mat)
  return(weight.mat.new)
}


#' @import Ckmeans.1d.dp
#' @importFrom parallel detectCores makeCluster stopCluster parLapply
#' @importFrom GenomeInfoDb seqlevels seqlengths
MD.Chr.zoning.Granges <- function(GOZ.ds){
  method <- GOZ.ds$input.data$ks.method
  ks <- GOZ.ds$input.data$ks
  no_cores <- GOZ.ds$input.data$ncores

  X.GRanges <- GOZ.ds$runtime.var$data.GRanges
  chr.all <- seqlevels(X.GRanges)
  # X.GRanges.new <- vector("list", length(chr.all))
  # names(X.GRanges.new) <- chr.all

  if(no_cores < 1) no_cores <- 1
  cl <- makeCluster(no_cores) #, type = "FORK"
  #for (chr in chr.all) {
  X.GRanges.new <- parLapply(cl, chr.all, function(chr){
    X.GRanges.chr <- X.GRanges[seqnames(X.GRanges) == chr]
    if(length(X.GRanges.chr) == 0) return(GRanges())

    X.GRanges.chr <- sort(X.GRanges.chr, by = ~ start + end)
    x <- GOZ.ds$runtime.var$Weight.matrix[names(X.GRanges.chr),]
    names(x) <- names(X.GRanges.chr)

    k.range <- NULL
    if(method == "optimal"){
      k.range <- c(1:400)
    }else if(method == "ks"){
      k.range <- as.numeric(ks[chr])
    }else{
      stop("ERROR: ks selecting bug!")

      # region.size <- as.numeric(method)
      # k.mean <- round(seqlengths(X.GRanges.chr)[chr] / region.size)
      # k.range <- as.numeric(k.mean) #c((max(c(k.mean-30, 1))):(k.mean+30))
    }

    x.output <- GRanges()
    if(!is.matrix(x)){
      x <- as.numeric(x)
      x.output <- x
      x.output <- Ckmeans.1d.dp(x = c(1:length(x)), y = x, k = k.range)$cluster
    }else{
      stop("ERROR: not yet supported!")

      # x.output <- x
      # x.output <- MDW_test(x = c(1:nrow(x)),
      #                      y = x,
      #                      Kmin = min(k.range), Kmax = max(k.range),
      #                      estimate_k = "BIC", method = "linear",
      #                      normalize_positive_weights = FALSE,
      #                      scale_weights = FALSE)$cluster
    }

    X.GRanges.chr.zone <- paste(chr, "_", x.output, sep="")

    X.GRanges.chr$zone <- factor(X.GRanges.chr.zone, levels = paste(chr, "_", sort(unique(x.output)), sep=""))
    #X.GRanges.new[[chr]] <- X.GRanges.chr
    return(X.GRanges.chr)
  })
  stopCluster(cl)
  names(X.GRanges.new) <- chr.all

  X.GRanges.new <- unlist(GRangesList(X.GRanges.new), use.names = FALSE)
  return(X.GRanges.new)
}


#' @importFrom GenomeInfoDb seqlevels seqlengths seqlevels<- seqlengths<-
#' @importFrom S4Vectors Rle
MD.Create.zone.GRanges <- function(GOZ.ds){
  X.GRanges <- GOZ.ds$runtime.var$data.GRanges
  chr.all <- seqlevels(X.GRanges)
  Zone.GRanges.list <- vector("list", length(chr.all))
  names(Zone.GRanges.list) <- chr.all

  for (chr in chr.all) {
    X.GRanges.chr <- X.GRanges[seqnames(X.GRanges) == chr]

    if(length(X.GRanges.chr) == 0){
      Zone.GRanges.list[[chr]] <- GRanges()
      next
    }

    zone.unique <- levels(X.GRanges.chr$zone)[levels(X.GRanges.chr$zone) %in% unique(X.GRanges.chr$zone)]

    # Get zone coordinates
    zone.start <- c()
    for (zone in zone.unique) {
      X.GRanges.chr.zone <- X.GRanges.chr[X.GRanges.chr$zone == zone]
      if(zone == zone.unique[1]){
        zone.start <- c(zone.start, 1)
      }else{
        zone.start <- c(zone.start, min(start(X.GRanges.chr.zone)))
      }
    }
    names(zone.start) <- zone.unique
    # zone.start <- sort(zone.start)
    # zone.unique <- names(zone.start)

    zone.end <- c(zone.start[-1] - 1, as.numeric(seqlengths(X.GRanges.chr)[chr]))
    names(zone.end) <- zone.unique

    Zone.GRanges.list[[chr]] <- GRanges(seqnames = Rle(chr),
                                        ranges = IRanges(start = zone.start, end = zone.end),
                                        zone = factor(names(zone.start), levels = names(zone.start)))
    names(Zone.GRanges.list[[chr]]) <- Zone.GRanges.list[[chr]]$zone
    ####
  }

  Zone.GRanges <- unlist(GRangesList(Zone.GRanges.list), use.names = FALSE)
  seqlevels(Zone.GRanges) <- seqlevels(X.GRanges)
  seqlengths(Zone.GRanges) <- seqlengths(X.GRanges)

  return(Zone.GRanges)
}


#' @importFrom stats aov p.adjust
#' @importFrom sjstats eta_sq
MD.rank.statistic <- function(GOZ.ds){
  Zone.GRanges <- GOZ.ds$runtime.var$zone.GRanges
  X.GRanges <- GOZ.ds$runtime.var$data.GRanges
  X <- GOZ.ds$input.data$data[names(GOZ.ds$runtime.var$data.GRanges),] #GOZ.ds$runtime.var$data.scaled
  colData <- GOZ.ds$input.data$colData

  Zone.stat <- list()

  p.value.adj <- NULL
  effect.size <- NULL

  Var.combination <- expand.grid(unique(colData[,as.character(all.vars(GOZ.ds$input.data$design)), drop = FALSE]))
  Var.combination.names <- apply(Var.combination, 1, paste, collapse = '.')
  Var.combination.all <- colData[,as.character(all.vars(GOZ.ds$input.data$design)), drop = FALSE]
  Var.combination.all <- apply(Var.combination.all, 1, paste, collapse = '.')
  colnames(Var.combination) <- as.character(all.vars(GOZ.ds$input.data$design))

  # Zone differential
  if(length(Var.combination.names) > 1){
    p.value.test <- GOZ.ds$input.data$p.value.test

    zone.stat <- lapply(names(Zone.GRanges), function(zone){
      test.res <- NULL

      X.GRanges.sub <- X.GRanges[X.GRanges$zone == zone]
      X.sub <- X[names(X.GRanges.sub),, drop = FALSE]

      acc.rank.tmp <- Obtain.gene.rank(X.sub)
      names(acc.rank.tmp) <- rownames(X.sub)
      acc.rank <- Accumulate.gene.rank(acc.rank.tmp)
      acc.rank.norm <- Normalize.statistic(x = acc.rank, Gene.num = nrow(X.sub), Cond.num = ncol(X.sub))


      # acc.rank.mat.tmp <- sapply(GOZ.ds$runtime.var$design.treatment.vals, function(x){
      #   return(sum(acc.rank[colData[,GOZ.ds$runtime.var$design.treatment.var] == x], na.rm = TRUE))
      # })
      #
      # if(p.value.test == "FunChisq"){
      #   acc.rank.mat <- matrix(c(acc.rank.mat.tmp, rep(sum(acc.rank.mat.tmp) / length(acc.rank.mat.tmp), length(acc.rank.mat.tmp))), ncol = 2)
      #   test.res <- fun.chisq.test(x = acc.rank.mat, method = "nfchisq")
      # }else if(p.value.test == "Chisq"){
      #   #acc.rank.mat <- matrix(c(acc.rank.mat.tmp, rep(sum(acc.rank.mat.tmp) / length(acc.rank.mat.tmp), length(acc.rank.mat.tmp))), ncol = 2)
      #   test.res <- chisq.test(x = acc.rank.mat.tmp)
      #   test.res$estimate <- sqrt(test.res$statistic / sum(acc.rank.mat.tmp))
      # }else if(p.value.test == "EFT"){
      #   acc.rank.mat <- matrix(c(acc.rank.mat.tmp, rep(sum(acc.rank.mat.tmp) / length(acc.rank.mat.tmp), length(acc.rank.mat.tmp))), ncol = 2)
      #   acc.rank.mat <- round(acc.rank.mat)
      #   test.res <- fun.chisq.test(x = acc.rank.mat, method = "exact")
      # }


      if(p.value.test == "ANOVA"){
        if(length(Var.combination.names) > 1){
          rank.mat.df <- data.frame(Type = colData[unlist(lapply(acc.rank.tmp, function(x){names(x)}), use.names = FALSE),
                                                   GOZ.ds$runtime.var$design.treatment.var],
                                    Gene = rep(names(acc.rank.tmp), sapply(acc.rank.tmp, function(x){length(x)})),
                                    Value = unlist(acc.rank.tmp, use.names = FALSE))
          res.aov.gene <- aov(Value ~ Type, data = rank.mat.df)
          test.res <- list(statistic = summary(res.aov.gene)[[1]][["F value"]][[1]],
                           estimate = eta_sq(res.aov.gene, partial = TRUE)$partial.etasq,
                           p.value = summary(res.aov.gene)[[1]][["Pr(>F)"]][[1]])
        }else if(length(Var.combination.names) == 1){
          # rank.mat.df <- data.frame(Sample = rep(colnames(X.sub), each = nrow(X.sub)),
          #                           Gene = rep(rownames(X.sub), ncol(X.sub)),
          #                           Value = log(c(X.sub) + 1))
          # res.aov.gene <- aov(Value ~ Gene + Sample, data = rank.mat.df)
          # test.res.eta.sq <- eta_sq(res.aov.gene, partial = TRUE)
          # test.res <- list(statistic = summary(res.aov.gene)[[1]]["Sample", "F value"],
          #                  estimate = test.res.eta.sq[test.res.eta.sq$term == "Sample", "partial.etasq"],
          #                  p.value = summary(res.aov.gene)[[1]]["Sample", "Pr(>F)"])
        }


        # if(F){
        #   acc.rank.mat.df <- data.frame(Type = as.character(colData[,GOZ.ds$runtime.var$design.treatment.var]),
        #                                 Value = as.numeric(acc.rank))
        #   res.aov <- aov(Value ~ Type, data = acc.rank.mat.df)
        #   test.res <- list(statistic = summary(res.aov)[[1]][["F value"]][[1]],
        #                    estimate = eta_sq(res.aov, partial = TRUE)$partial.etasq,
        #                    p.value = summary(res.aov)[[1]][["Pr(>F)"]][[1]])
        # }
      }
      return(list(zone.stat = acc.rank,
                  zone.stat.norm = acc.rank.norm,
                  zone.p.value = c(statistic = as.numeric(test.res$statistic),
                                   p.value = as.numeric(test.res$p.value),
                                   estimate = as.numeric(if(!is.null(test.res$estimate)) test.res$estimate else NA))))
    })

    zone.stat.norm <- lapply(zone.stat, function(x){x$zone.stat.norm})
    zone.p.value <- lapply(zone.stat, function(x){x$zone.p.value})
    zone.stat <- lapply(zone.stat, function(x){x$zone.stat})

    zone.stat <- do.call("rbind", zone.stat)
    rownames(zone.stat) <- names(Zone.GRanges)
    zone.stat.norm <- do.call("rbind", zone.stat.norm)
    rownames(zone.stat.norm) <- names(Zone.GRanges)

    zone.p.value <- do.call("rbind", zone.p.value)
    rownames(zone.p.value) <- names(Zone.GRanges)

    zone.p.value <- cbind(zone.p.value, p.adjust(zone.p.value[,"p.value"], "BH"))
    colnames(zone.p.value)[ncol(zone.p.value)] <- "p.value.adj"

    p.value.adj <- zone.p.value[,"p.value.adj"]
    effect.size <- zone.p.value[,"estimate"]

    Zone.stat$Zone.stat.per.sample <- zone.stat
    Zone.stat$Zone.stat.norm.per.sample <- zone.stat.norm
    Zone.stat$Zone.statistic <- zone.p.value
  }else{
    Zone.stat$Zone.stat.per.sample <- NULL
    Zone.stat$Zone.stat.norm.per.sample <- NULL
    Zone.stat$Zone.statistic <- NULL
  }

  if(!is.null(p.value.adj)){
    elementMetadata(Zone.GRanges) <- cbind(elementMetadata(Zone.GRanges),
                                           data.frame(p.value.adj = p.value.adj))
  }
  if(!is.null(effect.size)){
    elementMetadata(Zone.GRanges) <- cbind(elementMetadata(Zone.GRanges),
                                           data.frame(effect.size = effect.size))
  }

  return(list(Zone.GRanges = Zone.GRanges,
              Zone.stat = Zone.stat))
}


Obtain.gene.rank <- function(L){
  Gene.rank <- lapply(c(1:nrow(L)), function(x){
    Gene.rank.sub <- L[x,]
    Gene.rank.sub <- rank(Gene.rank.sub, na.last = NA, ties.method = "average")
    return(Gene.rank.sub)
  })
  return(Gene.rank)
}


Accumulate.gene.rank <- function(Gene.rank){
  Gene.rank.df <- do.call("rbind", Gene.rank)
  return(colSums(Gene.rank.df, na.rm = TRUE))
}


Normalize.statistic <- function(x, Gene.num, Cond.num){
  e <- Gene.num * (Cond.num + 1) / 2
  return((x - e) / (Gene.num*Cond.num - e))
}
