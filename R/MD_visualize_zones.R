# Package: GenomicOZone
# MD_visualize_zones.R
#   - Generate visualizations
# Created:
#   Hua Zhong
#   Dec 25, 2018


#' @importFrom stats median mad
#' @importFrom plyr ddply
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(Median = median(x[[col]], na.rm=TRUE),
      Mad = mad(x[[col]], na.rm=TRUE),
      Min = min(x[[col]], na.rm=TRUE),
      Max = max(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  return(data_sum)
}


#' @import ggplot2
get_legend<-function(p){
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom gridExtra grid.arrange
#' @importFrom utils setTxtProgressBar tail txtProgressBar
MD.Zone.Gene.path.plots <- function(GOZ.ds, plot.file, p.value.cutoff = 0.05, effect.size.rate = 0.05, log.exp = TRUE,
                                    plot.all.zones = FALSE, cancer.genes = NULL){
  colData <- GOZ.ds$input.data$colData
  Zone.GRanges <- GOZ.ds$runtime.var$zone.GRanges
  Gene.GRanges <- GOZ.ds$runtime.var$data.GRanges
  Gene.GRanges.names <- if(!is.null(Gene.GRanges$external_gene_name)){
    Gene.GRanges$external_gene_name
  }else if(!is.null(Gene.GRanges$Gene.name)){
    Gene.GRanges$Gene.name
  }else{
    names(Gene.GRanges)
  }

  data.mat <- NULL
  p.val.specific <- NULL
  effect.size.specific <- NULL

  if(!is.null(Zone.GRanges$p.value.adj)){
    X.value.var <- GOZ.ds$runtime.var$design.treatment.var
    X.value.unique <- GOZ.ds$runtime.var$design.treatment.vals
    Condition.label <- if(paste(X.value.var, ".label", sep='') %in% colnames(colData)){
      paste(X.value.var, ".label", sep='')
    }else{
      X.value.var
    }
    sample.rep.pick <- lapply(X.value.unique, function(y){
      return(which(colData[,X.value.var] == y))
    })
    names(sample.rep.pick) <- X.value.unique

    data.mat <- GOZ.ds$input.data$data[names(GOZ.ds$runtime.var$data.GRanges),]

    p.val.specific <- Zone.GRanges$p.value.adj
    effect.size.specific <- Zone.GRanges$effect.size
    names(p.val.specific) <- names(Zone.GRanges)
    names(effect.size.specific) <- names(Zone.GRanges)
  }else{
    stop("Gene pattern not ploted, because less than 2 differential conditions specified in colData!")
  }

  effect.size.specific.threshold <- sort(effect.size.specific, na.last = TRUE, decreasing = TRUE)[floor(length(effect.size.specific) * effect.size.rate)]

  zones.all <- names(Zone.GRanges)[order(effect.size.specific, na.last = TRUE, decreasing = TRUE)]
  zone.signif <- !is.na(p.val.specific) & !is.na(effect.size.specific) &
                  p.val.specific <= p.value.cutoff &
                  effect.size.specific >= effect.size.specific.threshold
  zone.signif <- zone.signif[zones.all]
  zone.color <- sapply(zone.signif, function(x){
    if(x){
      return("#F8766D")
    }else{
      return("gray")
    }
  })

  plots.all <- vector("list", length(zones.all))
  names(plots.all) <- zones.all

  if(!plot.all.zones && sum(zone.signif) == 0){
    message("No significant differential expressed zones!")
    return(NULL)
  }

  zones.all.process <- if(plot.all.zones) c(1:length(zones.all)) else c(1:length(zones.all))[zone.signif]
  #pb <- txtProgressBar(min = 1, max = length(zones.all.process), style = 3)
  for (pb.i in zones.all.process) {
    #setTxtProgressBar(pb, pb.i)

    zone <- zones.all[pb.i]

    #print(paste("[", which(zones.all == zone), "/", length(zones.all), "]: ", zone, sep=''))

    gene.pick <- Gene.GRanges$zone == zone
    data.mat.sub <- data.mat[names(Gene.GRanges)[gene.pick] ,, drop = FALSE]
    rownames(data.mat.sub) <- Gene.GRanges.names[gene.pick]

    plots.all.sub <- vector("list", 3)

    if(!is.null(cancer.genes)){
      rownames(data.mat.sub)[rownames(data.mat.sub) %in% cancer.genes] <- paste(rownames(data.mat.sub)[rownames(data.mat.sub) %in% cancer.genes], '*', sep='')
    }

    data.mat.sub.exp <- if(log.exp) log10(data.mat.sub + 1) else data.mat.sub
    data.mat.sub.rank <- do.call("rbind", Obtain.gene.rank(data.mat.sub))
    rownames(data.mat.sub.exp) <- rownames(data.mat.sub)
    colnames(data.mat.sub.exp) <- colnames(data.mat.sub)
    rownames(data.mat.sub.rank) <- rownames(data.mat.sub)
    colnames(data.mat.sub.rank) <- colnames(data.mat.sub)

    zone.gene.exp.df <- data.frame(Y.value = c(data.mat.sub.exp),
                                   X.value = factor(rep(colData[colnames(data.mat.sub.exp), Condition.label],
                                                        each = nrow(data.mat.sub.exp)), levels = X.value.unique),
                                   Gene = rep(rownames(data.mat.sub.exp), ncol(data.mat.sub.exp)))

    zone.gene.exp.df.sum <- data_summary(zone.gene.exp.df, varname="Y.value", groupnames=c("X.value", "Gene"))
    zone.gene.exp.df.sum.order <- zone.gene.exp.df.sum[zone.gene.exp.df.sum$X.value == tail(levels(zone.gene.exp.df.sum$X.value), 1),]
    zone.gene.exp.df.sum.order <- zone.gene.exp.df.sum.order[order(zone.gene.exp.df.sum.order$Median, decreasing = TRUE),]
    zone.gene.exp.df.sum <- zone.gene.exp.df.sum[order(match(zone.gene.exp.df.sum$Gene, unique(zone.gene.exp.df.sum.order$Gene))),]
    zone.gene.exp.df.sum$Gene <- factor(as.character(zone.gene.exp.df.sum$Gene), levels = as.character(unique(zone.gene.exp.df.sum$Gene)))
    gene.labels <- as.character(levels(zone.gene.exp.df.sum$Gene))

    zone.gene.rank.df <- data.frame(Y.value = c(data.mat.sub.rank),
                                    X.value = factor(rep(colData[colnames(data.mat.sub.rank), Condition.label],
                                                         each = nrow(data.mat.sub.rank)), levels = X.value.unique),
                                    Gene = rep(rownames(data.mat.sub.rank), ncol(data.mat.sub.rank)))

    legend.angle <- if(length(X.value.unique) >= 3) 30 else 0
    legend.angle.adjust <- if(length(X.value.unique) >= 3) 1 else 0.5
    col.num <- if(length(unique(zone.gene.exp.df.sum$Gene)) > 30) 2 else 1

    p <- ggplot(data = zone.gene.exp.df.sum, aes_string(x = "X.value", y = "Median", color = "Gene", group = "Gene")) +
      geom_hline(yintercept=0) +
      geom_errorbar(aes_string(ymin = "Min", ymax = "Max"), width=.1) +
      geom_line() +
      # expand_limits(x = c(1, length(levels(zone.gene.exp.df$X.value)))) +
      ggtitle(paste(zone, ": ", start(Zone.GRanges[zone]), " - ", end(Zone.GRanges[zone]), '\n',
                    "Adjusted p-value: ", signif(p.val.specific[zone], digits = 3), '\n',
                    "Effect size: ", signif(effect.size.specific[zone], digits = 3), sep="")) +
      labs(x = GOZ.ds$runtime.var$design.treatment.var, y = "Log10(exp+1)") +
      theme(plot.title = element_text(hjust = 0.5, color = zone.color[pb.i]),
            #axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            #panel.grid.major = element_blank(),
            #panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            legend.position = "right", #c(-Inf,Inf),
            legend.key = element_rect(size = 0.5),
            legend.key.size = unit(0.5, 'lines'),
            plot.margin = unit(c(1,0,0,1), "lines"),
            legend.margin=margin(t = 0, unit='cm')
      ) +
      guides(color=guide_legend(ncol=col.num, byrow=TRUE))

    p.legend <- get_legend(p)
    p <- p + theme(legend.position="none")

    p.box <- ggplot(data = zone.gene.rank.df, aes_string(x = "X.value", y = "Y.value", group = "X.value")) +
      geom_hline(yintercept=0) +
      geom_violin(fill = "#00BFC4", trim=TRUE) +
      geom_boxplot(width = 0.15, fill="white") + #color = "#00BFC4",
      stat_summary(data = zone.gene.rank.df, aes_string(x = "X.value", y = "Y.value", group = "1"),
                   fun.y = "median", geom = "line", color = "#F8766D") +
      # expand_limits(x = c(1, length(levels(zone.gene.rank.df$X.value)))) +
      # ggtitle(paste(zone, ": ", start(Zone.GRanges[zone]), " - ", end(Zone.GRanges[zone]), " (", length(gene.labels), " genes)", '\n',
      #               "Adjusted p-value: ", signif(p.val.specific[zone], digits = 3), sep="")) +
      labs(x = GOZ.ds$runtime.var$design.treatment.var, y = "Gene rank") +
      theme(plot.title = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x = element_text(size = 13, angle = legend.angle, hjust = legend.angle.adjust),
            #panel.grid.major = element_blank(),
            #panel.grid.minor = element_blank(),
            legend.title=element_blank(),
            legend.position="none",
            plot.margin = unit(c(1,0,1,1), "lines")
      )

    plots.all.sub[[1]] <- p
    plots.all.sub[[2]] <- p.box
    plots.all.sub[[3]] <- p.legend

    plots.all[[zone]] <- plots.all.sub
  }

  if(plot.all.zones){
    plot.file.tmp <- if(substr(plot.file, nchar(plot.file)-3, nchar(plot.file)) == ".pdf"){
      paste(substr(plot.file, 1, nchar(plot.file)-4), "_all_zones.pdf", sep='')
    }else{
      paste(plot.file, "_all_zones.pdf", sep='')
    }
    pdf(plot.file.tmp, width = length(X.value.unique)/2+3, height = 4.5)
    #pb <- txtProgressBar(min = 0, max = length(plots.all), style = 3)
    for (i in c(1:length(plots.all))) {
      #setTxtProgressBar(pb, i)
      if(length(X.value.unique) >= 3){
        grid.arrange(grobs = plots.all[[i]],
                     layout_matrix = rbind(c(1,3), c(2,3)),
                     widths = c(3.5, 1), heights = c(4.5,3.5))
      }else{
        grid.arrange(grobs = plots.all[[i]],
                     layout_matrix = rbind(c(1,3), c(2,3)),
                     widths = c(3.5, 1), heights = c(5,3))
      }
    }
    dev.off()
  }

  plots.all.signif <- plots.all[names(zone.signif[zone.signif])]
  if(length(plots.all.signif) > 0){
    pdf(plot.file, width = length(X.value.unique)/2+3, height = 4.5)
    #pb <- txtProgressBar(min = 0, max = length(plots.all.signif), style = 3)
    for (i in c(1:length(plots.all.signif))) {
      #setTxtProgressBar(pb, i)
      if(length(X.value.unique) >= 3){
        grid.arrange(grobs = plots.all.signif[[i]],
                     layout_matrix = rbind(c(1,3), c(2,3)),
                     widths = c(3.5, 1), heights = c(4.5,3.5))
      }else{
        grid.arrange(grobs = plots.all.signif[[i]],
                     layout_matrix = rbind(c(1,3), c(2,3)),
                     widths = c(3.5, 1), heights = c(5,3))
      }

    }
    dev.off()
  }else{
    message("No significant differential expressed zones!")
  }
}


#' @import ggplot2
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom grDevices dev.off pdf
#' @importFrom utils setTxtProgressBar txtProgressBar
MD.Chromosome.heatmap <- function(GOZ.ds, plot.file, p.value.cutoff = 0.05, effect.size.rate = 0.05,
                                  plot.width = NULL, plot.height = NULL){
  colData <- GOZ.ds$input.data$colData
  rownames(colData) <- colData[,1]
  Genes.GRanges <- GOZ.ds$runtime.var$data.GRanges
  Zone.GRanges <- GOZ.ds$runtime.var$zone.GRanges
  Exp <- GOZ.ds$input.data$data
  effect.size.threshold <- sort(Zone.GRanges$effect.size, na.last = TRUE, decreasing = TRUE)[floor(length(Zone.GRanges) * effect.size.rate)]

  Gene.npt.all.zero <- rownames(Exp)[apply(Exp, 1, function(x){!all(x == 0)})]

  plots.all <- vector("list", length(seqlevels(Zone.GRanges)))
  names(plots.all) <- seqlevels(Zone.GRanges)

  for (chr in seqlevels(Zone.GRanges)) {
    Zone.GRanges.chr <- Zone.GRanges[seqnames(Zone.GRanges) == chr]
    Genes.GRanges.chr <- Genes.GRanges[seqnames(Genes.GRanges) == chr]
    Zone.GRanges.chr.signif <- Zone.GRanges.chr[Zone.GRanges.chr$p.value.adj <= p.value.cutoff &
                                                Zone.GRanges.chr$effect.size >= effect.size.threshold]

    Gene.pick <- names(Genes.GRanges.chr)[names(Genes.GRanges.chr) %in% Gene.npt.all.zero]
    Genes.GRanges.chr <- Genes.GRanges.chr[Gene.pick]

    Rank.mat.chr <- GOZ.ds$runtime.var$zone.stat$Zone.stat.norm.per.sample[names(Zone.GRanges.chr),, drop = FALSE]

    chr.plot.df <- data.frame(Zone = factor(rep(names(Zone.GRanges.chr), ncol(Rank.mat.chr)), levels = names(Zone.GRanges.chr)),
                              #Gene = rep(names(Genes.GRanges.chr), ncol(Rank.mat.chr)),
                              Exp = c(Rank.mat.chr),
                              Sample = rep(colnames(Rank.mat.chr), each = nrow(Rank.mat.chr)),
                              Conition = factor(rep(colData[colnames(Rank.mat.chr), "Condition"], each = nrow(Rank.mat.chr)), levels = levels(colData$Condition)))

    p <- ggplot() +
          geom_tile(data = chr.plot.df, aes_string(x = "Zone", y = "Sample", fill = "Exp"), color = "white") +
          facet_grid(Conition ~ ., scales = "free_y", space = "free_y", switch="y") +
          scale_fill_gradient2(low = "blue", mid = "gray94", high = "red", limits=c(-1, 1), breaks=seq(-1,1,by=0.5)) +
          labs(title = paste(chr, " significant differentially expressed zones"))

    if(length(Zone.GRanges.chr.signif) > 0){
      Signif.rec.df <- data.frame(x.min = as.numeric(sapply(names(Zone.GRanges.chr.signif), function(x){which(names(Zone.GRanges.chr) == x)})) - 0.5,
                                  x.max = as.numeric(sapply(names(Zone.GRanges.chr.signif), function(x){which(names(Zone.GRanges.chr) == x)})) + 0.5,
                                  y.min = -Inf,
                                  y.max = Inf)
        p <- p +
          geom_rect(data = Signif.rec.df, aes_string(xmin="x.min", xmax="x.max", ymin="y.min", ymax="y.max"),
                    color = "red", fill = NA, size = 1, linetype = 5)
    }

    p <- p +
      theme(#text = element_text(size=20),
            plot.title = element_text(size = 20, hjust = 0.5),
            axis.title.x=element_blank(),
            axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5), #size = 10,
            axis.text.y = element_blank(),
            # panel.grid.major.y = element_line(colour = "darkgray", size = 0.1, linetype = "dashed"),
            # panel.grid.minor.y = element_line(colour = "darkgray", size = 0.1, linetype = "dashed"),
            # panel.grid.major.x = element_blank(),
            # panel.grid.minor.x = element_blank(),
            # axis.text.x=element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_blank(),
            #legend.position = "top",
            legend.position = "right", #c(0,1),
            # legend.justification = c(0, 0),
            # legend.direction = "horizontal",
            # legend.margin = margin(0,0,0,0),
            # legend.box.margin = margin(0,0,0,0),
            # panel.background = element_rect(fill = "white", colour = "black")
            strip.text.y = element_text(angle = 180)
      )

    plots.all[[chr]] <- p
  }

  pdf(plot.file, width = if(!is.null(plot.width)) plot.width else 15,
                 height = if(!is.null(plot.width)) plot.height else 6)
  for(chr in names(plots.all)){
    print(plots.all[[chr]])
  }
  dev.off()
}


#' @import ggplot2
#' @importFrom ggbio autoplot
#' @importFrom grDevices dev.off pdf
#' @importFrom utils setTxtProgressBar txtProgressBar
MD.Genome.plots <- function(GOZ.ds, plot.file, p.value.cutoff = 0.05, effect.size.rate = 0.05, plot.width = NULL, plot.height = NULL){
  Zone.GRanges <- GOZ.ds$runtime.var$zone.GRanges
  Zone.GRanges.sub <- NULL

  if(all(is.na(Zone.GRanges$effect.size))){
    Zone.GRanges.sub <- Zone.GRanges[!is.na(Zone.GRanges$p.value.adj) &
                                       Zone.GRanges$p.value.adj <= p.value.cutoff]
  }else{
    effect.size.threshold <- sort(Zone.GRanges$effect.size, na.last = TRUE, decreasing = TRUE)[floor(length(Zone.GRanges) * effect.size.rate)]
    Zone.GRanges.sub <- Zone.GRanges[!is.na(Zone.GRanges$p.value.adj) &
                                      Zone.GRanges$p.value.adj <= p.value.cutoff &
                                      Zone.GRanges$effect.size >= effect.size.threshold]
  }

  if(length(Zone.GRanges.sub) > 0){
    p <- autoplot(object = Zone.GRanges.sub, layout = "karyogram", fill = "#F8766D", color = "black") +
      labs(title=paste("Genome overview"), sep="") +
      theme(plot.title = element_text(size = 20, hjust = 0.5),
            legend.title = element_blank(),
            #legend.position = "right",
            text = element_text(size=20),
            legend.position = c(0,1),
            legend.justification = c(0, 0),
            legend.direction = "horizontal",
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(0,0,0,0))


    pdf(plot.file,
        width = if(!is.null(plot.width)) plot.width else 15,
        height = if(!is.null(plot.width)) plot.height else 10)
    print(p)
    dev.off()
  }else{
    message("No significant differential expressed zones!")
  }
  return(NULL)
}
