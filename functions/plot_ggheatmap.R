plot_ggheatmap <- function(obj, n=nrow(obj), norm=TRUE, log=TRUE, 
                           colnames.in.pdata="NewSampleID", 
                           col.labels=NULL,
                           row.labels=NULL,
                           facet.by=NULL,
                           annotate.cols=NULL,
                           dendrogram=FALSE,
                           col.annotation.offset=1,
                           col.annotation.width=4) {
  ##
  ## Better heatmaps with ggplot and grid.
  ##
  ## Args:
  ## colnames.in.pdata: the name of the column in pData(obj) that matches the 
  ##    colnames of the count matrix. Needs to match up for this to work.
  ## col.labels: what phenotype column to use for the columns (doesn't group,
  ##    just relabels) [default: colnames of count matrix]
  ## row.labels: what taxonomic rank or phenotype column to use for the rows 
  ##    (doesn't group, just relabels) [default: OTUs]
  ## facet.by: formula to pass to facet_grid. Set to NULL to prevent faceting.
  ## annotate.cols: {experimental} add annotations above the columns 
  ## col.annotation.offset: the number of y-axis units above the top to place
  ##    the annotations. Can be negative.
  ## col.annotation.width: the size of the annotation bar
  ## dendrogram: Use hierarchical clustering to produce a dendrogram of the OTUs
  ##    and reorder the heatmap accordingly. Dendrogram appears on right.
  ##
  if (!missing(annotate.cols) & !missing(facet.by)) {
    stop("Currently cannot facet plot and display column annotations at the same time.")
  }
  #-- Data munging
  counts <- MRcounts(obj, norm=norm, log=log)
  pdata <- pData(obj)
  taxa <- make_taxonomy_df(obj)
  # If we display fewer than the total number of OTUs, select the 
  # highest-variance OTUs first
  if (n != nrow(obj)) {
    otu.keepers <- which(rowSums(counts) > 0)
    otu.vars <- rowSds(counts[otu.keepers, ]) 
    otu.indices <- otu.keepers[order(otu.vars, decreasing=TRUE)[1:n]]
    counts <- counts[otu.indices, ]
    taxa <- taxa[otu.indices, ]
  }
  
  mcounts <- counts %>%
    # Melt the counts matrix for use with ggplot
    melt(varnames=c("OTU", colnames.in.pdata), value.name="count") %>%
    # Add in the phenotype data, referencing the column that matches the column 
    # names in the counts matrix
    merge(pdata, by=colnames.in.pdata) %>%
    # Add in the taxonomic ranks for each OTU
    merge(taxa, by.x="OTU", by.y="otu") %>%
    # Convert OTUs to factor
    mutate(OTU=factor(OTU)) %>%
    group_by(NewSampleID) %>%
    mutate(totals=sum(count))
  
  # Reorder the levels of the OTUs so that it will follow the dendrogram
  if (!missing(dendrogram)) {
    ddr <- as.dendrogram(hclust(dist(counts)))
    Rowv <- rowMeans(counts, na.rm = TRUE)
    hcr <- hclust(dist(counts))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    mcounts$OTU <- factor(mcounts$OTU, levels=levels(mcounts$OTU)[rowInd])  
  }
  
  #-- Creating the plot
  p <- ggplot(mcounts, aes_string(x=colnames.in.pdata, y="OTU", fill="count")) +
    geom_tile() +
    theme_classic() + 
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
    #     scale_fill_gradientn(colours=saturated_rainbow(10, saturation_limit = 0.9))
    scale_fill_gradientn(colours=c("white", "lightblue", "purple"))
  # Relabel the Y-axis if desired
  if (!is.null(row.labels)) {
    p <- p + scale_y_discrete(labels = function(x) mcounts[[row.labels]][mcounts$OTU == x]) +
      ylab(row.labels)
  }
  # Relabel the X-axis if desired
  if (!is.null(col.labels)) {
    p <- p + scale_x_discrete(labels = function(x) mcounts[[col.labels]][mcounts[[colnames.in.pdata]] == x]) +
      xlab(col.labels)
  }
  # Facet plot by formula in facet.by (passed as string)
  if (!is.null(facet.by)) {
    p <- p + facet_grid(facet.by, scales="free_x", space="free")
  }
  # Add column annotations (incompatible with faceting)
  if (!missing(annotate.cols)) {
    anno.bar.height <- length(levels(mcounts$OTU)) + col.annotation.offset
    # Have to call as.integer on a factor but aes_string evaluates the strings
    # too late, so this has to be done as a eval'd call
    call <- substitute(
      geom_segment(y=anno.bar.height, 
                   yend=anno.bar.height, 
                   size=col.annotation.width,
                   aes(x=as.integer(.)-0.5,
                       xend=as.integer(.)+0.5,
                       color=annotate.cols)),
      list(.=as.name(colnames.in.pdata), annotate.cols=as.name(annotate.cols)))
    p <- p + eval(call)
    # If we're plotting a dendrogram, we'll already be converting this to a gtable
    # so we'll take care of the clipping then
    if (missing(dendrogram))
      warning("Column annotations will be outside plotting area. Call noclip_plot to draw.")
  }
  # Add a dendrogram for the rows (currently dendrograms for columns are unsupported)
  if (!missing(dendrogram)) {
    require(ggdendro)
    require(grid)
    require(gridExtra)
    d <- ggdendrogram(ddr, rotate=TRUE) + theme(axis.text.x=element_blank(),
                                                axis.text.y=element_blank(),
                                                plot.margin=unit(c(1,0,0,0), 
                                                                 units="lines"))
    p <- p + theme(plot.margin=unit(c(1,0,0,0), units="lines"),
                   legend.position="left")
    gt.d <- ggplot_gtable(ggplot_build(d))
    gt.p <- ggplot_gtable(ggplot_build(p))
    gt.p$layout$clip[gt.p$layout$name=="panel"] <- "off"
    max.height <- unit.pmax(gt.p$heights[1:6], gt.d$heights[1:6])
    gt.d$heights[1:6] <- max.height
    gt.p$heights[1:6] <- max.height
    grid.arrange(gt.p, gt.d, ncol=2, widths=c(4,1))
    return(invisible(NULL))
  }
  return(p)
}

noclip_plot <- function(.ggplot, plot=TRUE) {
  ##
  ## Takes a ggplot2 plot object and turns off clipping, then plots the result.
  ##
  gt <- ggplot_gtable(ggplot_build(.ggplot))
  gt$layout$clip[gt$layout$name=="panel"] <- "off"
  if (plot)
    grid.draw(gt)
  return(gt)
}