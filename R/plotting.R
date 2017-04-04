#' Plotting functions for \pkg{\link{breakpointR}}
#' 
#' This function will produce several plots from a \code{\link{BreakPoint}} object.
#'
#' @param breakpoints A list of \code{\link{BreakPoint}} objects or a vector with files that contain such objects.
#' @param file Name of the file to plot to.
#' @return A list with \code{\link[ggplot2:ggplot]{ggplot}} objects.
#'
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "breakpointer", package="strandseqExampleData")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Run breakpointR
#'brkpts <- runBreakpointr(exampleFile, pairedEndReads=FALSE,
#'                          chromosomes=paste0('chr', c(1:22,'X')))
#'## Plot them
#'plotBreakpoints(brkpts)

plotBreakpoints <- function(files2plot, file=NULL) {
  
    if (class(files2plot) == class.breakpoint) {
        numplots <- 1
        files2plot <- list(files2plot)
    } else {
        numplots <- length(files2plot)
    } 
  
    plots <- list()
    for (i in 1:numplots) {
    
        filename <- basename(files2plot)[i]
        if (!is.null(filename)) {
          ptm <- startTimedMessage("Working on plot ", filename, " ...")
        } else {
          ptm <- startTimedMessage("Working on plot ", i, " ...")
        }
        data <- loadFromFiles(files2plot[[i]], check.class=class.breakpoint)[[1]]
        
        reads <- data$fragments
        chroms2plot <- seqlevels(reads)
        breaks <- data$breaks
        counts <- data$counts
        lib.metrics <- data$lib.metrics
        lib.metrics <- round(lib.metrics, digits = 5)
        lib.metrics <- paste(names(lib.metrics), lib.metrics, sep = '=')
        lib.metrics <- paste(lib.metrics, collapse = "  |  ")
        
        binned.data <- unlist(tileGenome(seqlengths(reads), tilewidth = 200000))
        
        #counts overlaps between bins and our reads
        Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='-']) 
        Crickreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='+'])
        bothreads <- Watsonreads + Crickreads
        
        mcols(binned.data)$bothreads <- bothreads
        mcols(binned.data)$Watsonreads <- Watsonreads
        mcols(binned.data)$Crickreads <- Crickreads
      
        #transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
        cum.seqlengths <- cumsum(as.numeric(seqlengths(binned.data)))
        cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
        names(cum.seqlengths.0) <- seqlevels(binned.data)
      
        #get positions of ends of each chromosome to plot lones between the chromosomes
        chr.lines <- data.frame( y=cum.seqlengths[-length(cum.seqlengths)] )
        #get positions of each chromosomes names
        chr.label.pos <- round( cum.seqlengths.0 + (0.5 * seqlengths(binned.data) ) )
      
        transCoord <- function(gr) {
            gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
            gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
            return(gr)
        }
      
        #binned.data <- unname(unlist(binned.data))
      
        trans.reads <- transCoord(binned.data)
        trans.breaks <- transCoord(breaks)
        trans.counts <- transCoord(counts)
      
        dfplot.reads <- as.data.frame(trans.reads)
        dfplot.breaks <- as.data.frame(trans.breaks)
        dfplot.counts <- as.data.frame(trans.counts)
      
        my_theme <- theme(
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()
        )
      
      
        ### PLOT READS
      
        #get midpoint values for each breakpoint
        dfplot.reads$midpoint <- dfplot.reads$start.genome + ( (dfplot.reads$start.genome - dfplot.reads$start.genome)/2 ) #get the midposition of each bin
      
        #filter bins with extremely high amount of reads
        Crickreads.outlier <- round(max(stats::runmed(dfplot.reads$Crickreads, 3)))
        Watsonreads.outlier <- round(max(stats::runmed(dfplot.reads$Watsonreads, 3)))
        limit <- min(Crickreads.outlier,Watsonreads.outlier) #take lower threshold in order to limit plotting of y axis
      
        #change red and blue color as you want
        dfplot.reads$mCrickreads <- -dfplot.reads$Crickreads
        ggplt1 <- ggplot(dfplot.reads) + geom_linerange(aes_string(ymin=0, ymax='mCrickreads', x='midpoint'), color="paleturquoise4", size=0.2)
        ggplt1 <- ggplt1 + geom_linerange(aes_string(ymin=0, ymax='Watsonreads', x='midpoint'), color="sandybrown", size=0.2) + scale_y_continuous(limits = c(-limit,limit))
        ggplt1 <- ggplt1 + geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black') + xlab(NULL) + ylab("Read counts") + scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + my_theme + ggtitle(filename, subtitle = lib.metrics)
        
        ### PLOT COUNTS
      
        #if (length(breaks)) {
            ## Scaling size of the rectangle to amout of reads in a given region
            scale <- (dfplot.counts[,c('Ws','Cs')] / dfplot.counts$width) * 1000000
          
            ## filter regions small regions => hard to see on the plot
            outlier.W <- round(max(stats::runmed(scale$Ws, 3)))
            outlier.C <- round(max(stats::runmed(scale$Cs, 3)))
        
            set.max.W <- round(max(scale$Ws[scale$Ws < outlier.W]))
            set.max.C <- round(max(scale$Cs[scale$Cs < outlier.C]))  
        
            scale$Ws[scale$Ws > outlier.W] <- set.max.W
            scale$Cs[scale$Cs > outlier.C] <- set.max.C
            
            ## putting scaled and filtered regions into a data frame
            names(scale) <- c('W.scaled', 'C.scaled')
            df.W <- cbind(dfplot.counts, scaled=scale$W.scaled, fill.strand=rep('W', length(scale$W.scaled)) )
            df.C <- cbind(dfplot.counts, scaled=-scale$C.scaled, fill.strand=rep('C', length(scale$W.scaled)) )
            dfplot.counts <- rbind(df.W, df.C)
        
            dfplot.counts <- dfplot.counts[dfplot.counts$width > 20000,]
        
            #get midpoint values for each breakpoint
            dfplot.breaks$midpoint <- dfplot.breaks$start.genome + ( (dfplot.breaks$end.genome - dfplot.breaks$start.genome) %/% 2)
        
            ggplt2 <- ggplot(dfplot.counts) + geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax='scaled', fill='fill.strand')) + scale_fill_manual(values=c("sandybrown","paleturquoise4"))
            ggplt2 <- ggplt2 + geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black')
            ggplt2 <- ggplt2 + geom_point(data=dfplot.breaks, aes_string(x='midpoint', y=0), size=5, color='red', shape=124, inherit.aes = FALSE) + xlab(NULL) + ylab("Breaks") + scale_x_continuous(expand = c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + my_theme
        
        
            ### PLOT STATES
        
            ggplt3 <- ggplot(dfplot.counts) + geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax=10, fill='states'))
            ggplt3 <- ggplt3 + geom_linerange(data=chr.lines, aes_string(ymin=0, ymax=10, x='y'), col='black') + xlab("Chromosomes") + ylab("States") +   scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) + scale_fill_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + my_theme
        #}
      
        ### PLOT ALL
      
        if (length(breaks)) {
            p <- cowplot::plot_grid(ggplt1, ggplt2, ggplt3, ncol=1, align="v", rel_heights = c(3,3,2))
        } else {
            p <- cowplot::plot_grid(ggplt1, ncol=1, align="v", rel_heights = 3)
        }
        plots[[length(plots)+1]] <- p
        stopTimedMessage(ptm)
    }
  
    ### PRINT TO FILE
  
    ## printing to a file or returning plot object
    if (!is.null(file)) {
        message("Printing to PDF ",file)
    
        grDevices::pdf(file, width=length(chroms2plot), height=5)
        bquiet = lapply(plots, print)
        d <- grDevices::dev.off()
    } else {
        return(plots)
    }
}



#' Genome wide heatmap of template inheritance states
#'
#' Plot a genome wide heatmap of template inheritance states.
#'
#' @param breakpoints A list of \code{\link{BreakPoint}} objects or a vector with files that contain such objects.
#' @param file Name of the file to plot to.
#' @param hotspots
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object of \code{NULL}, depending on option \code{file}.
#' 
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples 
#'## Get example BreakPoint objects to plot
#'data(example_BreakPoints)
#'## Plot the heatmap
#'plotHeatmap(example_BreakPoints)
#'
plotHeatmap <- function(files2plot, file=NULL, hotspots=NULL) {

    if (class(files2plot) == class.breakpoint) {
        stop("Cannot make heatmap from only one object.")
    } else {
        numlibs2plot <- length(files2plot)
    } 
  
    message("Preparing heatmap from ",numlibs2plot, " libraries")

    IDs <- list()
    grl <- GRangesList()
    breaks <- GRangesList()
    for (i in 1:numlibs2plot) {  
        data <- loadFromFiles(files2plot[[i]], check.class=class.breakpoint)
        IDs[[i]] <- data[[1]]$ID
        grl[[i]] <- data[[1]]$counts
        suppressWarnings( breaks[[i]] <- data[[1]]$breaks )
    }

    #transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
    cum.seqlengths <- cumsum(as.numeric(seqlengths(grl[[1]])))
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- seqlevels(grl[[1]])

    transCoord <- function(gr) {
        gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
        gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
        return(gr)
    }
    grl <- endoapply(grl, transCoord)
    
    seqlevels(breaks) <- seqlevels(grl)
    seqlengths(breaks) <- seqlengths(grl)

    # disjoin overlaping breaks
    breaks <- unlist(breaks, use.names=FALSE)
    ranges <- disjoin(breaks) # redefine ranges in df
    hits <- countOverlaps(ranges, breaks) # counts number of breaks overlapping at each range
    mcols(ranges)$hits <- hits

    disjoin.breaks <- transCoord(ranges)
    dfplot.disjoin.breaks <- as.data.frame(disjoin.breaks)

    # Chromosome lines for heatmap
    label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(grl[[1]]) )
    df.chroms <- data.frame(y=c(0,cum.seqlengths))
    
    my_theme <- theme(
      legend.position="none",
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()  
    )

    # Plot breaks summary
    ggplt1 <- ggplot(dfplot.disjoin.breaks) + geom_rect(aes_string(ymin=0, ymax='hits', xmin='start.genome', xmax='end.genome'), fill="red", color="red")
    ggplt1 <- ggplt1 + geom_vline(aes_string(xintercept='y'), data=df.chroms, col='black') + scale_y_continuous(expand = c(0,0)) +ylab("Break\ncounts") + my_theme

    # Data
    df <- list()
    for (i1 in 1:length(grl)) {
        df[[length(df)+1]] <- data.frame(start=grl[[i1]]$start.genome, end=grl[[i1]]$end.genome, seqnames=seqnames(grl[[i1]]), sample=IDs[[i1]], state=grl[[i1]]$states)
    }
    df <- do.call(rbind, df)

    ## PLOT
    if (is.null(hotspots)) {
      ggplt2 <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20))
      ggplt2 <- ggplt2 + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
    } else {
      hotspots.trans <- transCoord(hotspots)
      hotspot.midpoint <- hotspots.trans$start.genome + (hotspots.trans$end.genome - hotspots.trans$start.genome)/2
      hotspot.midpoint <- data.frame(y=hotspot.midpoint)
      ggplt2 <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20))
      ggplt2 <- ggplt2 + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black') + geom_hline(aes_string(yintercept='y'), data=hotspot.midpoint, col='red', alpha=0.5)
    }  
   
    #p <- align.plots(ggplt1, ggplt2, ncol = 1, heights = c(1,10), draw = F)
    
    ## PRINT TO FILE
    ## printing to a file or returning plot object
    if (!is.null(file)) {
        message("Printing to PDF ",file)
        height.cm <- length(files2plot) * 0.5
        width.cm <- 200
        grDevices::pdf(file, width=width.cm/5, height=height.cm/5)
        print(ggplt2)
        d <- grDevices::dev.off()
        #p <- arrangeGrob(p)
        #ggsave(file=file, p, width=width.cm/5, height=height.cm/5, limitsize = F, device = "pdf")
    } else {
        #return(p)
        return(ggplt2)
    }
}

#######################################################################################################################
## Helper function for plot alignment from egg package

align.plots <- function(..., plots = list(...)) 
{
  n <- length(plots)
  grobs <- lapply(plots, ggplotGrob)
  
  heights <- c(1,10)
  widths <- NULL
  ncol = 1
  
  if (is.numeric(heights)) 
    heights <- lapply(heights, unit, "null")
  if (grid::is.unit(heights)) 
    widths <- as.unit.list(heights)
  
  if (is.null(widths)) 
    widths <- lapply(rep(1, n), unit, "null")
  if (grid::is.unit(widths)) 
    widths <- as.unit.list(widths)

  splits <- rep(1, n)

  fg <- mapply(gtable_frame, g = grobs, width = widths, height = heights, SIMPLIFY = FALSE)
  spl <- split(fg, splits)

  cols <- lapply(spl, function(.c) do.call(gridExtra::rbind.gtable, .c))
  all <- do.call(gridExtra::cbind.gtable, cols)

  if (draw) {
    if (newpage) 
      grid.newpage()
    grid.draw(all)
  }
  invisible(all)
}


gtable_frame <- function(g, width = unit(1, "null"), height = unit(1, "null"), debug = FALSE) 
{
  panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]), 
                          ]
  ll <- unique(panels$l)
  tt <- unique(panels$t)
  fixed_ar <- g$respect
  if (fixed_ar) {
    ar <- as.numeric(g$heights[tt[1]])/as.numeric(g$widths[ll[1]])
    height <- width * ar
    g$respect <- FALSE
  }
  core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
  top <- g[seq(1, min(tt) - 1), seq(min(ll), max(ll))]
  bottom <- g[seq(max(tt) + 1, nrow(g)), seq(min(ll), max(ll))]
  left <- g[seq(min(tt), max(tt)), seq(1, min(ll) - 1)]
  right <- g[seq(min(tt), max(tt)), seq(max(ll) + 1, ncol(g))]
  fg <- grid::nullGrob()
  if (length(left)) {
    lg <- gtable::gtable_add_cols(left, unit(1, "null"), 
                                  0)
    lg <- gtable::gtable_add_grob(lg, fg, 1, l = 1)
  }
  else {
    lg <- fg
  }
  if (length(right)) {
    rg <- gtable::gtable_add_cols(right, unit(1, "null"))
    rg <- gtable::gtable_add_grob(rg, fg, 1, l = ncol(rg))
  }
  else {
    rg <- fg
  }
  if (length(top)) {
    tg <- gtable::gtable_add_rows(top, unit(1, "null"), 0)
    tg <- gtable::gtable_add_grob(tg, fg, t = 1, l = 1)
  }
  else {
    tg <- fg
  }
  if (length(bottom)) {
    bg <- gtable::gtable_add_rows(bottom, unit(1, "null"), 
                                  -1)
    bg <- gtable::gtable_add_grob(bg, fg, t = nrow(bg), l = 1)
  }
  else {
    bg <- fg
  }
  grobs = list(fg, tg, fg, lg, core, rg, fg, bg, fg)
  widths <- unit.c(sum(left$widths), width, sum(right$widths))
  heights <- unit.c(sum(top$heights), height, sum(bottom$heights))
  all <- gtable::gtable_matrix("all", grobs = matrix(grobs, ncol = 3, nrow = 3, byrow = TRUE), widths = widths, heights = heights)
  if (debug) {
    hints <- grid::rectGrob(gp = gpar(fill = NA, lty = 2, lwd = 0.2))
    tl <- expand.grid(t = 1:3, l = 1:3)
    all <- gtable::gtable_add_grob(all, replicate(9, hints, simplify = FALSE), t = tl$t, l = tl$l, z = Inf, name = "debug")
  }
  all[["layout"]][5, "name"] <- "panel"
  if (fixed_ar) 
    all$respect <- TRUE
  all
}
