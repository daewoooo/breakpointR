#' @import ggplot2
#' @import reshape2
#' @import grid
#' @importFrom cowplot plot_grid
NULL


#' Plotting function for BreakPointR
#' This function will take stored RData files and plot them into a single PDF file
#'
#' @param datapath location of RData files storing reads and breakpoint coordinates
#' @param plotLibraries subset of libaries to be plotted, otherwise all libraries will be plotted
#' @param file name of the file to store libraries in
#'
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @export

plotBreakpoints <- function(datapath, plotLibraries = NULL, file=NULL) {

files <- list.files(datapath, pattern=".RData$", full=T)

if (!is.null(plotLibraries)) {
	libs2plot <- plotLibraries
} else {
	libs2plot <- c(1:length(files))
}

plots <- list()
for (i in libs2plot) {
	message("Preparing plot for ",files[[i]])

	data <- get(load(files[[i]]))

	filename <- basename(files[i])

	reads <- data$fragments
	breaks <- data$breaks
	counts <- data$counts

	#if (length(reads) < 10000) next

	bin_size <- 200000
	binned.data <- GenomicRanges::GRangesList()
	chrom2use <- levels(seqnames(reads))
	seqlengths <- seqlengths(reads)

	for (j in chrom2use) {
		chrom_gr <- reads[reads@seqnames == j,]
		max_pos <- max(start(chrom_gr))
		numbins <- floor(max_pos/bin_size) #calculate number of bins
		ir <- successiveIRanges(rep(bin_size,numbins)) #create continuous ranges based on number of bins and the binsize
		lastBinWidth <- max_pos - (bin_size*numbins)
		lastBinStart <- (max_pos - lastBinWidth) + 1
		lastRange <- IRanges(start = lastBinStart, width = lastBinWidth) #calculate last range
		ir <- c(ir,lastRange)
		chr <- rep(j, length(ir))
		rows <- rep(0,length(ir))
		gr <- GRanges(seqnames=Rle(chr), ranges=ir, Watsonreads=rows, Crickreads=rows, bothreads=rows) #initialize GRanges object to store bin read counts

		#counts overlaps between bins and our reads
		Watsonreads <- GenomicRanges::countOverlaps(gr, chrom_gr[strand(chrom_gr)=='-']) 
		Crickreads <- GenomicRanges::countOverlaps(gr, chrom_gr[strand(chrom_gr)=='+'])
		bothreads <- Watsonreads + Crickreads
	
		mcols(gr)$bothreads <- bothreads
		mcols(gr)$Watsonreads <- Watsonreads
		mcols(gr)$Crickreads <- Crickreads
	
		binned.data[[j]] <- gr #store binned data for current chromosome in a GRanges list
	}

	seqlengths(binned.data) <- seqlengths

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

	binned.data <- unname(unlist(binned.data))

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
	ggplt1 <- ggplot(dfplot.reads) + geom_linerange(aes(ymin=0, ymax=-Crickreads, x=midpoint), color="paleturquoise4", size=0.2)
	ggplt1 <- ggplt1 + geom_linerange(aes(ymin=0, ymax=Watsonreads, x=midpoint), color="sandybrown", size=0.2) + scale_y_continuous(limits = c(-limit,limit))
	ggplt1 <- ggplt1 + geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black') + xlab("Chromosomes") + ylab("States") + scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + my_theme + ggtitle(filename)

### PLOT COUNTS

	if (length(breaks)) {
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
		ggplt2 <- ggplt2 + geom_point(data=dfplot.breaks, aes(x=midpoint, y=0), size=5, color='red', shape=124, inherit.aes = FALSE) + xlab(NULL) + ylab("Strands") + scale_x_continuous(expand = c(0,0)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + my_theme


### PLOT STATES

		ggplt3 <- ggplot(dfplot.counts) + geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax=10, fill='states'))
		ggplt3 <- ggplt3 + geom_linerange(data=chr.lines, aes_string(ymin=0, ymax=10, x='y'), col='black') + xlab("Chromosomes") + ylab("States") + 	scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) + scale_fill_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + my_theme
	}

### PLOT ALL

	if (length(breaks)) {
		p <- cowplot::plot_grid(ggplt1, ggplt2, ggplt3, ncol=1, align="v", rel_heights = c(3,3,2))
	} else {
		p <- cowplot::plot_grid(ggplt1, ncol=1, align="v", rel_heights = 3)
	}
	plots[[length(plots)+1]] <- p
}

### PRINT TO FILE

## printing to a file or returning plot object
	if (!is.null(file)) {
		message("Printing to PDF ",file)

		pdf(file, width=30, height=10)
		bquiet = lapply(plots, print)
		d <- dev.off()
	} else {
		return(plots)
	}
}



#' Genome wide heatmap of template inheritance states
#'
#' Plot a genome wide heatmap of template inheritance states.
#'
#' @param datapath location of RData files storing reads and breakpoint coordinates
#' @param plotLibraries subset of libaries to be plotted, otherwise all libraries will be plotted
#' @param file name of the file to store libraries in
#' 
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @export


plotHeatmap <- function(datapath, plotLibraries = NULL, file=NULL) {

	files <- list.files(datapath, pattern=".RData$", full=T)

	if (!is.null(plotLibraries)) {
		libs2plot <- plotLibraries
	} else {
		libs2plot <- c(1:length(files))
	}

	message("Preparing heatmap from ",length(libs2plot), " libraries")

	IDs <- list()
	grl <- GRangesList()
	for (i in libs2plot) {	
		data <- get(load(files[[i]]))
		IDs[[i]] <- basename(files[[i]])
		grl[[i]] <- data$counts
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

	# Data
	df <- list()
	for (i1 in 1:length(grl)) {
		df[[length(df)+1]] <- data.frame(start=grl[[i1]]$start.genome, end=grl[[i1]]$end.genome, seqnames=seqnames(grl[[i1]]), sample=IDs[[i1]], state=grl[[i1]]$states)
	}
	df <- do.call(rbind, df)
	# Chromosome lines
	label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(grl[[1]]) )
	df.chroms <- data.frame(y=c(0,cum.seqlengths))

	## PLOT
	ggplt <- ggplot(df) + geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) + scale_y_continuous(breaks=label.pos, labels=names(label.pos)) + coord_flip() + scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) + theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20))
	ggplt <- ggplt + geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')

	## PRINT TO FILE
	## printing to a file or returning plot object
	if (!is.null(file)) {
		message("Printing to PDF ",file)
		height.cm <- length(libs2plot) * 0.5
		width.cm <- 200
		pdf(file, width=width.cm/2.54, height=height.cm/2.54)
		print(ggplt)
		d <- dev.off()
	} else {
		return(ggplt)
	}
}


