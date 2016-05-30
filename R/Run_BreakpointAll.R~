#' Wrapper function for breakpointR
#'
#' This script will move through .bam files in a folder and perform several steps (see Details).
#'
#' 1. calculate deltaWs chromosome-by-chromsome
#' 2. localize breaks that pass zlim above the threshold
#' 3. genotype both sides of breaks to confirm whether strand state changes
#' 4. write a file of _reads, _deltaWs and _breaks in a chr fold -> can upload on to UCSC Genome browser
#' 5. write a file for each index with all chromosomes included -> can upload on to UCSC Genome browser
#'
#' @param bamfile A file with aligned reads in BAM format.
#' @param dataDirectory Output directory. If non-existent it will be created.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param chromosomes If only a subset of the chromosomes should be processed, specify them here.
#' @param windowsize The window size used to calculate deltaWs, either an integer or 'scale'.
#' @param trim The amount of outliers in deltaWs removed to calculate the stdev (10 will remove top 10\% and bottom 10\% of deltaWs).
#' @param peakTh The treshold that the peak deltaWs must pass to be considered a breakpoint.
#' @param zlim The number of stdev that the deltaW must pass the peakTh (ensures only significantly higher peaks are considered).
#' @param bg The amount of background introduced into the genotype test.
#' @param minReads The minimum number of reads required for genotyping.
#' @param writeBed If \code{TRUE}, will generate a bed of reads and deltaWs and breaks for upload onto the UCSC genome browser.
#' @param verbose Verbose messages if \code{TRUE}.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' @export

runBreakpointrALL <- function(datapath=NULL, dataDirectory='BreakPointR_analysis', pairedEndReads=TRUE, chromosomes=NULL, windowsize=100, scaleWindowSize=T, trim=10, peakTh=0.33, zlim=3.291, bg=0.02, minReads=10, writeBed=T, verbose=T, createCompositeFile=F, WC.cutoff=0.9) {

	options(warn=-1) # suppresses warnings in function
	files <- list.files(datapath, pattern=".bam$", full=T)

	if (!file.exists(dataDirectory)) {
		dir.create(dataDirectory)
	}

	deltas.all.files <- GenomicRanges::GRangesList()
	breaks.all.files <- GenomicRanges::GRangesList()
	counts.all.files <- GenomicRanges::GRangesList()

	if (createCompositeFile) {
		fragments <- createCompositeFile(file.list=files, chromosomes=chromosomes, pairedEndReads=pairedEndReads, WC.cutoff=WC.cutoff)
		deltas.breaks.counts.obj <- runBreakpointr(input.data=fragments, dataDirectory=dataDirectory, pairedEndReads=pairedEndReads, chromosomes=chromosomes, windowsize=windowsize, trim=trim, peakTh=peakTh, zlim=zlim, bg=bg, minReads=minReads, writeBed=writeBed, verbose=verbose)
		deltas <- unlist(deltas.breaks.counts.obj$deltas)
		breaks <- unlist(deltas.breaks.counts.obj$breaks)
		counts <- unlist(deltas.breaks.counts.obj$counts)
		parameters <- deltas.breaks.counts.obj$param

		names(deltas) <- NULL
		names(breaks) <- NULL
		names(counts) <- NULL

		deltas.all.files[['CompositeFile']] <- deltas
		counts.all.files[['CompositeFile']] <- counts

		if (length(breaks)) {
			breaks.all.files[['CompositeFile']] <- breaks
		}
	
	} else {

		for (bamfile in files) {
			message("Working on file ",bamfile)

			deltas.breaks.counts.obj <- runBreakpointr(input.data=bamfile, dataDirectory=dataDirectory, pairedEndReads=pairedEndReads, chromosomes=chromosomes, windowsize=windowsize, trim=trim, peakTh=peakTh, zlim=zlim, bg=bg, minReads=minReads, writeBed=writeBed, verbose=verbose)

			#deltas <- unlist(deltas.breaks.counts.obj$deltas)
			#breaks <- unlist(deltas.breaks.counts.obj$breaks)
			#counts <- unlist(deltas.breaks.counts.obj$counts)
			#parameters <- deltas.breaks.counts.obj$param

			#names(deltas) <- NULL
			#names(breaks) <- NULL
			#names(counts) <- NULL

			#deltas.all.files[[bamfile]] <- deltas
			#counts.all.files[[bamfile]] <- counts
			
			#if (length(breaks)) {
			#	breaks.all.files[[bamfile]] <- breaks
			#}
		}
	}
	
	if (createCompositeFile) {
		return(list(deltas=deltas.all.files, breaks=breaks.all.files, counts=counts.all.files, params=parameters))
	}
	
	options(warn=0) # turns warnings back on 
}

