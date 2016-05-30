#' Function to create composite Strand-seq file
#'
#' This script will move through .bam files in a folder
#' read in each individual file and go thru each chr
#' determine the chr is WW or CC based on WCcutoff
#' reverse compliment all reads in the WW file
#' append to a new composite file for that chr
#' order the composite file of each chr based on pos
#'
#' @param datapath Location where bam files to be processed are stored
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param keep.duplicate.reads A logical indicating whether or not duplicate reads should be kept.
#' @param WC.cutoff Percentage of WW or CC reads to consider chromosome being WW or CC
#' @importFrom gtools mixedsort
#' @import GenomicRanges
#' @author Ashley Sanders, David Porubsky
#' @export

#for testing 
#chromosomes=paste0('chr',c(1:22)) 
#pairedEndReads=TRUE
#datapath = '.'
#WC.cutoff=0.90

createCompositeFile <- function(file.list, chromosomes=NULL, pairedEndReads=TRUE, min.mapq=10, keep.duplicate.reads=FALSE, WC.cutoff=0.90) {
	message("Creating composite file from ", length(file.list), " bam files")

	composite.bam.grl <- GenomicRanges::GRangesList()
	for (bamfile in file.list) {
		message("Working on file ",bamfile)
		fragments <- suppressWarnings( bam2GRanges(bamfile, pairedEndReads=pairedEndReads, chromosomes=chromosomes, min.mapq=min.mapq, keep.duplicate.reads=keep.duplicate.reads) )
	        if (length(fragments) > 0) { mcols(fragments)$lib <- bamfile } # appends file name to reads
		composite.bam <- GenomicRanges::GRangesList()
		for (chr in unique(seqnames(fragments))) {
			#message("Working on chromosome ",chr)
			chr.fragments <- fragments[seqnames(fragments)==chr]
			tab.counts <- BiocGenerics::table(strand(chr.fragments))
			plus.count <- tab.counts[1]
			minus.count <- tab.counts[2]

			## Ashley's secret formula
			wcCall <- round( ( minus.count - plus.count ) / length(chr.fragments), digits=3 )

			## check if this is a pure (WW or CC) library
			if( wcCall <= -WC.cutoff | wcCall >= WC.cutoff ) {
				## if this is pure WW, reverse compliment all the reads
				if(wcCall <= -WC.cutoff) {
					chr.fragments.copy <- chr.fragments
					strand(chr.fragments)[strand(chr.fragments.copy) == '+'] <- '-'
					strand(chr.fragments)[strand(chr.fragments.copy) == '-'] <- '+'
				}
			composite.bam[[chr]] <- chr.fragments
			}
		}	
		composite.bam.grl[[bamfile]] <- unlist(composite.bam)
	}
	composite.data <- unlist(composite.bam.grl)
	names(composite.data) <- NULL
	seqlevels(composite.data) <- gtools::mixedsort( as.character(unique(seqnames(composite.data))) )
	return(sort(composite.data))
}

