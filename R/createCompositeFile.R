#' Create composite Strand-seq file
#'
#' This function will move through BAM files in a folder, read in each individual file and go through each chromosome, 
#' determine if the chromosome is WW or CC based on WCcutoff, reverse complement all reads in the WW file,
#' append to a new composite file for that chromosome, order the composite file of each chromosome based on position.
#'
#' @param file.list TODO:DAVID
#' @inheritParams readBamFileAsGRanges
#' @param WC.cutoff Percentage of WW or CC reads to consider chromosome being WW or CC
#' @param bg The amount of background introduced into the genotype test.
#' @return TODO:DAVID
#' @importFrom gtools mixedsort
#' @importFrom BiocGenerics table
#' @import GenomicRanges
#' @author Ashley Sanders, David Porubsky
#' @export
#' @examples
#'## Get some example files
#'inputfolder <- system.file("extdata", "breakpointer", package="strandseqExampleData")
#'files <- list.files(inputfolder, full.names=TRUE, pattern="bam$")
#'## Create the composite file
#'composite <- createCompositeFile(files, chromosomes=paste0('chr', c(1:22, 'X')),
#'                    pairedEndReads=FALSE)
#'                    
createCompositeFile <- function(file.list, chromosomes=NULL, pairedEndReads=TRUE, pair2frgm=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, WC.cutoff=0.90, bg=0.02) {

    message("Creating composite file from ", length(file.list), " bam files")

    composite.bam.grl <- GenomicRanges::GRangesList()
    for (bamfile in file.list) {
        message("Working on file ",bamfile)
        fragments <- suppressWarnings( readBamFileAsGRanges(bamfile, pairedEndReads=pairedEndReads, chromosomes=chromosomes, min.mapq=min.mapq, remove.duplicate.reads=remove.duplicate.reads, pair2frgm=pair2frgm) )
        if (length(fragments) > 0) { mcols(fragments)$lib <- bamfile } # appends file name to reads
        composite.bam <- GenomicRanges::GRangesList()
        for (chr in unique(seqnames(fragments))) {
            #message("Working on chromosome ",chr)
            chr.fragments <- fragments[seqnames(fragments)==chr]
            tab.counts <- BiocGenerics::table(strand(chr.fragments))
            plus.count <- tab.counts[1]
            minus.count <- tab.counts[2]

            # ## Ashley's secret formula
            # wcCall <- round( ( minus.count - plus.count ) / length(chr.fragments), digits=3 )
            # 
            # ## check if this is a pure (WW or CC) library
            # if( wcCall <= -WC.cutoff | wcCall >= WC.cutoff ) {
            #     ## if this is pure WW, reverse compliment all the reads
            #     if(wcCall <= -WC.cutoff) {
            #         chr.fragments.copy <- chr.fragments
            #         strand(chr.fragments)[strand(chr.fragments.copy) == '+'] <- '-'
            #         strand(chr.fragments)[strand(chr.fragments.copy) == '-'] <- '+'
            #     }
            #     composite.bam[[chr]] <- chr.fragments
            # }
            
            ## Alternatively use Fisher's test to genotype
            geno <- genotype.fisher(cReads=plus.count, wReads=minus.count, roiReads=plus.count+minus.count, backG=bg, minReads = 10)
            ## check if this is a pure (WW or CC) library
            if (geno$bestFit == 'ww' | geno$bestFit == 'cc') {
                ## if this is pure WW, reverse compliment all the reads
                if (geno$bestFit == 'ww') {
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
    composite.data <- sort(composite.data)
    return(composite.data)
}

