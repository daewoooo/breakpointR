#' Find breakpoints in Strand-seq data
#'
#' Find breakpoints in Strand-seq data. See section Details on how breakpoints are located.
#'
#' Breakpoints are located in the following way:
#' \enumerate{
#' \item calculate deltaWs chromosome-by-chromsome
#' \item localize breaks that pass zlim above the threshold
#' \item genotype both sides of breaks to confirm whether strand state changes
#' \item write a file of _reads, _deltaWs and _breaks in a chr fold -> can upload on to UCSC Genome browser
#' \item write a file for each index with all chromosomes included -> can upload on to UCSC Genome browser
#' }
#'
#' @param bamfile A file with aligned reads in BAM format.
#' @param ID A character string that will serve as identifier in downstream functions.
#' @param windowsize The window size used to calculate deltaWs, either number of reads or genomic size depending on binMethod
#' @inheritParams readBamFileAsGRanges
#' @param binMethod Method used to calculate optimal number of reads in the window ("size", "reads"). Default = "size"
#' @param trim The amount of outliers in deltaWs removed to calculate the stdev (10 will remove top 10\% and bottom 10\% of deltaWs).
#' @param peakTh The treshold that the peak deltaWs must pass to be considered a breakpoint.
#' @param zlim The number of stdev that the deltaW must pass the peakTh (ensures only significantly higher peaks are considered).
#' @param bg The amount of background introduced into the genotype test.
#' @param minReads The minimum number of reads required for genotyping.
#' @param maskRegions List of regions to be excluded from the analysis (format: chromosomes start end)
#' @return A \code{\link{BreakPoint}} object.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' @importFrom utils flush.console
#' @importFrom S4Vectors subjectHits
#' @export
#' @examples 
#'## Get an example file
#'exampleFolder <- system.file("extdata", "breakpointer", package="strandseqExampleData")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Run breakpointR
#'brkpts <- runBreakpointr(exampleFile, pairedEndReads=FALSE)
#'
runBreakpointr <- function(bamfile, ID=basename(bamfile), pairedEndReads=TRUE, chromosomes=NULL, windowsize=1e6, binMethod="size", trim=10, peakTh=0.33, zlim=3.291, bg=0.02, min.mapq=10, pair2frgm=FALSE, minReads=20, maskRegions=NULL, estimate.back=FALSE) {

    ## check the class of the bamfile, make GRanges object of file
    if ( class(bamfile) != "GRanges" ) {
        ptm <- startTimedMessage("Reading file ", bamfile, " ...")
        suppressWarnings( fragments <- readBamFileAsGRanges(bamfile, pairedEndReads=pairedEndReads, chromosomes=chromosomes, remove.duplicate.reads=TRUE, min.mapq=min.mapq, pair2frgm=pair2frgm) )
        stopTimedMessage(ptm)
    } else {
        fragments <- bamfile
        bamfile <- 'CompositeFile'
    }
  
    filename <- basename(bamfile)
    message("Working on ", filename)
    
    ## remove reads from maskRegions
    if (!is.null(maskRegions)) {
        mask <- findOverlaps(maskRegions, fragments)
        fragments <- fragments[-S4Vectors::subjectHits(mask)]        
    } 
  
    reads.all.chroms <- fragments
  
    deltas.all.chroms <- GenomicRanges::GRangesList()
    #seqlevels(deltas.all.chroms) <- seqlevels(fragments)
    breaks.all.chroms <- GenomicRanges::GRangesList()
    #seqlevels(breaks.all.chroms) <- seqlevels(fragments)
    counts.all.chroms <- GenomicRanges::GRangesList()
    #seqlevels(counts.all.chroms) <- seqlevels(fragments)
  
    for (chr in unique(seqnames(fragments))) {
        message("  Working on chromosome ",chr)
        fragments.chr <- fragments[seqnames(fragments)==chr]
        fragments.chr <- keepSeqlevels(fragments.chr, chr)
    
        ptm <- startTimedMessage("    calculating deltaWs ...")
        if (binMethod == "size") {
            ## normalize for size of chromosome one and for read counts of each chromosome
            if (windowsize > seqlengths(fragments)[chr]) {
                stopTimedMessage(ptm)
                next
            }
            tiles <- unlist(tileGenome(seqlengths(fragments)[chr], tilewidth = windowsize))
            counts <- countOverlaps(tiles, fragments.chr)
            reads.per.window <- max(10, round(mean(counts[counts>0], trim=0.05)))
          
            dw <- deltaWCalculator(fragments.chr, reads.per.window=reads.per.window)
        } else if (binMethod == "reads") {
            ## normalize only for size of the chromosome 1
            #reads.per.window <- max(10, round(windowsize/(seqlengths(fragments)[1]/seqlengths(fragments)[chr]))) # scales the bin to chr size, anchored to chr1 (249250621 bp)
            ## do not normalize to the size of the chromosome 1
            reads.per.window <- windowsize 
            dw <- deltaWCalculator(fragments.chr, reads.per.window=reads.per.window)
        }
        deltaWs <- dw[seqnames(dw)==chr]
        stopTimedMessage(ptm)
    
        ptm <- startTimedMessage("    finding breakpoints ...")
        breaks <- suppressMessages( breakSeekr(deltaWs, trim=trim, peakTh=peakTh, zlim=zlim) )
        stopTimedMessage(ptm)
    
        ## Genotyping
        if (length(breaks) >0 ) {
            maxiter <- 10
            iter <- 1
            ptm <- startTimedMessage("    genotyping ",iter, " ...")
            utils::flush.console()
            newBreaks <- GenotypeBreaks(breaks, fragments, backG=bg, minReads=minReads)
            prev.breaks <- breaks  
            breaks <- newBreaks
            stopTimedMessage(ptm)
            
            while (length(prev.breaks) > length(newBreaks) && !is.null(breaks) ) {
                utils::flush.console()
                iter <- iter + 1
                ptm <- startTimedMessage("    genotyping ",iter, " ...")
                newBreaks <- GenotypeBreaks(breaks, fragments, backG=bg, minReads=minReads)
                prev.breaks <- breaks  
                breaks <- newBreaks
                stopTimedMessage(ptm)
                if (iter == maxiter) { break }
            }
        } else {
            newBreaks<-NULL # assigns something to newBreaks if no peaks found
        }
    
        ### count reads between the breaks and write into GRanges
        if (is.null(newBreaks)) {
          
            Ws <- length(which(as.logical(strand(fragments.chr) == '-')))
            Cs <- length(which(as.logical(strand(fragments.chr) == '+')))
            chrRange <- GRanges(seqnames=chr, ranges=IRanges(start=1, end=seqlengths(fragments.chr)[chr]))
            counts <- cbind(Ws,Cs)
            mcols(chrRange) <- counts
            seqlengths(chrRange) <- seqlengths(fragments.chr)[chr]
      
            ## genotype entire chromosome
            WC.ratio <- (chrRange$Ws-chrRange$Cs)/sum(c(chrRange$Ws,chrRange$Cs))
            if (WC.ratio > 0.8) {
                state <- 'ww'
            } else if (WC.ratio < -0.8) {
                state <- 'cc'
            } else if (WC.ratio < 0.2 & WC.ratio > -0.2) {
                state <- 'wc'
            } else {
                state <- '?'
            }          
            chrRange$states <- state
            suppressWarnings( counts.all.chroms[[chr]] <- chrRange )
      
        } else {
          
            breaks.strand <- newBreaks
            strand(breaks.strand) <- '*'
            breakrange <- gaps(breaks.strand)
            breakrange <- breakrange[strand(breakrange)=='*']
      
            ## pull out reads of each line, and genotype in the fragments
            strand(breakrange) <- '-'
            Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
            strand(breakrange) <- '+'
            Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
            strand(breakrange) <- '*'
      
            ## pull out states for each region between the breakpoints
            concat.states <- paste(newBreaks$genoT, collapse="-")
            split.states <- unlist(strsplit(concat.states, "-"))
            states.idx <- seq(from=1,to=length(split.states), by=2)
            states.idx <- c(states.idx, length(split.states))
            states <- split.states[states.idx]
      
            counts <- cbind(Ws,Cs)
            mcols(breakrange) <- counts
            breakrange$states <- states
            suppressWarnings( counts.all.chroms[[chr]] <- breakrange )
        }
    
        ### write breaks and deltas into GRanges
        if (length(deltaWs) > 0) {
            suppressWarnings( deltas.all.chroms[[chr]] <- deltaWs[,'deltaW'] )  #select only deltaW metadata column to store
        }
        if (length(newBreaks) > 0) {
            suppressWarnings( breaks.all.chroms[[chr]] <- newBreaks )
        }
        
    }
    
    ## creating list of list where filename is first level list ID and deltas, breaks and counts are second list IDs
    reads.all.chroms <- unlist(reads.all.chroms, use.names=FALSE)  
    deltas.all.chroms <- unlist(deltas.all.chroms, use.names=FALSE)
    breaks.all.chroms <- unlist(breaks.all.chroms, use.names=FALSE)
    counts.all.chroms <- unlist(counts.all.chroms, use.names=FALSE)
    
    seqlevels(deltas.all.chroms) <- seqlevels(reads.all.chroms)
    seqlevels(breaks.all.chroms) <- seqlevels(reads.all.chroms)
    seqlevels(counts.all.chroms) <- seqlevels(reads.all.chroms)
    
    if (estimate.back) {
      counts.all.chroms[counts.all.chroms$states == 'ww'] -> ww
      counts.all.chroms[counts.all.chroms$states == 'cc'] -> cc
      mean(sum(ww$Cs) / sum(ww$Ws), sum(cc$Ws) / sum(cc$Cs)) -> bg
      #mean(sum(ww$Cs/width(ww)) / sum(ww$Ws/width(ww)), sum(cc$Ws/width(cc)) / sum(cc$Cs/width(cc))) -> bg
      #mean( c(ww$Cs / ww$Ws , cc$Ws / cc$Cs) )-> bg
      fisher <- sapply(1:length(counts.all.chroms), function(x) genotype.fisher(cReads=counts.all.chroms$Cs[x], wReads=counts.all.chroms$Ws[x], roiReads=counts.all.chroms$Cs[x]+counts.all.chroms$Ws[x], backG=bg, minReads=minReads))
      counts.all.chroms$new.states <- unlist(fisher[1,])
    }
    
    ## save set parameters for future reference
    parameters <- c(windowsize=windowsize, binMethod=binMethod, trim=trim, peakTh=peakTh, zlim=zlim, bg=bg, minReads=minReads)
  
    data.obj <- list(ID=ID, fragments=fragments, deltas=deltas.all.chroms, breaks=breaks.all.chroms, counts=counts.all.chroms, params=parameters)
    class(data.obj) <- class.breakpoint
  
    return(data.obj)
}

