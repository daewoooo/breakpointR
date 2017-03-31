#' Genotype the space between breakpoints
#'
#' Genotype the space between breakpoints with Fisher's exact test.
#'
#' @param breaks A \code{\link[GenomicRanges]{GRanges}} object with breakpoint coordinates.
#' @param fragments A \code{\link[GenomicRanges]{GRanges}} object with read fragments.
#' @param backG The percent (e.g. 0.02 = 2\%) of background reads allowed for WW or CC genotype calls.
#' @param minReads The minimal number of reads between two breaks required for genotyping.
#' @return A \code{\link[GenomicRanges]{GRanges}} object with genotyped breakpoint coordinates.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' @importFrom stats fisher.test
GenotypeBreaks <- function(breaks, fragments, backG=0.02, minReads=10)
{
    if (length(breaks)==0) {
        stop("argument 'breaks' is empty")
    }
    if (length(unique(as.character(seqnames(breaks))))>1) {
        stop("argument 'breaks' must contain only coordinates for one chromosome")
    }
    # Remove the non-used seqlevels
    breaks <- keepSeqlevels(breaks, value=as.character(unique(seqnames(breaks))))
  
    # create ranges between the breakpoints -> start and stops in a dataframe, use this to genotype between
    breaks.strand <- breaks
    strand(breaks.strand) <- '*'
    breakrange <- gaps(breaks.strand)
    breakrange <- breakrange[strand(breakrange)=='*']
    
    ## pull out reads of each line, and genotype in the fragments
    strand(breakrange) <- '-'
    breakrange$Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
    strand(breakrange) <- '+'
    breakrange$Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
    strand(breakrange) <- '*'
    breakrange$readNo <- breakrange$Ws + breakrange$Cs
    
    ## bestFit genotype each region by Fisher Exact test
    fisher <- sapply(1:length(breakrange), function(x) genotype.fisher(cReads=breakrange$Cs[x], wReads=breakrange$Ws[x], roiReads=breakrange$readNo[x], backG=backG, minReads=minReads))
  
    breakrange$genoT <- unlist(fisher[1,])
    breakrange$pVal <- unlist(fisher[2,])
    breakrange <- breakrange[!is.na(breakrange$genoT)]
    
    #genoTs <- cbind(breakrange$genoT[-length(breakrange)], breakrange$genoT[-1])
    #mask <- which(apply(genoTs, 1, function(x) all(is.na(x))))
    #breaks.strand <- breaks.strand[-mask]
    #breaks.strand$genoT <- genoTs[-mask,]
    
    ## remove break if genotype is the same on either side of it
    equal.on.either.side <- c(breakrange$genoT[-length(breakrange)] == breakrange$genoT[-1], TRUE)
    breakrange.new <- breakrange[!equal.on.either.side]
    start(breakrange.new) <- end(breakrange[!equal.on.either.side])
    end(breakrange.new) <- start(breakrange[which(!equal.on.either.side)+1])
    breakrange.new$genoT <- paste(breakrange$genoT[!equal.on.either.side], breakrange$genoT[which(!equal.on.either.side)+1], sep='-')
  
    ## Refine breakpoint regions to the highest deltaW in the given region
#    if (length(breakrange.new)) {  
#      hits <- findOverlaps(breakrange.new, breaks)
#      ToRefine <- split(breaks[subjectHits(hits)], queryHits(hits))
#      refined <- unlist(endoapply(ToRefine, RefineBreaks), use.names = F)
#      ranges(breakrange.new) <- ranges(refined)
#      breakrange.new$deltaW <- refined$deltaW
#    }  
    
    if (length(breakrange.new)) {
        ## Refine breakpoint regions to the highest deltaW in the given region
        hits <- GenomicRanges::findOverlaps(breakrange.new, breaks)
        ToRefine <- split(breaks[subjectHits(hits)], queryHits(hits))
        refined <- unlist(endoapply(ToRefine, RefineBreaks), use.names = F)
        ranges(breakrange.new) <- ranges(refined)
        breakrange.new$deltaW <- refined$deltaW
        return(breakrange.new)
    } else {
        return(breakrange.new<-NULL)
    }  
    
}


#' Genotype read numbers
#'
#' Genotype read numbers with Fisher's exact test.
#' 
#' @param cReads Number of Crick reads.
#' @param wReads Number of Watson reads.
#' @param roiReads Total number of reads.
#' @param backG Parameter for frequency of background reads.
#' @param minReads Minimal number of reads to perform the test.
#' @return A list with the $bestFit and $pval.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
genotype.fisher <- function(cReads, wReads, roiReads, backG=0.02, minReads=10)
{  ## FISHER EXACT TEST
    result <- list(bestFit=NA, pval=NA)
    if (length(roiReads)==0) {
        return(result)
    }
    if (is.na(roiReads)) {
        return(result)
    }
    if ( roiReads >= minReads ) {
        m <- matrix(c(cReads, wReads, round(roiReads*(1-backG)), ceiling(roiReads*backG)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
        CCpVal <- stats::fisher.test(m, alternative="greater")[[1]]
        m <- matrix(c(cReads, wReads, round(roiReads*0.5), ceiling(roiReads*0.5)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
        WCpVal <- 1 - stats::fisher.test(m, alternative="two.sided")[[1]]
        m <- matrix(c(wReads, cReads, round(roiReads*(1-backG)), ceiling(roiReads*backG)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Ws','Cs')))
        WWpVal <- stats::fisher.test(m, alternative="greater")[[1]]
        
        maxiter <- 5
        iter <- 0
        while (CCpVal == WCpVal & WCpVal == WWpVal) { #if pVals are equal, take 10% of reads and recalculate
            iter <- iter + 1
            cReads <-cReads*0.1
            wReads <- wReads*0.1
            roiReads <-roiReads*0.1
            m <- matrix(c(cReads, wReads, round(roiReads*(1-backG)), ceiling(roiReads*backG)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
            CCpVal <- stats::fisher.test(m, alternative="greater")[[1]]
            m <- matrix(c(cReads, wReads, round(roiReads*0.5), ceiling(roiReads*0.5)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
            WCpVal <- 1 - stats::fisher.test(m, alternative="two.sided")[[1]]
            m <- matrix(c(wReads, cReads, round(roiReads*(1-backG)), ceiling(roiReads*backG)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Ws','Cs')))
            WWpVal <- stats::fisher.test(m, alternative="greater")[[1]]
            if (iter == maxiter) { break }
        }
        
        pVal <- c(wc=WCpVal, cc=CCpVal, ww=WWpVal)
        result <- list(bestFit=names(pVal)[which.min(pVal)], pval=min(pVal))
        
        if (CCpVal == WCpVal & WCpVal == WWpVal) { 
          result$bestFit <- "?"
        }  
          
        return(result)
        
    } else { 
      
        return(result)
    }
}  

####################
# Helper functions #
####################

RefineBreaks <- function(gr) {
  maxDeltaW <- max(gr$deltaW)
  new.gr <- gr[gr$deltaW == maxDeltaW]
  new.gr <- mergeGR(new.gr)
  new.gr$deltaW <- maxDeltaW
  return(new.gr)
}

mergeGR <- function(gr) {
  gr <- sort(gr)
  #new.gr <- GRanges(seqnames=as.character(seqnames(gr))[1], ranges=IRanges(start=start(gr[1]), end=end(gr[length(gr)])))
  new.gr <- GRanges(seqnames=as.character(seqnames(gr))[1], ranges=IRanges(start=min(start(gr)), end=max(end(gr))))
  return(new.gr)
}


############################ old fisher

genotype <- function(cReads, wReads, roiReads, bg=0.02, minR=10, maxiter=10)
{  ## FISHER EXACT TEST
  if (length(roiReads)==0) {
    return(c(NA,NA))
  }
  if (is.na(roiReads)) {
    return(c(NA,NA))
  }
  if ( roiReads >= minR ) {
    CCpVal<- fisher.test(matrix(c(cReads, wReads, round(roiReads*(1-bg)), round(roiReads*bg)), ncol=2, byrow=T))[[1]]
    WCpVal<-fisher.test(matrix(c(cReads, wReads, round(roiReads*0.5), round(roiReads*0.5)), ncol=2, byrow=T))[[1]]
    WWpVal<- fisher.test(matrix(c(wReads, cReads, round(roiReads*(1-bg)), round(roiReads*bg)), ncol=2, byrow=T))[[1]]
    iter <- 1
    while (CCpVal == WCpVal & WCpVal == WWpVal) { #if pVals are equal, take 10% of reads and recalculate
      iter <- iter + 1
      cReads <-cReads*0.1
      wReads <- wReads*0.1
      roiReads <-roiReads*0.1
      CCpVal<- fisher.test(matrix(c(cReads, wReads, round(roiReads*(1-bg)), round(roiReads*bg)), ncol=2, byrow=T))[[1]]
      WCpVal<-fisher.test(matrix(c(cReads, wReads, round(roiReads*0.5), round(roiReads*0.5)), ncol=2, byrow=T))[[1]]
      WWpVal<- fisher.test(matrix(c(wReads, cReads, round(roiReads*(1-bg)), round(roiReads*bg)), ncol=2, byrow=T))[[1]]
      if (iter == maxiter) { break }
    }
    
    pVal<- cbind(CCpVal, WCpVal, WWpVal)
    maxP<- max(pVal)
    #if (pVal[which(pVal != maxP)][1] < 0.05 & pVal[which(pVal != maxP)][2] < 0.05) { signf <- '*'} else {signf <- 'ns'}
    if (maxP == CCpVal) { bestFit <- 'cc'} else if (maxP == WCpVal) { bestFit <- 'wc' } else{ bestFit <- 'ww'}
  } else { 
    return(c(NA,NA))
  }
  return(c(bestFit, maxP))
} 
