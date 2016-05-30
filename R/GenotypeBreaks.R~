#' Genotype the space between breakpoints
#'
#' Genotype the space between breakpoints with Fisher's exact test
#'
#' @param breaks A \code{\link{GRanges}} object with breakpoint coordinates.
#' @param fragments A \code{\link{GRanges}} object with read fragments.
#' @param backG The percent (e.g. 0.02 = 2\%) of background reads allowed for WW or CC genotype calls.
#' @param minReads The minimal number of reads between two breaks required for genotyping.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' @export
GenotypeBreaks<- function(breaks, fragments, backG=0.02, minReads=10)
{
	if (length(breaks)==0) {
		stop("argument 'breaks' is empty")
	}
	if (length(unique(as.character(seqnames(breaks))))>1) {
		stop("argument 'breaks' must contain only coordinates for one chromosome")
	}
	# Remove the non-used seqlevels
	breaks <- keepSeqlevels(breaks, value=as.character(unique(seqnames(breaks))))

  # function that genotypes a position based on number of + and - reads -> by Ashley Sanders
  genotype <- function(cReads, wReads, roiReads, bg=backG, minR=minReads, maxiter=10)
  {  ## FISHER EXACT TEST
		if (length(roiReads)==0) {
			return(c(NA,NA))
		}
		if (is.na(roiReads)) {
			return(c(NA,NA))
		}
    if( roiReads >= minR ) {
	    CCpVal<- fisher.test(matrix(c(cReads, wReads, round(roiReads*(1-bg)), round(roiReads*bg)), ncol=2, byrow=T))[[1]]
	    WCpVal<-fisher.test(matrix(c(cReads, wReads, round(roiReads*0.5), round(roiReads*0.5)), ncol=2, byrow=T))[[1]]
	    WWpVal<- fisher.test(matrix(c(wReads, cReads, round(roiReads*(1-bg)), round(roiReads*bg)), ncol=2, byrow=T))[[1]]
            iter <- 1
            while(CCpVal == WCpVal & WCpVal == WWpVal){ #if pVals are equal, take 10% of reads and recalculate
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
  ############################
  
  # create ranges between the breakpoints -> start and stops in a dataframe, use this to genotype between
	breaks.strand <- breaks
	strand(breaks.strand) <- '*'
	breakrange <- gaps(breaks.strand)
	breakrange <- breakrange[strand(breakrange)=='*']
	
  ## pull out reads of each line, and genotype in the fragments
	strand(breakrange) <- '-'
	Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
	strand(breakrange) <- '+'
	Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
	strand(breakrange) <- '*'
	readNo <- Ws + Cs
  
  ## bestFit genotype each region by Fisher Exact test
  fisher <- sapply(1:length(breakrange), function(x) genotype(Cs[x], Ws[x], readNo[x]))

	breakrange$readNo <- readNo
	breakrange$Ws <- Ws
	breakrange$Cs <- Cs
	breakrange$genoT <- fisher[1,]
	breakrange$pVal <- fisher[2,]
	breakrange <- breakrange[!is.na(breakrange$genoT)]

	## remove break if genotype is the same on either side of it
	equal.on.either.side <- c(breakrange$genoT[-length(breakrange)] == breakrange$genoT[-1], TRUE)
	breakrange.new <- breakrange[!equal.on.either.side]
	start(breakrange.new) <- end(breakrange[!equal.on.either.side])
	end(breakrange.new) <- start(breakrange[which(!equal.on.either.side)+1])
	breakrange.new$genoT <- paste(breakrange$genoT[!equal.on.either.side], breakrange$genoT[which(!equal.on.either.side)+1], sep='-')

  if (length(breakrange.new)) {	
  	return(breakrange.new)
  } else {
	return(breakrange.new<-NULL)
  }	
	
}
