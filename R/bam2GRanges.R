#' Import BAM file into GRanges
#'
#' Import aligned reads from a BAM file into a \code{\link{GRanges}} object.
#'
#' @param file Bamfile with aligned reads.
#' @param bamindex Bam-index file with or without the .bai ending. If this file does not exist it will be created and a warning is issued.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param min.mapq Minimum mapping quality when importing from BAM files.
#' @param keep.duplicate.reads A logical indicating whether or not duplicate reads should be kept.
#' @importFrom Rsamtools indexBam scanBamHeader ScanBamParam scanBamFlag
#' @importFrom GenomicAlignments readGAlignmentPairs readGAlignments first last
#' @author Aaron Taudt, David Porubsky, Ashley Sanders
#' @export

bam2GRanges <- function(file, bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10, keep.duplicate.reads=TRUE, pair2frgm=FALSE) {

	## Check if bamindex exists
	bamindex.raw <- sub('\\.bai$', '', bamindex)
	bamindex <- paste0(bamindex.raw,'.bai')
	if (!file.exists(bamindex)) {
		bamindex.own <- Rsamtools::indexBam(file)
		warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
		bamindex <- bamindex.own
	}
	file.header <- Rsamtools::scanBamHeader(file)[[1]]
	chrom.lengths <- file.header$targets
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}
	chroms2use <- intersect(chromosomes, chroms.in.data)
	if (length(chroms2use)==0) {
		chrstring <- paste0(chromosomes, collapse=', ')
		stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr',chromosomes), collapse=', '), ' instead.')
	}
	## Issue warning for non-existent chromosomes
	diff <- setdiff(chromosomes, chroms.in.data)
	if (length(diff)>0) {
		diffs <- paste0(diff, collapse=', ')
		warning(paste0('Not using chromosomes ', diffs, ' because they are not in the data.'))
	}
	## Import the file into GRanges
	gr <- GenomicRanges::GRanges(seqnames=Rle(chroms2use), ranges=IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))
	if (keep.duplicate.reads) {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq'))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq'))
		}
	} else {
		if (pairedEndReads) {
			data.raw <- GenomicAlignments::readGAlignmentPairs(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		} else {
			data.raw <- GenomicAlignments::readGAlignments(file, index=bamindex, param=Rsamtools::ScanBamParam(tag="XA", which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))
		}
	}

	if (pairedEndReads) {
		
		if (pair2frgm) {
			if (!is.null(min.mapq)) { min.mapq <- 0 }	
			data.prop.pairs <- data.raw[GenomicAlignments::isProperPair(data.raw)] #only proper pairs can be merged into single fragment (no negative ranges)

			data.first <- as(GenomicAlignments::first(data.prop.pairs), 'GRanges')
			data.last <- as(GenomicAlignments::last(data.prop.pairs), 'GRanges')

			## filter XA tag
			#data.first <- data.first[is.na(mcols(data.first)$XA)]
			#data.last <- data.last[is.na(mcols(data.last)$XA)]
			data.first.filt <- is.na(mcols(data.first)$XA)
			data.last.filt <- is.na(mcols(data.last)$XA)

			mask <- data.first.filt & data.last.filt

			data.first <- data.first[mask]
			data.last <- data.last[mask]

			## Filter by mapping quality
			if (!is.null(min.mapq)) {
				if (any(is.na(mcols(data.first)$mapq)) | any(is.na(mcols(data.last)$mapq))) {
					warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
					mcols(data.first)$mapq[is.na(mcols(data.first)$mapq)] <- -1
					mcols(data.last)$mapq[is.na(mcols(data.last)$mapq)] <- -1
				}
			
				data.first.filt <- mcols(data.first)$mapq >= min.mapq
				data.last.filt <- mcols(data.last)$mapq >= min.mapq
			
				mask <- data.first.filt & data.last.filt

				data.first <- data.first[mask]
				data.last <- data.last[mask]
			
				#take reads where not both mates have the expected mapping quality
				#data.first.single <- data.first[!mask]
				#data.last.single <- data.last[!mask]
			
				#data.singlets <- c(data.first.single, data.last.single)
				#data.singlets <- data.singlets[mcols(data.singlets)$mapq >= min.mapq] #filter singeltons by mapping quality
			}

			#split reads by directionality
			data.first.plus <- data.first[strand(data.first) == '+']
			data.first.minus <- data.first[strand(data.first) == '-']
			data.last.plus <- data.last[strand(data.last) == '+']
			data.last.minus <- data.last[strand(data.last) == '-']

			#merge pairs into a single range
			frag.plus.mapq <- data.first.plus$mapq + data.last.minus$mapq
			frag.minus.mapq <- data.first.minus$mapq + data.last.plus$mapq
		
			data.frag.plus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.plus), ranges=IRanges(start=start(data.first.plus), end=end(data.last.minus)), strand=strand(data.first.plus), mapq=frag.plus.mapq)
			seqlengths(data.frag.plus) <- seqlengths(data.first)
			data.frag.minus <- GenomicRanges::GRanges(seqnames=seqnames(data.first.minus), ranges=IRanges(start=start(data.last.plus), end=end(data.first.minus)), strand=strand(data.first.minus), mapq=frag.minus.mapq)
			seqlengths(data.frag.minus) <- seqlengths(data.first)
			#data.singlets <- GenomicRanges::GRanges(seqnames=seqnames(data.singlets), ranges=IRanges(start=start(data.singlets), end=end(data.singlets)), strand=strand(data.singlets), mapq=data.singlets$mapq)
			#seqlengths(data.singlets) <- seqlengths(data.first)
		
			#data <- sort(c(data.frag.plus, data.frag.minus, data.singlets))
			data <- sort(c(data.frag.plus, data.frag.minus))
		} else {
			#data.first <- as(GenomicAlignments::first(data.raw), 'GRanges')
			#data.last <- as(GenomicAlignments::last(data.raw), 'GRanges')
			#strand(data.last) <- strand(data.first)
			#data <- sort(c(data.first, data.last))
			data <- as(GenomicAlignments::first(data.raw), 'GRanges') #use only first mate of the read pair in subsequent analysis!!!

			## Filter by mapping quality
			if (!is.null(min.mapq)) {
				if (any(is.na(mcols(data)$mapq))) {
					warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
					mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
				}
			data <- data[mcols(data)$mapq >= min.mapq]
			}
			## Filter XA tag
			data <- data[is.na(mcols(data)$XA)]
		}

	} else {
		data <- as(data.raw, 'GRanges')

		## Filter by mapping quality
		if (!is.null(min.mapq)) {
			if (any(is.na(mcols(data)$mapq))) {
				warning(paste0(file,": Reads with mapping quality NA (=255 in BAM file) found and removed. Set 'min.mapq=NULL' to keep all reads."))
				mcols(data)$mapq[is.na(mcols(data)$mapq)] <- -1
			}
		data <- data[mcols(data)$mapq >= min.mapq]
		}
		## Filter XA tag
		data <- data[is.na(mcols(data)$XA)]
	}
	
	seqlevels(data) <- seqlevels(gr)
	return(data)
}


