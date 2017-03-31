#' Main function for the \pkg{\link{breakpointR}} package
#'
#' This function is an easy-to-use wrapper to find breakpoints with \code{\link[breakpointR]{runBreakpointr}} in parallel, write the results to file and find hotspots.
#'
#' @param inputfolder Folder with BAM files.
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU The numbers of CPUs that are used. Should not be more than available on your machine.
#' @param reuse.existing.files A logical indicating whether or not existing files in \code{outputfolder} should be reused.
#' @param createCompositeFile Merge single cell data into a single file to increase breakpoint mapping resolution (do not use for SCE mapping)
#' @param callHotSpots Search for regions of high abundance of breakpoints in single cells (set createCompositeFile=FALSE)

#' @inheritParams readBamFileAsGRanges
#' @inheritParams runBreakpointr
#' @inheritParams readBamFileAsGRanges
#' @inheritParams createCompositeFile

#' @return \code{NULL}
#' @author Aaron Taudt, David Porubsky, Ashley Sanders
#' @importFrom utils read.table
#' @import foreach
#' @import doParallel
#' @export
#' 
#' @examples
#'## Get some example files
#'inputfolder <- system.file("extdata", "breakpointer", package="strandseqExampleData")
#'outputfolder <- file.path(tempdir(), "breakpointer_example")
#'## Run breakpointer
#'breakpointer(inputfolder, outputfolder, chromosomes=paste0("chr", c(1:22, "X")),
#'             pairedEndReads=FALSE)
#'
breakpointer <- function(inputfolder, outputfolder, configfile=NULL, numCPU=1, reuse.existing.files=TRUE, windowsize=1e6, binMethod="size", pairedEndReads=FALSE, pair2frgm=FALSE, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, trim=10, peakTh=0.33, zlim=3.291, bg=0.02, minReads=10, createCompositeFile=FALSE, WC.cutoff=0.9, maskRegions=NULL, callHotSpots=FALSE, estimate.back=FALSE) {

#=======================
### Helper functions ###
#=======================
as.object <- function(x) {
    return(eval(parse(text=x)))
}

#========================
### General variables ###
#========================
conf <- NULL
if (is.character(configfile)) {
    ## Read config file ##
    errstring <- tryCatch({
        conf <- readConfig(configfile)
        errstring <- ''
    }, error = function(err) {
        errstring <- paste0("Could not read configuration file ",configfile)
    })
    if (errstring!='') {
        stop(errstring)
    }
}

### Mask Regions ###
if (!is.null(maskRegions)) {
    mask.df <- utils::read.table(maskRegions, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    if (any(mask.df$V1 %in% chromosomes)) {
        mask.df <- mask.df[mask.df$V1 %in% chromosomes,] #select only chromosomes submitted for analysis
        maskRegions <- GenomicRanges::GRanges(seqnames=mask.df$V1, IRanges(start=mask.df$V2, end=mask.df$V3))      
    } else {
        warning(paste0('Skipping maskRegions. Chromosomes names conflict: ', mask.df$V1[1], ' not equal ',chromosomes[1]))
        maskRegions <- NULL
    }  
}

## Put options into list and merge with conf
params <- list(numCPU=numCPU, reuse.existing.files=reuse.existing.files, windowsize=windowsize, binMethod=binMethod, pairedEndReads=pairedEndReads, pair2frgm=pair2frgm, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, trim=trim, peakTh=peakTh, zlim=zlim, bg=bg, minReads=minReads, createCompositeFile=createCompositeFile, WC.cutoff=WC.cutoff, maskRegions=maskRegions, callHotSpots=callHotSpots, estimate.back=estimate.back)
conf <- c(conf, params[setdiff(names(params),names(conf))])

## Input checks
# None at the moment

## Helpers
windowsize <- conf[['windowsize']]

## Set up the directory structure ##
datapath <- file.path(outputfolder,'data')
browserpath <- file.path(outputfolder,'browserfiles')

## Delete old directory if desired ##
if (conf[['reuse.existing.files']]==FALSE) {
    if (file.exists(outputfolder)) {
        message("Deleting old directory ",outputfolder)
        unlink(outputfolder, recursive=TRUE)
    }
}
if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
}
if (!file.exists(datapath)) { dir.create(datapath) }
if (!file.exists(browserpath)) { dir.create(browserpath) }
## Make a copy of the conf file
writeConfig(conf, configfile=file.path(outputfolder, 'breakpointR.config'))


#=====================================
# Find breakpoints and write BED files
#=====================================
files <- list.files(inputfolder, full.names=TRUE, recursive=FALSE, pattern='.bam$')

### Binning ###
if (createCompositeFile) {
    conf[['numCPU']] <- 1 #always use only one CPU for composite file analysis
    savename <- file.path(datapath, paste0('compositeFile', '.RData'))
    fragments <- createCompositeFile(file.list=files, chromosomes=conf[['chromosomes']], pairedEndReads=conf[['pairedEndReads']], pair2frgm=conf[['pair2frgm']], min.mapq=conf[['min.mapq']], WC.cutoff=conf[['WC.cutoff']])    
    ## Find breakpoints
    if (!file.exists(savename)) {  
        breakpoints <- runBreakpointr(bamfile=fragments, ID='compositeFile', pairedEndReads=conf[['pairedEndReads']], chromosomes=conf[['chromosomes']], windowsize=conf[['windowsize']], binMethod=conf[['binMethod']], trim=conf[['trim']], peakTh=conf[['peakTh']], zlim=conf[['zlim']], bg=conf[['bg']], minReads=conf[['minReads']], maskRegions=conf[['maskRegions']], estimate.back=conf[['estimate.back']])
        save(breakpoints, file=savename)
    } else {
        breakpoints <- get(load(savename))
    }
  
    ## Write BED file
    savename.breakpoints <- file.path(browserpath, paste0('compositeFile', '_breakPoints.bed.gz'))
    if (!file.exists(savename.breakpoints)) {
        writeBedFile(index='compositeFile', outputDirectory=browserpath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks)
    }  
  

} else if (numCPU > 1) {
    ## Parallelization ##
    message("Using ",conf[['numCPU']]," CPUs")
    cl <- parallel::makeCluster(conf[['numCPU']])
    doParallel::registerDoParallel(cl)
  
    message("Finding breakpoints ...", appendLF=FALSE); ptm <- proc.time()
    temp <- foreach (file = files, .packages=c('breakpointR')) %dopar% {
        savename <- file.path(datapath,paste0(basename(file),'.RData'))
        ## Find breakpoints
        if (!file.exists(savename)) {
            tC <- tryCatch({
                breakpoints <- runBreakpointr(bamfile=file, ID=basename(file), pairedEndReads=conf[['pairedEndReads']], pair2frgm=conf[['pair2frgm']], min.mapq=conf[['min.mapq']], chromosomes=conf[['chromosomes']], windowsize=conf[['windowsize']], binMethod=conf[['binMethod']], trim=conf[['trim']], peakTh=conf[['peakTh']], zlim=conf[['zlim']], bg=conf[['bg']], minReads=conf[['minReads']], maskRegions=conf[['maskRegions']], estimate.back=conf[['estimate.back']])
            }, error = function(err) {
                stop(file,'\n',err)
            })
            save(breakpoints, file=savename)
      
        } else {
            breakpoints <- get(load(savename))
        }
        ## Write BED file
        savename.breakpoints <- file.path(browserpath,paste0(basename(file), '_breakPoints.bed.gz'))
        if (!file.exists(savename.breakpoints)) {
            writeBedFile(index=basename(file), outputDirectory=browserpath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks)
        }
    
    }
    parallel::stopCluster(cl)
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

} else {
    conf[['numCPU']] <- 1 #if to use only one CPU or CPU argument not defined
    temp <- foreach (file = files, .packages=c('breakpointR')) %do% {
        savename <- file.path(datapath,paste0(basename(file),'.RData'))
        ## Find breakpoints
        if (!file.exists(savename)) {
            breakpoints <- runBreakpointr(bamfile=file, ID=basename(file), pairedEndReads=conf[['pairedEndReads']], pair2frgm=conf[['pair2frgm']], min.mapq=conf[['min.mapq']], chromosomes=conf[['chromosomes']], windowsize=conf[['windowsize']], binMethod=conf[['binMethod']], trim=conf[['trim']], peakTh=conf[['peakTh']], zlim=conf[['zlim']], bg=conf[['bg']], minReads=conf[['minReads']], maskRegions=conf[['maskRegions']], estimate.back=conf[['estimate.back']])
            save(breakpoints, file=savename)
        } else {
            breakpoints <- get(load(savename))
        }
        ## Write BED file
        savename.breakpoints <- file.path(browserpath,paste0(basename(file), '_breakPoints.bed.gz'))
        if (!file.exists(savename.breakpoints)) {
            writeBedFile(index=basename(file), outputDirectory=browserpath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks)
        }
    
    }
}

## Write masked regions to BED file
if (!is.null(maskRegions)) {
    writeBedFile(index='MaskedRegions', outputDirectory=browserpath, maskedRegions=maskRegions)
}

## Compile all breaks using disjoin function
if (createCompositeFile==FALSE) {
    files <- list.files(datapath, pattern=".RData$", full.names=TRUE)
  
    breaks.all.files <- GenomicRanges::GRangesList()
    for (file in files) {
        breakpoints <- get(load(file))$breaks
        if (length(breakpoints)) {
            suppressWarnings( breaks.all.files[[file]] <- breakpoints ) #TODO check if this can be done without warnings
        }  
    }
    
    breaks <- unlist(breaks.all.files, use.names=FALSE)
    ranges <- disjoin(breaks) # redefine ranges in df
    hits <- countOverlaps(ranges, breaks) # counts number of breaks overlapping at each range
    mcols(ranges)$hits <- hits # appends hits as a metacolumn in ranges
  
    ## write a bedgraph of data (overlapping breaks)
    writeBedFile(index='BreakpointSummary', outputDirectory=browserpath, breaksGraph=ranges)
}

## Search for SCE hotspots
if (callHotSpots) {
    hotspots <- hotspotter(gr.list=breaks.all.files, bw = 1000000)
    savename <- file.path(datapath,'SCE_HotSpots.RData')
    if (length(hotspots)) {
        save(hotspots, file=savename)
        writeBedFile(index='HotSpots', outputDirectory=browserpath, hotspots=hotspots)
    }
}

}
