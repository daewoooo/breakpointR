#' BreakPoint object
#'
#' The \code{BreakPoint} object is output of the function \code{\link{runBreakpointr}} and is basically a list with various entries. The class() attribute of this list was set to "BreakPoint". Entries can be accessed with the list operators '[[]]' and '$'.
#'
#' @name BreakPoint
#' @return
#' \item{fragments}{A \code{\link[GenomicRanges]{GRanges}} object with read fragments.}
#' \item{deltas}{A \code{\link[GenomicRanges]{GRanges}} object with deltaWs.}
#' \item{breaks}{A \code{\link[GenomicRanges]{GRanges}} object containing the breakpoint coordinates.}
#' \item{counts}{A \code{\link[GenomicRanges]{GRanges}} object with the regions between breakpoints.}
#' \item{params}{A vector with parameters that were used to obtain the results.}
#'
#' @seealso runBreakpointr
NULL
