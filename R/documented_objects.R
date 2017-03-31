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


#' BreakPoint objects for demonstration purposes
#'
#' A list of \code{\link{BreakPoint}} objects for demonstration purposes in examples of package \pkg{\link{breakpointR}}.
#'
#' @docType data
#' @name example_BreakPoints
#' @format A list with \code{\link{BreakPoint}} objects.
#' @examples
#'data(example_BreakPoints)
#'length(example_BreakPoints)
#'## Have a look at the first object
#'example_BreakPoints[[1]]
NULL
