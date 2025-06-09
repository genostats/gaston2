#' @export
ind.stats <- function(x) {
if (class(x) == "snp.matrix")
    getIndStats(x@ptr)
else getIndStatsDosage(x@ptr)
}