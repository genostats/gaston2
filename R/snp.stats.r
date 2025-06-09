#' @export
snp.stats <- function(x) {
if (class(x) == "snp.matrix")
    getSNPStats(x@ptr)
else getSNPStatsDosage(x@ptr)
}