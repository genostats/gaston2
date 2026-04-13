#' @export
ind.stats <- function(x, use_datastruct = FALSE) {
    if (use_datastruct) {
        if (class(x) == "snp.matrix") {
            getIndStats_DataStruct(x@ptr)
        } else {
            getIndStatsDosage_DataStruct(x@ptr)
        }
    } else {
        if (class(x) == "snp.matrix") {
            getIndStats(x@ptr)
        } else {
            getIndStatsDosage(x@ptr)
        }
    }
}
