#' @export
snp.stats <- function(x, use_datastruct = FALSE) {
    if (use_datastruct) {
        if (class(x) == "snp.matrix") {
            getSNPStats_DataStruct(x@ptr)
        } else {
            getSNPStatsDosage_DataStruct(x@ptr)
        }
    } else {
        if (class(x) == "snp.matrix") {
            getSNPStats(x@ptr)
        } else {
            getSNPStatsDosage(x@ptr)
        }
    }
}
