#' Converting custom datastruct to R data.frame
#' 
#' @name as.data.frame.data.struct
#'
#' @exportS3Method as.data.frame data.struct
as.data.frame.data.struct <- function(x, row.names = NULL, optional = FALSE, ..., cols) {
  if(missing(cols)) {
    cols <- 1:ncolDataStruct(x)
  } else if(is.character(cols)) {
    ds_names <- colNamesDataStruct(x)
    cols <- match(cols, ds_names)
    cols <- cols[ !is.na(cols) ]
  } else {
    cols <- as.integer(cols)
    cols <- cols[ cols > 0L & cols <= ncolDataStruct(x) ]
  }
  DataStructToSEXP(x, cols - 1L)
}


