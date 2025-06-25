
# rbind
rbind2_dose_memory <- function(x, y, ...) {
  print("Calling bindIndstoDosagematrixMemory_")
  new_ptr <- bindIndstoDosagematrixMemory_(x@ptr, y@ptr)
  new("dose.matrix", ptr = new_ptr, file = NULL, type = "memory")
}

#' @export
rbind2_dose_disk <- function(x, y, file, ...) {
  print("Calling bindIndstoSNPmatrixDisk_")
  # add a check that file correct ?or let c++ do it ?
  new_ptr <- bindIndstoDosagematrixDisk_(x@ptr, y@ptr, file)
  new("dose.matrix", ptr = new_ptr, file = file, type = "disk")
}


#' @export
setMethod("rbind2", c(x = "dose.matrix", y = "dose.matrix"),
  function(x, y, ..., file = NULL) {
    print("using my rbinf dose.mat")
    if(...length() > 0)
      type <- match.arg(..1, c("disk", "memory"))
    else
      type <- x@type
      # sinon je prends le type de la première

    #get file from ... if it was passed ? 
    #else do an error
    # also modify the path with all prefix from r
    if (type == "disk") {
      if(type == "disk") {
        if(...length() > 1) 
          file <- ..2
        else
          file <- tempfile("gaston2")
      } else {
        file <- NULL
      }
    }
    if(type == "disk")
      rbind2_dose_disk(x, y, file = file, ...)
    else
      rbind2_dose_memory(x, y, ...)
  }
)

# for cbind,
# because the new matrix only exists as an "empty" object
# and all its content is simply pointed to 
# and stored in other matrix effectively handling it
# the type is a bit weird. 
# Could be defined as in memory without having any SNPs loaded
# but seems the most logical thing to do ?
# or else add a way to keep track of all files
# that it's refering to but seems a bit tedious...

#' @export
setMethod("cbind2", c(x = "dose.matrix", y = "dose.matrix"), 
  function(x, y, ...) {
  print("Calling cbind_Dosagematrix")
  new_ptr <- cbind_Dosagematrix(x@ptr, y@ptr)
  new("dose.matrix", ptr = new_ptr, type = "memory")
  }
)