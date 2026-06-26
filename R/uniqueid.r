uniqueid <- function(famid, id) {
  if(!anyDuplicated(id)) return(id)

  # try to build unique ids as famid:id
  r <- paste(famid, id, sep = ":")
  if(anyDuplicated(r)) {
    warning("Duplicated famid:id")
  }

  return(r)
}
