
#' @export
read.snp.matrix <- function(basename, memory = TRUE) {

  bed <- path.expand(paste0(basename, ".bed"))
  bim <- path.expand(paste0(basename, ".bim"))
  fam <- path.expand(paste0(basename, ".fam"))

  if(!file.exists(bed)) { # peut-être on a donné le .bed pour basename
    if(length(grep("\\.bed$", basename)) > 0) {
      basename <- sub("\\.bed$", "", basename)
      bed <- path.expand(paste0(basename, ".bed"))
      bim <- path.expand(paste0(basename, ".bim"))
      fam <- path.expand(paste0(basename, ".fam"))
    }
  }

  if(!file.exists(bed)) stop("file ", bed, " not found")
  if(!file.exists(bim)) stop("file ", bim, " not found")
  if(!file.exists(fam)) stop("file ", fam, " not found")

  if(memory) {
    ptr <- readBedFileMemory_(bed, bim, fam)
  } else {
    ptr <- readBedFileDisk_(bed, bim, fam)
  }
  
  new("snp.matrix", ptr = ptr, file = bed, type = ifelse(memory, "memory", "disk"))
}
