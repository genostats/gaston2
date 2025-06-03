
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

  nbInds <- R.utils::countLines(fam)
  nbSNPs <- R.utils::countLines(bim)

  if(memory) {
    ptr <- readBedFileMemory_(bed, nbInds, nbSNPs)
  } else {
    ptr <- readBedFileDisk_(bed, nbInds, nbSNPs)
  }

  new("snp.matrix", ptr = ptr, file = bed, type = ifelse(memory, "memory", "disk"))
}
