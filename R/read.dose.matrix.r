#' @export
read.dose.matrix <- function(basename, memory = TRUE) {

  dosf <- path.expand(paste0(basename, ".dosf"))
  bim <- path.expand(paste0(basename, ".bim"))
  fam <- path.expand(paste0(basename, ".fam"))

  if(!file.exists(dosf)) { # peut-être on a donné le .dosf pour basename
    if(length(grep("\\.dosf$", basename)) > 0) {
      basename <- sub("\\.dosf$", "", basename)
      dosf <- path.expand(paste0(basename, ".dosf"))
      bim <- path.expand(paste0(basename, ".bim"))
      fam <- path.expand(paste0(basename, ".fam"))
    }
  }

  if(!file.exists(dosf)) stop("file ", dosf, " not found")
  if(!file.exists(bim)) stop("file ", bim, " not found")
  if(!file.exists(fam)) stop("file ", fam, " not found")

  if(memory) {
    ptr <- readDosageFileMemory_(dosf, bim, fam)
  } else {
    ptr <- readDosageFileDisk_(dosf , bim, fam)
  }
  
  new("dose.matrix", ptr = ptr, file = dosf, type = ifelse(memory, "memory", "disk"))
}