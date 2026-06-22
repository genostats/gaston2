# fonction qui vérifie si "dist" est set
check.dist <- function(x) {
  md <- mean(snp.stats(x)$dist == 0)
  if(is.na(md)) stop("NA values in 'dist'")
  if(md == 1) stop("'dist' doesn't seem set, use 'set.dist'")
  if(md > 0.1) warning("Suprising amount of 0s in 'dist', might be incorrectly set?")
}
