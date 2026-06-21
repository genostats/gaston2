
# to eval an expression in an environment generated from a data struct
# used in select.inds and select.snps

eval_ds <- function(expr, ds, enclos = parent.frame()) {
  # get a data frame with only the relevant columns
  df <- as.data.frame(ds, cols = all.names(expr))
  eval(expr, df, enclos)
}
