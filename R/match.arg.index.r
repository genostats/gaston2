
# d'après le code match.arg
# fait comme match.arg mais au lieu de renvoyer la valeur de l'argument, renvoie son indice 
# pour le moment utilisé dans set.hwe 
match.arg.index <- function (arg, c.index = TRUE) {
  formal.args <- formals(sys.function(sysP <- sys.parent()))
  choices <- eval(formal.args[[as.character(substitute(arg))]], envir = sys.frame(sysP))
  if (is.null(arg)) 
    return (1L - c.index)
  else if (!is.character(arg)) 
    stop("'arg' must be NULL or a character vector")

  if (identical(arg, choices)) 
    return (1L - c.index)

  if (length(arg) != 1L) 
     stop(gettextf("'%s' must be of length 1", "arg"), domain = NA)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L)) 
     stop(sprintf("'arg' should be one of %s"), paste(choices, collapse = ", "))

  i <- i[i > 0L]
  if(length(i) > 1) stop("there is more than one match in 'match.arg'")
  i - c.index
}

