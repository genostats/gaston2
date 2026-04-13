#' @name data.struct
#' @rdname data.struct
#' @title Data Structure
#' 
#' @description A container for columns with dtaas
#' @exportClass data.struct
setClass("data.struct", slots = c(ptr = "externalptr"))

setMethod("show", "data.struct", 
  function(object) {
    if(isnullptr(object@ptr)) {
      cat("A data.struct with a broken external ptr\n")
    } else {
      cat("A data.struct with", ncolDataStruct(object@ptr), "columns\n")
      # list all Columns present
      vecNames <- colNamesDataStruct(object@ptr)
      for (name in vecNames) {
        cat("$ ", name, ": ")
        cat(tolower(getcolTypeDataStruct(object@ptr, name)), " ")
        # also cat beginning of columns... with a threshold of end of line 
        # je voulais aller que aux 3/4 mais peut etre que c'est nul en fait
        maxcharleft <- (getOption("width") * 0.75) - nchar(name) - 6 #approx de nchar(type) + ": " ...
        # brouillon pour un getter
        # for (i in )
        # val <- getcolValueDatastruct(object@ptr, name, i)
        # if (nchar(val) + 3 < maxcharleft) cat(val, " ")
        # else cat("...") break en gros
        cat("\n") 
      }
    }
  } 
)

# # créer un getter de column par nom ou par indice POUR LE DERNIER PAS OUBLIER LE -1 !!!
# setMethod("[", c(x = "data.struct", column_name),
#   function(x, column_name) {
#   }
# )


### THIS IS THE STR METHOD FOR DATAFRAME §§
# getS3method("str", "data.frame")
#function (object, ...) 
# {
#     if (!is.data.frame(object)) {
#         warning("str.data.frame() called with non-data.frame -- coercing to one.")
#         object <- data.frame(object)
#     }
#     cl <- oldClass(object)
#     cl <- cl[cl != "data.frame"]
#     if (0 < length(cl)) 
#         cat("Classes", paste(sQuote(cl), collapse = ", "), "and ")
#     cat("'data.frame':\t", nrow(object), " obs. of  ", (p <- length(object)), 
#         " variable", if (p != 1) 
#             "s", if (p > 0) 
#             ":", "\n", sep = "")
#     if (length(l <- list(...)) && any("give.length" == names(l))) 
#         invisible(NextMethod("str", ...))
#     else invisible(NextMethod("str", give.length = structure(FALSE, 
#         from = "data.frame"), ...))
# }
# <bytecode: 0x571b143ce640>
# <environment: namespace:utils>

# TO THINK : est-ce que je veux un lien avec houba ? Parce qu'a priori non ms 
# peut etre si les stats sont très grosses ? Utiliser le threshold de houba ?

# TO THINK : de quoi est-ce que j'ai besoin
# le nombre de Columns
# leurs noms
# leurs types
# leurs premiers éléments? il va falloir des templates, pas sûre que Rcpp aime
# 
# TODO add une fonction interne qui permets d'ajouter une Column simplement
# TODO faire un lien avec DataFrameToDataStruct que j'ai
# TODO ajouter un attribut avec un vecteur de Columns
