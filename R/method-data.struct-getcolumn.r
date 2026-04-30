#' Getters for a column in a datastruct, by name

#' @export
setMethod("$", "data.struct",
# should be able to fetch one column by name and send it in a new data.struct
    function(x, name) {
        if (isnullptr(x@ptr)) stop("This data.struct has a broken external ptr")
        if (!is.character(name)) stop("Please give a name")

        ds_names = colNamesDataStruct(x@ptr)
        if (!(name %in% ds_names))
            stop("column ", name, " does not exist")
        # do I want to stop ou send back a NULL comme data.frame ? 
        else
            new("data.struct", ptr = extractcolDataStruct_(x@ptr, name), matrixptr = x@matrixptr)
    }
)

#' @export
setMethod("[", "data.struct",
# should be able to fetch a list of columns, and send them in a new data.struct
    function(x, i, j, ..., drop = FALSE) { # je peux pas vraiment faire un drop 
    # mais sinon R rale
        if (isnullptr(x@ptr)) stop("This data.struct has a broken external ptr")
        if (!missing(j) || drop) stop("This is not a real extract yet, just give a list of column names")
        # i doit avoir un vecteur de nom ou un nom unique
        if (!is.character(i)) stop("Please give column name(s)")

        ds_names = colNamesDataStruct(x@ptr)
        for (name in i) {
        if (!(name %in% ds_names))
            stop("column ", name, " does not exist")
        }
        new("data.struct", ptr = extractcolDataStruct_(x@ptr, i), matrixptr = x@matrixptr)
    }
)


# This one will only handle one at a time
#' @export
setMethod("[[", "data.struct",
    function(x, i, j, ...) { # R is making a mess, I is what gets teh name
        if (isnullptr(x@ptr)) stop("This data.struct has a broken external ptr")

        if (length(i) != 1) stop("Please only give one column")
        if (!is.character(i)) stop("Please give a column name")
        ds_names = colNamesDataStruct(x@ptr)
        if (!(i %in% ds_names))
            stop("column ", i, " does not exist")
        # do I want to stop ou send back a NULL comme data.frame ? 
        else extractSEXPDataStruct(x@ptr, i)
    }
)