require(gaston2)

filename <- system.file("extdata", "LCT.bed", package="gaston2")

# -------- getting a data.struct for stats
a <- read.snp.matrix(filename)
ds_snp <- snp.stats(a, TRUE)
ds_ind <- ind.stats(a, TRUE)

#--------converting to dataframes and compare with gaston
df_snp <- as.data.frame(ds_snp)
df_ind <- as.data.frame(ds_ind)

stopifnot(all(df_snp == snp.stats(a)))
stopifnot(all(df_ind == ind.stats(a)))

# ------ compares str(df) and show(ds)
df_str <- capture.output(str(snp.stats(a)))
ds_show <- capture.output(show(ds_snp))

data("LCT", package = "gaston")
x <- gaston::as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
# TODO ADD gatson compare
stopifnot(all(table(x@snps["A2"])==table(df_snp["A2"])))

#--------- checking c° from dataframe
test_ind <- data.struct(df_ind)

#--------- checking le addcolumn
#adanced R
df <- data.frame(x = 1:3, y = c("a", "b", "c"))
to_add <- data.frame(z = 4:6, t = c("D", "E", "F"))

ds_add <- data.struct(df)

# checking the addcol with a full dataframe
addcolumn(ds_add, to_add)

# checking the addcol with a vector
v <- c("THIS", "IS", "A", "NEW", "COLUMN")
addcolumn(ds_add, "COL", vec=v)
#will now modify in place
addcolumn(ds_add, "COL", vec=c("SORRY", "THIS", "IS","A", "NEW", "COLUMN" ))

ds_add

# n'identifie que le dataframe et va pas chopper le nom 
addcolumn(ds_add, "NEWCOLUMN", x@ped["N0"])



#-------- checking it does modify the snp.stats / ind.stats
# This does work but will give a warning and downgrade to list
addcolumn(ds_snp, "COL", vec=v)

## THIS IS COMMENTED BCOS tryCatch never frees the DataFrame => leak
#warning msg : Column sizes are not equal in DataFrame::push_back, object degrading to List
# tryCatch({snp.stats(a)}, 
#          warning = function(w) {
#            return(TRUE)
#          })
warn <- snp.stats(a)

addcolumn(ds_snp, "COL", vec=rep("NEW!", 607))
#no warnings now
sil <- snp.stats(a)




# ---------- getcolumn
d1 <- ds_snp$COL
#should fail
tryCatch({ds_snp$WRONG}, 
  error = function(e) {
  return(TRUE)
})

# #All of this should work
d2 <- ds_snp['COL']
d3 <- ds_snp[['COL']]
d4 <- ds_snp[c('chr', 'COL')]

#should fail
tryCatch({ds_snp[c('chr', 'COL', 'WRONG')]}, 
         error = function(e) {
           return(TRUE)
         })

# #---------- Booléen
addcolumn(ds_add, "TEST_BOOL", vec=c(TRUE, FALSE, FALSE))

df_bool <- data.frame(x = c(TRUE, FALSE, TRUE, TRUE, FALSE), y = c(FALSE, FALSE, FALSE, FALSE, FALSE))
data.struct(df_bool)

