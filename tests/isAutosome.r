require(gaston2)
gaston2:::isAutosome_( 1:30 )
gaston2:::setGastonOptions( list(autosomes = 1:12, x = 13, y = 15, mt = 14))
gaston2:::isAutosome_( 1:30 )

