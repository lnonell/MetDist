find.dml.n <- function(x,alpha=0.05) length(x[!is.na(x) & x<alpha])
