res.chr <-function(res=rnb.meth2comp.models.adj.p.s){
  annot <- strsplit(rownames(res),split="_")
  chr <-unlist(lapply(annot, function(l) l[[1]]))
  return(as.data.frame(table(chr)))
}