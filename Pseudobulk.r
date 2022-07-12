getPseudobulk <- function(mat, celltype) {
   mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
     cells <- names(celltype)[celltype==ct]
     pseudobulk <- rowSums(mat[, cells])
     return(pseudobulk)
   }))
   colnames(mat.summary) <- levels(celltype)
   return(mat.summary)
}
