# Patch for SeuratDisk by shinichinamba :
#  https://github.com/mojaveazure/seurat-disk/issues/147#issuecomment-1721031207
###### WORKAROUND ######
# assigning the previous version of the `[[` function
# for the Assay class to the SeuratDisk package environment

"[[.Assay" <- function(x, i, ..., drop = FALSE) {
    if (missing(x = i)) {
        i <- colnames(x = slot(object = x, name = "meta.features"))
    }
    data.return <- slot(object = x, name = "meta.features")[, i, drop = FALSE, ...]
    if (drop) {
        data.return <- unlist(x = data.return, use.names = FALSE)
        names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
    }
    return(data.return)
}
environment(`[[.Assay`) <- asNamespace("SeuratObject")
rlang::env_unlock(asNamespace("SeuratDisk"))
assign("[[.Assay", `[[.Assay`, asNamespace("SeuratDisk"))
lockEnvironment(asNamespace("SeuratDisk"), bindings = TRUE)
rm(`[[.Assay`)
###### WORKAROUND ######
