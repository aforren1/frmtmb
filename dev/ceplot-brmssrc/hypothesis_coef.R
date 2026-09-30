hypothesis_coef <- 
function (x, hypothesis, alpha, ...) 
{
    stopifnot(is.array(x), length(dim(x)) == 3)
    levels <- dimnames(x)[[2]]
    coefs <- dimnames(x)[[3]]
    x <- lapply(seq_along(levels), function(l) structure(as.data.frame(x[, 
        l, ]), names = coefs))
    out <- vector("list", length(levels))
    for (l in seq_along(levels)) {
        out[[l]] <- .hypothesis(x[[l]], hypothesis, class = "", 
            alpha = alpha, combine = FALSE, ...)
        for (i in seq_along(out[[l]])) {
            out[[l]][[i]]$summary$Group <- levels[l]
        }
    }
    out <- unlist(out, recursive = FALSE)
    out <- as.list(matrix(out, ncol = length(hypothesis), byrow = TRUE))
    out <- combine_hlist(out, class = "", alpha = alpha)
    out$hypothesis$Group <- factor(out$hypothesis$Group, levels)
    out$hypothesis <- move2start(out$hypothesis, "Group")
    out
}
