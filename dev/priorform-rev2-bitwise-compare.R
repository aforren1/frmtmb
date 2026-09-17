# Reviewer, lane wt-priorform: compare the ref and lane digests from
# dev/priorform-rev-bitwise.R with identical().
#   Rscript dev/priorform-rev-bitwise-compare.R A|B
part <- commandArgs(trailingOnly = TRUE)[1]
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform/dev"
r <- readRDS(file.path(root, sprintf("priorform-rev-bitwise-%s-ref.rds", part)))
l <- readRDS(file.path(root, sprintf("priorform-rev-bitwise-%s-laner1.rds", part)))
cmp <- function(a, b) {
  if (is.character(a) || is.character(b)) {
    return(if (identical(a, b)) "same error" else
      paste("DIFFERENT:", substr(paste(a), 1, 60), "|", substr(paste(b), 1, 60)))
  }
  f <- c("ll", "par", "obj", "cov", "est")
  bad <- f[!vapply(f, function(k) identical(a[[k]], b[[k]]), TRUE)]
  if (length(bad)) paste("DIFFER in", paste(bad, collapse = ",")) else "identical"
}
n_fit <- 0L; n_same <- 0L
if (part == "A") {
  for (nm in union(names(r), names(l))) {
    v <- cmp(r[[nm]], l[[nm]])
    n_fit <- n_fit + 1L; n_same <- n_same + (v == "identical")
    cat(sprintf("%-26s %s\n", nm, v))
  }
} else {
  for (nm in union(names(r), names(l))) {
    a <- r[[nm]]; b <- l[[nm]]
    if (!identical(a$err, b$err)) {
      cat(sprintf("%-50s ERR ref='%s' lane='%s'\n", nm, a$err, b$err))
    }
    fn <- union(names(a$fits), names(b$fits))
    for (f in fn) {
      v <- if (is.null(a$fits[[f]]) || is.null(b$fits[[f]])) "MISSING in one" else
        cmp(a$fits[[f]], b$fits[[f]])
      n_fit <- n_fit + 1L; n_same <- n_same + (v == "identical")
      if (v != "identical") cat(sprintf("%-50s %-12s %s\n", nm, f, v))
    }
  }
}
cat(sprintf("part %s: %d fits, %d identical in ll, par, objective, cov.fixed and estimates\n",
            part, n_fit, n_same))
