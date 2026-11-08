# Short check of the 0.80.0 release fixes: bernoulli sig4 55 of the nlstart
# reviewer's fresh2 rows. Is the default fit at a maximum? Its gradient,
# the Newton step's gain from the outer Hessian at the stop, the Hessian's
# eigenvalues, and the best of the review's 16 random starts plus a start
# at the default fit's own estimates.
#
#   Rscript dev/review-rel-sig4.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(REV_ARM = strsplit(lib, ";")[[1]][1])
src <- readLines("dev/review-nlstart-fresh2.R")
eval(parse(text = src[1:(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
fam <- "bernoulli"; body <- "sig4"; seed <- 55
d <- gen(fam, body, seed * 1000 + match(body, names(B)))
fm <- form(fam, body)
r <- tt(frm(fm, family = F[[fam]], data = d))
cat("default ll", sprintf("%.6f", llv(r)), "conv", r$opt$convergence,
    r$opt$message, "warnings", length(attr(r, "w")), "\n")
cat("  messages:", paste(attr(r, "mc"), collapse = ", "), "\n")
p <- r$opt$par
g <- r$obj$gr(p)
H <- stats::optimHess(p, r$obj$fn, r$obj$gr)
ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
step <- tryCatch(solve(H, as.numeric(g)), error = function(e) rep(NA, length(p)))
cat("  max |grad|", format(max(abs(g)), digits = 3),
    " Newton gain", format(0.5 * sum(as.numeric(g) * step), digits = 3),
    " Hessian eigen min/max", format(min(ev), digits = 3), "/",
    format(max(ev), digits = 3), "\n")
cat("  estimates", paste(sprintf("%s=%.5g", names(r$estimates$beta),
                                 r$estimates$beta), collapse = "; "), "\n")
nm <- names(r$estimates$beta)
np <- sub("_.*$", "", nm)
sc <- B[[body]]$sc[np]
sc[is.na(sc)] <- 0.5
set.seed(seed * 13 + 7)
starts <- lapply(1:16, function(i) {
  setNames(rnorm(length(nm)) * sc * if (i > 8) 3 else 1, nm)
})
if (!is.null(B[[body]]$tr)) {
  tr <- setNames(rep(0, length(nm)), nm)
  ic <- grepl("_\\(Intercept\\)$", nm)
  tr[ic] <- B[[body]]$tr[np[ic]]
  starts <- c(starts, list(tr))
}
starts <- c(starts, list(setNames(unname(r$estimates$beta), nm)))
best <- -Inf; bf <- NULL
for (i in seq_along(starts)) {
  f <- tt(frm(fm, family = F[[fam]], data = d, se = FALSE,
              start = list(beta = starts[[i]])))
  v <- llv(f)
  if (!is.na(v) && v > best) { best <- v; bf <- f; bi <- i }
}
cat("best of", length(starts), "starts", sprintf("%.6f", best), "(start",
    bi, ") short", format(best - llv(r), digits = 3), "\n")
cat("  best estimates", paste(sprintf("%s=%.5g", names(bf$estimates$beta),
                                      bf$estimates$beta), collapse = "; "), "\n")
cat("  max |estimate difference|",
    format(max(abs(bf$estimates$beta - r$estimates$beta)), digits = 3), "\n")
cat("SIG4 DONE\n")
