# Merge check 7 (parity2 x refits): influence() on an s(exp(x)) fit runs
# on the merged build and its deletions keep the fit's smooth basis. The
# check of tests/testthat/test-refit-plumbing.R ("an influence() deletion
# keeps the fit's smooth basis") with the smooth reading exp(x), and the
# same fit with exp(x) as a data column.
#   Rscript dev/rel083-kt/refits/mc7-influence-sexp.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")),
    as.character(packageVersion("frmtmb")), "\n")
set.seed(8)
d <- data.frame(x = stats::runif(60, 0, 3), z = stats::rnorm(60),
                blk = factor(rep(1:5, 12)))
d$y <- sin(d$x) + 0.3 * d$z + stats::rnorm(60, 0, 0.3)
d$ex <- exp(d$x)
fe <- frm(y ~ z + s(exp(x), k = 5), data = d)
fc <- frm(y ~ z + s(ex, k = 5), data = d)
cat("fit s(exp(x)) logLik", format(as.numeric(logLik(fe)), digits = 15),
    "| s(ex) logLik", format(as.numeric(logLik(fc)), digits = 15), "\n")
cat("fixef names (s(exp(x))):", rownames(fixef(fe)), "\n")
ie <- tryCatch(influence(fe, groups = "blk", data = d), error = function(e) e)
if (inherits(ie, "error")) {
  cat("influence(s(exp(x))) ERROR:", conditionMessage(ie), "\n")
  quit(status = 1)
}
ic <- influence(fc, groups = "blk", data = d)
cat("influence(s(exp(x))) ran; fixed dims", dim(ie$fixed), "\n")
fxe <- grep("fx1$", colnames(ie$fixed), value = TRUE)
fxc <- grep("fx1$", colnames(ic$fixed), value = TRUE)
cat("null-space columns:", fxe, "|", fxc, "\n")
a <- ie$fixed[, c("z", fxe)]
b <- ic$fixed[, c("z", fxc)]
cat("max |influence s(exp(x)) - influence s(ex)| over 5 deletions:",
    format(max(abs(unname(a) - unname(b)))), "\n")
# the deletion of block 1 with the fit's basis written out
si <- fe$frame$linpreds[[1L]]$smooths[[1L]]
keep <- d[d$blk != "1", ]
M <- frmtmb:::smooth_basis_at(si, keep)
nr <- sum(si$nr)
for (k in seq_len(nr)) keep[[paste0("r", k)]] <- M[, k]
keep$fx <- M[, nr + 1L]
keep$id <- factor(1)
fo <- stats::as.formula(paste0("y ~ z + fx + homdiag(0 + ",
                               paste0("r", seq_len(nr), collapse = " + "),
                               " | id)"))
ref <- fixef(frm(fo, data = keep))[, 1]
got <- ie$fixed["1", c("z", fxe)]
moved <- abs(got - ie$fixed_full[c("z", fxe)])
r <- max(abs(got - ref[c("z", "fx")]) / moved)
cat("block 1: influence", format(got, digits = 10), "| fit's basis written out",
    format(ref[c("z", "fx")], digits = 10), "| moved", format(moved, digits = 4),
    "\nratio (test bound 1e-3):", format(r), if (r < 1e-3) "PASS" else "FAIL",
    "\n")
cat("DONE\n")
