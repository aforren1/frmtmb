# The measurements the findings quote that the reproduction script does
# not carry: what storing `fit$data` costs, which refusal an
# autocorrelation term reaches first, and whether the D5 item is a live
# defect on the base build.
#
#   Rscript dev/adefects-evidence.R > dev/adefects-log/evidence.txt 2>&1
#
# ADEFECTS_LIB picks the build (the lane's by default); the reference
# build and the user library follow it, so frmtmb.latent and the other
# extensions come from the round's reference build either way.
lib <- Sys.getenv("ADEFECTS_LIB", "C:/Users/adf44/source/r/adefects-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", lib, "\n")
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.latent",
    format(packageVersion("frmtmb.latent")), "brms",
    format(packageVersion("brms")), "\n")

addr <- function(x) {
  # the SEXP address, which is the only way to tell a SHARED object from
  # an equal copy; utils::object.size() counts a shared object once per
  # path it is reached by, so it cannot answer this
  out <- utils::capture.output(.Internal(inspect(x)))[1L]
  # the first whitespace-separated token; the rest of the line carries a
  # REFERENCE COUNT, which differs between two reads of the same object
  # and made a first spelling of this report FALSE on a shared SEXP
  strsplit(trimws(out), "[[:space:]]+")[[1L]][1L]
}

cat("\n== D4 what the `data` element costs\n")
set.seed(20260917)
A <- diag(6)
A[A == 0] <- 0.3
dimnames(A) <- list(1:6, 1:6)
dd <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
dd$y <- dd$x + rnorm(6)[dd$g] + rnorm(30)
fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))
a1 <- addr(fa$data)
a2 <- addr(fa$frame[["data_frame"]])
cat("  address of fit$data              ", a1, "\n")
cat("  address of fit$frame$data_frame  ", a2, "\n")
cat("  SAME OBJECT:", identical(a1, a2), "\n")
# the inverse case, so that TRUE above is not true of any two equal
# frames: an identical COPY has a different address
cat("  control, an equal copy of the same frame:",
    identical(a1, addr(data.frame(fa$data))), "\n")
mf <- fa$frame[["data_frame"]]
cat("  control, object.size double counts a shared object:\n")
cat("    object.size(list(mf))      ",
    as.numeric(utils::object.size(list(mf))), "bytes\n")
cat("    object.size(list(mf, mf))  ",
    as.numeric(utils::object.size(list(mf, mf))), "bytes\n")
cat("    object.size(mf)            ",
    as.numeric(utils::object.size(mf)), "bytes\n")
cat("  so the +", as.numeric(utils::object.size(mf)),
    " bytes object.size() reports on the fit is the double count, ",
    "not an allocation\n", sep = "")
if (requireNamespace("lobstr", quietly = TRUE)) {
  cat("  lobstr::obj_size(fit), which does NOT double count:",
      as.numeric(lobstr::obj_size(fa)), "\n")
} else {
  cat("  lobstr is not installed; the address identity is the evidence\n")
}
cat("  names(fit):", paste(names(fa), collapse = " "), "\n")
cat("  brms keeps the model frame there too: dim(ex1$data) =",
    paste(dim(brms:::rename_pars(
      get("brmsfit_example1", envir = asNamespace("brms")))$data),
      collapse = " x "), "\n")

cat("\n== D1/D2 which refusal an autocorrelation term reaches first\n")
# brm:106 and brm:108 are written WITHOUT cov = TRUE, so the ledger row
# turns on the grammar check running before the cov check
set.seed(20260917)
d <- data.frame(y = rnorm(24), x = rnorm(24), g1 = rep(1:2, 12),
                g2 = rep(1:2, each = 12), g = rep(1:4, 6))
one <- function(e) {
  tryCatch({
    frm(e, d)
    "NO ERROR"
  }, error = function(err) substr(conditionMessage(err), 1, 70))
}
cat("  ar(x + y, g)      :", one(y ~ ar(x + y, g)), "\n")
cat("  ar(gr = g1/g2)    :", one(y ~ ar(gr = g1/g2)), "\n")
cat("  ar(x, g)          :", one(y ~ ar(x, g)), "\n")

cat("\n== D5 frmtmb.latent's hmm/lca readers on an atomic argument\n")
cat("  every exported reader, given 1, a string, NULL and an empty list\n")
ns <- asNamespace("frmtmb.latent")
readers <- grep("^(hmm|lca)_", getNamespaceExports(ns), value = TRUE)
readers <- readers[vapply(readers, function(n) {
  f <- get(n, envir = ns)
  is.function(f) && length(formals(f)) >= 1L
}, NA)]
bad <- character(0)
for (nm in sort(readers)) {
  f <- get(nm, envir = ns)
  for (arg in list(1, "a", NULL, list())) {
    cls <- tryCatch({
      f(arg)
      "NO ERROR"
    }, error = function(e) paste(class(e), collapse = ","))
    if (!grepl("frmtmb_error", cls)) {
      bad <- c(bad, sprintf("%s(%s) -> %s", nm,
                            if (is.null(arg)) "NULL" else
                              paste(deparse(arg), collapse = ""), cls))
    }
  }
}
cat("  readers scanned:", length(readers), "\n")
cat("  arguments that do NOT get a classed frmtmb refusal:",
    length(bad), "\n")
for (b in bad) cat("    ", b, "\n", sep = "")

cat("\nDONE\n")
