# Short check of fix round 2: the release fixes, as exact edits on a source copy, installed
# in rev-predread-lib/release with frmtmb.sample, to show they do what they say and pass
# the tests they reach. Each pattern must match exactly once.
src <- "C:/Users/adf44/source/r/frmtmb-wt-predread"
d <- "C:/Users/adf44/source/r/rev-predread-src/release"
lib <- "C:/Users/adf44/source/r/rev-predread-lib/release"
unlink(d, recursive = TRUE)
dir.create(d, recursive = TRUE)
for (p in c("DESCRIPTION", "NAMESPACE", ".Rbuildignore")) file.copy(file.path(src, p), file.path(d, p))
for (p in c("R", "inst", "data", "man")) file.copy(file.path(src, p), d, recursive = TRUE)
ed <- function(f, old, new) {
  p <- file.path(d, f)
  x <- paste(readLines(p, warn = FALSE), collapse = "\n")
  parts <- strsplit(x, old, fixed = TRUE)[[1L]]
  if (length(parts) != 2L) stop("pattern matched ", length(parts) - 1L, " times: ", old)
  writeLines(paste0(parts[1L], new, parts[2L]), p, useBytes = TRUE)
  cat("edited", f, "\n")
}
# F1: a category on a new me() row reads the draw's own quantile, on 99 equal-weight
# points of the draw; the logit carries the coefficient's share
ed("R/predict.R",
   "  at <- if (!is.null(mm)) matrix(mm$plug, dim(out)[1L], dim(out)[3L])\n",
   paste0(
     "  if (!is.null(mm)) return(prob_band_draw(out, probs, mm, up))\n",
     "  at <- NULL\n"))
ed("R/predict.R",
   "prob_band_quantiles <- function(out, probs, mm = NULL, upper = 1) {",
   paste0(
     "prob_band_draw <- function(out, probs, mm, up) {\n",
     "  q <- setdiff(seq_len(dim(out)[2L]), 1:2)\n",
     "  # a middle category has a hump in the draw, and no interval about one point\n",
     "  # reads that, so each row reads the draw's own quantile, on 99 points of\n",
     "  # equal weight; the logit carries the coefficients' share of the variance\n",
     "  vals <- lapply(stats::qnorm(stats::ppoints(99L)), mm$val)\n",
     "  vr <- matrix(mm$var, dim(out)[1L], dim(out)[3L])\n",
     "  for (k in seq_len(dim(out)[3L])) {\n",
     "    V <- vapply(vals, function(v) as.matrix(v)[, k], numeric(dim(out)[1L]))\n",
     "    V <- matrix(V, dim(out)[1L], length(vals))\n",
     "    # each row is sorted once and read at its two ends, where quantile() sorts per end\n",
     "    V <- matrix(t(apply(V, 1L, sort, na.last = TRUE)), nrow(V), ncol(V))\n",
     "    m <- as.numeric(out[, 1L, k]) / up\n",
     "    se <- as.numeric(out[, 2L, k]) / up\n",
     "    v2 <- pmax(vr[, k] / up^2, .Machine$double.eps)\n",
     "    infl <- sqrt(1 + pmax(se^2 - v2, 0) / v2)\n",
     "    inside <- m > 0 & m < 1\n",
     "    for (i in seq_along(q)) {\n",
     "      h <- (ncol(V) - 1L) * probs[i] + 1\n",
     "      v <- (V[, floor(h)] + (h - floor(h)) * (V[, ceiling(h)] - V[, floor(h)])) / up\n",
     "      r <- stats::plogis(\n",
     "        stats::qlogis(m) + (stats::qlogis(v) - stats::qlogis(m)) * infl)\n",
     "      # a probability at 0 or 1 has no logit, and no spread either\n",
     "      r[which(!inside)] <- m[which(!inside)]\n",
     "      out[, q[i], k] <- r * up\n",
     "    }\n",
     "  }\n",
     "  out\n",
     "}\n\n",
     "prob_band_quantiles <- function(out, probs, mm = NULL, upper = 1) {"))
# F2: a row of NA in a band and in a difference curve
ed("R/curve-cov.R",
   "    if (!length(E@x) || !any(E@x != 0)) return(S)\n",
   "    if (!length(E@x) || !any(E@x != 0, na.rm = TRUE)) return(S)\n")
ed("R/curve-cov.R",
   "    if (!length(E@x) || !any(E@x != 0)) return(NULL)\n",
   "    if (!length(E@x) || !any(E@x != 0, na.rm = TRUE)) return(NULL)\n")
ed("R/curve-cov.R",
   "  if (!any(E != 0)) NULL else unname(E)",
   "  if (!any(E != 0, na.rm = TRUE)) NULL else unname(E)")
ed("R/curve-cov.R",
   "  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)\n  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)",
   paste0(
     "  # a row of NA in the grid has no place in the band; the others do\n",
     "  gone <- !is.finite(diag(S))\n",
     "  if (any(gone)) {\n",
     "    S <- S[!gone, !gone, drop = FALSE]\n",
     "    keep <- keep[!gone]\n",
     "    div <- div[!gone]\n",
     "    m <- nrow(S)\n",
     "  }\n",
     "  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)\n  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = m)"))
# F3: the warning says what is known. The pin in helper-warnings.R stays as it is.
ed("R/interop.R",
   paste0(
     "  frm_warning(\"marginaleffects forms its standard errors from the population-\",\n",
     "              \"level coefficients alone, so the standard error of a level \",\n",
     "              \"from predictions() holds no variance of \",\n",
     "              paste(omits, collapse = \", or of \"),\n",
     "              \". A slope, contrast or test that reads those coefficients \",\n",
     "              \"alone is not affected. frm_linpred(newdata = , se.fit = TRUE) \",\n",
     "              \"and fitted(newdata = ) carry it. vcov = FALSE gives the \",\n",
     "              \"estimates alone, with no standard error\", call. = FALSE,\n"),
   paste0(
     "  frm_warning(\"marginaleffects forms its standard errors from the population-\",\n",
     "              \"level coefficients alone. It holds fixed \",\n",
     "              paste(omits, collapse = \", and \"),\n",
     "              \". The standard error of a level from predictions(), or of an \",\n",
     "              \"average of levels, can be too wide or too narrow, often by a \",\n",
     "              \"factor of 2 or more. On a log or logit link, or with a \",\n",
     "              \"group-level slope, a slope or contrast can be off too. A test \",\n",
     "              \"of a coefficient is not affected. frm_linpred(newdata = , \",\n",
     "              \"se.fit = TRUE) and fitted(newdata = ) give a level's standard \",\n",
     "              \"error. vcov = FALSE gives the estimates alone, with no \",\n",
     "              \"standard error\", call. = FALSE,\n"))
ed("R/interop.R",
   "\"the group-level effects (re_formula = NA drops them, and then it agrees)\"",
   "\"the group-level effects (re_formula = NA drops them, and then a level agrees)\"")
ed("R/compat.R",
   "every row reads as new, and get_vcov() gets no rows. A slope, contrast or test of a coefficient is right. frm_linpred",
   "every row reads as new, and get_vcov() gets no rows. A test of a coefficient is right, and so is a slope or contrast on an identity link. frm_linpred")
# F4: a node on a grid of about 500 rows makes a key past an environment name's 10000 bytes
ed("R/me.R",
   paste0(
     "  memo <- new.env(parent = emptyenv())
",
     "  function(z, d = 0) {
",
     "    key <- paste(c(z, d), collapse = \" \")
",
     "    old <- memo[[key]]
",
     "    if (!is.null(old)) return(old)
"),
   paste0(
     "  # a list: z has one value per row, and an environment's name holds 10000
",
     "  # bytes, which about 500 rows pass. A fit reads a few hundred nodes at most
",
     "  memo <- list()
",
     "  function(z, d = 0) {
",
     "    for (m in memo) if (m$d == d && identical(m$z, z)) return(m$value)
"))
ed("R/me.R",
   paste0(
     "    memo[[key]] <- what(sh)
",
     "    memo[[key]]
"),
   paste0(
     "    value <- what(sh)
",
     "    memo[[length(memo) + 1L]] <<- list(z = z, d = d, value = value)
",
     "    value
"))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_LIBS = paste(lib, "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6", sep = ";"))
r <- system2(file.path(R.home("bin"), "R.exe"), c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load", paste0("--library=", lib), shQuote(d)), stdout = TRUE, stderr = TRUE)
cat(tail(r, 1), "\n")
r <- system2(file.path(R.home("bin"), "R.exe"), c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load", paste0("--library=", lib), shQuote(file.path(src, "extensions/frmtmb.sample"))), stdout = TRUE, stderr = TRUE)
cat(tail(r, 1), "\n")
