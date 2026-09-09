# Reviewer checks, round two.
#
# 1. Does probe 09's YAML extraction reproduce the step verbatim?
# 2. Are the R/ edits examples-only, i.e. is the installed code at the
#    base commit still identical to what the tree parses to?
# 3. Has the new gradient assertion been SEEN to fail? Run its exact
#    expression with the defect's gradient substituted for rstan's.
LIB <- "C:/Users/adf44/source/r/lanelib-tmbstan"
REF <- "C:/Users/adf44/source/r/reflib-r2"
.libPaths(c(LIB, REF, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-tmbstan"

cat("== 1. probe 09's extraction against the YAML ==\n")
yaml <- file.path(root, ".github/workflows/check-frmtmb-sample.yaml")
src <- readLines(yaml)
i <- grep("^      - name: Pin a tmbstan", src)
cat("line after the name is 'run: |':",
    identical(trimws(src[[i + 1L]]), "run: |"), "\n")
j <- grep("^        shell: Rscript \\{0\\}", src)
j <- min(j[j > i])
body <- src[(i + 2L):(j - 1L)]
cat("every body line carries the 10-space block indent or is blank:",
    all(grepl("^          ", body) | !nzchar(trimws(body))), "\n")
ded <- sub("^          ", "", body)
cat("dedent is lossless (re-indent restores the YAML):",
    identical(paste0("          ", ded), body), "\n")
cat("the extracted body parses:",
    !inherits(try(parse(text = paste(ded, collapse = "\n")),
                  silent = TRUE), "try-error"), "\n")
cut <- grep("^install.packages", ded)[1L]
cat("the probe evaluates lines 1 to", cut - 1L, "of", length(ded),
    "; the stop() branches after the verdict are at lines",
    paste(grep("^\\s*stop\\(", ded), collapse = ", "), "\n")

cat("\n== 2. are the R/ edits examples-only ==\n")
# roxygen comments do not survive into an installed package, so a
# difference between the base commit's INSTALLED code and what the
# tree parses to is a real code change, not a documentation one.
ns <- asNamespace("frmtmb.sample")
detach_ok <- TRUE
refns <- new.env(parent = emptyenv())
lazyLoad(file.path(REF, "frmtmb.sample", "R", "frmtmb.sample"),
         envir = refns)
treeenv <- new.env(parent = globalenv())
for (f in list.files(file.path(root, "extensions/frmtmb.sample/R"),
                     pattern = "[.]R$", full.names = TRUE)) {
  eval(parse(f, keep.source = FALSE), envir = treeenv)
}
nm <- sort(ls(refns, all.names = TRUE))
diffs <- character(0)
for (k in nm) {
  if (!exists(k, envir = treeenv, inherits = FALSE)) {
    diffs <- c(diffs, paste0(k, " (absent from the tree)")); next
  }
  a <- get(k, envir = refns)
  b <- get(k, envir = treeenv)
  if (is.function(a) && is.function(b)) {
    if (!identical(deparse(a), deparse(b))) diffs <- c(diffs, k)
  } else if (!identical(a, b)) {
    diffs <- c(diffs, paste0(k, " (non-function)"))
  }
}
cat("objects in the base install:", length(nm), "\n")
cat("objects whose code differs from the tree:", length(diffs), "\n")
if (length(diffs)) cat(paste0("  ", diffs), sep = "\n")
extra <- setdiff(ls(treeenv, all.names = TRUE), nm)
cat("objects the tree defines that the base install does not:",
    if (length(extra)) paste(extra, collapse = ", ") else "none", "\n")

cat("\n== 3. the gradient assertion, seen to fail ==\n")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
library(testthat)
# the block's own model and seed
set.seed(4021)
dd <- data.frame(x = stats::rnorm(80))
dd$y <- stats::rnorm(80, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
sf <- suppressWarnings(suppressMessages(
  as_tmbstan(fit, chains = 1L, iter = 20L, warmup = 10L,
             refresh = 0, seed = 11)))
np <- rstan::get_num_upars(sf)
set.seed(4030)
pts <- list(rep(0, np), fit$obj$par + 0.3, fit$obj$par - 0.7,
            stats::rnorm(np))
arm <- function(subst) {
  fired <- 0L
  for (p in pts) {
    u <- as.numeric(p)
    g_real <- as.numeric(rstan::grad_log_prob(sf, u))
    g_obj <- -as.numeric(fit$obj$gr(u))
    g_bad <- -u
    g_stan <- if (subst) g_bad else g_real
    scale <- max(abs(g_obj))
    lhs <- max(abs(g_stan - g_obj)) / scale
    rhs <- max(abs(g_bad - g_obj)) / scale / 1e6
    cat(sprintf("  lhs=%.3e  rhs=%.3e  passes=%s  ulp=%s\n", lhs, rhs,
                lhs < rhs,
                if (subst) "n/a" else
                  max(abs(as.numeric(sapply(seq_along(g_stan),
                    function(i) (g_stan[i] - g_obj[i]))))) == 0))
    if (!(lhs < rhs)) fired <- fired + 1L
  }
  fired
}
cat(" healthy arm (rstan's own gradient):\n")
f1 <- arm(FALSE)
cat(" defect arm (gradient replaced by -u):\n")
f2 <- arm(TRUE)
cat("blocks failing: healthy", f1, " defect", f2, "of", length(pts),
    "\n")
cat("bitwise identical on the healthy arm:",
    identical(as.numeric(rstan::grad_log_prob(sf, fit$obj$par + 0.3)),
              -as.numeric(fit$obj$gr(fit$obj$par + 0.3))), "\n")
