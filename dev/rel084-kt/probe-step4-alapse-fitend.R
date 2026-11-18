# rel084-kt, merge step 4 (alapse x gpselimit at the fit end): on alapse's
# D (y ~ x, sigma ~ z) and E (y ~ x, quantile ~ z) asym_laplace() fits of
# dev/alapse-rq.R (same data: seeds 4000 + s and 5000 + s, s = 1..25),
# record every condition class each fit raises (full class vectors), and
# whether any is a gp or boundary message saying the fit holds:
# frmtmb_gp_basis, frmtmb_boundary_fit, frmtmb_limit_fit,
# frmtmb_converged_at_max.
# Usage: Rscript dev/rel084-kt/probe-step4-alapse-fitend.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
bad <- c("frmtmb_gp_basis", "frmtmb_boundary_fit", "frmtmb_limit_fit",
         "frmtmb_converged_at_max")
conds_of <- function(expr) {
  said <- list()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      said[[length(said) + 1L]] <<- c("w", class(w))
      invokeRestart("muffleWarning")
    }, message = function(m) {
      said[[length(said) + 1L]] <<- c("m", class(m))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      said[[length(said) + 1L]] <<- c("E", class(e))
      NULL
    })
  list(v = v, said = said)
}
tot <- c(n = 0, local = 0, bad = 0)
allcls <- character(0)
one <- function(id, s, d, form) {
  r <- conds_of(frm(form, family = asym_laplace(), data = d))
  cl <- unique(unlist(lapply(r$said, `[`, -1L)))
  str <- paste(vapply(r$said, function(z) paste0(z[1], ":", paste(z[-1], collapse = "/")), ""),
               collapse = " ; ")
  hit <- intersect(cl, bad)
  code <- if (is.null(r$v)) NA else r$v$opt$convergence
  cat(sprintf("%s %2d code %s | local %s | gp/boundary hit [%s] | %s\n", id, s,
              code, "frmtmb_ald_local" %in% cl, paste(hit, collapse = ","), str))
  tot[["n"]] <<- tot[["n"]] + 1
  tot[["local"]] <<- tot[["local"]] + ("frmtmb_ald_local" %in% cl)
  tot[["bad"]] <<- tot[["bad"]] + (length(hit) > 0)
  allcls <<- union(allcls, cl)
}
for (s in 1:25) {
  set.seed(4000 + s)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- 1 + 0.5 * d$x + exp(0.4 * d$z) * rnorm(200)
  one("D", s, d, bf(y ~ x, sigma ~ z))
  set.seed(5000 + s)
  d <- data.frame(x = rnorm(200), z = rnorm(200))
  d$y <- 1 + 0.5 * d$x + rnorm(200)
  one("E", s, d, bf(y ~ x, quantile ~ z))
}
cat(sprintf("TOTAL fits %d | ald_local %d | with a gp/boundary/limit/converged_at_max class %d\n",
            tot[["n"]], tot[["local"]], tot[["bad"]]))
cat("ALL CLASSES SEEN:", paste(sort(allcls), collapse = " "), "\n")
