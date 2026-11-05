# Shared set-up of the review's scripts. Sourced from the lane's worktree.
#   REV_LIB   the library holding the frmtmb under test
#             (default the lane's; "before" means rellib-r16 alone)
rev_root <- "C:/Users/adf44/source/r/frmtmb-wt-predread"
setwd(rev_root)
rev_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
rev_r16 <- "C:/Users/adf44/source/r/rellib-r16"
rev_lib <- Sys.getenv("REV_LIB", "C:/Users/adf44/source/r/wt-predread-lib")
if (identical(rev_lib, "before")) {
  .libPaths(c(rev_r16, rev_user))
} else {
  .libPaths(c(rev_lib, rev_r16, rev_user))
}
Sys.setenv(FRMTMB_STAN_CACHE = file.path(rev_root, "dev/stan-cache"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
rev_E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = rev_E)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
rev_ev <- function(txt, pkg) {
  eval(parse(text = txt)[[1L]], envir = list2env(list(),
       parent = if (pkg == "brms") asNamespace("brms") else
         asNamespace("frmtmb")))
}
# a brmsfit at the frm() estimate, the same spelling evaluated by each package
rev_pair <- function(ftxt, famtxt, d, ndraws = 20, ..., fit_args = list()) {
  ff <- rev_ev(ftxt, "frmtmb")
  fb <- rev_ev(ftxt, "brms")
  af <- rev_ev(famtxt, "frmtmb")
  ab <- rev_ev(famtxt, "brms")
  fit <- quiet(do.call(frm, c(list(ff, family = af, data = d), fit_args)))
  fx <- quiet(rev_E$brms_fixed_fit(fb, ab, d, fit, ndraws = ndraws, ...))
  list(fit = fit, fx = fx)
}
rev_fake <- function(fit, n = 4L) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  structure(list(stanfit = NULL,
                 draws = matrix(0, n, length(lab),
                                dimnames = list(NULL, lab)),
                 fit = fit),
            class = "frmtmb_draws")
}
try_num <- function(e) {
  tryCatch(quiet(e), error = function(err) {
    structure(NA_real_, err = substr(gsub("[\r\n\t]+", " ",
                                          conditionMessage(err)), 1, 120))
  })
}
rd <- function(a, b) {
  if (!is.null(attr(a, "err")) || !is.null(attr(b, "err"))) return(NA_real_)
  rev_E$eq_rel(a, b)
}
# frmtmb_draws whose every draw IS the frm() estimate, for a comparison with a
# Fixed_param brmsfit at the same values
rev_draws_at <- function(fit, n = 4L) {
  ns <- asNamespace("frmtmb.sample")
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  idx <- ns$draws_par_index(fit)
  row <- numeric(length(lab))
  for (cp in names(idx)) {
    v <- fit$estimates[[cp]]
    if (cp == "betad" && length(fx <- fit$frame[["betad_fixed_idx"]])) {
      v <- v[setdiff(seq_along(v), fx)]
    }
    row[idx[[cp]]] <- v
  }
  m <- matrix(rep(row, each = n), n, length(lab), dimnames = list(NULL, lab))
  structure(list(stanfit = NULL, draws = ns$draws_to_natural(m, fit), fit = fit),
            class = "frmtmb_draws")
}
