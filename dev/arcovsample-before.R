# Lane wt-arcovsample, step 1: the refusals as the REFERENCE build
# (frmtmb 0.64.0) gives them, for brms's default cov = FALSE ARMA.
#
# Run against C:/Users/adf44/source/r/rellib-r3 only. Nothing here
# installs anything.
#
#   Rscript dev/arcovsample-before.R > dev/arcovsample-log/before.txt

LIB <- "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
cat("frmtmb        ", format(packageVersion("frmtmb")), "\n")
cat("frmtmb.sample ", format(packageVersion("frmtmb.sample")), "\n")
cat("tmbstan       ", format(packageVersion("tmbstan")), "\n")

suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))

# ---- data: ragged groups with a gap, shuffled rows -------------------
# seed 4021; the same shape the arcov lane used, small enough to sample
arcov_data <- function(seed = 4021L, ng = 8L, nt = 10L) {
  set.seed(seed)
  dd <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
  dd <- dd[-c(3L, 17L, 41L, 55L), ]          # interior gaps
  dd$x <- stats::rnorm(nrow(dd))
  u <- stats::rnorm(ng, 0, 0.6)
  e <- stats::rnorm(nrow(dd))
  dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$g)] + e
  dd[sample(nrow(dd)), ]
}

dd <- arcov_data()
cat("N =", nrow(dd), " groups =", length(unique(dd$g)), "\n")

say <- function(lab, expr) {
  v <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(v, "condition")) {
    cat("\n[", lab, "] ", class(v)[1L], ": ",
        conditionMessage(v), "\n", sep = "")
  } else {
    sz <- if (is.null(dim(v))) length(v) else dim(v)
    cat("\n[", lab, "] OK ",
        paste(class(v), collapse = "/"), " size=",
        paste(sz, collapse = "x"), "\n", sep = "")
  }
  invisible(v)
}

cases <- list(
  gauss_ar1 = list(
    form = bf(y ~ x + ar(t, g) + (1 | g)), fam = gaussian()),
  student_arma11 = list(
    form = bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)),
    fam = student())
)

for (nm in names(cases)) {
  cat("\n================ ", nm, " ================\n")
  cs <- cases[[nm]]
  ds <- suppressWarnings(suppressMessages(
    frm_sample(cs$form, family = cs$fam, data = dd,
               chains = 1, iter = 400, refresh = 0, seed = 11)))
  cat("draws: ", nrow(ds$draws), " x ", ncol(ds$draws), "\n", sep = "")

  say("log_lik", log_lik(ds))
  say("loo", loo(ds))
  say("waic", waic(ds))
  say("psis", psis(ds))
  say("loo_compare", loo_compare(ds, ds))
  say("loo_moment_match", loo_moment_match(ds))
  say("bayes_R2", bayes_R2(ds))
  # loo_R2() is not in frmtmb.sample's NAMESPACE; recorded as absent
  say("loo_subsample", loo_subsample(ds))
  say("pp_check-dens_overlay", pp_check(ds, type = "dens_overlay"))
  say("pp_check-loo_pit", pp_check(ds, type = "loo_pit_overlay"))
  say("posterior_predict", posterior_predict(ds))
  say("posterior_epred", posterior_epred(ds))
  say("kfold", kfold(ds))
}

cat("\nDONE\n")
