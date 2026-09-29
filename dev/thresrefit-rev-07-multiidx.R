## REVIEW claim 5, second attempt: reach a multi-index prior entry that is
## NOT the ordinal threshold vector, so that the new refusal in
## draw_prior_entry() can be checked against what the base build did
## there. Coverage is supplied by priors alone, so newparams naming
## cannot stop the call before the draw.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 140)

say <- function(tag, expr) {
  cat("==", tag, "\n")
  out <- tryCatch(expr, error = function(e) e, warning = function(w) w)
  if (inherits(out, "condition")) {
    cat("   ", class(out)[1L], ":",
        gsub("\n", " ", conditionMessage(out)), "\n\n")
  } else {
    cat("    OK\n"); print(out); cat("\n")
  }
  invisible(out)
}

## the entries the resolver builds, read directly, so the index lengths
## are a measurement and not a guess
entries_of <- function(form, family, data, prior) {
  bform <- frmtmb:::resolve_deferred_families(
    frmtmb:::as_bform(form, family), data)
  spec <- frmtmb:::parse_spec(bform)
  frame <- frmtmb:::assemble_frame(spec, data)
  shim <- list(spec = spec, frame = frame,
               estimates = frame[["par_template"]])
  e <- frmtmb:::resolve_prior_input(shim, frmtmb:::as_priorlist(prior))$entries
  for (i in seq_along(e)) {
    cat("    entry", i, ": comp =", e[[i]]$comp,
        " idx =", paste(e[[i]]$idx, collapse = ","),
        " (length", length(e[[i]]$idx), ") scale =",
        e[[i]]$scale %||% "NULL", "\n")
  }
  invisible(e)
}

## ---- class "cor": an LKJ density over a 3 x 3 block's correlations ---
set.seed(801)
n <- 120
d <- data.frame(g = factor(rep(1:20, each = 6)), x = stats::rnorm(n),
                z = stats::rnorm(n))
d$y <- stats::rnorm(n)
pr_cor <- set_prior("normal(0, 1)", class = "b") +
  set_prior("normal(0, 1)", class = "Intercept") +
  set_prior("normal(0, 1)", class = "sd") +
  set_prior("normal(0, 1)", class = "sigma") +
  set_prior("lkj(2)", class = "cor")
cat("entries for the (1 + x + z | g) model:\n")
try(entries_of(bf(y ~ x + z + (1 + x + z | g)), gaussian(), d, pr_cor),
    silent = FALSE)
say("frm_simulate with a prior on class cor",
    {
      s <- frm_simulate(bf(y ~ x + z + (1 + x + z | g)), family = gaussian(),
                        data = d, prior = pr_cor, nsim = 2, seed = 820)
      attr(s, "pars")
    })

## ---- class "ar" at order 2: two thetaac parameters -------------------
set.seed(802)
d2 <- data.frame(t = 1:80, gg = factor(rep(1, 80)), x = stats::rnorm(80))
d2$y <- stats::rnorm(80)
pr_ar <- set_prior("normal(0, 1)", class = "b") +
  set_prior("normal(0, 1)", class = "Intercept") +
  set_prior("normal(0, 1)", class = "sigma") +
  set_prior("normal(0, 0.3)", class = "ar")
cat("entries for the ar(p = 2) model:\n")
try(entries_of(bf(y ~ x + ar(time = t, gr = gg, p = 2)), gaussian(), d2,
               pr_ar), silent = FALSE)
say("frm_simulate with a prior on class ar at order 2",
    {
      s <- frm_simulate(bf(y ~ x + ar(time = t, gr = gg, p = 2)),
                        family = gaussian(), data = d2, prior = pr_ar,
                        nsim = 2, seed = 821)
      attr(s, "pars")
    })

## ---- the SAME two models in the FIT path, which must be untouched ---
cat("#### the fit path on the same two models\n")
f1 <- suppressWarnings(frm(bf(y ~ x + z + (1 + x + z | g)),
                          family = gaussian(), data = d, prior = pr_cor))
cat("cor model: logLik =", sprintf("%.14f", logLik(f1)),
    " objective =", sprintf("%.14f", f1$opt$objective), "\n")
cat("   theta =", paste(sprintf("%.12g", f1$estimates[["theta"]]),
                        collapse = " "), "\n")
f2 <- suppressWarnings(frm(bf(y ~ x + ar(time = t, gr = gg, p = 2)),
                          family = gaussian(), data = d2, prior = pr_ar))
cat("ar model: logLik =", sprintf("%.14f", logLik(f2)),
    " objective =", sprintf("%.14f", f2$opt$objective), "\n")
cat("   thetaac =", paste(sprintf("%.12g", f2$estimates[["thetaac"]]),
                          collapse = " "), "\n")
cat("DONE rev-07\n")
