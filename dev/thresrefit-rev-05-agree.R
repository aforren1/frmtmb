## REVIEW claim 4: the worker says influence rows equal a hand-pinned fit
## to 2.4e-5 relative at worst. Reproduce on the worker's seeds and on
## fresh ones, and settle whether the residual is optimizer noise on a
## flat ridge or a different model, by comparing LOG-LIKELIHOODS and by
## evaluating one fit's objective at the other's parameter vector.
lib <- Sys.getenv("FRMTMB_LIB")
lib <- if (identical(lib, "base")) "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-thresrefit-lib"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("## lib =", lib, " ver =",
    as.character(utils::packageVersion("frmtmb")), "\n\n")
options(width = 150)

mk <- function(seed, n = 50) {
  set.seed(seed)
  x <- stats::rnorm(n)
  cp <- cbind(stats::plogis(-0.7 - 0.5 * x), stats::plogis(0.6 - 0.5 * x),
              stats::plogis(2.3 - 0.5 * x))
  y <- 1L + rowSums(stats::runif(n) > cp)
  top <- which(y == 4L)
  if (!length(top)) return(NULL)
  y[top[-1L]] <- 3L
  data.frame(x = x, y = y)
}

## the leave-one-out refit influence() itself performs, reproduced so its
## log-likelihood can be read: same spec, same pin, same template rule
loo_refit <- function(fit, data, drop) {
  ctl <- fit$control %||% frmtmb_control()
  ctl$verbose <- FALSE
  pin <- if (exists("thres_pin_of_fit", envir = asNamespace("frmtmb"),
                    inherits = FALSE)) {
    frmtmb:::thres_pin_of_fit(fit)
  } else NULL
  fr <- if (is.null(pin)) {
    frmtmb:::assemble_frame(fit$spec, data[-drop, , drop = FALSE],
                            sparse_x = isTRUE(ctl$sparse_x),
                            data2 = fit$data2 %||% list())
  } else {
    frmtmb:::assemble_frame(fit$spec, data[-drop, , drop = FALSE],
                            sparse_x = isTRUE(ctl$sparse_x),
                            data2 = fit$data2 %||% list(), thres_pin = pin)
  }
  tpl <- fr$par_template
  for (cp in setdiff(names(tpl), "b")) {
    if (length(fit$estimates[[cp]]) == length(tpl[[cp]])) {
      tpl[[cp]][] <- fit$estimates[[cp]]
    }
  }
  frmtmb:::fit_assembled(fit$spec, fr, fit$bform, fit$call,
                         REML = fit$REML, start = NULL, control = ctl,
                         se = FALSE, lower = fit$lower, upper = fit$upper,
                         prior = fit$prior,
                         quadrature = isTRUE(fit$quadrature),
                         importance = fit$importance$draws %||% 0L,
                         template = tpl, data2 = fit$data2 %||% list())
}

seeds <- c(501, 1501, 1502, 1503, 1504)
for (sd0 in seeds) {
  dd <- mk(sd0)
  if (is.null(dd)) { cat("seed", sd0, ": no top category, skipped\n"); next }
  itop <- which(dd$y == 4L)
  if (length(itop) != 1L) { cat("seed", sd0, ": bad build\n"); next }
  for (fn in c("cumulative", "sratio")) {
    fam <- switch(fn, cumulative = cumulative(), sratio = sratio())
    fit <- suppressWarnings(frm(bf(y ~ x), family = fam, data = dd))
    inf <- suppressWarnings(influence(fit, force = TRUE))
    got <- inf$fixed[itop, ]
    ref <- suppressWarnings(frm(bf(y | thres(3) ~ x), family = fam,
                               data = dd[-itop, , drop = FALSE]))
    want <- frmtmb:::get_coef.frmtmb_fit(ref)
    rr <- try(loo_refit(fit, dd, itop), silent = TRUE)
    cat("== seed", sd0, fn, " table(y) =",
        paste(table(dd$y), collapse = "/"), " deleted row =", itop, "\n")
    for (nm in names(want)) {
      cat(sprintf("   %-10s got %.15g  want %.15g  rel %.3g\n", nm,
                  got[[nm]], want[[nm]],
                  abs(got[[nm]] - want[[nm]]) / abs(want[[nm]])))
    }
    sp <- apply(inf$fixed, 2, stats::sd)
    ident <- setdiff(names(want), "tau_raw_3")
    cat("   max rel over identified:",
        signif(max(abs(got[ident] - want[ident]) / abs(want[ident])), 4),
        "  max |diff|/sd(col):",
        signif(max(abs(got[ident] - want[ident]) / sp[ident]), 4), "\n")
    if (!inherits(rr, "try-error")) {
      l1 <- as.numeric(logLik(rr)); l2 <- as.numeric(logLik(ref))
      cat(sprintf("   logLik loo-refit = %.14f\n", l1))
      cat(sprintf("   logLik reference = %.14f   rel = %.4g  abs = %.3g\n",
                  l2, abs(l1 - l2) / abs(l2), abs(l1 - l2)))
      cat("   n_tau loo-refit =", length(rr$estimates[["tau_raw"]]),
          " reference =", length(ref$estimates[["tau_raw"]]),
          "  coef names identical =",
          identical(names(frmtmb:::get_coef.frmtmb_fit(rr)), names(want)),
          "\n")
      cat("   loo-refit coefs match the influence row bitwise:",
          identical(unname(frmtmb:::get_coef.frmtmb_fit(rr)),
                    unname(got)), "\n")
      ## the ridge test: the reference objective at the loo-refit's own
      ## parameter vector, against the reference's own optimum
      o <- ref$obj
      if (!is.null(o) && is.function(o$fn)) {
        p_ref <- o$par
        p_loo <- rr$obj$par
        if (identical(names(p_ref), names(p_loo))) {
          f_ref <- o$fn(p_ref); f_loo <- o$fn(p_loo)
          cat(sprintf(paste0("   reference objective at its own optimum ",
                             "%.14f, at the loo point %.14f, gap %.3g\n"),
                      f_ref, f_loo, f_loo - f_ref))
        } else {
          cat("   objective par names differ, ridge test skipped\n")
        }
      }
    } else {
      cat("   loo refit reproduction FAILED:",
          conditionMessage(attr(rr, "condition")), "\n")
    }
    cat("\n")
  }
}
cat("DONE rev-05\n")
