## Reviewer, claim 2, the two shapes the first identity run could not
## build: a nonlinear fit and an importance-corrected fit.
## usage: Rscript gradcheck-rev-07-ident2.R <core-lib> <out.rds>
a <- commandArgs(TRUE)
LIB <- a[1]; OUT <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                             warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w)
}
row_of <- function(tag, r) {
  f <- r$fit
  if (inherits(f, "condition")) {
    return(list(tag = tag, err = conditionMessage(f)))
  }
  g <- try(drop(f$obj$gr(f$opt$par)), silent = TRUE)
  list(tag = tag, err = NA_character_, coef = coef(f), par = f$opt$par,
       logLik = as.numeric(logLik(f)),
       objective = as.numeric(f$obj$fn(f$opt$par)),
       opt_obj = f$opt$objective, convergence = f$opt$convergence,
       gmax = if (inherits(g, "try-error")) NA_real_ else max(abs(g)),
       grad = if (inherits(g, "try-error")) NA_real_ else g,
       warns = sort(r$w))
}

models <- list()
models$nonlinear <- function() {
  set.seed(7003)
  n <- 700
  dd <- data.frame(x = runif(n, 0, 5))
  dd$y <- rnorm(n, 3 * exp(-0.8 * dd$x), 0.2)
  frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian(), data = dd)
}
models$nonlinear_short <- function() {
  set.seed(7003)
  n <- 700
  dd <- data.frame(x = runif(n, 0, 5))
  dd$y <- rnorm(n, 3 * exp(-0.8 * dd$x), 0.2)
  frm(bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) + gaussian(), data = dd,
      control = frmtmb_control(restarts = 0,
                               optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                              iter.max = 1000,
                                              eval.max = 1000)))
}
models$importance <- function() {
  set.seed(7007)
  ng <- 45
  dd <- data.frame(g = factor(rep(seq_len(ng), 20)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.6)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd,
      importance = 1000L)
}
models$importance_big <- function() {
  set.seed(7107)
  ng <- 300
  dd <- data.frame(g = factor(rep(seq_len(ng), 20)))
  dd$x <- rnorm(nrow(dd))
  re <- rnorm(ng, 0, 0.6)
  dd$y <- rpois(nrow(dd), exp(0.3 + 0.5 * dd$x + re[dd$g]))
  frm(bf(y ~ x + (1 | g)), family = poisson(), data = dd,
      importance = 500L)
}
out <- list()
for (nm in names(models)) {
  r <- gw(models[[nm]]())
  out[[nm]] <- row_of(nm, r)
  o <- out[[nm]]
  cat(sprintf("%-18s err=%s conv=%s logLik=%s gmax=%s nwarn=%d\n", nm,
              o$err %||% NA, format(o$convergence %||% NA),
              format(o$logLik %||% NA, digits = 15),
              format(o$gmax %||% NA, digits = 6),
              length(o$warns %||% character(0))))
  for (w in o$warns %||% character(0)) cat("     WARN:", w, "\n")
  if (!inherits(r$fit, "condition")) {
    d <- tryCatch(diagnose(r$fit, quiet = TRUE), error = function(e) NULL)
    if (!is.null(d)) {
      cat("     max_grad=", format(d$max_grad, digits = 5),
          " grad_proj=", format(d$grad_proj %||% NA, digits = 5),
          " headroom=", format(d$grad_headroom %||% NA, digits = 5), "\n",
          sep = "")
    }
  }
}
saveRDS(out, OUT)
cat("DONE ident2\n")
