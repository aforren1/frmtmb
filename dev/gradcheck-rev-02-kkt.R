## Reviewer, claim 1(a): the KKT sign test.
## usage: Rscript gradcheck-rev-02-kkt.R <core-lib>
## Every row prints the fit's own numbers, whether the gradient warning
## fired, and (on the lane build) the three new diagnose() fields.
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")

gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w,
       grad = any(grepl("Large maximum absolute gradient", w, fixed = TRUE)))
}
report <- function(tag, r) {
  f <- r$fit
  d <- tryCatch(diagnose(f, quiet = TRUE), error = function(e) NULL)
  nm <- frmtmb:::outer_par_names(f)
  cat("### ", tag, "\n", sep = "")
  cat("  conv=", f$opt$convergence,
      "  logLik=", format(as.numeric(logLik(f)), digits = 12),
      "  GRADWARN=", r$grad, "\n", sep = "")
  cat("  par: ", paste(sprintf("%s=%.8g", nm, f$opt$par),
                       collapse = "  "), "\n", sep = "")
  g <- drop(f$obj$gr(f$opt$par))
  cat("  grad:", paste(sprintf("%s=%.4g", nm, g), collapse = "  "), "\n")
  bx <- tryCatch(frmtmb:::fit_outer_box(f), error = function(e) NULL)
  if (!is.null(bx)) {
    cat("  box lower:", paste(format(bx$lower, digits = 5),
                              collapse = " "), "\n")
    cat("  box upper:", paste(format(bx$upper, digits = 5),
                              collapse = " "), "\n")
  }
  if (!is.null(d)) {
    cat("  max_grad=", format(d$max_grad, digits = 5),
        "  grad_proj=", format(d$grad_proj %||% NA, digits = 5),
        "  headroom=", format(d$grad_headroom %||% NA, digits = 5),
        "  held=[", paste(d$grad_bound_held %||% character(0),
                          collapse = ","), "]", sep = "")
    cat("  pdHess=", d$pdHess, "\n", sep = "")
  }
  if (length(r$w)) for (x in r$w) cat("  WARN: ", x, "\n", sep = "")
  cat("\n")
  invisible(d)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- A1: a parameter ON a bound whose gradient points INTO the
## feasible region. The bound is NOT binding and the fit is genuinely
## short, so the check must still warn. Constructed by starting the
## optimizer at the bound and capping it at one iteration.
set.seed(9101)
n <- 400
dA <- data.frame(x = rnorm(n))
dA$y <- rnorm(n, 1 + 2 * dA$x, 1)
rA <- gw(frm(bf(y ~ x), family = gaussian(), data = dA,
             prior = set_prior("", class = "b", lb = 0.1),
             start = list(beta = c(0, 0.1)),
             control = frmtmb_control(
               restarts = 0,
               optCtrl = list(iter.max = 1, eval.max = 3))))
dd <- report("A1 lb=0.1 on b, started AT the bound, iter.max=1", rA)
refA <- frm(bf(y ~ x), family = gaussian(), data = dA,
            prior = set_prior("", class = "b", lb = 0.1))
cat("A1 shortfall vs converged bounded fit:",
    format(as.numeric(logLik(refA)) - as.numeric(logLik(rA$fit)),
           digits = 8), "log-lik units\n")
cat("A1 converged bounded fit x =",
    format(refA$opt$par[[2]], digits = 8), "(bound 0.1, so NOT binding)\n\n")

## ---- A2: an interior parameter and a bound-held parameter in the same
## fit. ub = 0.1 with a true slope of 2 binds; sigma is interior.
set.seed(101)
n <- 300
dB <- data.frame(x = rnorm(n))
dB$y <- rnorm(n, 1 + 2 * dB$x, 1)
rB <- gw(frm(bf(y ~ x), family = gaussian(), data = dB,
             prior = set_prior("", class = "b", ub = 0.1)))
report("A2 ub=0.1 on b (binds), sigma interior", rB)

## ---- A3: an lb AND a ub active in the same fit, on two coefficients.
set.seed(9103)
n <- 600
dC <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
dC$y <- rnorm(n, 1 + 2 * dC$x1 - 2 * dC$x2, 1)
rC <- gw(frm(bf(y ~ x1 + x2), family = gaussian(), data = dC,
             prior = c(set_prior("", class = "b", coef = "x1", ub = 0.1),
                       set_prior("", class = "b", coef = "x2", lb = -0.1))))
report("A3 ub on x1 AND lb on x2, both bind", rC)

## ---- A4: a bound on a standard deviation, which lives on the log
## scale internally. lb = 2 on an sd whose truth is 0.2 binds from below.
set.seed(9104)
ng <- 40
dD <- data.frame(g = factor(rep(seq_len(ng), 25)))
dD$x <- rnorm(nrow(dD))
re <- rnorm(ng, 0, 0.2)
dD$y <- rnorm(nrow(dD), 1 + 0.5 * dD$x + re[dD$g], 1)
rD <- gw(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dD,
             prior = set_prior("", class = "sd", lb = 2)))
report("A4 lb=2 on class sd (log scale internally), binds from below",
       rD)

## ---- A5: par_units. An autoscaled fit, so the check's natural-unit
## gradient and diagnose()$max_grad are on different scales.
set.seed(9105)
n <- 4000
dE <- data.frame(xs = rnorm(n) * 1e-4 + 5e-4)
dE$y <- rnorm(n, 1 + 2000 * dE$xs, 1)
rE <- gw(frm(bf(y ~ xs), family = gaussian(), data = dE))
dE1 <- report("A5 autoscaled column (spread 1e-4), autoscale default", rE)
cat("A5 par_units:", paste(format(rE$fit$par_units, digits = 6),
                           collapse = " "), "\n")
g <- drop(rE$fit$obj$gr(rE$fit$opt$par))
cat("A5 raw max|g| =", format(max(abs(g)), digits = 6),
    " natural-unit max|g| =",
    format(max(abs(g * (rE$fit$par_units %||% 1))), digits = 6), "\n\n")

## ---- A6: the same autoscaled shape WITH a binding bound, so the KKT
## test and the natural-unit projection meet in one fit.
rF <- gw(frm(bf(y ~ xs), family = gaussian(), data = dE,
             prior = set_prior("", class = "b", ub = 100)))
report("A6 autoscaled column AND ub=100 on b", rF)
cat("A6 par_units:", paste(format(rF$fit$par_units, digits = 6),
                           collapse = " "), "\n")

cat("DONE kkt\n")
