## Reviewer, claim 1: the pieces the reachable fits in -02- could not
## reach, plus the diagnose() clean-line channel.
## usage: Rscript gradcheck-rev-03-probe.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w,
       grad = any(grepl("Large maximum absolute gradient", w, fixed = TRUE)))
}
has_new <- exists("grad_verdict", envir = asNamespace("frmtmb"))
cat("build has grad_verdict():", has_new, "\n\n")

## ---- B1: a parameter ON a bound whose gradient points INTO the
## feasible set. Not reachable by an ordinary fit (the optimizer leaves
## the bound), so the fit object is poked onto the bound and the verdict
## recomputed. This is the case that must NOT be excluded.
set.seed(9201)
n <- 400
d1 <- data.frame(x = rnorm(n))
d1$y <- rnorm(n, 1 - 2 * d1$x, 1)          # truth well BELOW the bound
f1 <- frm(bf(y ~ x), family = gaussian(), data = d1,
          prior = set_prior("", class = "b", ub = 0.1))
nm <- frmtmb:::outer_par_names(f1)
j <- match("x", nm)
cat("B1 unconstrained-in-box optimum x =", format(f1$opt$par[[j]], digits = 8),
    " (ub 0.1 does NOT bind)\n")
if (has_new) {
  f2 <- f1
  f2$cache <- new.env(parent = emptyenv())
  f2$cache$bounds <- frmtmb:::fit_outer_box(f1)
  f2$opt$par[[j]] <- 0.1                    # put it ON the bound
  g2 <- drop(f2$obj$gr(f2$opt$par))
  act <- frmtmb:::grad_bound_active(f2$opt$par, g2, f2$cache$bounds)
  v2 <- frmtmb:::grad_verdict(f2, g2, f1$control)
  cat("B1 at x = 0.1 the gradient is", format(g2[[j]], digits = 6),
      "(positive = the minimizer wants to move DOWN, into the box)\n")
  cat("B1 grad_bound_active:", paste(act, collapse = " "), "\n")
  cat("B1 verdict: proj =", format(v2$proj, digits = 6),
      " headroom =", format(v2$headroom, digits = 6),
      " warn =", v2$warn, " held = [",
      paste(v2$bound_held, collapse = ","), "]\n")
  sh <- as.numeric(f2$obj$fn(f2$opt$par)) - as.numeric(f1$obj$fn(f1$opt$par))
  cat("B1 true shortfall at that point:", format(sh, digits = 8),
      "log-lik units\n\n")
}

## ---- B2: can grad_proj EXCEED max_grad? diagnose() prints the two on
## adjacent lines. max_grad is the RAW largest component; grad_proj is
## the largest component in the natural units par_units works in.
set.seed(9202)
n <- 3000
d2 <- data.frame(xs = rnorm(n) * 1e-5, z = rnorm(n))
d2$y <- rnorm(n, 1 + 2e5 * d2$xs + 0.5 * d2$z, 1)
r2 <- gw(frm(bf(y ~ xs + z), family = gaussian(), data = d2,
             prior = set_prior("", class = "b", coef = "z", ub = 0.1),
             control = frmtmb_control(
               restarts = 0,
               optCtrl = list(rel.tol = 1e-3, x.tol = 1e-3,
                              iter.max = 1000, eval.max = 1000))))
d <- diagnose(r2$fit, quiet = TRUE)
cat("B2 par_units:", paste(format(r2$fit$par_units %||% 1, digits = 6),
                           collapse = " "), "\n")
cat("B2 max_grad =", format(d$max_grad, digits = 6),
    " grad_proj =", format(d$grad_proj %||% NA, digits = 6),
    " held = [", paste(d$grad_bound_held %||% character(0),
                       collapse = ","), "]\n")
cat("B2 grad_proj > max_grad :", isTRUE((d$grad_proj %||% -Inf) > d$max_grad),
    "\n")
cat("B2 printed block:\n")
writeLines(paste0("   ", capture.output(diagnose(r2$fit))))
cat("\n")

## ---- B3: the clean line. Old: hardcoded max_grad < 1e-3 on the RAW
## gradient. New: the fit's own verdict. Two fits, one where grad_tol is
## tightened and one where it is loosened.
set.seed(9203)
n <- 20000
d3 <- data.frame(x = rnorm(n))
d3$y <- rnorm(n, 1 + 2 * d3$x, 1)
for (gt in c(1e-3, 1e-8, 1e-1)) {
  r <- gw(frm(bf(y ~ x), family = gaussian(), data = d3,
              control = frmtmb_control(grad_tol = gt)))
  out <- capture.output(diagnose(r$fit))
  dd <- diagnose(r$fit, quiet = TRUE)
  cat("B3 grad_tol =", format(gt), " max_grad =",
      format(dd$max_grad, digits = 4),
      " headroom =", format(dd$grad_headroom %||% NA, digits = 4),
      " GRADWARN =", r$grad,
      " clean line =", any(grepl("No convergence problems", out)), "\n")
}
cat("\n")

## ---- B4: a benign indefinite Hessian. A variance component truly at
## zero, and a student-t whose nu runs off the top of its link. Both are
## correct fits whose unconstrained curvature can be unusable.
set.seed(9204)
ng <- 60
d4 <- data.frame(g = factor(rep(seq_len(ng), 40)))
d4$x <- rnorm(nrow(d4))
d4$y <- rnorm(nrow(d4), 1 + 0.7 * d4$x, 1)      # NO group effect at all
r4 <- gw(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d4))
d <- diagnose(r4$fit, quiet = TRUE)
cat("B4a zero variance component, n =", nrow(d4),
    " conv =", r4$fit$opt$convergence,
    " max_grad =", format(d$max_grad, digits = 4),
    " headroom =", format(d$grad_headroom %||% NA, digits = 4),
    " pdHess =", d$pdHess, " GRADWARN =", r4$grad, "\n")

set.seed(9205)
n <- 20000
d5 <- data.frame(x = rnorm(n))
d5$y <- rnorm(n, 1 + 2 * d5$x, 1)               # no heavy tails
r5 <- gw(frm(bf(y ~ x), family = student(), data = d5))
d <- diagnose(r5$fit, quiet = TRUE)
cat("B4b student() on gaussian data, nu runs up. conv =",
    r5$fit$opt$convergence,
    " max_grad =", format(d$max_grad, digits = 4),
    " headroom =", format(d$grad_headroom %||% NA, digits = 4),
    " pdHess =", d$pdHess, " GRADWARN =", r5$grad, "\n")
if (length(r5$w)) for (x in r5$w) cat("   WARN:", x, "\n")
cat("\n")

## ---- B5: is the lane's warning set a SUBSET of the base build's? The
## code says stage 1 is the old test unchanged, so print the three
## quantities that decide it on one fit that trips.
set.seed(9206)
n <- 20000
d6 <- data.frame(x = rnorm(n), z = rnorm(n))
d6$y <- rpois(n, exp(0.3 + 0.5 * d6$x - 0.2 * d6$z))
r6 <- gw(frm(bf(y ~ x + z), family = poisson(), data = d6))
d <- diagnose(r6$fit, quiet = TRUE)
cat("B5 poisson n = 20000: max_grad =", format(d$max_grad, digits = 6),
    " grad_proj =", format(d$grad_proj %||% NA, digits = 6),
    " headroom =", format(d$grad_headroom %||% NA, digits = 6),
    " GRADWARN =", r6$grad, "\n")
cat("B5 logLik =", format(as.numeric(logLik(r6$fit)), digits = 15), "\n")
cat("B5 coef   =", paste(format(coef(r6$fit), digits = 15),
                         collapse = " "), "\n")
cat("DONE probe\n")
