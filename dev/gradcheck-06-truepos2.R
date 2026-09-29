# True positives the OPTIMIZER'S OWN CODE does not catch.
#
# Every fit in gradcheck-02-truepos.R stops at its iteration cap and so
# already raises "Optimizer did not report convergence". A gradient check
# earns its place only on fits nlminb calls converged, so those are
# constructed here.
#
#   Rscript dev/gradcheck-06-truepos2.R base|lane

args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[[1L]] else "base"
source("dev/gradcheck-helpers.R")
gc_libpaths(which_lib)
library(frmtmb)
cat("library:", which_lib, " frmtmb", format(packageVersion("frmtmb")),
    "\n\n")

rows <- list()
add <- function(r) {
  rows[[length(rows) + 1L]] <<- r
  gc_print(r)
}
hr <- function(s) cat("\n== ", s, "\n", sep = "")

## ------------------------------------------------------------------
hr("U1. a badly scaled predictor, autoscale off")
# The case R/autoscale.R exists for. nlminb reports code 0 and the fit
# is far from the optimum: the shortfall against the standardized fit is
# the evidence.
set.seed(401)
n <- 2000
dd <- data.frame(x = rnorm(n))
dd$xbig <- dd$x * 1e6
dd$y <- rnorm(n, 1 + 2e-6 * dd$xbig, 1)
good <- frm(bf(y ~ x), family = gaussian(), data = dd)
cp <- gc_catch(frm(bf(y ~ xbig), family = gaussian(), data = dd,
                   control = frmtmb_control(autoscale = FALSE)))
f <- cp$value
add(gc_row("scaled 1e6 predictor", 401, f, cp$warnings))
cat("  logLik shortfall against the well-scaled fit:",
    format(as.numeric(logLik(good)) - as.numeric(logLik(f)),
           digits = 6), "\n")

## ------------------------------------------------------------------
hr("U2. nlminb's own tolerances loosened: code 0, not at the optimum")
set.seed(402)
n <- 1500
d2 <- data.frame(x = rnorm(n), z = rnorm(n))
d2$y <- rpois(n, exp(0.4 + 0.7 * d2$x - 0.4 * d2$z))
ref2 <- frm(bf(y ~ x + z), family = poisson(), data = d2)
for (rt in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  cp <- gc_catch(frm(bf(y ~ x + z), family = poisson(), data = d2,
                     control = frmtmb_control(
                       restarts = 0,
                       optCtrl = list(rel.tol = rt, x.tol = rt,
                                      iter.max = 1000, eval.max = 1000))))
  f <- cp$value
  add(gc_row(paste0("poisson, rel.tol = ", rt), 402, f, cp$warnings))
  cat("  logLik shortfall:",
      format(as.numeric(logLik(ref2)) - as.numeric(logLik(f)),
             digits = 6), "\n")
}

## ------------------------------------------------------------------
hr("U3. the same, on a GLMM (Laplace, so the Hessian is numeric)")
set.seed(403)
n <- 1200
d3 <- data.frame(x = rnorm(n), g = factor(rep(1:40, 30)))
d3$y <- rpois(n, exp(0.4 + 0.6 * d3$x + rnorm(40, 0, 0.7)[d3$g]))
ref3 <- frm(bf(y ~ x + (1 | g)), family = poisson(), data = d3)
for (rt in c(1e-2, 1e-3, 1e-4)) {
  cp <- gc_catch(frm(bf(y ~ x + (1 | g)), family = poisson(), data = d3,
                     control = frmtmb_control(
                       restarts = 0,
                       optCtrl = list(rel.tol = rt, x.tol = rt,
                                      iter.max = 1000, eval.max = 1000))))
  f <- cp$value
  add(gc_row(paste0("poisson GLMM, rel.tol = ", rt), 403, f, cp$warnings))
  cat("  logLik shortfall:",
      format(as.numeric(logLik(ref3)) - as.numeric(logLik(f)),
             digits = 6), "\n")
}

## ------------------------------------------------------------------
hr("U4. a near-collinear design: a flat ridge nlminb stops on")
set.seed(404)
n <- 1000
x1 <- rnorm(n)
d4 <- data.frame(x1 = x1, x2 = x1 + rnorm(n, 0, 1e-5))
d4$y <- rnorm(n, 1 + 2 * d4$x1, 1)
cp <- gc_catch(frm(bf(y ~ x1 + x2), family = gaussian(), data = d4,
                   control = frmtmb_control(restarts = 0)))
f <- cp$value
add(gc_row("near-collinear ridge", 404, f, cp$warnings))
cat("  condition number of the design:",
    format(kappa(cbind(1, d4$x1, d4$x2)), digits = 4), "\n")

## ------------------------------------------------------------------
hr("U5. a genuinely flat direction: a bump off its own support")
set.seed(405)
nw <- 45
pd <- data.frame(w = seq_len(nw))
pd$y <- rnorm(nw, 2, 0.5)
cp <- gc_catch(frm(
  bf(y ~ exp(lamp) * exp(-0.5 * ((w - pk) / exp(lsig))^2),
     lamp + pk + lsig ~ 1, nl = TRUE),
  family = gaussian(), data = pd))
f <- cp$value
if (inherits(f, "gc_error")) {
  cat("  ERROR:", f, "\n")
} else {
  add(gc_row("flat nonlinear peak", 405, f, cp$warnings))
  for (w in cp$warnings) cat("   -", substr(w, 1, 110), "\n")
}

## ------------------------------------------------------------------
hr("U6. complete separation in a bernoulli fit")
set.seed(406)
n <- 200
d6 <- data.frame(x = rnorm(n))
d6$y <- as.integer(d6$x > 0)
cp <- gc_catch(frm(bf(y ~ x), family = bernoulli(), data = d6,
                   control = frmtmb_control(restarts = 0)))
f <- cp$value
add(gc_row("complete separation", 406, f, cp$warnings))

## ------------------------------------------------------------------
hr("U7. how far can a CORRECT fit's headroom go? bigger n")
for (spec in list(list("gaussian", 1e5), list("gaussian", 4e5),
                  list("nearzero", 1e5), list("nearzero", 4e5))) {
  nn <- as.integer(spec[[2L]])
  set.seed(407)
  if (identical(spec[[1L]], "gaussian")) {
    d <- data.frame(x = rnorm(nn), z = rnorm(nn))
    d$y <- rnorm(nn, 1 + 0.8 * d$x - 0.4 * d$z, 1.5)
    cp <- gc_catch(frm(bf(y ~ x + z), family = gaussian(), data = d))
  } else {
    q <- nn %/% 25L
    d <- data.frame(x = rnorm(nn), g = factor(rep_len(seq_len(q), nn)))
    d$y <- rnorm(nn, 1 + 0.6 * d$x, 1)
    cp <- gc_catch(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d))
  }
  f <- cp$value
  if (inherits(f, "gc_error")) {
    cat("  ERROR:", f, "\n"); next
  }
  add(gc_row(paste0(spec[[1L]], " n = ", nn), 407, f, cp$warnings))
}

saveRDS(rows, file.path("dev", paste0("gradcheck-06-", which_lib, ".rds")))
cat("\nrows:", length(rows), "  warned:",
    sum(vapply(rows, function(r) isTRUE(r$warned), TRUE)), "\n")
