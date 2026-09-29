## Reviewer re-check (b), (c), (d).
## usage: Rscript gradcheck-rev-19-recheck2.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-gradcheck"
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w)
}

## ================= (b) the fuzz grep, on REAL warnings ===============
## The worker's dev/gradcheck-15-fuzzgrep.R feeds hand-built lists to
## grad_warning_msg(). That tests the formatter, not the package: if
## grad_verdict() ever produced a shape the formatter never sees, the
## assertion would still pass. Here the messages come out of actual fits
## and the pattern is read from the harness source.
env <- new.env()
sys.source(file.path(WT, "tests/testthat/helper-fuzz.R"), envir = env)
pat <- get("FUZZ_NONCONVERGENCE", envir = env)
cat("FUZZ_NONCONVERGENCE, read from helper-fuzz.R:\n  ", pat, "\n\n")

branch_fits <- list()

## B-i plain, short, finite headroom
branch_fits$plain_finite <- function() {
  set.seed(402); n <- 1500
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rpois(n, exp(0.4 + 0.7 * dd$x - 0.4 * dd$z))
  frm(bf(y ~ x + z), family = poisson(), data = dd,
      control = frmtmb_control(restarts = 0,
                               optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                              iter.max = 1000,
                                              eval.max = 1000)))
}
## B-ii one bound held, short free set, finite headroom
branch_fits$one_bound_finite <- function() {
  set.seed(9301); n <- 3000
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rpois(n, exp(0.4 + 0.9 * dd$x - 0.6 * dd$z))
  frm(bf(y ~ x + z), family = poisson(), data = dd,
      prior = set_prior("", class = "b", coef = "x", ub = 0.1),
      control = frmtmb_control(restarts = 0,
                               optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                              iter.max = 2000,
                                              eval.max = 2000)))
}
## B-iii plain, curvature unusable
branch_fits$plain_unusable <- function() {
  set.seed(201); n <- 500
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rnorm(n, 3 + 2 * dd$x - 1.5 * dd$z, 1)
  frm(bf(y ~ x + z), family = gaussian(), data = dd,
      start = list(beta = c(50, -40, 30)),
      control = frmtmb_control(restarts = 0,
                               optCtrl = list(iter.max = 1, eval.max = 3)))
}
## B-iv two bounds held, short free set
branch_fits$two_bounds <- function() {
  set.seed(9401); n <- 2500
  dd <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z = rnorm(n))
  dd$y <- rpois(n, exp(0.3 + 0.9 * dd$x1 - 0.9 * dd$x2 + 0.5 * dd$z))
  frm(bf(y ~ x1 + x2 + z), family = poisson(), data = dd,
      prior = c(set_prior("", class = "b", coef = "x1", ub = 0.1),
                set_prior("", class = "b", coef = "x2", lb = -0.1)),
      control = frmtmb_control(restarts = 0,
                               optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                                              iter.max = 2000,
                                              eval.max = 2000)))
}
## B-v a bound held AND the curvature unusable
branch_fits$bound_unusable <- function() {
  set.seed(9402); n <- 800
  dd <- data.frame(x = rnorm(n), z = rnorm(n))
  dd$y <- rnorm(n, 3 + 2 * dd$x - 1.5 * dd$z, 1)
  frm(bf(y ~ x + z), family = gaussian(), data = dd,
      prior = set_prior("", class = "b", coef = "x", ub = 0.1),
      start = list(beta = c(50, 0.1, 30)),
      control = frmtmb_control(restarts = 0,
                               optCtrl = list(iter.max = 1, eval.max = 3)))
}
allok <- TRUE
seen <- character(0)
for (nm in names(branch_fits)) {
  r <- gw(branch_fits[[nm]]())
  gm <- grep("Large maximum absolute gradient", r$w, value = TRUE)
  if (!length(gm)) {
    cat(sprintf("%-18s NO GRADIENT WARNING RAISED (branch not reached)\n",
                nm))
    next
  }
  hit <- grepl(pat, gm)
  allok <- allok && all(hit)
  d <- diagnose(r$fit, quiet = TRUE)
  shape <- paste0("held=", length(d$grad_bound_held %||% character(0)),
                  " headroom=",
                  if (is.finite(d$grad_headroom %||% NA)) "finite" else "NA")
  seen <- c(seen, shape)
  cat(sprintf("%-18s grep %-5s  [%s]\n", nm, all(hit), shape))
  cat("   ", gm, "\n\n")
}
cat("every REAL warning matched FUZZ_NONCONVERGENCE:", allok, "\n")
cat("distinct message shapes reached from real fits:",
    length(unique(seen)), "of 4 (held 0/1+ x headroom finite/NA)\n")
## the ABSENT case: a pattern that must NOT match, so the grep above is
## not passing on something that matches everything
cat("control, pattern against an unrelated sentence:",
    grepl(pat, "Model fitted successfully with no problems"), "\n\n")

## ================= (c) importance withholding =========================
set.seed(506)
ng <- 60
d4 <- data.frame(g = factor(rep(seq_len(ng), 25)))
d4$x <- rnorm(nrow(d4))
re <- rnorm(ng, 0, 0.6)
d4$y <- rpois(nrow(d4), exp(0.3 + 0.5 * d4$x + re[d4$g]))
for (rt in c(1e-3, 1e-2, 1e-1)) {
  r <- gw(frm(bf(y ~ x + (1 | g)), family = poisson(), data = d4,
              importance = 500L,
              control = frmtmb_control(
                restarts = 0,
                optCtrl = list(rel.tol = rt, x.tol = rt,
                               iter.max = 2000, eval.max = 2000))))
  d <- diagnose(r$fit, quiet = TRUE)
  out <- capture.output(diagnose(r$fit))
  cat(sprintf(paste0("(c) importance rel.tol=%.0e  max_grad=%-10s",
                     " grad_proj=%-6s proj_par=%-6s headroom=%-6s",
                     " held=%-6s clean=%-5s nwarn=%d\n"),
              rt, format(d$max_grad, digits = 4),
              format(d$grad_proj %||% NA),
              format(d$grad_proj_par %||% NA),
              format(d$grad_headroom %||% NA),
              format(length(d$grad_bound_held %||% character(0))),
              any(grepl("No convergence problems", out)), length(r$w)))
  cat("     any 'Newton step' line printed:",
      any(grepl("Newton step", out)), "\n")
  cat("     any 'held by a bound' line printed:",
      any(grepl("held by a bound", out)), "\n")
  cat("     verdict cached on the fit:",
      !is.null(r$fit$cache$grad_verdict), "\n")
  for (w in r$w) cat("     WARN:", substr(w, 1, 100), "\n")
}
## and the same model WITHOUT importance, so the withholding is shown to
## be about the importance branch and not about this design
r0 <- gw(frm(bf(y ~ x + (1 | g)), family = poisson(), data = d4,
             control = frmtmb_control(
               restarts = 0,
               optCtrl = list(rel.tol = 1e-1, x.tol = 1e-1,
                              iter.max = 2000, eval.max = 2000))))
d0 <- diagnose(r0$fit, quiet = TRUE)
out0 <- capture.output(diagnose(r0$fit))
cat(sprintf(paste0("(c) SAME design, no importance: max_grad=%s",
                   " headroom=%s clean=%s nwarn=%d newton_line=%s\n"),
            format(d0$max_grad, digits = 4),
            format(d0$grad_headroom %||% NA, digits = 4),
            any(grepl("No convergence problems", out0)), length(r0$w),
            any(grepl("Newton step", out0))))
cat("\n")

## ================= (d) the quadraticity yardstick =====================
set.seed(402)
n <- 1500
dq <- data.frame(x = rnorm(n), z = rnorm(n))
dq$y <- rpois(n, exp(0.4 + 0.7 * dq$x - 0.4 * dq$z))
rq <- gw(frm(bf(y ~ x + z), family = poisson(), data = dq,
             control = frmtmb_control(
               restarts = 0,
               optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                              iter.max = 1000, eval.max = 1000))))
d <- diagnose(rq$fit, quiet = TRUE)
p <- rq$fit$opt$par
g <- drop(rq$fit$obj$gr(p))
H <- rq$fit$obj$he(p)
ch <- chol((H + t(H)) / 2)
step <- -backsolve(ch, backsolve(ch, g, transpose = TRUE))
f0 <- as.numeric(rq$fit$obj$fn(p))
realized <- f0 - as.numeric(rq$fit$obj$fn(p + step))
half <- f0 - as.numeric(rq$fit$obj$fn(p + 0.5 * step))
nonquad <- abs(half / realized / 0.75 - 1)
yard <- max(nonquad, sqrt(.Machine$double.eps))
err <- abs(d$grad_headroom / realized - 1)
cat("(d) realized =", format(realized, digits = 10), "\n")
cat("(d) half-step drop =", format(half, digits = 10),
    " half/realized =", format(half / realized, digits = 10),
    " (0.75 exactly for a quadratic)\n")
cat("(d) nonquad (MEASURED from the fit) =", format(nonquad, digits = 6),
    "\n")
cat("(d) yard =", format(yard, digits = 6),
    " (floor sqrt(eps) =", format(sqrt(.Machine$double.eps), digits = 4),
    "active:", yard > nonquad, ")\n")
cat("(d) headroom error =", format(err, digits = 6),
    " error/yard =", format(err / yard, digits = 6),
    " assertion (< 10) passes:", err / yard < 10, "\n")
## is it discriminating? how wrong could the headroom be and still pass?
cat("(d) the assertion admits a headroom wrong by up to a factor of",
    format(1 + 10 * yard, digits = 6), "\n")
## a third measurement of the same quantity, at t = 0.25, so the
## quadratic law is checked rather than assumed
q25 <- f0 - as.numeric(rq$fit$obj$fn(p + 0.25 * step))
cat("(d) quarter-step drop / realized =", format(q25 / realized, digits = 8),
    " (0.4375 exactly for a quadratic)\n")
cat("DONE recheck-bcd\n")
