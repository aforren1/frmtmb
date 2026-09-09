# rev-ndt: does gddm() really read every dpar at the FIRST ROW of its
# condition, and is that a silent wrong answer for dpars other than
# `ndt`?
#
# The lane filed this as a defect it found and did not fix. It is
# checked here BY CONSTRUCTION rather than by reading gd_densities():
# the data of every row but the first of a condition is changed and the
# objective is evaluated at the SAME parameter vector. A likelihood that
# does not move is a likelihood that never read those rows.
#
# The probe is taken at a FITTED parameter vector, not at the start: at
# the start every regression coefficient is zero, so the linear
# predictor does not depend on the covariate and the probe would pass
# for the wrong reason. That is the first version of this script.
#
# Seed 5, worktree build.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(5)
d <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$cond <- rep(1:2, length.out = nrow(d))
d$x <- rnorm(nrow(d))
first <- vapply(split(seq_len(nrow(d)), d$cond), min, integer(1))

obj_of <- function(dat, form) {
  frm(form, family = gddm(control = ctl), data = dat,
      dry_run = "objective")$obj
}

probe <- function(form, label) {
  fit <- frm(form, family = gddm(control = ctl), data = d)
  p <- fit$opt$par
  cat("\n== ", label, "\n", sep = "")
  cat("  fitted coefficients:",
      paste(names(unlist(fixef(fit))),
            sprintf("%.6g", unlist(fixef(fit))), sep = "=",
            collapse = "  "), "\n")
  f0 <- as.numeric(obj_of(d, form)$fn(p))
  d2 <- d
  j <- setdiff(seq_len(nrow(d)), first)
  d2$x[j] <- d2$x[j] + 100
  f1 <- as.numeric(obj_of(d2, form)$fn(p))
  d3 <- d
  d3$x[first] <- d3$x[first] + 0.01
  f2 <- as.numeric(obj_of(d3, form)$fn(p))
  cat(sprintf("  fn, data as drawn                      %.17g\n", f0))
  cat(sprintf("  fn, x + 100 on %d of %d rows (not first) %.17g  ",
              length(j), nrow(d), f1))
  cat("identical:", identical(f0, f1), "\n")
  cat(sprintf("  fn, x + 0.01 on the %d FIRST rows only   %.17g  ",
              length(first), f2))
  cat("moved:", !identical(f0, f2), "\n")
  invisible(fit)
}

probe(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1, bias = 0.5),
      "mu ~ x")
probe(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ x, bias = 0.5),
      "ndt ~ x")
probe(bf(rt | vint(upper, cond) ~ 1, bs ~ x, ndt ~ 1, bias = 0.5),
      "bs ~ x")

# Does anything refuse or warn?
w <- NULL
form <- bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1, bias = 0.5)
m <- withCallingHandlers(
  tryCatch(frm(form, family = gddm(control = ctl), data = d),
           error = function(e) paste("ERR:", conditionMessage(e))),
  warning = function(z) { w <<- c(w, conditionMessage(z))
                          invokeRestart("muffleWarning") })
cat("\nfit of a within-condition-varying dpar:",
    if (is.character(m)) m else "FITTED, no error\n")
cat("warnings:", if (is.null(w)) "NONE" else paste(w, collapse = " / "),
    "\n")

# The ridge: move mu.x and compensate the intercept at the first rows'
# x alone. If the solver only ever sees the first row of a condition,
# the likelihood cannot tell the two apart.
ob <- obj_of(d, form)
p <- unname(m$opt$par)
nm <- names(m$opt$par)
ib <- which(nm == "beta")
cat("\nridge probe, beta slots:", paste(ib, collapse = " "), "\n")
x1 <- d$x[first[1L]]
x2 <- d$x[first[2L]]
cat(sprintf("  x at the two first rows: %.6f  %.6f\n", x1, x2))
q <- p
q[ib[2L]] <- p[ib[2L]] + 0.5
q[ib[1L]] <- p[ib[1L]] - 0.5 * x1
cat(sprintf("  fn(p) %.17g\n  fn(q) %.17g\n  gap %.3g\n",
            as.numeric(ob$fn(p)), as.numeric(ob$fn(q)),
            abs(as.numeric(ob$fn(p)) - as.numeric(ob$fn(q)))))
cat("  (exact compensation at condition 1's first row only; a solver",
    "that read every row could not match this)\n")
