# lane gddm: reproduce the first-row defect on the BEFORE arm.
#
# gd_densities() reads every distributional parameter at the FIRST ROW
# of its condition, so a dpar that varies WITHIN a condition never
# reaches the density. Checked by construction rather than by reading
# the source: the data of every row but the first of a condition is
# changed and the objective is evaluated at the SAME parameter vector.
# A likelihood that does not move never read those rows.
#
# The probe is taken at a FITTED parameter vector. At the starting
# values every regression coefficient is zero, so the linear predictor
# does not depend on the covariate and the probe would pass for the
# wrong reason.
#
# Seed 5. Arm chosen by GDDM_LIB; default is the shared reference build
# of the base commit.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rellib-0552")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "\n")
cat("frmtmb", format(packageVersion("frmtmb")),
    " frmtmb.eam", format(packageVersion("frmtmb.eam")), "\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(5)
d <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$cond <- rep(1:2, length.out = nrow(d))
d$x <- rnorm(nrow(d))
d$s <- factor(rep(1:4, length.out = nrow(d)))
first <- vapply(split(seq_len(nrow(d)), d$cond), min, integer(1))

obj_of <- function(dat, form) {
  frm(form, family = gddm(control = ctl), data = dat,
      dry_run = "objective")$obj
}

probe <- function(form, label) {
  fit <- try(frm(form, family = gddm(control = ctl), data = d),
             silent = TRUE)
  cat("\n== ", label, "\n", sep = "")
  if (inherits(fit, "try-error")) {
    cat("  REFUSED: ", conditionMessage(attr(fit, "condition")), "\n",
        sep = "")
    return(invisible(NULL))
  }
  p <- fit$opt$par
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
probe(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1, bias ~ x),
      "bias ~ x")
probe(bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
         bias = 0.5),
      "mu ~ 1 + (1 | s), a random effect crossing conditions")

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

# The control: the guarded thing PRESENT. A vreal() covariate that
# varies within a condition is already refused, by the check this one
# is modelled on.
dv <- d
# NOT on row parity: the condition index alternates on row parity
# too, so a covariate that did would be constant WITHIN each
# condition and the check would be right not to fire. The first
# version of this script made that mistake and read the silence as
# the check failing.
dv$coh <- rep(c(0.3, 0.3, 0.7, 0.7), length.out = nrow(dv))
r <- tryCatch(
  frm(bf(rt | vint(upper, cond) + vreal(coh) ~ 1, bs ~ 1, ndt ~ 1,
         bias = 0.5),
      family = gddm(drift = gddm_drift_coherence(), control = ctl),
      data = dv, dry_run = "objective"),
  error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("\nvreal() varying within a condition:",
    if (is.character(r)) r else "FITTED, no error", "\n")

# And the same vreal() held constant within each condition, which must
# NOT be refused.
dv2 <- d
dv2$coh <- ifelse(dv2$cond == 1L, 0.3, 0.7)
r2 <- tryCatch({
  frm(bf(rt | vint(upper, cond) + vreal(coh) ~ 1, bs ~ 1, ndt ~ 1,
         bias = 0.5),
      family = gddm(drift = gddm_drift_coherence(), control = ctl),
      data = dv2, dry_run = "objective")
  "accepted"
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("vreal() constant within a condition:", r2, "\n")
