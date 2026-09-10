# lane gddm: two probes the first script got wrong, redone.
#
# 1. The random-effect case. The first script perturbed `x` in a model
#    that has no `x` in it, so "no move" was the wrong reason. Here the
#    subject codes themselves are permuted on every row but the first
#    of each condition.
# 2. The vreal() control. The first script's coherence alternated on
#    row parity and the condition index also alternates on row parity,
#    so the covariate WAS constant within condition and the existing
#    check was right not to fire. Here it varies inside a condition.
#
# Seed 5, the same design as gddm-firstrow.R. Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rellib-0552")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(5)
d <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$cond <- rep(1:2, length.out = nrow(d))
d$x <- rnorm(nrow(d))
d$s <- factor(rep(1:4, length.out = nrow(d)))
first <- vapply(split(seq_len(nrow(d)), d$cond), min, integer(1))

# A random effect on mu, grouped by a factor that crosses conditions.
form <- bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
           bias = 0.5)
fit <- frm(form, family = gddm(control = ctl), data = d)
p <- fit$opt$par
obj_of <- function(dat) frm(form, family = gddm(control = ctl),
                            data = dat, dry_run = "objective")$obj
f0 <- as.numeric(obj_of(d)$fn(p))
d2 <- d
j <- setdiff(seq_len(nrow(d)), first)
# every non-first row moved to one single subject level
d2$s[j] <- factor("1", levels = levels(d$s))
f1 <- as.numeric(obj_of(d2)$fn(p))
d3 <- d
d3$s[first] <- factor(c("2", "3"), levels = levels(d$s))
f2 <- as.numeric(obj_of(d3)$fn(p))
cat("\n== mu ~ 1 + (1 | s), s crossing conditions\n")
cat(sprintf("  fn, data as drawn                        %.17g\n", f0))
cat(sprintf("  fn, s collapsed on %d of %d rows (not first) %.17g  ",
            length(j), nrow(d), f1))
cat("identical:", identical(f0, f1), "\n")
cat(sprintf("  fn, s relabelled on the %d FIRST rows only  %.17g  ",
            length(first), f2))
cat("moved:", !identical(f0, f2), "\n")
cat("  random-effect sd:", format(unlist(VarCorr(fit))), "\n")

# The vreal() control, done properly this time: coherence varies WITHIN
# each condition, so the existing tapply check must refuse.
dv <- d
dv$coh <- rep(c(0.3, 0.3, 0.7, 0.7), length.out = nrow(dv))
cat("\ncoherence by condition:\n")
print(table(dv$cond, dv$coh))
r <- tryCatch(
  frm(bf(rt | vint(upper, cond) + vreal(coh) ~ 1, bs ~ 1, ndt ~ 1,
         bias = 0.5),
      family = gddm(drift = gddm_drift_coherence(), control = ctl),
      data = dv, dry_run = "objective"),
  error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("vreal() varying within a condition:",
    if (is.character(r)) r else "FITTED, no error", "\n")
