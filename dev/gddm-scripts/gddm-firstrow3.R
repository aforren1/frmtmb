# lane gddm: the random-effect (Z) half of the first-row defect, at a
# parameter vector where the random effect is NOT zero.
#
# The fitted vector is useless for this probe: a subject deviation that
# only ever enters through one row per condition is unidentified, so the
# variance component collapses to 7.6e-11 and every b is zero. Probing
# there would report "no move" for the wrong reason, the same way
# probing at the start values does for a fixed effect. So the b vector
# and its log standard deviation are SET, and the objective is evaluated
# at that vector.
#
# Seed 5. Arm chosen by GDDM_LIB.

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
d$s <- factor(rep(1:4, length.out = nrow(d)))
first <- vapply(split(seq_len(nrow(d)), d$cond), min, integer(1))

form <- bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
           bias = 0.5)
obj_of <- function(dat) frm(form, family = gddm(control = ctl),
                            data = dat, dry_run = "objective")$obj
ob <- obj_of(d)
p <- ob$par
cat("parameter slots:", paste(names(p), collapse = " "), "\n")
p[names(p) == "beta"] <- c(2.0, 1.0, 0.3)[seq_len(sum(names(p) == "beta"))]
p[names(p) == "b"] <- c(-0.8, 0.4, 0.9, -0.5)[seq_len(sum(names(p) == "b"))]
p[names(p) == "theta"] <- log(0.7)
cat("probe vector:", paste(sprintf("%s=%.3g", names(p), p),
                           collapse = " "), "\n")

f0 <- as.numeric(ob$fn(p))
d2 <- d
j <- setdiff(seq_len(nrow(d)), first)
d2$s[j] <- factor("1", levels = levels(d$s))
f1 <- as.numeric(obj_of(d2)$fn(p))
d3 <- d
d3$s[first] <- factor(c("2", "3"), levels = levels(d$s))
f2 <- as.numeric(obj_of(d3)$fn(p))

cat("\n== mu ~ 1 + (1 | s), s crossing conditions, b set nonzero\n")
cat(sprintf("  fn, data as drawn                          %.17g\n", f0))
cat(sprintf("  fn, s collapsed on %d of %d rows (not first) %.17g  ",
            length(j), nrow(d), f1))
cat("identical:", identical(f0, f1), "\n")
cat(sprintf("  fn, s relabelled on the %d FIRST rows only   %.17g  ",
            length(first), f2))
cat("moved:", !identical(f0, f2), "\n")
