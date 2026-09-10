# lane gddm: the Z half, probed on the JOINT objective.
#
# The marginal (Laplace) objective is the wrong instrument for this: the
# subject deviations are integrated out and are exchangeable, so
# permuting which subject sits on a condition's first row leaves the
# marginal likelihood unchanged whether or not the other rows are read.
# The joint objective takes b as an argument, so a b that is set rather
# than integrated shows whether rows other than the first ever see it.
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
lp <- ob$env$last.par
cat("joint parameter slots:", paste(names(lp), collapse = " "), "\n")
lp[names(lp) == "beta"] <- 2.0
lp[names(lp) == "theta"] <- log(0.7)
bi <- which(names(lp) == "b")
lp[bi] <- c(-0.8, 0.4, 0.9, -0.5)[seq_along(bi)]
cat("b entries:", length(bi), "\n")

joint <- function(o, v) as.numeric(o$env$f(v, order = 0))
f0 <- joint(ob, lp)
d2 <- d
j <- setdiff(seq_len(nrow(d)), first)
d2$s[j] <- factor("1", levels = levels(d$s))
f1 <- joint(obj_of(d2), lp)
d3 <- d
d3$s[first] <- factor(c("2", "3"), levels = levels(d$s))
f2 <- joint(obj_of(d3), lp)

cat("\n== mu ~ 1 + (1 | s), joint objective, b set\n")
cat(sprintf("  joint, data as drawn                       %.17g\n", f0))
cat(sprintf("  joint, s collapsed on %d of %d rows (not first) %.17g  ",
            length(j), nrow(d), f1))
cat("identical:", identical(f0, f1), "\n")
cat(sprintf("  joint, s relabelled on the %d FIRST rows only  %.17g  ",
            length(first), f2))
cat("moved:", !identical(f0, f2), "\n")
