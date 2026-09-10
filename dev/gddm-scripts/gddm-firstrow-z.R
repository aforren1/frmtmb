# lane gddm: the random-effect (Z) half of the first-row defect.
#
# The marginal objective is the wrong instrument on its own: the subject
# deviations are exchangeable, so permuting which subject sits on a
# condition's first row leaves the Laplace likelihood unchanged whether
# or not the other rows are read, and the variance component collapses
# so the fitted vector has every b at zero for the wrong reason.
#
# What the marginal objective DOES expose is the inner solution. A
# deviation the likelihood never reads is pinned at exactly zero by its
# own prior, whatever the outer parameters are. So the probe is: put the
# outer vector somewhere the data pull hard against, evaluate the
# marginal likelihood, and read the inner b vector. A b that is bitwise
# 0 is a subject the density never saw.
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
d$s <- factor(rep(1:4, length.out = nrow(d)))
d$blk <- rep(1:2, length.out = nrow(d))
# crossing: two conditions, every subject in both
d$cond_cross <- d$blk
# nested: the condition index names the subject as well, which is what
# ?gddm asks for and what gddm_conditions() builds
d$cond_nest <- gddm_conditions(d, blk, s)

form <- bf(rt | vint(upper, cond) ~ 1 + (1 | s), bs ~ 1, ndt ~ 1,
           bias = 0.5)

probe <- function(cndcol, label) {
  dd <- d
  dd$cond <- dd[[cndcol]]
  ob <- frm(form, family = gddm(control = ctl), data = dd,
            dry_run = "objective")$obj
  p <- ob$par
  # somewhere the data pull against, so the inner problem has work to do
  p[names(p) == "beta"] <- 3.0
  p[names(p) == "theta"] <- log(0.8)
  invisible(ob$fn(p))
  b <- ob$env$last.par[names(ob$env$last.par) == "b"]
  cat("\n== ", label, ": ", length(unique(dd$cond)), " conditions\n",
      sep = "")
  cat("  subjects on a condition's first row:",
      paste(sort(unique(as.character(
        dd$s[match(sort(unique(dd$cond)), dd$cond)]))), collapse = ","),
      "\n")
  cat("  inner b:", paste(sprintf("%.17g", b), collapse = "  "), "\n")
  cat("  b bitwise zero:",
      sum(vapply(b, function(z) identical(z, 0), logical(1))), "of",
      length(b), "\n")
}

probe("cond_cross", "s crosses the condition index")
probe("cond_nest", "the condition index names s")
