# lane gddm: the check runs at FRAME ASSEMBLY, and frame assembly runs
# again for every post-fit path that takes new rows. A refusal that
# fires on predict() or on conditional_effects() would be a false alarm
# in the worst place, so each path is run on a CORRECT model.
#
# Seed 11. Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

set.seed(11)
ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
n <- 200L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$coh <- factor(rep(c("lo", "hi"), length.out = n))
d$blk <- factor(rep(1:2, each = 2L, length.out = n))
d$cond <- gddm_conditions(d, coh, blk)
cat("conditions:", length(unique(d$cond)), "\n")

fit <- frm(bf(rt | vint(upper, cond) ~ coh, bs ~ blk, ndt ~ 1,
              bias = 0.5),
           family = gddm(control = ctl), data = d)
cat("logLik:", format(as.numeric(logLik(fit))), "\n")

try_path <- function(label, expr) {
  r <- tryCatch({
    v <- force(expr)
    paste0("ok (", paste(utils::head(format(as.numeric(unlist(v))), 2L),
                         collapse = ", "), " ...)")
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat("  ", label, ": ", substr(r, 1, 160), "\n", sep = "")
}

cat("\n== post-fit paths on a correct model\n")
try_path("fitted()", fitted(fit))
try_path("residuals()", residuals(fit))
try_path("predict(type = link)", predict(fit, type = "link"))
try_path("simulate(nsim = 2)", simulate(fit, nsim = 2))

nd <- unique(d[, c("coh", "blk", "cond")])
nd$rt <- 0.6
nd$upper <- 1L
try_path("predict(newdata = one row per condition)",
         predict(fit, newdata = nd, type = "link"))

# The shape that could bite: newdata whose condition index no longer
# separates the covariate. A user who builds newdata by hand and
# forgets to rebuild the index gets the same refusal, which is right.
nd2 <- nd
nd2$cond <- 1L
try_path("predict(newdata = index collapsed by hand)",
         predict(fit, newdata = nd2, type = "link"))

try_path("conditional_effects()", conditional_effects(fit))
try_path("influence()", tryCatch(influence(fit),
                                 error = function(e) stop(e)))
