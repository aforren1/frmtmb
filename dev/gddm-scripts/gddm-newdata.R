# lane gddm: does the guard reach the NEWDATA path, and does it need to?
#
# The refusal is registered as a frame check, and frmtmb runs frame
# checks in assemble_frame(). predict(newdata = ) rebuilds a frame, so
# the question is whether the check runs there. It does not, and this
# says what that costs: a link-scale prediction never touches the
# density, but a response-scale one goes through post$mean_fn into
# gd_densities(), which reads the condition's first row.
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
fit <- frm(bf(rt | vint(upper, cond) ~ coh, bs ~ blk, ndt ~ 1,
              bias = 0.5),
           family = gddm(control = ctl), data = d)

nd <- unique(d[, c("coh", "blk", "cond")])
nd$rt <- 0.6
nd$upper <- 1L
nd <- nd[order(nd$cond), ]
bad <- nd
bad$cond <- 1L

show <- function(label, e) {
  r <- tryCatch(paste(format(as.numeric(unlist(e)), digits = 10),
                      collapse = "  "),
                error = function(z) paste("REFUSED:",
                                          conditionMessage(z)))
  cat("  ", label, ": ", substr(r, 1, 140), "\n", sep = "")
}

cat("\n== newdata, four conditions correctly indexed\n")
show("predict link", predict(fit, newdata = nd, type = "link"))
show("predict response", predict(fit, newdata = nd, type = "response"))

cat("\n== the same four rows with the index collapsed to one\n")
show("predict link", predict(fit, newdata = bad, type = "link"))
show("predict response", predict(fit, newdata = bad, type = "response"))
cat("  (a link prediction never touches the density; a response one",
    "does, so the collapsed index gives every row the first row's",
    "mean)\n")
