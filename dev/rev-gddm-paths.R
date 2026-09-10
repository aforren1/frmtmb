# rev-gddm: paths that reach the first-row read without going through
# assemble_frame(), and the mixture() reachability claim.
#
# Seed 31. Arm from GDDM_LIB; default is this review's install of the
# worktree, so anything found here is found ON the fixed arm.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

cat("== 1. gddm_simulate(), which assembles no frame at all\n")
# The parameters are rep_len()'d to n, so a per-trial vector is the
# documented shape. gd_densities() then reads each one at the first row
# of its condition, and the condition here comes from `coh` alone.
set.seed(31)
n <- 400L
mu_vec <- rep(c(-2.5, 2.5), each = n / 2L)
a <- gddm_simulate(n, mu = mu_vec, bs = 1.5, ndt = 0.2, coh = 0,
                   control = gddm_control(t_max = 2))
cat("   mu = c(-2.5 x 200, +2.5 x 200), coh = 0 (one condition)\n")
cat("     upper-response rate, first half", mean(a$upper[1:200]),
    " second half", mean(a$upper[201:400]), "\n")
b <- gddm_simulate(n, mu = mu_vec, bs = 1.5, ndt = 0.2,
                   coh = rep(0:1, each = n / 2L),
                   control = gddm_control(t_max = 2))
cat("   the same call with coh separating the two halves\n")
cat("     upper-response rate, first half", mean(b$upper[1:200]),
    " second half", mean(b$upper[201:400]), "\n")
cat("   any error or warning from the first call: ")
w <- NULL
withCallingHandlers(
  gddm_simulate(20, mu = c(rep(-2.5, 10), rep(2.5, 10)), coh = 0,
                control = gddm_control(t_max = 2)),
  warning = function(z) { w <<- conditionMessage(z)
                          invokeRestart("muffleWarning") })
cat(if (is.null(w)) "none\n" else paste(w, "\n"))

cat("\n== 2. mixture(gddm(), gddm()), with the dpar names it wants\n")
d <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$cond <- rep(1:2, length.out = nrow(d))
d$x <- rnorm(nrow(d))
ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
sp <- bf(rt | vint(upper, cond) ~ x,
         bs1 ~ 1, ndt1 ~ 1, bias1 = 0.5,
         mu2 ~ 1, bs2 ~ 1, ndt2 ~ 1, bias2 = 0.5)
for (stage in c("frame", "objective")) {
  r <- tryCatch({
    frm(sp, family = mixture(gddm(control = ctl), gddm(control = ctl)),
        data = d, dry_run = stage)
    "ACCEPTED"
  }, error = function(e) paste("refused:",
                               substr(conditionMessage(e), 1, 100)))
  cat("   dry_run =", stage, ":", r, "\n")
}
# and the same with mu1 named explicitly
sp2 <- bf(rt | vint(upper, cond) ~ 1,
          mu1 ~ x, bs1 ~ 1, ndt1 ~ 1, bias1 = 0.5,
          mu2 ~ 1, bs2 ~ 1, ndt2 ~ 1, bias2 = 0.5)
r <- tryCatch({
  frm(sp2, family = mixture(gddm(control = ctl), gddm(control = ctl)),
      data = d, dry_run = "frame")
  "ACCEPTED"
}, error = function(e) paste("refused:",
                             substr(conditionMessage(e), 1, 100)))
cat("   mu1 named explicitly, dry_run = frame :", r, "\n")

cat("\n== 3. the check on the paths that DO assemble a frame\n")
d2 <- d
for (what in c("get_prior", "frm_simulate")) {
  r <- tryCatch({
    if (what == "get_prior") {
      get_prior(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                   bias = 0.5),
                family = gddm(control = ctl), data = d2)
    } else {
      frm_simulate(bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1,
                      bias = 0.5),
                   family = gddm(control = ctl), data = d2,
                   pars = list(), nsim = 1)
    }
    "ACCEPTED, no refusal"
  }, error = function(e) paste("refused:",
                               substr(conditionMessage(e), 1, 60)))
  cat("   ", what, ": ", r, "\n", sep = "")
}
