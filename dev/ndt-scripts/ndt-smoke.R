# ndt lane, control: with NO ndt_group() the new parameterization is the
# OLD one written differently, so a small fit must reproduce the
# reference build's log-likelihood and estimates.
#
# Run against BOTH libraries:
#   Rscript --vanilla dev/ndt-scripts/ndt-smoke.R new
#   Rscript --vanilla dev/ndt-scripts/ndt-smoke.R ref
# Seed 4242.

a <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(a)) a[1L] else "new"
.libPaths(if (identical(which_lib, "ref")) {
  c("C:/Users/adf44/source/r/reflib-r2",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
} else {
  c("C:/Users/adf44/source/r/ndt-lib",
    "C:/Users/adf44/AppData/Local/R/win-library/4.6")
})
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("arm:", which_lib, " frmtmb.eam",
    as.character(utils::packageVersion("frmtmb.eam")),
    "from", dirname(dirname(getNamespaceInfo("frmtmb.eam", "path"))), "\n")

set.seed(4242)
ns <- 6L
nt <- 60L
u <- rnorm(ns, 0, 0.12)
s <- rep(seq_len(ns), each = nt)
cond <- rep(rep(0:1, each = nt / 2L), times = ns)
d <- ddm_simulate(ns * nt, mu = 0.4 + 0.9 * cond, bs = 1.4,
                  ndt = 0.25 * exp(u[s]), bias = 0.5)
d$s <- factor(s)
d$cond <- factor(cond, labels = c("a", "b"))

fit <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
           family = wiener(), data = d)
cat(sprintf("plain      logLik %.10f  conv %d\n",
            as.numeric(logLik(fit)), fit$opt$convergence))
b <- unlist(fixef(fit))
cat("plain      fixef  ",
    paste(names(b), formatC(b, digits = 10, format = "g"), sep = "=",
          collapse = "  "), "\n")
cat("plain      ndt link name:", family(fit)$links$ndt$name, "\n")
cat("plain      predict ndt   :",
    formatC(predict(fit, dpar = "ndt", type = "response")[1],
            digits = 10, format = "g"), "\n")
cat("plain      fitted[1]     :",
    formatC(fitted(fit)[1], digits = 10, format = "g"), "\n")
cat("plain      min(rt)       :", formatC(min(d$rt), digits = 10,
                                          format = "g"), "\n")

# the sv arm, which is the one that converged to the wall
fit2 <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
            family = wiener(variability = "sv"), data = d)
cat(sprintf("sv         logLik %.10f  conv %d\n",
            as.numeric(logLik(fit2)), fit2$opt$convergence))

# max_ndt, which under the new scheme still sets one bound for every row
fit3 <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
            family = wiener(max_ndt = 0.15), data = d)
cat(sprintf("max_ndt    logLik %.10f  conv %d\n",
            as.numeric(logLik(fit3)), fit3$opt$convergence))

# a race family, to show the shared helper did not move rdm or lba
dr <- rdm_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
fr <- frm(bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
          family = rdm(2), data = dr)
cat(sprintf("rdm        logLik %.10f  conv %d\n",
            as.numeric(logLik(fr)), fr$opt$convergence))
dl <- lba_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.4, ndt = 0.2)
fl <- frm(bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
          family = lba(2), data = dl)
cat(sprintf("lba        logLik %.10f  conv %d\n",
            as.numeric(logLik(fl)), fl$opt$convergence))
dg <- wiener_gng_simulate(400, mu = 1.2, bs = 1.4, ndt = 0.25,
                          deadline = 1.5)
fg <- frm(bf(rt | dec(responded) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = wiener_gng(deadline = 1.5), data = dg)
cat(sprintf("wiener_gng logLik %.10f  conv %d\n",
            as.numeric(logLik(fg)), fg$opt$convergence))
