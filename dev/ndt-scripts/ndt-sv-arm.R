# ndt lane: the eam-sv arm's condition effect, with and without the
# per-subject bound, on the SAME data. The tier's sv row now misses the
# truth on `mu.condb` and the question is whether the bound moved it.
#
# Seed 20260908, the tier's own. 30 x 400, sv = 0.4.
# Run: Rscript --vanilla dev/ndt-scripts/ndt-sv-arm.R
.libPaths(c("C:/Users/adf44/source/r/ndt-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
           sd_mu = 0.35, sd_lbs = 0.20, sd_lndt = 0.12, sv = 0.4)
ns <- 30L
nt <- 400L
set.seed(20260908L)
u_mu <- rnorm(ns, 0, tr$sd_mu)
u_bs <- rnorm(ns, 0, tr$sd_lbs)
u_nd <- rnorm(ns, 0, tr$sd_lndt)
s <- rep(seq_len(ns), each = nt)
cond <- rep(rep(0:1, each = nt / 2L), times = ns)
d <- ddm_simulate(ns * nt, mu = tr$mu0 + tr$mu_cond * cond + u_mu[s],
                  bs = tr$bs * exp(u_bs[s]), ndt = tr$ndt * exp(u_nd[s]),
                  bias = 0.5, sv = tr$sv)
d$s <- factor(s)
d$cond <- factor(cond, labels = c("a", "b"))
t0 <- tr$ndt * exp(u_nd)
key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]

report <- function(tag, fit) {
  ci <- suppressWarnings(stats::confint(fit))
  j <- grep("condb", rownames(ci), fixed = TRUE)[1L]
  b <- unlist(fixef(fit))
  hat <- as.numeric(ndt_time(fit, newdata = key))
  cat(tag, "logLik", format(as.numeric(logLik(fit)), digits = 10),
      "| conv", fit$opt$convergence,
      "| maxgrad", format(max(abs(fit$obj$gr(fit$opt$par))), digits = 3),
      "\n")
  cat(tag, "mu.condb", format(b[["mu.condb"]], digits = 6),
      "(", format(ci[j, 1], digits = 6), ",",
      format(ci[j, 2], digits = 6), ") covers 0.9:",
      ci[j, 1] <= 0.9 && ci[j, 2] >= 0.9, "\n")
  cat(tag, "mu0", format(b[["mu.(Intercept)"]], digits = 6),
      "| sv", if (is.null(b[["sv.(Intercept)"]])) NA else
        format(exp(b[["sv.(Intercept)"]]), digits = 6),
      "| sd(ndt)", format(sd(hat), digits = 5),
      "truth", format(sd(t0), digits = 5),
      "| ndt bias ms", format(1000 * mean(hat - t0), digits = 4), "\n\n")
}

fg <- frm(bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s),
             bs ~ 1 + (1 | s), ndt ~ 1 + (1 | s), bias = 0.5),
          family = wiener(variability = "sv"), data = d, se = TRUE)
report("sv-per-group ", fg)

fo <- frm(bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
             ndt ~ 1 + (1 | s), bias = 0.5),
          family = wiener(variability = "sv"), data = d, se = TRUE)
report("sv-global    ", fo)

pg <- frm(bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s),
             bs ~ 1 + (1 | s), ndt ~ 1 + (1 | s), bias = 0.5),
          family = wiener(), data = d, se = TRUE)
report("plain-pg(sv) ", pg)
