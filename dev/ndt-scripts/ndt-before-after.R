# ndt lane: the BEHAVIORAL failure, on the same data, on both builds.
#
# Two groups, non-decision times 0.18 and 0.40. The slow group's true
# ndt is above the FAST group's fastest response, so a single global
# bound cannot express it at any coefficient value. Seed 909.
#
#   Rscript --vanilla dev/ndt-scripts/ndt-before-after.R ref
#   Rscript --vanilla dev/ndt-scripts/ndt-before-after.R new
a <- commandArgs(trailingOnly = TRUE)
arm <- if (length(a)) a[1L] else "new"
.libPaths(if (identical(arm, "ref")) {
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
cat("arm:", arm, "frmtmb.eam",
    as.character(utils::packageVersion("frmtmb.eam")), "\n")

set.seed(909)
nt <- 250L
t0 <- c(0.18, 0.40)
d <- ddm_simulate(2L * nt, mu = 1.4, bs = 1.3, ndt = rep(t0, each = nt),
                  bias = 0.5)
d$s <- factor(rep(c("a", "b"), each = nt))
own <- c(tapply(d$rt, d$s, min))
cat("global min(rt)      :", format(min(d$rt), digits = 8), "\n")
cat("per-group min(rt)   :", format(own, digits = 8), "\n")
cat("true ndt            :", format(t0, digits = 8), "\n")
cat("slow group's truth above the global bound:", t0[2L] > min(d$rt),
    "\n")

# The GLOBAL bound, which is what both builds do without ndt_group().
old <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 0 + s, bias = 0.5),
           family = wiener(), data = d, se = TRUE)
lk <- family(old)$links$ndt
t_old <- lk$linkinv(unlist(fixef(old))[c("ndt.sa", "ndt.sb")])
if (identical(arm, "new")) t_old <- t_old * min(d$rt)
cat("GLOBAL bound  logLik", format(as.numeric(logLik(old)), digits = 10),
    " conv", old$opt$convergence, " pdHess",
    paste(old$sdr$pdHess, collapse = ""), " maxgrad",
    format(max(abs(old$obj$gr(old$opt$par))), digits = 3), "\n")
cat("GLOBAL bound  ndt per group:", format(t_old, digits = 8), "\n")
cat("GLOBAL bound  slow group relative error:",
    format(abs(t_old[2L] - t0[2L]) / t0[2L], digits = 4), "\n")
cat("GLOBAL bound  slow group ndt is at the wall (within 1e-6 of it):",
    abs(t_old[2L] - min(d$rt)) < 1e-6, "\n")
ci <- suppressWarnings(confint(old))
j <- grep("ndt", rownames(ci))
cat("GLOBAL bound  link-scale ndt rows, Wald interval and se:\n")
print(cbind(ci[j, 1:2, drop = FALSE],
            se = (ci[j, 2] - ci[j, 1]) / (2 * 1.959964)))

if (!identical(arm, "new")) {
  cat("PER-GROUP bound: not available on this build\n")
  quit(save = "no")
}

new <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1, ndt ~ 0 + s,
              bias = 0.5), family = wiener(), data = d, se = TRUE)
one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
t_new <- as.numeric(ndt_time(new, newdata = one))
cat("PER-GROUP     logLik", format(as.numeric(logLik(new)), digits = 10),
    " conv", new$opt$convergence, " pdHess",
    paste(new$sdr$pdHess, collapse = ""), " maxgrad",
    format(max(abs(new$obj$gr(new$opt$par))), digits = 3), "\n")
cat("PER-GROUP     ndt per group:", format(t_new, digits = 8), "\n")
cat("PER-GROUP     slow group relative error:",
    format(abs(t_new[2L] - t0[2L]) / t0[2L], digits = 4), "\n")
cat("PER-GROUP     fast group relative error:",
    format(abs(t_new[1L] - t0[1L]) / t0[1L], digits = 4), "\n")
cat("log-likelihood gain, same parameter count:",
    format(as.numeric(logLik(new)) - as.numeric(logLik(old)),
           digits = 6), "\n")
