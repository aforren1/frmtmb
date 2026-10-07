# Reviewer of lane optima, final check, n1 and n2: summary() of a
# weight at a face, and the order of the simo_ columns on the review's
# two-term design against brms 2.23.0's order (recorded in the re-check,
# dev/optima-rev2-out/readers.txt).
#   Rscript dev/optima-rev3-n1n2.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(1)
dw <- data.frame(x = sample(0:3, 100, TRUE))
dw$y <- 0.05 * dw$x + rnorm(100)
print(summary(frm(bf(y ~ mo(x)), data = dw, family = gaussian())))
set.seed(5)
n <- 300
d <- data.frame(x1 = sample(0:3, n, TRUE), x2 = sample(0:4, n, TRUE),
                z = rnorm(n), g = factor(sample(1:15, n, TRUE)))
d$y <- c(0, 1, 1, 2)[d$x1 + 1] + c(0, 0.3, 0.6, 0.6, 1)[d$x2 + 1] +
  0.3 * d$z + rnorm(15, sd = 0.5)[d$g] +
  rnorm(n, sd = exp(c(0, 0.2, 0.2, 0.4)[d$x1 + 1]))
ft <- frm(bf(y ~ mo(x1) * z + mo(x2)), data = d, family = gaussian())
s <- suppressWarnings(suppressMessages(
  frm_sample(ft, chains = 1, iter = 300, seed = 3, cores = 1,
             refresh = 0)))
v <- grep("^simo", variables(s), value = TRUE)
brms_order <- c(paste0("simo_mox11[", 1:3, "]"), paste0("simo_mox21[", 1:4, "]"),
                paste0("simo_mox1:z1[", 1:3, "]"))
cat("frmtmb:", v, "\nbrms:  ", brms_order, "\nidentical order:",
    identical(v, brms_order), "\n")
cat("all variables:", variables(s), "\n")
# the inverse by name: a reader still reproduces the chart
im <- frmtmb.sample:::draws_internal_matrix(s)
m <- as.matrix(s)
for (tm in frmtmb::mo_frame_terms(ft)) {
  zc <- im[, intersect(tm$names, colnames(im)), drop = FALSE]
  cat(tm$zeta, "max |mo_simplex(chart) - simo_|",
      format(max(abs(t(apply(zc, 1, frmtmb::mo_simplex)) -
                     m[, tm$names])), digits = 3), "\n")
}
pe <- posterior_epred(s)
cat("posterior_epred finite:", all(is.finite(pe)), "\n")
