# Fit agreement of frmtmb's fcor() and frmtmb.spatial's sar() with brms
# 2.23.0 fits on the same data. Seeds below; 4 chains x 2000 iterations.
# Run: Rscript dev/sarfcor-brms-fit.R > dev/sarfcor-out/brms-fit.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.spatial"),
    "| brms", format(packageVersion("brms")), "\n")

lattice <- function(k) {
  grid <- expand.grid(r = 1:k, c = 1:k)
  n <- nrow(grid)
  outer(seq_len(n), seq_len(n), function(i, j) {
    as.numeric(abs(grid$r[i] - grid$r[j]) + abs(grid$c[i] - grid$c[j]) == 1)
  })
}
A <- lattice(8)
n <- nrow(A)
Wrow <- A / rowSums(A)
# symmetric normalization: real eigenvalues by construction, which brms
# needs (it passes eigen()'s values to Stan as they come)
dg <- 1 / sqrt(rowSums(A))
Wsym <- A * outer(dg, dg)

set.seed(11)
d <- data.frame(x = rnorm(n))
IA <- diag(n) - 0.6 * Wsym
d$ylag <- as.vector(solve(IA, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
d$yerr <- 1 + 0.5 * d$x + as.vector(solve(IA, rnorm(n, 0, 0.7)))
IAr <- diag(n) - 0.6 * Wrow
d$ylagr <- as.vector(solve(IAr, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
V <- 0.7^abs(outer(1:n, 1:n, "-"))
d$yfc <- 1 + 0.5 * d$x + as.vector(t(chol(V)) %*% rnorm(n)) * 0.8
dd2 <- list(Wsym = Wsym, Wrow = Wrow, V = V)

cmp <- function(lab, fit, bfit, map) {
  dr <- as.data.frame(bfit)
  for (nm in names(map)) {
    est <- map[[nm]][1]
    se <- map[[nm]][2]
    pm <- median(dr[[nm]])
    psd <- sd(dr[[nm]])
    cat(sprintf("%-10s %-12s ML %9.5f SE %8.5f | brms median %9.5f sd %8.5f | (ML - median)/SE %7.3f  sd/SE %6.3f\n",
                lab, nm, est, se, pm, psd, (est - pm) / se, psd / se))
  }
  rh <- brms::rhat(bfit)
  cat(sprintf("%-10s max Rhat %.4f\n", lab, max(rh, na.rm = TRUE)))
}
fe_map <- function(fit, extra = list()) {
  cf <- summary(fit)$fixed
  out <- list(b_Intercept = c(cf["Intercept", 1], cf["Intercept", 2]),
              b_x = c(cf["x", 1], cf["x", 2]))
  sp <- summary(fit)$spec_pars
  out$sigma <- c(sp["sigma", 1], sp["sigma", 2])
  c(out, extra)
}
ac_map <- function(fit, nm) {
  cv <- confint_varcorr(fit)
  r <- cv[cv$term == nm, ]
  # the delta-method SE of rho on its own scale
  se <- (r$upr - r$lwr) / (2 * qnorm(0.975))
  stats::setNames(list(c(r$estimate, se)), nm)
}
run <- function(lab, form, ac = NULL, seed) {
  f <- frm(form, data = d, data2 = dd2)
  b <- brms::brm(form, data = d, data2 = dd2, chains = 4, iter = 2000,
                 seed = seed, refresh = 0, silent = 2,
                 backend = "rstan")
  cmp(lab, f, b, fe_map(f, if (!is.null(ac)) ac_map(f, ac)))
  invisible(list(f = f, b = b))
}
t0 <- proc.time()
r1 <- run("lag", ylag ~ x + sar(Wsym, type = "lag"), "lagsar", 101)
r2 <- run("error", yerr ~ x + sar(Wsym, type = "error"), "errorsar", 102)
r3 <- run("fcor", yfc ~ x + fcor(V), NULL, 103)
# brms on a row-standardized W, whose eigen() values come back complex
r4 <- tryCatch(run("lag-row", ylagr ~ x + sar(Wrow, type = "lag"),
                   "lagsar", 104),
               error = function(e) {
                 cat("lag-row brms ERROR:", conditionMessage(e), "\n")
               })
cat("elapsed", (proc.time() - t0)[["elapsed"]], "s\n")
