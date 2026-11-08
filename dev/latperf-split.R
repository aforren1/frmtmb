# Lane latperf, fix round (review M1), from dev/review-latperf-pattern.R:
# positions that its likelihood does not couple, so the pair can fall outside
# the pattern of the factor and the row goes to the direct solve. A
# zero-inflated or hurdle family separates eta_mu from eta_zi on y > 0 rows.
# Counts the rows sent to jcs_rowvar_solve(), the columns jcs_W_pairs()
# solves, and times fitted().
# Usage: Rscript dev/latperf-split.R <n> <design>   (seed 1)
#   designs: zipzi   zero_inflated_poisson, zi ~ (1 | id), latent on mu
#            hurdle  hurdle_poisson, hu ~ (1 | id), latent on mu
#            zipmu   zero_inflated_poisson, (1 | id) on mu (control)
a <- commandArgs(TRUE)
n <- as.integer(a[1])
design <- a[2]
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat(sprintf("lib: %s frmtmb %s\n", dirname(find.package("frmtmb")),
            format(packageVersion("frmtmb"))))
cat("n", n, "design", design, "\n")
nsolve <- 0L; ncols <- 0L
old_solve <- get("jcs_rowvar_solve", envir = ns)
unlockBinding("jcs_rowvar_solve", ns)
assign("jcs_rowvar_solve", function(S, Ar, ri) {
  nsolve <<- nsolve + nrow(Ar)
  old_solve(S, Ar, ri)
}, envir = ns)
ncols <- 0L
old_pairs <- get("jcs_W_pairs", envir = ns)
unlockBinding("jcs_W_pairs", ns)
assign("jcs_W_pairs", function(S, i, j) {
  out <- old_pairs(S, i, j)
  ncols <<- ncols + attr(out, "columns")
  out
}, envir = ns)
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n),
                id = factor(sample.int(50, n, TRUE)))
e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
u <- rnorm(50, 0, 0.5)
mu <- exp(.5 + .3 * d$x + e + if (design == "zipmu") u[d$id] else 0)
zi <- plogis(-1 + if (design == "zipmu") 0 else u[d$id])
y <- rpois(n, mu)
d$y <- ifelse(runif(n) < zi, 0L, y)
if (design == "hurdle") {
  yp <- y
  while (any(yp == 0)) yp[yp == 0] <- rpois(sum(yp == 0), mu[yp == 0])
  d$y <- ifelse(runif(n) < zi, 0L, yp)
}
f <- switch(design,
  zipzi = bf(y ~ x + ar(t, g, cov = TRUE), zi ~ (1 | id)),
  hurdle = bf(y ~ x + ar(t, g, cov = TRUE), hu ~ (1 | id)),
  zipmu = bf(y ~ x + (1 | id) + ar(t, g, cov = TRUE)))
fam <- if (design == "hurdle") hurdle_poisson() else zero_inflated_poisson()
tf <- system.time(fit <- suppressWarnings(frm(f, family = fam, data = d)))
cat(sprintf("RESULT fit time %.2f code %s\n", tf[["elapsed"]],
            fit$opt$convergence))
cat(sprintf("RESULT rows with y > 0: %d of %d\n", sum(d$y > 0), n))
nsolve <- 0L; ncols <- 0L
tt <- system.time(v <- suppressWarnings(fitted(fit)))[["elapsed"]]
cat(sprintf("RESULT fitted time %.2f direct-solve rows %d columns %d se_na %d\n", tt,
            nsolve, ncols, sum(is.na(v[, "Est.Error"]))))
nsolve <- 0L; ncols <- 0L
tl <- system.time(l <- frm_linpred(fit, se.fit = TRUE))[["elapsed"]]
cat(sprintf("RESULT linpred(mu) time %.2f direct-solve rows %d\n", tl, nsolve))
if (n <= 2000) {
  # against the dense route on the same fit
  f2 <- fit
  f2$cache <- list2env(mget(setdiff(ls(fit$cache), "Vjoint"),
                            envir = fit$cache))
  old <- get("jcs_joint_cov", envir = ns)
  unlockBinding("jcs_joint_cov", ns)
  assign("jcs_joint_cov", function(fit) NULL, envir = ns)
  w <- suppressWarnings(fitted(f2))
  assign("jcs_joint_cov", old, envir = ns)
  r <- v[, "Est.Error"] / w[, "Est.Error"] - 1
  cat(sprintf("RESULT fitted vs dense max |rel| %.2e\n", max(abs(r))))
}
