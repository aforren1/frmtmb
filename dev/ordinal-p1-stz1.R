# B3: sum-to-zero with one threshold per vector (every thres(gr = ) level
# at two categories, and the ungrouped two-category case). Fit, the post-
# fit methods, and the log-likelihood against brms's R-side density
# with every threshold at 0. Seed 20261001.
# Output: dev/ordinal-p1-log-stz1.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261001)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
d$y2 <- 1L + (stats::rlogis(n) + 0.8 * d$x > 0.3)
step <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, ": ", paste(format(r), collapse = " "), "\n", sep = "")
  invisible(r)
}
for (fam in c("cumulative", "sratio", "acat")) {
  cat("\n==", fam, "thres(gr = g), every level one threshold, disc ~ 0 + z\n")
  f <- tryCatch(suppressMessages(
    frm(bf(y2 | thres(gr = g) ~ x, disc ~ 0 + z),
        family = get(fam)(threshold = "sum_to_zero"), data = d)),
    error = function(e) e)
  if (inherits(f, "error")) { cat("ERROR:", conditionMessage(f), "\n"); next }
  step("length(tau_raw)", length(f$estimates$tau_raw))
  step("variables", variables(f))
  step("fixef rows", rownames(fixef(f)))
  eta <- d$x * fixef(f)["x", "Estimate"]
  disc <- exp(d$z * fixef(f)["disc_z", "Estimate"])
  P <- get(paste0("d", fam), asNamespace("brms"))(1:2, eta = eta,
    thres = matrix(0, n, 1), disc = disc, link = "logit")
  ll <- sum(log(P[cbind(seq_len(n), d$y2)]))
  step("logLik rel diff to brms's density at thresholds 0",
       abs(as.numeric(logLik(f)) - ll) / abs(ll))
  step("max |fitted - brms|", max(abs(fitted(f)[, "Estimate", ] - P)))
  step("summary", capture.output(print(summary(f)))[8:14])
  step("simulate", table(simulate(f, nsim = 1, seed = 1)[[1]]))
  step("predict newdata", dim(predict(f, newdata = d[1:3, ])))
}
cat("\n== ungrouped two categories, cumulative sum_to_zero\n")
f <- suppressMessages(frm(y2 ~ x, family = cumulative(threshold = "sum_to_zero"),
                          data = d))
step("variables", variables(f))
g <- stats::glm(I(y2 == 1) ~ 0 + I(-x), family = binomial, data = d)
step("logLik rel diff to glm(y == 1 ~ 0 + -x)",
     abs(as.numeric(logLik(f)) - as.numeric(logLik(g))) /
       abs(as.numeric(logLik(g))))
