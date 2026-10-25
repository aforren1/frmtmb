# Review: frm_sample() on cross-lane models: exp(x) + sar(W) with a
# lagsar prior and a coef = "expx" prior; cox + bhaz() + exp(x) with a
# dirichlet() sbhaz prior; a lagsar and a theta/thetaac_1 prior together.
#   Rscript dev/review-rel075-sample.R > dev/review-rel075-log/sample.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.survival)
  library(frmtmb.sample)
})
show <- function(label, expr) {
  cat("\n=====", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("MESSAGE:", conditionMessage(m))
    invokeRestart("muffleMessage")
  }), error = function(e) {
    cat("ERROR [", paste(class(e), collapse = ","), "]:",
        conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) print(r)
  invisible(r)
}
k <- 6
g <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A)
n <- nrow(W)
set.seed(20261024)
d <- data.frame(x = runif(n, 0, 1.5), z = rnorm(n))
IA <- diag(n) - 0.4 * W
d$y <- as.vector(solve(IA, 1 + 0.5 * exp(d$x) - 0.3 * d$z +
                          rnorm(n, 0, 0.5)))
fit <- frm(y ~ exp(x) + z + sar(W), data = d, data2 = list(W = W))
pr <- c(set_prior("normal(-0.5, 0.02)", class = "lagsar"),
        set_prior("normal(3, 0.01)", class = "b", coef = "expx"))
ds <- show("S1 frm_sample(sar + exp(x)) with lagsar and expx priors",
           {
             x <- frm_sample(fit, chains = 1, iter = 600, warmup = 300,
                             refresh = 0, seed = 7, prior = pr)
             m <- as.matrix(x)
             list(cols = colnames(m),
                  means = round(colMeans(m[, intersect(colnames(m),
                    c("b_expx", "lagsar", "b_z", "b_Intercept"))]), 3),
                  fixef_rows = rownames(fixef(x)))
           })
show("S2 default prior announced (no prior given)", {
  x <- frm_sample(fit, chains = 1, iter = 300, warmup = 150, refresh = 0,
                  seed = 8)
  m <- as.matrix(x)
  range(m[, "lagsar"])
})
show("S3 lagsar and theta thetaac_1 priors together (one parameter)",
     validate_prior(c(set_prior("normal(0, 0.5)", class = "lagsar"),
                      set_prior("normal(0, 1)", class = "theta",
                                coef = "thetaac_1")),
                    y ~ exp(x) + z + sar(W), data = d, data2 = list(W = W)))
show("S3b frm with both: MAP rho",
     {
       f <- frm(y ~ exp(x) + z + sar(W), data = d, data2 = list(W = W),
                prior = c(set_prior("normal(0, 0.5)", class = "lagsar"),
                          set_prior("normal(0, 1)", class = "theta",
                                    coef = "thetaac_1")))
       prior_summary(f)
     })
set.seed(20261024)
dsv <- data.frame(x = runif(120, 0, 1.5), z = rnorm(120))
dsv$time <- rexp(120, exp(-1 + 0.5 * exp(dsv$x) + 0.3 * dsv$z))
dsv$cens <- as.numeric(runif(120) < 0.2)
fc <- frm(time | cens(cens) + bhaz(df = 4) ~ exp(x) + z, data = dsv,
          family = cox())
show("S4 frm_sample(cox + bhaz + exp(x)) with dirichlet(2) on sbhaz", {
  x <- frm_sample(fc, chains = 1, iter = 400, warmup = 200, refresh = 0,
                  seed = 9, prior = set_prior("dirichlet(2)",
                                              class = "sbhaz"))
  m <- as.matrix(x)
  list(cols = colnames(m), fixef_rows = rownames(fixef(x)),
       simplex_sums = range(rowSums(m[, grep("^sbhaz", colnames(m)),
                                       drop = FALSE])))
})
show("S5 frm_sample with dirichlet on lagsar (refusal timing)", {
  frm_sample(fit, chains = 1, iter = 100, warmup = 50, refresh = 0,
             prior = set_prior("dirichlet(1)", class = "lagsar"))
})
cat("\nDONE\n")
