## Reviewer, claim 4, BASE arm: what ranef()/coef() on draws returned for
## the rr and animal-model blocks before the lane. Same data and seeds as
## dev/brmsnames-rev-excluded.R.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
set.seed(21)
G <- 12; n <- 240
d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), g = factor(rep(1:G, 20)),
                id = factor(rep(1:G, 20)))
A <- diag(G); A[cbind(1:(G - 1), 2:G)] <- 0.25; A[cbind(2:G, 1:(G - 1))] <- 0.25
rownames(A) <- colnames(A) <- levels(d$id)
d$y <- 1 + d$x1 + rnorm(G)[d$g] * d$x1 + rnorm(G)[d$g] * 0.5 * d$x2 + rnorm(n)
t1 <- function(e) tryCatch(q(e), error = function(err) paste("ERROR:", conditionMessage(err)))
ms <- list(
  rr = list(bf(y ~ x1 + rr(x1 + x2 | g, d = 1)), NULL),
  animal = list(bf(y ~ 1 + (1 | gr(id, cov = A)) + (1 | id)), list(A = A))
)
for (nm in names(ms)) {
  cat("\n== base", nm, "==\n")
  fit <- q(frm(ms[[nm]][[1]], family = gaussian(), data = d, data2 = ms[[nm]][[2]]))
  ds <- q(frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 4))
  re <- t1(ranef(ds))
  cat("ranef(ds) class:", class(re), "\n")
  if (!is.character(re)) { str(unclass(re)[1:min(2, length(re))], max.level = 1)
    print(utils::head(unclass(re)[[length(re)]], 3)) } else cat(re, "\n")
  cf <- t1(coef(ds)); cat("coef(ds):\n"); if (is.character(cf)) cat(cf, "\n") else str(cf, max.level = 2)
}
