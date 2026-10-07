# Lane setier, punch round 1 (B2): the cost of separation_check() on
# large bernoulli fits, dense and sparse_x designs, at n = 1e5 and 1e6.
# The check alone, minimum of 3 calls, against two controls measured the
# same way on the same fit: one objective evaluation, and the product
# X %*% beta the check needs anyway. Memory is gc()'s max used.
#   Rscript dev/setier-sepcost.R <lib> <n> [case ...]
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
n <- as.numeric(args[2])
which_cases <- if (length(args) > 2) args[-(1:2)] else
  c("dense_p20", "sparse_f500")
cat("lib", find.package("frmtmb"), " n", n, "\n")
tm <- function(expr) {
  gc(reset = TRUE)
  t <- system.time(v <- expr)[["elapsed"]]
  g <- gc()
  list(v = v, t = t, mb = sum(g[, ncol(g)]))
}
best3 <- function(f) {
  r <- lapply(1:3, function(i) tm(f()))
  list(t = min(vapply(r, `[[`, 0, "t")), mb = max(vapply(r, `[[`, 0, "mb")),
       v = r[[1]]$v)
}
cases <- list(
  dense_p20 = function() {
    set.seed(1)
    X <- matrix(rnorm(n * 19), n)
    d <- data.frame(X)
    d$y <- rbinom(n, 1, plogis(-2 + X %*% c(3, rep(0.2, 18))))
    list(d = d, f = as.formula(paste("y ~", paste(names(d)[1:19],
                                                  collapse = " + "))),
         ctl = frmtmb_control())
  },
  sparse_f500 = function() {
    set.seed(2)
    d <- data.frame(f = factor(sample(500, n, TRUE)), x = rnorm(n))
    d$y <- rbinom(n, 1, plogis(-2 + rnorm(500, 0, 1.5)[d$f] + 2 * d$x))
    list(d = d, f = y ~ f + x, ctl = frmtmb_control(sparse_x = TRUE))
  })
for (nm in which_cases) {
  cs <- cases[[nm]]()
  r <- tm(suppressWarnings(frm(cs$f, family = bernoulli(), data = cs$d,
                                control = cs$ctl)))
  fit <- r$v
  lp <- fit$frame$linpreds[[1]]
  cat(sprintf("%s n=%g: fit %.2f s, max mem %.0f MB, code %d, X %s\n", nm,
              n, r$t, r$mb, fit$opt$convergence, class(lp$X)[1]))
  s <- best3(function() ns$separation_check(fit))
  c1 <- best3(function() fit$obj$fn(fit$opt$par))
  c2 <- best3(function() as.numeric(lp$X %*% fit$estimates$beta[lp$idx]))
  cat(sprintf(paste0("   separation_check %.3f s (%.0f MB, named %s) | ",
                     "one objective %.3f s | X %%*%% beta %.3f s\n"),
              s$t, s$mb, !is.null(s$v), c1$t, c2$t))
}
