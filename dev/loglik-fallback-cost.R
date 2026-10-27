# Lane loglik, fix round 1, item 3: the cost of the dense fallback of
# autocor_row_loglik() (a registered term with no jacobian slot), by n,
# to set the size above which log_lik() warns. The term is sar()'s own
# lag map registered again without its slot, over a row-standardized
# path graph. Per n: the minimum of 5 calls, and the peak memory of one
# call (gc()'s max used, reset before it). Seed 7.
# Output: dev/loglik-fallback-cost.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spatial)})
cat("lib:", find.package("frmtmb"), "\n")
frmtmb_local_registries()
e <- frmtmb:::ac_entry("sar")
frmtmb_register_autocor("sarfb", args = c("M", "type"), parse = e$parse,
                        build = e$build, whiten = e$whiten,
                        unwhiten = e$unwhiten, npar = 1L, start = e$start,
                        natural = e$natural, sigma_varies = FALSE,
                        bounds = e$bounds)
set.seed(7)
for (n in c(500, 1000, 1500, 2000, 3000)) {
  P <- Matrix::sparseMatrix(i = c(1:(n - 1), 2:n), j = c(2:n, 1:(n - 1)),
                            x = 1, dims = c(n, n))
  W <- as.matrix(P / Matrix::rowSums(P))
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  fit <- frm(y ~ x + sarfb(W), data = d, data2 = list(W = W),
             dry_run = "objective")
  class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
  fit$estimates$thetaac[] <- 0.5
  dp <- eval_dpars(fit)
  run <- function() suppressWarnings(autocor_row_loglik(fit, dp, "y"))
  tt <- min(vapply(1:5, function(k) system.time(run())[["elapsed"]], 0))
  invisible(gc(reset = TRUE))
  m0 <- sum(gc()[, 2])
  run()
  m1 <- sum(gc()[, 6])
  w <- tryCatch(autocor_row_loglik(fit, dp, "y"),
                frmtmb_autocor_dense_fallback = function(w) "warns")
  cat(sprintf("n = %5d  call %.3f s  peak Mb %.0f  4000 draws %.2f h  %s\n",
              n, tt, m1 - m0, 4000 * tt / 3600,
              if (identical(w, "warns")) "warns" else "silent"))
}
