# What brms 2.23.0's sigma is under each residual correlation term:
# a marginal sd or an innovation sd. Reads brms's own Stan code
# (stancode(), standata(); nothing is compiled) for every term frmtmb
# has: ar/ma/arma with cov = TRUE (order 1 and 2), with cov = FALSE
# (order 1 and 2), cosy(), unstr(), fcor(), sar() lag and error, car(),
# and ar(cov = TRUE) under sigma ~ x and under student().
# Usage: Rscript dev/sigmainn-brms-code.R
# Log: dev/sigmainn-out/brms-code.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
ng <- 6; k <- 5
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d))
n <- nrow(d)
# a path graph over the rows, row-standardized for sar()
W <- matrix(0, n, n)
W[cbind(1:(n - 1), 2:n)] <- 1
W[cbind(2:n, 1:(n - 1))] <- 1
Ws <- W / rowSums(W)
M <- 0.5^abs(outer(1:n, 1:n, "-"))
# car() over the 6 groups, a ring
A <- matrix(0, ng, ng)
for (i in 1:ng) A[i, (i %% ng) + 1] <- A[(i %% ng) + 1, i] <- 1
rownames(A) <- colnames(A) <- seq_len(ng)

show <- function(lbl, f, data = d, data2 = NULL, family = gaussian(),
                 pat = NULL) {
  cat("\n==========", lbl, "\n")
  code <- tryCatch(stancode(f, data = data, data2 = data2, family = family),
                   error = function(e) paste("ERROR:", conditionMessage(e)))
  if (length(code) == 1L && startsWith(code, "ERROR:")) {
    cat(code, "\n")
    return(invisible())
  }
  lines <- strsplit(code, "\n")[[1]]
  keep <- grep(paste(c("sigma", "chol_cor", "cholesky_cor", "_lpdf\\(",
                       "lpdf\\(Y", "target \\+=", "Lcortime", "cosy",
                       "lagsar", "errorsar", "Msar", "sdcar", "rcar",
                       "Mi\\b", "Ar\\b", "sdcar", "fcor", "Lfcor",
                       "ar\\[", "ma\\[", "Err", "err", pat),
                     collapse = "|"), lines)
  cat(paste(sprintf("%4d: %s", keep, lines[keep]), collapse = "\n"), "\n")
  # the functions block's definition of each covariance helper, whole
  for (fn in c("cholesky_cor_ar1", "cholesky_cor_ma1", "cholesky_cor_arma1",
               "cholesky_cor_cosy", "normal_time_het_lpdf",
               "normal_time_hom_lpdf", "normal_lagsar_lpdf",
               "normal_errorsar_lpdf", "normal_fcor_hom_lpdf",
               "student_t_time_hom_lpdf", "normal_time_hom_flex_lpdf")) {
    i0 <- grep(paste0("^\\s*(matrix|real) ", fn, "\\("), lines)
    if (!length(i0)) next
    i1 <- i0 + which(grepl("^  \\}", lines[(i0 + 1):length(lines)]))[1]
    cat("  ---", fn, "\n")
    cat(paste("   ", lines[i0:i1], collapse = "\n"), "\n")
  }
}

show("ar(cov = TRUE), p = 1", y ~ x + ar(t, g, cov = TRUE))
show("ma(cov = TRUE), q = 1", y ~ x + ma(t, g, cov = TRUE))
show("arma(cov = TRUE), 1, 1", y ~ x + arma(t, g, cov = TRUE))
show("ar(cov = TRUE), p = 2", y ~ x + ar(t, g, p = 2, cov = TRUE))
show("ma(cov = TRUE), q = 2", y ~ x + ma(t, g, q = 2, cov = TRUE))
show("ar(cov = FALSE), p = 1", y ~ x + ar(t, g))
show("ar(cov = FALSE), p = 2", y ~ x + ar(t, g, p = 2))
show("arma(cov = FALSE), 1, 1", y ~ x + arma(t, g))
show("cosy()", y ~ x + cosy(t, g))
show("unstr()", y ~ x + unstr(t, g))
show("fcor()", y ~ x + fcor(M), data2 = list(M = M))
show("sar(lag)", y ~ x + sar(W, type = "lag"), data2 = list(W = Ws))
show("sar(error)", y ~ x + sar(W, type = "error"), data2 = list(W = Ws))
show("car()", y ~ x + car(A, gr = g), data2 = list(A = A))
show("ar(cov = TRUE) under sigma ~ x",
     bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x))
show("ar(cov = TRUE) under student()", y ~ x + ar(t, g, cov = TRUE),
     family = student())
show("ma(cov = TRUE) under sigma ~ x",
     bf(y ~ x + ma(t, g, cov = TRUE), sigma ~ x))
show("cosy() under sigma ~ x", bf(y ~ x + cosy(t, g), sigma ~ x))
show("sar(lag) under sigma ~ x",
     bf(y ~ x + sar(W, type = "lag"), sigma ~ x), data2 = list(W = Ws))
show("fcor() under sigma ~ x", bf(y ~ x + fcor(M), sigma ~ x),
     data2 = list(M = M))
cat("\n========== default_prior(ar cov = TRUE, sigma ~ x)\n")
print(default_prior(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), data = d))
cat("\n========== default_prior(ar cov = TRUE)\n")
print(default_prior(y ~ x + ar(t, g, cov = TRUE), data = d))
