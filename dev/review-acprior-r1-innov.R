# Reviewer, fix round 1, item 1: the innovation factor c against brms's
# OWN Stan functions (rstan::expose_stan_functions() on the functions
# block of brms's stancode), not a hand transcription. brms's residual
# covariance is sigma_i^2 L L' with L = cholesky_cor_*(ar, ma, n);
# frmtmb's is sigma_m^2 R. With sigma_i = sigma_m * c the two agree
# iff R = c^2 L L', which checks the sign convention of ma too.
# c = exp(innov_logc(th)) from the lane library.
# Usage: Rscript dev/review-acprior-r1-innov.R
# Log: dev/review-acprior-log/r1-innov.txt
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "; brms",
    format(packageVersion("brms")), "\n")
set.seed(2)
d <- data.frame(t = rep(1:6, 10), g = rep(1:10, each = 6), y = rnorm(60))
fns <- character(0)
for (fo in list(y ~ ar(t, g, cov = TRUE), y ~ ma(t, g, cov = TRUE),
                y ~ arma(t, g, cov = TRUE))) {
  sc <- stancode(fo, data = d)
  a <- regexpr("functions [{]", sc)
  b <- regexpr("\ndata [{]", sc)
  fns <- c(fns, substr(sc, a + attr(a, "match.length"), b - 1L))
}
# each block ends with the closing brace of `functions {`; drop it
body <- vapply(fns, function(s) sub("[}][[:space:]]*$", "", s), "")
# keep one copy of each function by splitting on the brms function names
prog <- paste0("functions {\n", paste(body, collapse = "\n"), "\n}\n")
# duplicate definitions across the three blocks: keep the arma block,
# which brms writes with every cholesky_cor_* it needs, plus ar1 and ma1
keep <- function(s, nm) {
  i <- regexpr(paste0("matrix ", nm, "[(]"), s)
  if (i < 0) return(NULL)
  ch <- strsplit(substr(s, i, nchar(s)), "")[[1]]
  # the function ends where its first brace closes
  depth <- cumsum((ch == "{") - (ch == "}"))
  open <- which(ch == "{")[1]
  j <- open - 1L + which(depth[open:length(ch)] == 0)[1]
  paste(ch[seq_len(j)], collapse = "")
}
src <- paste(fns, collapse = "\n")
parts <- unlist(lapply(c("cholesky_cor_ar1", "cholesky_cor_ma1",
                         "cholesky_cor_arma1"), function(nm) keep(src, nm)))
stopifnot(length(parts) == 3L)
prog <- paste0("functions {\n", paste(parts, collapse = "\n"), "\n}\n")
cat(prog)
env <- new.env()
rstan::expose_stan_functions(rstan::stanc(model_code = prog), env = env)
for (kind in c("ar", "ma", "arma")) {
  fo <- switch(kind, ar = y ~ ar(t, g, cov = TRUE),
               ma = y ~ ma(t, g, cov = TRUE),
               arma = y ~ arma(t, g, cov = TRUE))
  fr <- frm(fo, data = d, dry_run = "frame")
  ac <- fr$autocor[[1]]
  inn <- frmtmb:::autocor_innov(ac)
  for (rep in 1:4) {
    th <- rnorm(ac$npar, 0, 1.3)
    nat <- as.numeric(frmtmb:::autocor_natural(th, ac))
    ar <- if (kind != "ma") nat[1] else 0
    ma <- if (kind == "ma") nat[1] else if (kind == "arma") nat[2] else 0
    R <- as.matrix(frmtmb:::autocor_cor(th, ac))
    n <- nrow(R)
    L <- switch(kind, ar = env$cholesky_cor_ar1(ar, n),
                ma = env$cholesky_cor_ma1(ma, n),
                arma = env$cholesky_cor_arma1(ar, ma, n))
    B <- L %*% t(L)
    c2 <- exp(2 * frmtmb:::innov_logc(th, inn))
    cat(sprintf("%-4s ar %+.4f ma %+.4f | max|R - c^2 LL'| %.2e | c^2 %.6f vs 1/B11 %.6f\n",
                kind, ar, ma, max(abs(R - c2 * B)), c2, 1 / B[1, 1]))
  }
}
