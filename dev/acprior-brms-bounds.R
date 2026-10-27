# brms 2.23.0: can a user's lb/ub replace the declared bound of ar (cov =
# TRUE), ma and cosy? And brms's refusals for p > 1 under cov = TRUE.
# Log: dev/acprior-out/brms-bounds.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
n <- 40
d <- data.frame(y = rnorm(n), x = rnorm(n), t = rep(1:10, 4),
                g = rep(1:4, each = 10))
decl <- function(f, pr) {
  sc <- tryCatch(stancode(f, data = d, prior = pr),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  ln <- strsplit(sc, "\n")[[1]]
  cat(paste0("  | ", grep("ERROR|ar;|ma;|cosy;|lprior \\+=.*(ar|ma|cosy)",
                          ln, value = TRUE)), sep = "\n")
}
cosy <- y ~ x + cosy(time = t, gr = g)
ar1c <- y ~ x + ar(time = t, gr = g, p = 1, cov = TRUE)
cat("--- cosy, lb = -0.1\n")
decl(cosy, set_prior("normal(0, 0.5)", class = "cosy", lb = -0.1))
cat("--- cosy, ub = 0.8 only\n")
decl(cosy, set_prior("normal(0, 0.5)", class = "cosy", ub = 0.8))
cat("--- cosy, bounds only, no density\n")
decl(cosy, set_prior("", class = "cosy", lb = 0.2))
cat("--- ar cov, lb = 0\n")
decl(ar1c, set_prior("normal(0, 0.5)", class = "ar", lb = 0))
cat("--- ar cov, lb = -2\n")
decl(ar1c, set_prior("normal(0, 0.5)", class = "ar", lb = -2))
cat("--- ar cov, p = 2\n")
decl(y ~ x + ar(time = t, gr = g, p = 2, cov = TRUE), NULL)
cat("--- ma cov, q = 2\n")
decl(y ~ x + ma(time = t, gr = g, q = 2, cov = TRUE), NULL)
cat("--- ar cond, coef-level prior\n")
decl(y ~ x + ar(time = t, gr = g, p = 2),
     set_prior("normal(0, 0.5)", class = "ar", coef = "1"))
cat("--- default_prior row for ar p = 2 cond with lb set\n")
print(validate_prior(set_prior("normal(0, 0.5)", class = "ar", lb = -1,
                               ub = 1),
                     y ~ x + ar(time = t, gr = g, p = 2), data = d))
