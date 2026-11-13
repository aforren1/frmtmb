# Short check: a dpar-class support box through a log link, at an lb whose
# exp(log(lb)) rounds below lb, against one that round-trips exactly.
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(digits = 12)
w <- c(seq(0.31, 0.99, by = 0.01), seq(1.1, 3, by = 0.1))
below <- w[exp(log(w)) < w]; exact <- w[exp(log(w)) == w]
cat("lb with exp(log(lb)) < lb:", head(below, 5), "\n")
cat("lb exact:", head(exact, 5), "\n")
set.seed(1)
n <- 120
d <- data.frame(z = rnorm(n))
d$y <- 1 + 0.3 * d$z + rnorm(n, sd = 0.3)
for (lb in c(below[1], exact[1])) {
  for (bnd in c(FALSE, TRUE)) {
    pr <- if (bnd) set_prior(sprintf("uniform(%.17g, 5)", lb), class = "sigma",
                             lb = lb, ub = 5)
          else set_prior(sprintf("uniform(%.17g, 5)", lb), class = "sigma")
    f <- tryCatch(suppressWarnings(frm(y ~ z, data = d, prior = pr)),
                  error = function(e) conditionMessage(e))
    if (is.character(f)) {
      cat(sprintf("lb %.2f bounds %-5s ERROR %s\n", lb, bnd, substr(f, 1, 90)))
    } else {
      cat(sprintf("lb %.2f bounds %-5s conv %d obj %.10f sigma %.17g msg %s\n",
                  lb, bnd, f$opt$convergence, f$opt$objective,
                  exp(f$opt$par[["sigma"]] %||% NA), f$opt$message))
    }
  }
}
