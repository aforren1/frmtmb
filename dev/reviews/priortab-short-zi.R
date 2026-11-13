# Short check: zi uniform(0, 0.3) boxed by the support (lane) against the
# same box written as lb/ub, both arms.
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(digits = 10)
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 0; invisible(rnorm(10)); invisible(rnorm(n))  # keep the RNG path of the guard script
dd <- d; dd$cnt <- rpois(n, 2) * rbinom(n, 1, 0.5)
run <- function(lbl, pr) {
  f <- tryCatch(suppressWarnings(frm(cnt ~ z, data = dd,
                                     family = zero_inflated_poisson(),
                                     prior = pr)),
                error = function(e) conditionMessage(e))
  if (is.character(f)) { cat(lbl, ": ERROR", substr(f, 1, 120), "\n"); return() }
  cat(sprintf("%-28s conv %d obj %.10f par %s\n", lbl, f$opt$convergence,
              f$opt$objective, paste(signif(f$opt$par, 8), collapse = " ")))
}
run("ML", NULL)
run("uniform, no bound", set_prior("uniform(0, 0.3)", class = "zi"))
run("uniform, lb 0 ub 0.3", set_prior("uniform(0, 0.3)", class = "zi",
                                      lb = 0, ub = 0.3))
run("flat, ub 0.3", set_prior("", class = "zi", ub = 0.3))
