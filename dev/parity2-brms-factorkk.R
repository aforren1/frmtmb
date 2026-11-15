# Lane parity2: y ~ factor(kk) * x with integer kk, no random effects,
# data seed 21, n 200. What brms 2.23.0's conditional_effects() does
# (fixed_param fit), and what frmtmb does on the lane library.
#   Rscript dev/parity2-brms-factorkk.R > dev/parity2-out/factorkk.txt 2>&1
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(brms); library(frmtmb)})
cat("brms", as.character(packageVersion("brms")), " frmtmb lib",
    find.package("frmtmb"), "\n")
set.seed(21)
n <- 200
d <- data.frame(x = rnorm(n), kk = sample(1:6, n, TRUE))
d$y <- d$x + d$kk * 0.3 + rnorm(n)
show <- function(label, expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(label, ": ERROR ", conditionMessage(r), "\n", sep = "")
  } else {
    cat(label, ": ", paste(names(r), collapse = " "), "\n", sep = "")
    for (e in names(r)) {
      v <- r[[e]][["kk"]]
      cat("   ", e, " kk values: ", paste(head(unique(v), 6), collapse = " "),
          " class ", class(v)[1], "\n", sep = "")
    }
  }
}
fb <- suppressMessages(suppressWarnings(
  brm(y ~ factor(kk) * x, data = d, algorithm = "fixed_param", chains = 1,
      iter = 20, warmup = 0, init = 0, refresh = 0, seed = 1, silent = 2)))
show("brms default", conditional_effects(fb))
show("brms effects = 'x:kk'", conditional_effects(fb, effects = "x:kk"))
ff <- frm(bf(y ~ factor(kk) * x), family = gaussian(), data = d)
show("frmtmb default", conditional_effects(ff))
show("frmtmb effects = 'x:kk'", conditional_effects(ff, effects = "x:kk"))
