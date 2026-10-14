# Lane optstart: test-se-check.R's gr(g, by = f) fixture (seed 11), from
# the lane's start and from 0.70.0's (given as `start`): objective,
# code, iterations and theta, with restarts = 0 and 1.
#   Rscript dev/optstart-grby.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(11)
dd <- data.frame(x = stats::rnorm(160), g = factor(rep(1:16, 10)))
dd$f <- factor(ifelse(as.integer(dd$g) <= 8, "a", "b"))
uu <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(160, 1 + 0.5 * dd$x + uu[dd$g, 1] +
                       uu[dd$g, 2] * dd$x, 1)
fg <- bf(y ~ x + (1 + x | gr(g, by = f)))
old <- list(beta = c(mean(dd$y), 0), betad = log(stats::sd(dd$y)))
print(unlist(par_template(fg, data = dd, family = gaussian())))
for (st in list(NULL, old)) for (rs in 0:1) {
  f <- suppressWarnings(suppressMessages(frm(fg, family = gaussian(),
    data = dd, start = st, control = frmtmb_control(restarts = rs))))
  cat(sprintf("%s restarts %d: objective %.10f code %d %s iter %d\n",
              if (is.null(st)) "lane start" else "old start ", rs,
              f$opt$objective, f$opt$convergence, f$opt$message,
              f$opt$iterations))
  print(round(f$opt$par[names(f$opt$par) == "theta"], 4))
}
