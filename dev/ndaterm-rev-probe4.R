# Reviewer of lane ndaterm: simulate(newdata, censored = TRUE), and
# whether CE holds a se() and a dec() variable. Seed 7.
# Usage: Rscript dev/ndaterm-rev-probe4.R <arm>   (mine | r9)
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n))
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
d$cc <- ifelse(d$cnt > 2, "right", "none")
d$cnt2 <- pmin(d$cnt, 2)
f <- suppressMessages(frm(bf(cnt2 | cens(cc) ~ x), family = poisson(),
                          data = d))
nd <- data.frame(x = c(3, 3, 3), cc = c("none", "none", "none"))
r <- tryCatch(as.matrix(simulate(f, nsim = 2000, newdata = nd, seed = 1,
                                 censored = TRUE)),
              error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r)) cat(r, "\n") else {
  cat("simulate(newdata x = 3, cc = none, censored = TRUE) row max:",
      apply(r, 1, max), "\n")
  cat("row means:", rowMeans(r), " (latent mean exp(0.5 + 0.9) = ",
      exp(1.4), ")\n")
  cat("training cc of rows 1:3:", d$cc[1:3], "\n")
}
r2 <- as.matrix(simulate(f, nsim = 2000, newdata = nd, seed = 1))
cat("simulate(newdata, censored = FALSE) row means:", rowMeans(r2), "\n")
# se(): the CE grid value
d$s <- runif(n, 0.2, 1)
d$g <- 1 + d$x + rnorm(n, 0, d$s)
fg <- frm(bf(g | se(s) ~ x), family = gaussian(), data = d)
ce <- conditional_effects(fg, "x", method = "predict", ndraws = 401)[[1]]
cat("CE se(s) grid s:", unique(ce$s), " mean(s) =", mean(d$s), "\n")
cat("done\n")
