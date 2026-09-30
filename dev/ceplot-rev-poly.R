# Reviewer re-check (lane ceplot, punch round 1): effects read through a
# transform, y ~ poly(x, 2) + log(abs(z) + 1), against brms 2.23.0 run
# at ONE draw equal to frmtmb's estimate (fixed_param), so brms's
# estimate__ is the curve at those parameters. Grid, conditions and
# estimate are compared. Data seed 13.
#   Rscript dev/ceplot-rev-poly.R > dev/ceplot-rev-log/p1/poly.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
set.seed(13)
d <- data.frame(x = rnorm(150), z = rnorm(150), w = rnorm(150))
d$y <- rnorm(150, 1 + 0.5 * d$x - 0.3 * d$x^2 + 0.8 * log(abs(d$z) + 1))
f <- frm(bf(y ~ poly(x, 2) + log(abs(z) + 1)), family = gaussian(), data = d)
b0 <- fixef(f)[, "Estimate"]
print(b0)
sd <- standata(y ~ poly(x, 2) + log(abs(z) + 1), data = d)
cat("brms X columns:", colnames(sd$X), "\n")
beta <- unname(b0[-1])
xm <- colMeans(sd$X[, -1, drop = FALSE])
init <- list(list(b = array(beta, length(beta)),
                  Intercept = unname(b0[1]) + sum(xm * beta),
                  sigma = sigma(f)))
b <- suppressMessages(suppressWarnings(
  brm(y ~ poly(x, 2) + log(abs(z) + 1), data = d, algorithm = "fixed_param",
      chains = 1, iter = 1, warmup = 0, init = init, refresh = 0, seed = 1,
      silent = 2)))
cat("brms b at the draw vs frmtmb:",
    max(abs(as.numeric(as_draws_matrix(b)[1, c("b_Intercept",
      grep("^b_poly|^b_log", variables(b), value = TRUE))]) - b0)), "\n")
for (eff in c("x", "z", "x:z", "w")) {
  a <- tryCatch(suppressWarnings(conditional_effects(b, effects = eff,
                                                     resolution = 5)),
                error = function(e) conditionMessage(e))
  r <- tryCatch(suppressWarnings(conditional_effects(f, effects = eff,
                                                     resolution = 5)),
                error = function(e) conditionMessage(e))
  if (is.character(a) || is.character(r)) {
    cat(eff, ": brms", if (is.character(a)) substr(a, 1, 90) else "answers",
        "| frmtmb", if (is.character(r)) substr(r, 1, 90) else "answers",
        "\n")
    next
  }
  a <- a[[1]]
  r <- r[[1]]
  cols <- intersect(c("x", "z"), names(a))
  gd <- max(vapply(cols, function(v) max(abs(a[[v]] - r[[v]])), 1))
  cat(sprintf("%-4s grid max diff %.3g; estimate max diff %.3g (scale %.3g)\n",
              eff, gd, max(abs(a$estimate__ - r$estimate__)),
              max(abs(a$estimate__))))
}
