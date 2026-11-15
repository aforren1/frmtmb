# Review of lane silent2, item 5235: an addition term that reads a grid
# column AND a global vector of the fit's length, on a CE grid of that
# length (resolution = n), and on newdata of that length.
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib", "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
scaled <- frmtmb_family(
  "scaled_normal", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- dpars[["mu"]] * aterms[["vreal1"]]
    -0.5 * ((y - m) / dpars[["sigma"]])^2 - log(dpars[["sigma"]]) - 0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y / aterms[["vreal1"]]),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) dpars[["mu"]] * aterms[["vreal1"]]),
  sim = function(dpars, aterms, n) stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]]))
set.seed(3)
n <- 120
d <- data.frame(x = rnorm(n), s = 1)
gv <- rep_len(c(1, 3), n)
d$y <- (1 + 0.5 * d$x) * gv + rnorm(n, 0, 0.3)
f <- suppressMessages(frm(bf(y | vreal(s * gv) ~ x), family = scaled, data = d))
r <- tryCatch(conditional_effects(f, effects = "x", resolution = n)[[1]],
              error = function(e) conditionMessage(e))
if (is.character(r)) cat("CE res 120:", r, "\n") else {
  cat("CE res 120 ran; estimate rows 1:4:", format(r$estimate__[1:4], digits = 4),
      "; x rows 1:4:", format(r$x[1:4], digits = 4), "\n")
}
nd <- data.frame(x = seq(-2, 2, length.out = n), s = 1)
p <- tryCatch(fitted(f, newdata = nd)[1:4, "Estimate"], error = function(e) conditionMessage(e))
cat("newdata 120:", if (is.character(p)) p else format(p, digits = 4), "\n")
