# Item 720: a vreal() payload conditional_effects() cannot evaluate on
# its grid. Usage: Rscript dev/silent-ce-probe.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")

# a gaussian whose mean is mu times a known per-row scale in vreal1
scaled <- frmtmb_family(
  "scaled_normal", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- dpars[["mu"]] * aterms[["vreal1"]]
    -0.5 * ((y - m) / dpars[["sigma"]])^2 - log(dpars[["sigma"]]) -
      0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y / aterms[["vreal1"]]),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vreal1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]])
  })
set.seed(3)
n <- 120
d <- data.frame(x = rnorm(n))
# the scale alternates by row; written as a data column it is held at
# its mean on the grid, written as an expression of the row count it
# cannot be evaluated on a grid of another length
d$s <- rep_len(c(1, 3), n)
d$y <- (1 + 0.5 * d$x) * d$s + rnorm(n, 0, 0.3)
f_col <- frm(bf(y | vreal(s) ~ x), family = scaled, data = d)
f_expr <- frm(bf(y | vreal(rep_len(c(1, 3), 120)) ~ x), family = scaled,
              data = d)
cat("logLik col", format(as.numeric(logLik(f_col)), digits = 10),
    " expr", format(as.numeric(logLik(f_expr)), digits = 10), "\n")
show <- function(lab, expr) {
  cat("==", lab, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) cat(" ", r, "\n") else {
    df <- r[[1L]]
    cat("  rows", nrow(df), " estimate__ range",
        paste(format(range(df$estimate__), digits = 6), collapse = " to "),
        "\n")
  }
}
show("col mean", conditional_effects(f_col, effects = "x"))
show("expr mean", conditional_effects(f_expr, effects = "x"))
show("col predict", conditional_effects(f_col, effects = "x",
                                        method = "predict"))
show("expr predict", conditional_effects(f_expr, effects = "x",
                                         method = "predict"))
show("expr predict, resolution 120",
     conditional_effects(f_expr, effects = "x", method = "predict",
                         resolution = 120))
cat("predict(newdata) on 5 rows:\n")
r <- tryCatch(predict(f_expr, newdata = d[1:5, ]),
              error = function(e) paste("ERROR:", conditionMessage(e)),
              warning = function(w) paste("WARNING:", conditionMessage(w)))
print(r)
