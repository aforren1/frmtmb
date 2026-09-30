# Smoke test of the post-fit paths on disc and threshold-structure fits:
# conditional_effects(categorical = TRUE), emmeans, VarCorr, influence,
# frm_bootstrap, refit through update(), residuals, confint, loo-free
# methods. Seed 20260930. Output: dev/ordinal-log-smoke2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(20260930)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                f = factor(sample(c("p", "q", "r"), n, TRUE)))
re <- rnorm(8, 0, 0.5)
u <- stats::rlogis(n) / exp(0.4 * d$z) + 0.8 * d$x + re[d$g]
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
try1 <- function(label, expr) {
  cat("\n==", label, "==\n")
  tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  })
}
fits <- list(
  equi = frm(bf(y ~ x + f + (1 | g), disc ~ 0 + z),
             family = cumulative(threshold = "equidistant"), data = d),
  stz = frm(bf(y ~ x + f + (1 | g), disc ~ 0 + z),
            family = sratio(threshold = "sum_to_zero"), data = d),
  acat_disc = frm(bf(y ~ x + f + (1 | g), disc ~ 0 + z),
                  family = acat(), data = d)
)
for (nm in names(fits)) {
  f <- fits[[nm]]
  cat("\n######", nm, "######\n")
  try1("VarCorr", print(VarCorr(f)))
  try1("ce categorical", {
    ce <- conditional_effects(f, "x", categorical = TRUE)
    print(dim(ce[[1]])); print(head(ce[[1]][, c("x", "cats__", "estimate__")], 3))
  })
  try1("ce default", {
    ce <- conditional_effects(f, "z")
    print(head(ce[[1]][, c("z", "estimate__")], 3))
  })
  try1("emmeans", print(emmeans::emmeans(f, ~ f)))
  try1("residuals", print(head(residuals(f), 3)))
  try1("confint", print(head(confint(f), 8)))
  try1("influence", {
    inf <- influence(f, groups = "g")
    print(inf)
  })
  try1("bootstrap", {
    b <- frm_bootstrap(f, nsim = 3, seed = 1)
    print(b)
  })
  try1("update refit", {
    f2 <- update(f, data = d[d$y < 5, ])
    print(variables(f2))
  })
  try1("predict newdata", print(predict(f, newdata = d[1:3, ])))
  try1("fitted newdata", print(dim(fitted(f, newdata = d[1:3, ]))))
  try1("hypothesis", print(hypothesis(f, "x > 0")$hypothesis))
  try1("validate_prior", print(validate_prior(
    set_prior("normal(0, 2)", class = "b"), bf(y ~ x + f + (1 | g), disc ~ 0 + z),
    family = f$family, data = d)))
}
