# brms 2.23.0 reference for items 1 and 3, on the item-1 data.
#   R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win \
#   Rscript dev/predfix-brms.R > dev/predfix-log/brms.txt
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "StanHeaders",
    format(packageVersion("StanHeaders")), "\n")
set.seed(20260921)
n <- 150
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
u <- rnorm(15, 0, 0.8)
e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2))
d$y1 <- 1 + 0.5 * d$x + e[, 1]
d$y2 <- -0.3 * d$x + e[, 2]
d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(n)
shp <- function(lbl, v) {
  cat(lbl, ": dim ", paste(dim(v), collapse = " x "), "\n", sep = "")
  dn <- dimnames(v)
  for (i in seq_along(dn)) {
    cat("  dimnames[[", i, "]]: ",
        if (is.null(dn[[i]])) "NULL" else paste(head(dn[[i]], 5),
                                                collapse = ", "),
        "\n", sep = "")
  }
}

cat("\n== item 1: fitted() on a multivariate fit\n")
bm <- brm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = d,
          chains = 1, iter = 600, refresh = 0, seed = 1)
shp("fitted(mv)", fitted(bm))
shp("fitted(mv, resp = 'y2')", fitted(bm, resp = "y2"))
shp("fitted(mv, resp = c('y1', 'y2'))", fitted(bm, resp = c("y1", "y2")))
shp("fitted(mv, newdata = 3 rows)", fitted(bm, newdata = d[1:3, ]))
shp("fitted(mv, dpar = 'sigma')", fitted(bm, dpar = "sigma"))
shp("fitted(mv, scale = 'linear')", fitted(bm, scale = "linear"))
shp("predict(mv)", predict(bm))
print(fitted(bm, newdata = d[1:3, ]))

cat("\n== item 3: an unseen level without allow_new_levels\n")
bg <- brm(y ~ x + (1 | g), data = d, chains = 1, iter = 600, refresh = 0,
          seed = 1)
nd <- data.frame(x = 0, g = factor("new1"))
for (call in c("fitted", "predict", "posterior_epred")) {
  r <- tryCatch({
    do.call(call, list(bg, newdata = nd))
    "answered"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(call, "(no flag): ", r, "\n", sep = "")
}
r <- tryCatch({
  fitted(bg, newdata = nd, sample_new_levels = "gaussian")
  "answered"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("fitted(sample_new_levels = 'gaussian' alone): ", r, "\n", sep = "")
r <- tryCatch({
  fitted(bg, newdata = nd, re_formula = NA)
  "answered"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("fitted(re_formula = NA, no flag): ", r, "\n", sep = "")
r <- tryCatch({
  fitted(bg, newdata = nd, re_formula = NA, allow_new_levels = TRUE)
  "answered"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("fitted(re_formula = NA, allow_new_levels = TRUE): ", r, "\n", sep = "")
