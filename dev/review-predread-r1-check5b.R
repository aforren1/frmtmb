# Short check, item 5: the callers of the multivariate shape change against brms.
# Seed 8301; also the contrast's pairing rule in frm_curve(), item 2.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
.libPaths(c(.libPaths()[1L], "C:/Users/adf44/source/r/rev-predread-lib/extra", .libPaths()[-1L]))
set.seed(8301)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 1 + d$x + rnorm(n)
d$y2 <- -d$x + 0.5 * d$z + rnorm(n)
nd <- data.frame(x = c(-1, 0, 1), z = 0.2)
mv <- rev_pair("mvbf(bf(y ~ x), bf(y2 ~ x + z), rescor = FALSE)", "gaussian()", d)
sh <- function(lab, e) {
  r <- tryCatch(quiet(e), error = function(err) paste("ERROR:", substr(conditionMessage(err), 1, 120)))
  if (is.character(r)) cat(lab, ":", r, "\n") else {
    cat(lab, ": class", class(r)[1], " dim", paste(dim(r) %||% length(r), collapse = "x"), " names", paste(names(r), collapse = ","), "\n")
    if (is.data.frame(r)) print(utils::head(as.data.frame(r)[, intersect(c("rowid", "group", "estimate", "Predicted", "Response", "Row"), names(r)), drop = FALSE], 6))
    else print(signif(utils::head(as.numeric(r), 6), 4))
  }
}
sh("frmtmb insight get_predicted", insight::get_predicted(mv$fit, data = nd))
sh("brms insight get_predicted", insight::get_predicted(mv$fx, data = nd))
sh("frmtmb predictions", marginaleffects::predictions(mv$fit, newdata = nd))
sh("brms predictions", marginaleffects::predictions(mv$fx, newdata = nd))
sh("frmtmb get_predict link", marginaleffects::get_predict(mv$fit, newdata = nd, type = "link"))
sh("brms get_predict link", marginaleffects::get_predict(mv$fx, newdata = nd, type = "link"))
cat("--- the contrast's pairing rule, seed 7402\n")
set.seed(7402)
n <- 200
e <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
e$xo <- e$lat + rnorm(n, 0, 0.4)
e$y <- 1 + 1.2 * e$lat - 0.5 * e$z + rnorm(n, 0, 0.5)
fit <- quiet(frm(bf(y ~ me(xo, sx) + z), family = gaussian(), data = e))
a <- data.frame(xo = c(1, 1, 2), sx = 0.6, z = c(1, 2, 1))
b <- data.frame(xo = c(1, 2, 2), sx = 0.6, z = c(0, 0, 0))
ct <- frm_curve(fit, newdata = a, contrast = b, simultaneous = FALSE)
cat("pairs (x 1,1) (x 1,2) (x 2,2): se", signif(ct$.se, 4), "\n")
S <- attr(ct, "Sigma")
print(signif(S, 3))
