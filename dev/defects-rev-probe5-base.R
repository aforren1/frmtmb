# Reviewer of lane defects: dev/defects-rev-probe5.R on the base build.
# with an ordinal cs() response, and residuals() of that fit.
.libPaths(c(
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
m <- function(e) tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
set.seed(5)
n <- 150
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$yo <- cut(d$x + rnorm(n), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
d$yg <- d$z + rnorm(n)
f <- m(frm(bf(yg ~ x) + bf(yo ~ x + cs(z), family = sratio()) +
             set_rescor(FALSE), data = d))
if (is.character(f)) cat(f, "\n") else {
  a <- m(fitted(f, scale = "linear"))
  cat("mv resp=NULL:", if (is.character(a)) a else paste(dim(a), collapse = "x"), "\n")
  b <- m(fitted(f, scale = "linear", resp = "yo"))
  cat("mv resp='yo':", if (is.character(b)) b else paste(dim(b), collapse = "x"), "\n")
  u <- frm(yo ~ x + cs(z), family = sratio(), data = d)
  cat("equal to univariate layers:", isTRUE(all.equal(unclass(b)[, "Estimate", ], unclass(fitted(u, scale = "linear"))[, "Estimate", ], tolerance = 1e-4)), "\n")
  g <- m(fitted(f, scale = "linear", resp = "yg"))
  cat("mv resp='yg':", if (is.character(g)) g else paste(dim(g), collapse = "x"), "\n")
  cat("residuals mv (ordinal inside):", substr(m(residuals(f)), 1, 150), "\n")
  cat("residuals resp='yg':", paste(dim(m(residuals(f, resp = "yg"))), collapse = "x"), "\n")
}
# order of responses: yo first, gaussian second
f2 <- m(frm(bf(yo ~ x + cs(z), family = sratio()) + bf(yg ~ x) +
              set_rescor(FALSE), data = d))
if (!is.character(f2)) {
  a <- m(fitted(f2, scale = "linear"))
  cat("mv (yo first) resp=NULL:", if (is.character(a)) a else paste(dim(a), collapse = "x"), "\n")
  a2 <- m(fitted(f2, scale = "linear", resp = "yg"))
  cat("mv (yo first) resp='yg':", if (is.character(a2)) a2 else paste(dim(a2), collapse = "x"), "\n")
}
