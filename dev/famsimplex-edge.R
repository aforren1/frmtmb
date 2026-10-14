# Edge cases of the simplex families on the lane build, beside brms
# 2.23.0's standata() on the same data where brms has an answer.
# Output: dev/famsimplex-edge.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cat("MESSAGE:", conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  if (length(w)) cat("WARNINGS:", w, sep = "\n  ")
  if (is.character(r) && length(r) == 1) cat(r, "\n") else print(r)
  invisible(r)
}
brms_sd <- function(Y, dd, fam) {
  dd$Y <- Y
  tryCatch({
    s <- brms::standata(brms::bf(Y ~ x), data = dd, family = fam)
    paste("brms accepts; N =", s$N)
  }, error = function(e) paste("brms ERROR:", conditionMessage(e)),
  warning = function(w) paste("brms WARNING:", conditionMessage(w)))
}
dd <- fs_data_dirichlet(seed = 11)
Y <- dd$Y
dd$Y <- NULL
fitY <- function(Y, fam = dirichlet()) {
  d <- dd
  d$Y <- Y
  fit <- frm(bf(Y ~ x), family = fam, data = d)
  paste("frm fits; nobs =", nobs(fit))
}
# 1. a row off the simplex by 1e-9, 1e-7 (one row of 300) and 1e-5
for (eps in c(1e-9, 1e-7, 1e-5)) {
  Yb <- Y
  Yb[1, 1] <- Yb[1, 1] + eps
  cat("\n=== row 1 off by", eps, "\n")
  cat(brms_sd(Yb, dd, brms::dirichlet()), "\n")
  show("frm", fitY(Yb))
}
# 2. the same 1e-7 departure on every row: brms's mean test refuses it
Yb <- Y
Yb[, 1] <- Yb[, 1] + 1e-7
cat("\n=== every row off by 1e-7\n")
cat(brms_sd(Yb, dd, brms::dirichlet()), "\n")
show("frm", fitY(Yb))
# 3. a zero, a negative, an NA row
for (case in c("zero", "negative", "NA")) {
  Yb <- Y
  Yb[1, ] <- switch(case, zero = c(0.5, 0.5, 0), negative = c(1.2, -0.2, 0),
                    "NA" = c(NA, 0.5, 0.5))
  cat("\n===", case, "\n")
  cat(brms_sd(Yb, dd, brms::dirichlet()), "\n")
  show("frm", fitY(Yb))
  show("frm logistic_normal", fitY(Yb, logistic_normal()))
}
# 4. one column, a vector, an unnamed matrix, cbind()
cat("\n=== shapes\n")
show("one column", fitY(Y[, 1, drop = FALSE]))
d2 <- dd
d2$y1 <- Y[, 1]; d2$y2 <- Y[, 2]; d2$y3 <- Y[, 3]
show("vector", frm(bf(y1 ~ x), family = dirichlet(), data = d2))
show("cbind", fixef(frm(bf(cbind(y1, y2, y3) ~ x), family = dirichlet(),
                        data = d2)))
d3 <- dd
d3$Y <- unname(Y)
show("unnamed", fixef(frm(bf(Y ~ x), family = dirichlet(), data = d3)))
# 5. refcat: a named one, NA, an unknown one
d4 <- dd
d4$Y <- Y
show("refcat b", fixef(frm(bf(Y ~ x), family = dirichlet(refcat = "b"),
                           data = d4)))
show("refcat b logLik", c(logLik(frm(bf(Y ~ x), family = dirichlet(), data = d4)),
                          logLik(frm(bf(Y ~ x),
                                     family = dirichlet(refcat = "b"),
                                     data = d4))))
show("refcat NA", {
  f <- frm(bf(Y ~ x), family = dirichlet(refcat = NA), data = d4)
  print(summary(f))
  fixef(f)
})
show("refcat zz", frm(bf(Y ~ x), family = dirichlet(refcat = "zz"), data = d4))
show("ln refcat NA", logistic_normal(refcat = NA))
show("bad link", dirichlet(link = "probit"))
show("bad link_phi", dirichlet(link_phi = "logit"))
show("link fields", c(dirichlet()$link, dirichlet()$link_phi,
                      logistic_normal()$link, logistic_normal()$link_sigma,
                      brmsfamily("dirichlet2")$link))
show("brms family object", fixef(frm(bf(Y ~ x),
                                     family = brms::dirichlet(refcat = "c"),
                                     data = d4)))
# 6. K = 2 and K = 4
set.seed(5)
d5 <- data.frame(x = rnorm(200))
P <- fs_softmax(cbind(0, 0.3 + 0.5 * d5$x))
d5$Y <- fs_rdirichlet(P * 10)
colnames(d5$Y) <- c("lo", "hi")
show("K = 2 dirichlet", fixef(frm(bf(Y ~ x), family = dirichlet(), data = d5)))
show("K = 2 ln", summary(frm(bf(Y ~ x), family = logistic_normal(),
                             data = d5)))
set.seed(6)
d6 <- data.frame(x = rnorm(300))
M <- cbind(0.2 + 0.4 * d6$x, -0.1 * d6$x, 0.5 + 0 * d6$x)
C <- matrix(c(1, 0.3, -0.2, 0.3, 1, 0.4, -0.2, 0.4, 1), 3)
d6$Y <- fs_rlogistic_normal(M, diag(c(0.5, 0.7, 0.9)) %*% C %*%
                              diag(c(0.5, 0.7, 0.9)))
colnames(d6$Y) <- c("a", "b", "c", "d")
f6 <- frm(bf(Y ~ x), family = logistic_normal(), data = d6)
show("K = 4 ln summary", summary(f6))
# 7. invariance: logistic_normal's maximum under another reference
show("ln refcat logLik", c(logLik(f6),
                           logLik(frm(bf(Y ~ x),
                                      family = logistic_normal(refcat = "c"),
                                      data = d6))))
# 8. a category's own formula and weights
show("muc ~ 1", fixef(frm(bf(Y ~ x, muc ~ 1), family = dirichlet(),
                          data = d4)))
d4$w <- rep(c(1, 2), 150)
show("weights", fixef(frm(bf(Y | weights(w) ~ x), family = dirichlet(),
                          data = d4)))
show("trials refused", frm(bf(Y | trials(w) ~ x), family = dirichlet(),
                           data = d4))
show("mixture refused", frm(bf(Y ~ x), family = mixture(dirichlet(),
                                                        dirichlet()),
                            data = d4))
