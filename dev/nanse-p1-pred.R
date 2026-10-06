# Punch round 1, B2: predictions along a lost direction on the
# reviewer's spread ridge (k = 10), every prediction route, with the
# warnings each raises. The identified mu = a_k + b keeps its SE.
#   Rscript dev/nanse-p1-pred.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(1)
k <- 10
f <- factor(rep(seq_len(k), each = 5))
dd <- data.frame(f = f, y = rnorm(k)[f] + rnorm(k * 5, 0, 0.5))
fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE),
                            data = dd))
nd <- dd[c(1, 6), , drop = FALSE]
show <- function(lab, expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    paste("ERROR:", conditionMessage(e))
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  cat("\n--", lab, "\n")
  if (is.character(v)) cat(v, "\n") else print(v)
  cat("   warnings:", length(w), if (length(w)) substr(w[1], 1, 120), "\n")
}
show("fitted(dpar = 'a')", fitted(fit, newdata = nd, dpar = "a"))
show("fitted(dpar = 'b')", fitted(fit, newdata = nd, dpar = "b"))
show("fitted() mu", fitted(fit, newdata = nd))
show("predict(se.fit, dpar a)", frm_linpred(fit, newdata = nd, dpar = "a",
                                            se.fit = TRUE))
show("frm_lp_basis(dpar a) se_nonest", frm_lp_basis(fit, newdata = nd,
                                                   dpar = "a")$se_nonest)
show("frm_joint_cov V[1:3,1:3]", frm_joint_cov(fit)$V[1:3, 1:3])
# a linear ridge: an exactly duplicated column is dropped by the frame,
# so the lost direction comes from a bound instead; the control is a
# mo() fit whose lost coordinate has an exactly zero gradient
set.seed(71)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
dm <- data.frame(income, ls)
dm$age <- rnorm(100, mean = 40, sd = 10)
fm <- suppressWarnings(frm(ls ~ mo(income) * age, data = dm))
show("mo seed 71 fitted() Est.Error", fitted(fm, newdata = dm[1:3, ]))
ce <- withCallingHandlers(conditional_effects(fm, "income:age"),
                          warning = function(x) invokeRestart("muffleWarning"))
cat("mo seed 71 conditional_effects se__ finite:",
    sum(is.finite(ce[[1]]$se__)), "of", nrow(ce[[1]]), "\n")
