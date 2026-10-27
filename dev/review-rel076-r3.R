# Review: does release fix R3's bound in test-mo-coding.R fail when
# marginaleffects loses the simplex? Run the test block's own
# computation as written, then again with a mutant interop_extra_comps()
# that always leaves the mo() simplex out (0.75.0's behavior).
#   Rscript dev/review-rel076-r3.R > dev/review-rel076-out/r3.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
rel <- "C:/Users/adf44/source/r/frmtmb-wt-rel076"
# the file's data helpers
src <- readLines(file.path(rel, "tests/testthat/test-mo-coding.R"))
env <- new.env()
for (h in list.files(file.path(rel, "tests/testthat"), "^helper.*[.]R$",
                     full.names = TRUE)) sys.source(h, env)
ex <- parse(text = src)
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("mo_code_data", "mo_gap_fit")) {
    eval(e, env)
  }
}
run <- function(lab) {
  d <- env$mo_code_data()
  fit <- frm(y ~ mo(k) + z, data = d)
  m <- stats::lm(y ~ factor(k) + z, data = d)
  nd <- data.frame(k = 1:4, z = 0)
  ex <- stats::predict(m, newdata = nd, se.fit = TRUE)$se.fit *
    sqrt(m$df.residual / nrow(d))
  se <- marginaleffects::predictions(fit, newdata = nd)$std.error
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
  jc <- get_joint_cov(fit)
  A <- as.matrix(lb$A)
  Vb <- as.matrix(jc$V)[lb$coef_pos, lb$coef_pos, drop = FALSE]
  zc <- which(jc$names[lb$coef_pos] == "zeta1")
  A0 <- A[, -zc, drop = FALSE]
  se0 <- sqrt(rowSums((A0 %*% Vb[-zc, -zc, drop = FALSE]) * A0))
  lhs <- max(abs(se / ex - 1)); rhs <- 1e-4 * max(abs(se0 / ex - 1))
  cat(lab, ": se", signif(se, 6), "| exact", signif(ex, 6), "| se0",
      signif(se0, 6), "\n  max|se/ex - 1| =", signif(lhs, 4),
      " bound =", signif(rhs, 4), " test", if (lhs < rhs) "PASSES" else
        "FAILS", "\n")
}
environment(run) <- list2env(list(env = env), parent = asNamespace("frmtmb"))
run("as built")
ns <- asNamespace("frmtmb")
orig <- get("interop_extra_comps", ns)
mut <- function(model, simplex = TRUE) orig(model, simplex = FALSE)
unlockBinding("interop_extra_comps", ns)
assign("interop_extra_comps", mut, envir = ns)
run("mutant: simplex left out")
assign("interop_extra_comps", orig, envir = ns)
