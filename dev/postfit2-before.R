# Lane postfit2, the "before" arm on the base build rellib-r3 (0.65.0):
# the claim that conditional_effects() "also covers what brms calls
# conditional_smooths()", measured. Data as in
# dev/postfit2-brms-compare.R (seed 1, n = 150).
#
#   Rscript dev/postfit2-before.R > dev/postfit2-log/before.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("frmtmb %s from %s", format(packageVersion("frmtmb")),
    find.package("frmtmb"))
say("exists conditional_smooths: %s",
    exists("conditional_smooths", envir = asNamespace("frmtmb")))
set.seed(1)
n <- 150
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(2 * pi * d$x) + d$z * (d$f == "b") + rnorm(n, 0, 0.3)
d$y2 <- 2 * exp(0.8 * d$z) * (1 + 0.3 * (d$f == "b")) + rnorm(n, 0, 0.3)
f1 <- frm(bf(y ~ f + s(x) + s(z, by = f)), family = gaussian(), data = d)
ce <- conditional_effects(f1, "x")$x
# the s(x) term alone, from the public delta-method seam: the columns of
# the design that belong to s(x), times their estimates
lb <- frm_lp_basis(f1, newdata = ce, re_formula = NA)
pm <- frmtmb:::joint_pos_map(frmtmb:::get_joint_cov(f1))
chat <- vapply(lb$coef_pos, function(k) {
  f1$estimates[[pm$comp[k]]][pm$idx[k]]
}, 0)
sx <- grepl("s(x)", lb$coef_names, fixed = TRUE)
term <- as.vector(lb$A[, sx, drop = FALSE] %*% chat[sx])
se_term <- sqrt(rowSums((lb$A[, sx] %*% lb$V[sx, sx]) * lb$A[, sx]))
say("conditional_effects(fit, 'x') estimate__ minus the s(x) term: min %.4f, max %.4f",
    min(ce$estimate__ - term), max(ce$estimate__ - term))
say("its se__ over the term's own delta-method se: min %.3f, max %.3f",
    min(ce$se__ / se_term), max(ce$se__ / se_term))
say("so the display is the expected response, not the term: no argument of conditional_effects() returns the term")
f4 <- frm(bf(y2 ~ a * exp(b * z), a ~ 1 + f, b ~ 1, nl = TRUE),
          family = gaussian(), data = d)
say("nonlinear fit, default band: %s",
    tryCatch({
      conditional_effects(f4, "z")
      "answered"
    }, error = conditionMessage))
for (a in c("spaghetti", "select_points", "too_far")) {
  args <- list(f1, "x", 0.1)
  names(args) <- c("", "", a)
  if (a == "spaghetti") args[[3]] <- TRUE
  say("%s: %s", a, tryCatch({
    do.call(conditional_effects, args)
    "answered"
  }, error = conditionMessage))
}
say("surface = TRUE: %s", tryCatch({
  conditional_effects(f1, "x:z", surface = TRUE)
  "answered"
}, error = conditionMessage))
for (fn in c("make_conditions", "update_adterms")) {
  say("exists %s: %s", fn, exists(fn, envir = asNamespace("frmtmb")))
}
