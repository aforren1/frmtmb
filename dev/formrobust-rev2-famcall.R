# Re-check item 5: a pooled update() of a family with a non-link option
# (huber(k =), categorical(levels =), hurdle_cumulative(threshold =)).
# Data seed 95. REVLIB="" for the base arm.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(95)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 0.5 * d$x + rt(n, 3)
d$ycat <- factor(sample(c("p", "q", "r"), n, TRUE))
print(names(formals(huber)))
tryf <- function(e) tryCatch(e, error = function(err) {
  cat("ERROR:", conditionMessage(err), "\n"); NULL })
run <- function(label, f0, newf) {
  cat("\n==", label, "==\n")
  if (is.null(f0)) return()
  fam0 <- f0$spec$responses[[1]]$family
  u <- tryf(suppressMessages(update(f0, newf)))
  if (is.null(u)) return()
  famu <- u$spec$responses[[1]]$family
  cat("call family:", deparse1(u$call$formula$family %||%
                                 u$call$formula[["family"]]), "\n")
  diffs <- setdiff(union(names(fam0), names(famu)), "links")
  for (fld in diffs) {
    a <- fam0[[fld]]; b <- famu[[fld]]
    if (!is.function(a) && !is.function(b) && !identical(a, b)) {
      cat("family field", fld, "differs: original", format(a), " updated",
          format(b), "\n")
    }
  }
  ref <- tryf(frm(bf(newf, sigma ~ z), data = d, family = fam0))
  if (!is.null(ref)) cat("updated logLik == direct fit with original family:",
                         identical(logLik(u), logLik(ref)), "\n")
}
f_h <- tryf(frm(bf(y ~ x, sigma ~ z, family = huber(k = 3)), data = d))
run("huber(k = 3) in bf()", f_h, y ~ x + z)
f_h2 <- tryf(frm(bf(y ~ x, sigma ~ z), family = huber(k = 3), data = d))
run("huber(k = 3) as frm(family =)", f_h2, y ~ x + z)
cat("\n== the updated huber fit against a default-k direct fit ==\n")
u <- suppressMessages(update(f_h, y ~ x + z))
r_def <- frm(bf(y ~ x + z, sigma ~ z, family = huber()), data = d)
r_k3 <- frm(bf(y ~ x + z, sigma ~ z, family = huber(k = 3)), data = d)
cat("update logLik", format(logLik(u), digits = 12), "; direct k default",
    format(logLik(r_def), digits = 12), "; direct k = 3",
    format(logLik(r_k3), digits = 12), "\n")
cat("identical to default-k fit:", identical(logLik(u), logLik(r_def)), "\n")
print(formals(huber)$k)
cat("\n== categorical(levels =) ==\n")
f_c <- tryf(frm(bf(ycat ~ x, muq ~ z, family = categorical(levels = c("r", "p", "q"))),
                data = d))
if (!is.null(f_c)) {
  print(rownames(fixef(f_c)))
  uc <- tryf(suppressMessages(update(f_c, ycat ~ x + z)))
  if (!is.null(uc)) {
    print(rownames(fixef(uc)))
    cat("call family:", deparse1(uc$call$formula[["family"]]), "\n")
  }
}
