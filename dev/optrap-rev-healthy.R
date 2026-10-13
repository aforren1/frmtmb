# Reviewer of lane optrap, item 8: the two "healthy fit" tests the lane
# moved to a response of unit spread (test-se-tier.R: sleepstudy / 64;
# test-se-check.R: ls / 16). On the response in its own units: does the
# fit stay silent, and how far is vcov() from sdreport()'s cov.fixed?
# Also: does par_template() show the new starts (NEWS)?
#   Rscript dev/optrap-rev-healthy.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
for (e in parse("tests/testthat/test-se-check.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(deparse(e[[2]]), "se_mo_data")) eval(e, env)
}
said <- function(expr) {
  w <- character()
  m <- character()
  v <- withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")
    })
  list(value = v, n = length(w) + length(m), what = c(w, m))
}
cmp <- function(lab, r) {
  f <- r$value
  V <- vcov(f, full = TRUE)
  ref <- RTMB::sdreport(f$obj)$cov.fixed
  rel <- max(abs(V - ref) / sqrt(outer(diag(ref), diag(ref))))
  cat(sprintf("%-22s said %d | identical to sdreport %s | max |V - ref| / (se se) %.3g | units %s\n",
              lab, r$n, identical(unname(V), unname(ref)), rel,
              if (is.null(f$par_units)) "none" else
                paste(signif(unique(f$par_units), 3), collapse = ",")))
  for (x in r$what) cat("   ", substr(x, 1, 120), "\n")
}
ss <- sleepstudy
cmp("sleepstudy ms", said(frm(Reaction ~ Days + (Days | Subject), data = ss)))
ss$Reaction <- ss$Reaction / 64
cmp("sleepstudy / 64", said(frm(Reaction ~ Days + (Days | Subject),
                                data = ss)))
d <- env$se_mo_data(7)
cmp("mo ls", said(frm(ls ~ mo(income), data = d)))
d$ls <- d$ls / 16
cmp("mo ls / 16", said(frm(ls ~ mo(income), data = d)))
pt <- par_template(Reaction ~ Days + (Days | Subject), data = sleepstudy)
cat("par_template theta start:\n")
print(pt[grepl("theta|sd", pt$class %||% pt[[1]]) | TRUE, ][1:6, ])
