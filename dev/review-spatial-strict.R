# Reviewer: how strict is the seam? Slips a registrant can make, and
# what the fit does with each.
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
try1 <- function(label, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 300)))
  if (is.character(r) && length(r) == 1L && startsWith(r, "ERROR")) {
    cat(label, "->", r, "\n")
  } else if (inherits(r, "frmtmb_fit")) {
    cat(label, "-> FIT logLik", format(as.numeric(logLik(r)), digits = 10),
        "\n")
  } else cat(label, "-> ok\n")
  invisible(r)
}
lp_iid <- function(b, theta, blk) {
  sum(RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE))
}
vc_iid <- function(theta, blk) {
  d <- blk[["dim"]]
  V <- diag(exp(2 * theta[1]), d)
  dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
  V
}
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
d$y <- d$x + rnorm(20)[d$g] + rnorm(100)
ref <- frm(y ~ x + (1 | g), data = d)
cat("reference (1 | g) logLik", format(as.numeric(logLik(ref)), digits = 10),
    "\n")
reg <- function(nm, ...) {
  try1(paste("register", nm), frmtmb_register_covstruct(nm, ...))
}
reg("ok1", lpdf = lp_iid, vcov = vc_iid, npar = 1L, sd_idx = 1L)
try1("fit ok1", frm(y ~ x + ok1(1 | g), data = d))
reg("badstart", lpdf = lp_iid, vcov = vc_iid, npar = 1L, start = c(0, 0, 0))
try1("fit badstart (start length 3, npar 1)",
     frm(y ~ x + badstart(1 | g), data = d))
reg("badsd", lpdf = lp_iid, vcov = vc_iid, npar = 1L, sd_idx = 4L)
try1("fit badsd (sd_idx 4, npar 1)", frm(y ~ x + badsd(1 | g), data = d))
reg("badvc", lpdf = lp_iid, vcov = function(theta, blk) matrix(1, 3, 3),
    npar = 1L, sd_idx = 1L)
f <- try1("fit badvc (vcov 3 x 3 for dim 1)",
          frm(y ~ x + badvc(1 | g), data = d))
if (inherits(f, "frmtmb_fit")) {
  try1("VarCorr badvc", print(VarCorr(f)))
  try1("predict new level badvc", print(predict(
    f, newdata = data.frame(g = "new", x = 0), allow_new_levels = TRUE)))
}
reg("veclp", lpdf = function(b, theta, blk) {
  RTMB::dnorm(b, 0, exp(theta[1]), log = TRUE)
}, vcov = vc_iid, npar = 1L, sd_idx = 1L)
try1("fit veclp (lpdf returns a vector)", frm(y ~ x + veclp(1 | g),
                                              data = d))
reg("signflip", lpdf = function(b, theta, blk) -lp_iid(b, theta, blk),
    vcov = vc_iid, npar = 1L, sd_idx = 1L)
try1("fit signflip (lpdf returns the NEGATIVE log density)",
     frm(y ~ x + signflip(1 | g), data = d))
reg("np_fun_bad", lpdf = lp_iid, vcov = vc_iid,
    npar = function(blk) -1L)
try1("fit np_fun_bad (npar(blk) = -1)",
     frm(y ~ x + np_fun_bad(1 | g), data = d))
# names
reg("log", lpdf = lp_iid, vcov = vc_iid, npar = 1L, sd_idx = 1L)
try1("y ~ log(x + 3) after registering log",
     frm(y ~ log(x + 3) + (1 | g), data = d))
try1("y ~ x + log(1 | g)", frm(y ~ x + log(1 | g), data = d))
reg("poly", lpdf = lp_iid, vcov = vc_iid, npar = 1L, sd_idx = 1L)
try1("y ~ poly(x, 2) after registering poly",
     frm(y ~ poly(x, 2) + (1 | g), data = d))
reg("cor_ar", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("ar", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("Intercept", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("sd", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("me2", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("trials", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("dec", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("weibull", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("wiener", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("cox", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("sar", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("fcor", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("gp", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
reg("car", lpdf = lp_iid, vcov = vc_iid, npar = 1L)
try1("y ~ x + car(1 | g) after a bar car registered from outside",
     frm(y ~ x + car(1 | g), data = d))
# same name from two places: the second silently replaces the first?
e <- new.env()
reg("ok1", lpdf = lp_iid, vcov = function(theta, blk) {
  matrix(99, 1, 1, dimnames = list(blk[["cnms"]], blk[["cnms"]]))
}, npar = 1L, sd_idx = 1L)
f <- try1("fit ok1 after re-registration", frm(y ~ x + ok1(1 | g), data = d))
if (inherits(f, "frmtmb_fit")) print(VarCorr(f))
# a package namespace registering a name another package holds
ns1 <- new.env(); ns2 <- new.env()
cat("entries' package slot:", frmtmb:::cs_entry("ok1")[["package"]] %||% "NULL",
    "\n")
cat("cs_check_entry over every registered entry:\n")
for (nm in frmtmb:::cs_names()) {
  m <- frmtmb:::cs_check_entry(nm, frmtmb:::cs_entry(nm))
  if (!is.null(m)) cat("  ", m, "\n")
}
cat("done\n")
