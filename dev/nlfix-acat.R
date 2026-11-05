# Ruling (b): the acat ridge with the return-code arm gone. Run under R's BLAS and
# under OpenBLAS (the second by that build's Rscript). Data of nl_saddle_data().
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "BLAS:",
    La_library(), "\n")
set.seed(9202)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
d$y <- ordered(cut(1.5 * exp(0.6 * d$x) + rlogis(n),
                   c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
body <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
show <- function(lab, expr) {
  msg <- character(0); warn <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(r, "error")) { cat(lab, "ERROR", conditionMessage(r), "\n"); return(invisible(NULL)) }
  cat(sprintf("%-28s ll %.4f code %d | %s\n   msg: %s\n", lab, as.numeric(logLik(r)),
              r$opt$convergence, paste(round(r$estimates$beta, 3), collapse = ","),
              substr(paste(msg, collapse = " // "), 1, 260)))
  invisible(r)
}
show("acat default", frm(body, family = acat(), data = d))
show("acat se = FALSE", frm(body, family = acat(), data = d, se = FALSE))
show("acat start", frm(body, family = acat(), data = d,
                       start = list(beta = c("b1_(Intercept)" = 1.5,
                                             "b2_(Intercept)" = 0.6))))
for (fam in list(cumulative(), sratio(), cratio())) {
  show(fam$family, frm(body, family = fam, data = d))
}
ctl <- frmtmb_control(optCtrl = list(iter.max = 0))
r <- show("cumulative iter.max = 0", frm(body, family = cumulative(), data = d, control = ctl))
cat("start read:", r$estimates$beta, "\n")
ctl <- frmtmb_control(optCtrl = list(iter.max = 3))
show("cumulative iter.max = 3", frm(body, family = cumulative(), data = d, control = ctl))
show("acat check_se = ignore", frm(body, family = acat(), data = d,
     control = frmtmb_control(check_se = "ignore")))
