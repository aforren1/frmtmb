# Reviewer: subset(), K = 12, dirichlet2 newdata and links, multinomial
# parser arms, confint aliases. Output: dev/famsimplex-rev-misc.log
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
  if (length(w)) cat("WARNINGS:", unique(w), sep = "\n  ")
  if (is.character(r) && length(r) == 1) cat(r, "\n") else print(r)
  invisible(r)
}
# 1. subset() with an off-simplex row the subset removes
dd <- fs_data_dirichlet(seed = 21)
d <- dd
d$keep <- TRUE
d$keep[2] <- FALSE
f_all <- frm(bf(Y ~ x), family = dirichlet(), data = dd[-2, ])
show("fit on rows 1, 3:300", c(nobs(f_all), logLik(f_all)))
show("subset(keep), row 2 valid", {
  f <- frm(bf(Y | subset(keep) ~ x), family = dirichlet(), data = d)
  c(nobs(f), logLik(f))
})
d$Y[2, ] <- c(0.9, 0.9, 0.9)
show("subset(keep), row 2 off the simplex", {
  f <- frm(bf(Y | subset(keep) ~ x), family = dirichlet(), data = d)
  c(nobs(f), logLik(f))
})
show("gaussian subset(keep) nobs", {
  d$yy <- rnorm(300)
  nobs(frm(bf(yy | subset(keep) ~ x), family = gaussian(), data = d))
})
# 2. K = 12, logLik against the hand density
set.seed(22)
K <- 12
d12 <- data.frame(x = rnorm(400))
sl <- seq(-0.5, 0.6, length.out = K - 1)
ic <- seq(-1, 1, length.out = K - 1)
E <- cbind(0, outer(d12$x, sl) + matrix(ic, 400, K - 1, byrow = TRUE))
d12$Y <- fs_rdirichlet(fs_softmax(E) * 30)
f12 <- frm(bf(Y ~ x), family = dirichlet(), data = d12)
b <- fixef(f12)[, "Estimate"]
Eh <- cbind(0, vapply(2:K, function(k) {
  b[[paste0("mu", k, "_Intercept")]] + b[[paste0("mu", k, "_x")]] * d12$x
}, numeric(400)))
phi <- summary(f12)$spec_pars["phi", "Estimate"]
A <- fs_softmax(Eh) * phi
ref <- sum(lgamma(rowSums(A)) - rowSums(lgamma(A)) +
             rowSums((A - 1) * log(d12$Y)))
cat(sprintf("K=12 logLik %.10f hand %.10f rel %.3e phi %.3f\n",
            logLik(f12), ref, abs(logLik(f12) - ref) / abs(ref), phi))
cat("max |slope - truth| / se:",
    max(abs(b[paste0("mu", 2:K, "_x")] - sl) /
          fixef(f12)[paste0("mu", 2:K, "_x"), "Est.Error"]), "\n")
show("K=12 fitted dimnames", dimnames(fitted(f12, newdata = d12[1:2, ,
                                                        drop = FALSE]))[[3]])
show("K=12 refcat \"10\"", rownames(fixef(frm(bf(Y ~ x),
     family = dirichlet(refcat = "10"), data = d12)))[1:10])
# 3. dirichlet2 newdata and links
d2 <- fs_data_dirichlet2(seed = 15)
f2 <- frm(bf(Y ~ x), family = brmsfamily("dirichlet2"), data = d2)
nd <- d2[c(5, 1, 3), "x", drop = FALSE]
show("dir2 fitted newdata", fitted(f2, newdata = nd))
show("dir2 predict newdata dims", dimnames(predict(f2, newdata = nd)))
for (lk in c("log", "softplus", "squareplus", "identity", "logm1")) {
  show(paste("dirichlet2 link", lk), {
    f <- frm(bf(Y ~ x), family = brmsfamily("dirichlet2", lk), data = d2)
    c(logLik = as.numeric(logLik(f)), conv = f$opt$convergence %||% NA)
  })
}
for (lk in c("log", "identity", "softplus", "squareplus")) {
  show(paste("dirichlet link_phi", lk), {
    f <- frm(bf(Y ~ x), family = dirichlet(link_phi = lk), data = dd)
    as.numeric(logLik(f))
  })
}
show("ln link_sigma identity", as.numeric(logLik(frm(bf(Y ~ x),
     family = logistic_normal(link_sigma = "identity"),
     data = fs_data_logistic_normal(seed = 23)))))
# 4. confint alias for the one correlation
fl <- frm(bf(Y ~ x), family = logistic_normal(),
          data = fs_data_logistic_normal(seed = 13))
show("confint parm lncor__b__c", confint(fl, parm = "lncor__b__c"))
show("confint parm lncor_raw_1 profile", confint(fl, parm = "lncor_raw_1",
                                                 method = "profile"))
show("confint_varcorr", confint_varcorr(fl))
show("vcov full names", rownames(vcov(fl, full = TRUE)))
# 5. names with special characters
d <- dd
colnames(d$Y) <- c("a.b", "ab", "c_d")
show("names a.b, ab, c_d", frm(bf(Y ~ x), family = dirichlet(), data = d))
colnames(d$Y) <- c("a b", "x-1", "c")
show("names 'a b', 'x-1', c", rownames(fixef(frm(bf(Y ~ x),
                                                 family = dirichlet(),
                                                 data = d))))
# 6. CSV-style rounding
d <- dd
d$Y <- round(dd$Y, 4)
show("rounded to 4 decimals", nobs(frm(bf(Y ~ x), family = dirichlet(),
                                       data = d)))
# 7. a single-row newdata and K = 2 shapes
show("dir fitted 1-row newdata dims", dim(fitted(frm(bf(Y ~ x),
     family = dirichlet(), data = dd), newdata = dd[1, , drop = FALSE])))
