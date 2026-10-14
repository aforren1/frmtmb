# Reviewer: response validation and refcat on the lane build beside brms
# 2.23.0's standata(). Output: dev/famsimplex-rev-valid.log
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
bsd <- function(f, d, fam) {
  tryCatch({
    s <- brms::standata(f, data = d, family = fam)
    paste0("brms accepts; N = ", s$N, "; ncat = ", s$ncat,
           "; refcat = ", s$refcat %||% "-")
  }, error = function(e) paste("brms ERROR:", conditionMessage(e)),
  warning = function(w) paste("brms WARNING:", conditionMessage(w)))
}
dd <- fs_data_dirichlet(seed = 21)
Y <- dd$Y
# 1. refcat by index, on named and on unnamed columns
d1 <- dd
d1$Y <- unname(Y)
cat("\n=== refcat = 2, unnamed columns\n")
cat(bsd(brms::bf(Y ~ x), d1, brms::dirichlet(refcat = 2)), "\n")
show("frm dirichlet(refcat = 2)", fixef(frm(bf(Y ~ x),
                                            family = dirichlet(refcat = 2),
                                            data = d1)))
show("frm dirichlet(refcat = \"2\")",
     fixef(frm(bf(Y ~ x), family = dirichlet(refcat = "2"), data = d1)))
show("frm brms::dirichlet(refcat = 2)",
     fixef(frm(bf(Y ~ x), family = brms::dirichlet(refcat = 2), data = d1)))
cat("\n=== refcat = 2, named columns\n")
cat(bsd(brms::bf(Y ~ x), dd, brms::dirichlet(refcat = 2)), "\n")
show("frm", frm(bf(Y ~ x), family = dirichlet(refcat = 2), data = dd))
cat("\n=== logistic_normal(refcat = NA)\n")
show("brms", brms::logistic_normal(refcat = NA))
show("frm", logistic_normal(refcat = NA))
cat("\n=== dirichlet(refcat = NA_character_)\n")
show("brms", brms::dirichlet(refcat = NA_character_)$refcat)
show("frm", dirichlet(refcat = NA_character_)$simplex)
# 2. a value of exactly 1 with zeros; a row of NA; a row partly NA
for (case in c("one", "allNA", "rowsum2")) {
  Yb <- Y
  Yb[1, ] <- switch(case, one = c(1, 0, 0), allNA = c(NA, NA, NA),
                    rowsum2 = c(1, 0.5, 0.5))
  d <- dd
  d$Y <- Yb
  cat("\n===", case, "\n")
  cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
  show("frm dirichlet", nobs(frm(bf(Y ~ x), family = dirichlet(), data = d)))
  show("frm ln", nobs(frm(bf(Y ~ x), family = logistic_normal(), data = d)))
}
# 3. an NA in a predictor only, and an NA row that breaks the sum
d <- dd
d$x[2] <- NA
d$Y[2, ] <- c(0.9, 0.9, 0.9)
cat("\n=== off-simplex row whose predictor is NA\n")
cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
show("frm", nobs(frm(bf(Y ~ x), family = dirichlet(), data = d)))
# 4. subset() dropping the bad row
d <- dd
d$keep <- TRUE
d$keep[2] <- FALSE
d$Y[2, ] <- c(0.9, 0.9, 0.9)
cat("\n=== subset() drops the off-simplex row\n")
cat(bsd(brms::bf(Y | subset(keep) ~ x), d, brms::dirichlet()), "\n")
show("frm", nobs(frm(bf(Y | subset(keep) ~ x), family = dirichlet(),
                     data = d)))
# 5. a data.frame column holding a data.frame; a tibble
d <- dd
d$Y <- NULL
d$Y <- as.data.frame(Y)
cat("\n=== data.frame column\n")
cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
show("frm", fixef(frm(bf(Y ~ x), family = dirichlet(), data = d)))
# 6. K = 12 categories, unnamed: the order of categories and dpars
set.seed(22)
K <- 12
d12 <- data.frame(x = rnorm(400))
E <- cbind(0, outer(d12$x, seq(-0.5, 0.6, length.out = K - 1)) +
             matrix(seq(-1, 1, length.out = K - 1), 400, K - 1, byrow = TRUE))
d12$Y <- fs_rdirichlet(fs_softmax(E) * 30)
cat("\n=== K = 12 unnamed\n")
s12 <- brms::standata(brms::bf(Y ~ x), data = d12, family = brms::dirichlet())
cat("brms dpars:", names(brms::brmsterms(brms::bf(Y ~ x),
                                         family = brms::dirichlet())$dpars),
    "\n")
cat("brms ncat:", s12$ncat, "\n")
cat("brms X names:", grep("^X_", names(s12), value = TRUE), "\n")
f12 <- frm(bf(Y ~ x), family = dirichlet(), data = d12)
cat("frm coef:", rownames(fixef(f12)), "\n")
b <- fixef(f12)[, "Estimate"]
Eh <- cbind(0, vapply(2:K, function(k) {
  b[[paste0("mu", k, "_Intercept")]] + b[[paste0("mu", k, "_x")]] * d12$x
}, numeric(400)))
phi <- exp(b[["phi_Intercept"]] %||% NA)
if (is.na(phi)) phi <- summary(f12)$spec_pars["phi", "Estimate"]
A <- fs_softmax(Eh) * phi
ref <- sum(lgamma(rowSums(A)) - rowSums(lgamma(A)) +
             rowSums((A - 1) * log(d12$Y)))
cat(sprintf("K=12 logLik %.10f hand %.10f rel %.3e\n", logLik(f12), ref,
            abs(logLik(f12) - ref) / abs(ref)))
cat("truth slopes:", round(seq(-0.5, 0.6, length.out = K - 1), 3), "\n")
cat("est slopes:  ", round(b[paste0("mu", 2:K, "_x")], 3), "\n")
show("K=12 refcat = \"10\"", rownames(fixef(frm(bf(Y ~ x),
                              family = dirichlet(refcat = "10"),
                              data = d12)))[1:3])
# 7. special characters in category names
d <- dd
colnames(d$Y) <- c("a.b", "ab", "c_d")
cat("\n=== names a.b, ab, c_d\n")
cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
show("frm", frm(bf(Y ~ x), family = dirichlet(), data = d))
colnames(d$Y) <- c("a b", "x-1", "c")
cat("\n=== names 'a b', 'x-1', c\n")
show("brms dpars", names(brms::brmsterms(brms::bf(Y ~ x),
     family = brms::validate_family(brms::dirichlet()))$dpars))
cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
show("frm", rownames(fixef(frm(bf(Y ~ x), family = dirichlet(), data = d))))
# 8. ln row off by 1e-7: brms fits; frm refuses (deliberate)
d <- fs_data_logistic_normal(seed = 23)
d$Y[1, 1] <- d$Y[1, 1] + 1e-7
cat("\n=== ln, one row off by 1e-7\n")
cat(bsd(brms::bf(Y ~ x), d, brms::logistic_normal()), "\n")
show("frm", frm(bf(Y ~ x), family = logistic_normal(), data = d))
# 9. typical CSV data: proportions rounded to 4 decimals, renormalized
d <- dd
d$Y <- round(Y, 4)
cat("\n=== rounded to 4 decimals\n")
cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
show("frm", nobs(frm(bf(Y ~ x), family = dirichlet(), data = d)))
d$Y <- d$Y / rowSums(d$Y)
cat(bsd(brms::bf(Y ~ x), d, brms::dirichlet()), "\n")
show("frm renormalized", nobs(frm(bf(Y ~ x), family = dirichlet(), data = d)))
# 10. dirichlet2 under each link brms lists
d2 <- fs_data_dirichlet2(seed = 24)
for (lk in c("log", "softplus", "squareplus", "identity", "logm1")) {
  show(paste("dirichlet2 link", lk), {
    f <- frm(bf(Y ~ x), family = brmsfamily("dirichlet2", lk), data = d2)
    c(logLik = as.numeric(logLik(f)), conv = f$opt$convergence %||% NA)
  })
}
# 11. dirichlet link_phi each
for (lk in c("log", "identity", "softplus", "squareplus")) {
  show(paste("dirichlet link_phi", lk), {
    f <- frm(bf(Y ~ x), family = dirichlet(link_phi = lk), data = dd)
    as.numeric(logLik(f))
  })
}
show("ln link_sigma identity", as.numeric(logLik(frm(bf(Y ~ x),
     family = logistic_normal(link_sigma = "identity"),
     data = fs_data_logistic_normal(seed = 23)))))
