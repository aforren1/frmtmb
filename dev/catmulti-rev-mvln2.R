# Reviewer: lncor names for a response whose name is not a Stan name.
# Log dev/catmulti-rev-out/mvln2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(46)
n <- 300
dd <- data.frame(x = stats::rnorm(n))
W <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .4, .4, 1), 2))
E <- cbind(0, 0.4 + 0.5 * dd$x + 0.6 * W[, 1], -0.3 - 0.4 * dd$x + 0.9 * W[, 2])
Y <- exp(E) / rowSums(exp(E)); colnames(Y) <- c("a", "b", "c")
dd$Y_a <- Y; dd$a <- Y[, 1]; dd$b <- Y[, 2]; dd$c <- Y[, 3]
dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
TRY <- function(nm, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(nm, ":\n"); print(r); invisible(r)
}
f2 <- frm(mvbf(bf(Y_a ~ x, family = logistic_normal()), bf(w ~ x)), data = dd)
dp <- default_prior(f2)
print(as.data.frame(dp)[, c("class", "coef", "resp", "dpar")])
TRY("hypothesis lncor_Y_a__b__c", hypothesis(f2, "lncor_Y_a__b__c = 0", class = NULL)$hypothesis[, 1:3])
TRY("hypothesis lncor_Ya__b__c", hypothesis(f2, "lncor_Ya__b__c = 0", class = NULL)$hypothesis[, 1:3])
TRY("set_prior resp Ya", summary(frm(mvbf(bf(Y_a ~ x, family = logistic_normal()), bf(w ~ x)), data = dd,
        prior = set_prior("lkj(2)", class = "lncor", resp = "Ya")))$spec_pars)
TRY("set_prior sigma resp Ya", rownames(summary(frm(mvbf(bf(Y_a ~ x, family = logistic_normal()), bf(w ~ x)),
        data = dd, prior = set_prior("normal(0, 1)", class = "b", resp = "Ya")))$spec_pars))
TRY("confint alias lncor_Y_a__b__c", confint(f2, parm = "lncor_Y_a__b__c"))
f3 <- frm(mvbf(bf(cbind(a, b, c) ~ x, family = logistic_normal()), bf(w ~ x)), data = dd)
TRY("cbind variables", variables(f3))
TRY("cbind hypothesis", hypothesis(f3, "lncor_cbind(a, b, c)__b__c = 0", class = NULL)$hypothesis[, 1:3])
TRY("cbind summary spec rows", rownames(summary(f3)$spec_pars))
cat("brms names (stancode and default_prior, no compile):\n")
bf2 <- brms::bf(Y_a ~ x, family = brms::logistic_normal()) + brms::bf(w ~ x) + brms::set_rescor(FALSE)
sc <- strsplit(brms::make_stancode(bf2, data = dd), "\n")[[1]]
cat(grep("ncor_", sc, value = TRUE), sep = "\n")
bdp <- as.data.frame(brms::default_prior(bf2, data = dd))
print(bdp[, c("prior", "class", "resp", "dpar")])
bf3 <- brms::bf(cbind(a, b, c) ~ x, family = brms::logistic_normal()) + brms::bf(w ~ x) + brms::set_rescor(FALSE)
r <- tryCatch(unique(as.data.frame(brms::default_prior(bf3, data = dd))$resp),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("brms cbind resp names:", r, "\n")
cat("DONE\n")
