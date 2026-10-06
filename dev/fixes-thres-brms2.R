# Lane fixes, item 4a: brms's draw labels under draw_ids and ndraws for
# posterior_linpred(incl_thres = TRUE), and its answer with dpar given.
#   Rscript dev/fixes-thres-brms2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
suppressMessages({library(testthat); library(frmtmb)})
henv <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = henv)
}
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
fit <- frm(y ~ x, family = cumulative(), data = d)
bb <- henv$brms_fixed_fit(brms::bf(y ~ x), brms::cumulative(), d, fit,
                          ndraws = 6)
a <- brms::posterior_linpred(bb, incl_thres = TRUE, draw_ids = c(2, 5))
cat("draw_ids = c(2, 5): dim", dim(a), "\n"); str(dimnames(a))
a <- brms::posterior_linpred(bb, incl_thres = TRUE, ndraws = 3)
cat("ndraws = 3: dim", dim(a), "\n"); str(dimnames(a))
a <- brms::posterior_linpred(bb, incl_thres = TRUE, dpar = "mu")
cat("dpar = 'mu': dim", dim(a), "\n")
a <- brms::posterior_linpred(bb, incl_thres = TRUE, transform = TRUE)
cat("transform = TRUE: dim", dim(a), "\n")
a <- tryCatch(brms::posterior_linpred(bb, incl_thres = NA),
              error = function(e) conditionMessage(e))
cat("incl_thres = NA:", if (is.character(a)) a else dim(a), "\n")
a <- brms::posterior_linpred(bb, incl_thres = TRUE,
                             newdata = data.frame(x = c(0, 1)))
cat("newdata 2 rows: dim", dim(a), "\n"); print(a[1, , ])
cat("frmtmb at the same estimates:\n")
print(frmtmb:::ord_thres_linpred(fit, newdata = data.frame(x = c(0, 1))))
