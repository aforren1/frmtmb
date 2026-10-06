# Reviewer of lane fixes, claim 4a follow-up: NA pattern of grouped
# thresholds, what brms's hurdle_cumulative columns are, and
# frmtmb.sample's array under draw_ids and newdata.
#   Rscript dev/fixes-rev-thres2.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = normalizePath("dev/fixes-rev-stan-cache"))
suppressMessages({library(testthat); library(frmtmb); library(frmtmb.sample)})
cat("LIB", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
henv <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path("tests/testthat", h), envir = henv)
}
brms_fixed_fit <- henv$brms_fixed_fit
set.seed(8080)
n <- 220
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
fit <- frm(yg | thres(gr = g) ~ x, family = cumulative(), data = d)
bb <- brms_fixed_fit(brms::bf(yg | thres(gr = g) ~ x), brms::cumulative(), d,
                     fit, ndraws = 3)
a <- frmtmb:::ord_thres_linpred(fit)
b <- unname(brms::posterior_linpred(bb, incl_thres = TRUE)[1, , ])
cat("grouped: NA pattern identical", identical(is.na(a), is.na(b)),
    " NA count", sum(is.na(a)), " rows of level b", sum(d$g == "b"),
    " max abs/rel diff", max(abs(a - b) / pmax(1, abs(b)), na.rm = TRUE), "\n")

# hurdle_cumulative: brms's columns 1..K against (1 - hu) * (tau - eta)
fh <- frm(yh ~ x, family = hurdle_cumulative(), data = d)
bh <- brms_fixed_fit(brms::bf(yh ~ x), brms::hurdle_cumulative(), d, fh,
                     ndraws = 3)
ph <- brms::posterior_linpred(bh, incl_thres = TRUE)[1, , ]
hu <- as.numeric(frm_linpred(fh, dpar = "hu", type = "response"))
eta <- as.numeric(frm_linpred(fh, dpar = "mu", type = "link"))
fe <- fixef(fh)[, "Estimate"]
tau <- fe[grep("^Intercept", names(fe))]
mine <- (1 - hu) * outer(-eta, tau, "+")
cat("hurdle: col 0 equals hu:", isTRUE(all.equal(unname(ph[, 1]), hu)),
    " cols 1..K equal (1 - hu) * (tau - eta):",
    isTRUE(all.equal(unname(ph[, -1]), unname(mine))), "\n")
# brms's own linpred without thresholds, for reference
cat("hurdle incl_thres = FALSE dim:",
    dim(brms::posterior_linpred(bh, incl_thres = FALSE)), "\n")

# frmtmb.sample: draws, draw_ids, newdata
fs <- frm(bf(y ~ x + cs(z), disc ~ 0 + z), family = cratio(), data = d)
ds <- suppressMessages(suppressWarnings(
  frm_sample(fs, chains = 1, iter = 200, refresh = 0, seed = 5)))
pl <- posterior_linpred(ds, incl_thres = TRUE, draw_ids = c(2, 5))
cat("sample: dim", dim(pl), " dimnames:", deparse1(dimnames(pl)), "\n")
idx <- frmtmb.sample:::draws_par_index(ds$fit)
f5 <- frmtmb.sample:::draws_fit_at(ds, 5, idx)
cat("draw 5 equals ord_thres_linpred at draw 5:",
    identical(unname(pl[2, , ]), frmtmb:::ord_thres_linpred(f5)), "\n")
nd <- data.frame(x = c(-1, 2), z = c(0.3, -0.4), y = 1L)
pn <- posterior_linpred(ds, incl_thres = TRUE, newdata = nd, ndraws = 3)
cat("newdata dim", dim(pn), "\n")
r <- tryCatch(posterior_linpred(ds, incl_thres = NA),
              error = function(e) conditionMessage(e))
cat("incl_thres = NA:", if (is.character(r)) r else "returned", "\n")
cat("transform = TRUE ignores it, dim:",
    dim(posterior_linpred(ds, incl_thres = TRUE, transform = TRUE,
                          ndraws = 2)), "\n")
