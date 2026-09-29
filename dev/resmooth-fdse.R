# Lane wt-resmooth. A defect the fix EXPOSES and does not close: the
# finite-difference standard-error route (a category probability, and an
# autocor cov = FALSE quantity) differences no smooth coefficients at
# re_formula = NA, so the Est.Error it reports there leaves out the
# uncertainty of a curve that is now in the estimate. The analytic route
# carries it (lp_delta_A() adds a b column per smooth part), so the two
# routes disagree. Through 0.64.0 the group-indexed smooth was dropped
# from the estimate as well, which was at least self-consistent; a
# POPULATION smooth has had this gap all along.
#   Rscript dev/resmooth-fdse.R > dev/resmooth-fdse.txt
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")
set.seed(21)
ng <- 6
d <- data.frame(g = factor(rep(seq_len(ng), each = 25)),
                x = stats::runif(ng * 25, -2, 2))
lat <- stats::rnorm(ng, 0, 1)[d$g] * sin(d$x) + 0.5 * d$x +
  stats::rlogis(nrow(d))
d$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
              ordered = TRUE)
fit <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)) + cumulative(),
                            data = d))
nd <- d[c(3, 40, 90), c("x", "g")]
bl <- fit$frame[["re_blocks"]]
sm_b <- sort(unique(unlist(lapply(bl[vapply(bl, function(b) {
  b[["covstruct"]] %in% c("smooth", "gp", "hsgp")
}, NA)], `[[`, "b_idx"))))
f <- function(x) frmtmb:::fitted_point(x, nd, re_formula = NA)
got <- fitted(fit, newdata = nd, re_formula = NA)[, "Est.Error", ]
ref <- frmtmb:::fit_fd_se(fit, f, b_idx = sm_b)
cat("smooth b positions differenced by the reference:", length(sm_b), "\n")
cat("shipped Est.Error at re_formula = NA (rows x categories):\n")
print(round(got, 6))
cat("reference that also differences the smooth's coefficients:\n")
print(round(matrix(as.vector(ref), nrow = nrow(nd)), 6))
rel <- abs(as.vector(ref) - as.vector(got)) / abs(as.vector(ref))
cat(sprintf("max relative gap %.3e, median %.3e\n", max(rel),
            stats::median(rel)))
# the value the route gave before: no b differenced at all. This is the
# shortfall the fix removes, and the margin the test asserts.
bare <- as.vector(frmtmb:::fit_fd_se(fit, f, b_idx = NULL))
sh <- (as.vector(ref) - bare) / as.vector(ref)
cat("share of the SE that the smooth's coefficients supply: ")
cat(sprintf("median %.4f, min %.4f, max %.4f\n", stats::median(sh),
            min(sh), max(sh)))
# the analytic route, for contrast: it DOES carry the smooth's b
se_na <- frm_linpred(fit, newdata = nd, re_formula = NA, se.fit = TRUE)
se_null <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
cat("analytic se(eta) at NA :", sprintf("%.6f", se_na$se.fit), "\n")
cat("analytic se(eta) at NULL:", sprintf("%.6f", se_null$se.fit), "\n")
