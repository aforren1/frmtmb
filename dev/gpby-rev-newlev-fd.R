# Reviewer: fitted()'s category Est.Error at an unseen level of (1 | g)
# beside an exact gp() past the data, against the closed form
# sqrt(plain^2 + (dP/deta)^2 * extra_var) with dP/deta of the cumulative
# logit written out and extra_var from frm_lp_basis().
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(11)
n <- 120
d <- data.frame(x = round(runif(n, 0, 5), 1),
                g = factor(sample(letters[1:8], n, TRUE)))
lat <- sin(d$x) + rnorm(8)[as.integer(d$g)] + rlogis(n)
d$y <- factor(cut(lat, quantile(lat, c(0, 0.3, 0.6, 1)),
                  include.lowest = TRUE, labels = FALSE), ordered = TRUE)
fit <- suppressWarnings(frm(bf(y ~ gp(x) + (1 | g)), data = d,
                            family = cumulative()))
nd <- data.frame(x = c(5.6, 6.3, 2.55), g = factor(c("zz", "zz", "yy")))
fv <- fitted(fit, newdata = nd, allow_new_levels = TRUE)[, "Est.Error", ]
f <- function(fx) ns$fitted_point(fx, nd, NULL, "response", NULL, NULL, TRUE)
b_idx <- ns$smooth_b_idx(fit)
plain <- ns$fit_fd_se(fit, f, b_idx = b_idx,
                      b_batch = ns$re_b_batches(fit, nd, NULL, TRUE, b_idx))
ev <- frm_lp_basis(fit, newdata = nd, allow_new_levels = TRUE)$extra_var
eta <- as.numeric(frm_linpred(fit, newdata = nd, type = "link",
                              allow_new_levels = TRUE))
fx <- fixef(fit)
tau <- fx[grepl("^Intercept", rownames(fx)), "Estimate"]
gg <- cbind(-dlogis(tau[1] - eta), dlogis(tau[1] - eta) - dlogis(tau[2] - eta),
            dlogis(tau[2] - eta))
want <- sqrt(plain^2 + gg^2 * ev)
sdg <- frmtmb::VarCorr(fit)$g$sd[1, "Estimate"]
cat(sprintf("extra_var %s | sd_g^2 %.6f\n", paste(sprintf("%.6f", ev), collapse = " "), sdg^2))
cat(sprintf("Est.Error / closed form: max |ratio - 1| %.3e\n", max(abs(fv / want - 1))))
