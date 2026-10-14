# Why eam's density at brms draws differs from brms::log_lik() by 2.9 in
# dev/famuni-wiener-fit2.R. Reads the cached brmsfit.
# Output: dev/famuni-out/wiener-lldbg.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
suppressMessages(library(brms))
bfit <- readRDS(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out/wiener-brmsfit.rds")
dat <- bfit$data
dr <- as_draws_df(bfit)
k <- 7L
ll1 <- brms::log_lik(bfit, draw_ids = k)
mu <- dr$b_Intercept[k] + dr$b_x[k] * dat$x
bs <- exp(dr$b_bs_Intercept[k] + dr$b_bs_x[k] * dat$x)
e1 <- ddm_lpdf_both(dat$rt - dr$ndt[k], mu, bs, dr$bias[k], dat$d)
r1 <- brms::dwiener(dat$rt, alpha = bs, tau = dr$ndt[k], beta = dr$bias[k],
                    delta = mu, resp = dat$d, log = TRUE)
cat("draw", k, ": eam vs brms log_lik", max(abs(e1 - ll1)),
    "; eam vs brms::dwiener", max(abs(e1 - r1)),
    "; brms::dwiener vs brms log_lik", max(abs(r1 - ll1)), "\n")
pe <- brms::posterior_epred(bfit, draw_ids = k, dpar = "mu")
cat("brms mu from prepare_predictions vs hand mu", max(abs(pe - mu)), "\n")
pb <- brms::posterior_epred(bfit, draw_ids = k, dpar = "bs")
cat("brms bs vs hand bs", max(abs(pb - bs)), "\n")
cat(".draw of row k:", dr$.draw[k], "\n")
print(head(cbind(ll1 = as.numeric(ll1), e1, r1)))
