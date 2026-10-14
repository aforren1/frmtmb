# Reviewer probe: can frmtmb's exported ordinal family objects answer
# brms's log_lik_<ordinal>(i, prep) and posterior_predict_<ordinal>(i,
# prep) with a thin adapter (no internal function), as the survey says
# they cannot? Data are brms 2.23.0's tests.log_lik.R:446 block.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(446)
ns <- 50; nobs <- 8; nthres <- 3; ncat <- nthres + 1
mu <- array(rnorm(ns * nobs), dim = c(ns, nobs))
disc <- rexp(ns)
# brms's block uses all-zero thresholds; use increasing ones as well so
# the comparison against a closed form is informative
thres <- t(replicate(ns, sort(rnorm(nthres))))
Y <- rep(1:ncat, 2)
ad_ll <- function(fun, i, link = "logit") {
  fam <- frmtmb::brmsfamily(fun, link = link)
  lp <- fam$ord_lpdf_make(identity)
  vapply(seq_len(ns), function(s) {
    as.numeric(lp(Y[i], list(mu = mu[s, i], disc = disc[s]), list(),
                  list(tau_raw = thres[s, ])))
  }, numeric(1))
}
for (f in c("cumulative", "sratio", "cratio", "acat")) {
  ll <- tryCatch(sapply(1:nobs, ad_ll, fun = f),
                 error = function(e) conditionMessage(e))
  cat(f, ": ", if (is.character(ll)) ll else
        paste("dim", paste(dim(ll), collapse = " x "),
              "finite", all(is.finite(ll))), "\n", sep = "")
}
# cumulative logit closed form: P(Y = k) = F(disc (t_k - mu)) - F(...)
cf <- sapply(1:nobs, function(i) {
  vapply(seq_len(ns), function(s) {
    tt <- c(-Inf, thres[s, ], Inf)
    k <- Y[i]
    log(plogis(disc[s] * (tt[k + 1] - mu[s, i])) -
          plogis(disc[s] * (tt[k] - mu[s, i])))
  }, numeric(1))
})
ll <- sapply(1:nobs, ad_ll, fun = "cumulative")
cat("cumulative against the closed form: max abs diff",
    format(max(abs(ll - cf)), digits = 3), "\n")
# brms's own value on the same prep
prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
prep$dpars <- list(mu = mu, disc = disc)
prep$thres$thres <- thres
prep$data <- list(Y = Y, ncat = ncat)
prep$family <- brms::cumulative("logit")
prep$family$link <- "logit"
bl <- sapply(1:nobs, brms:::log_lik_cumulative, prep = prep)
cat("cumulative against brms:::log_lik_cumulative: max abs diff",
    format(max(abs(ll - bl)), digits = 3), "\n")
fam <- frmtmb::brmsfamily("cumulative")
cat("ord_sim_make present:", is.function(fam$ord_sim_make), "\n")
print(args(fam$sim))
# posterior_predict side: one draw per posterior draw at observation i
for (f in c("cumulative", "sratio", "cratio", "acat")) {
  fm <- frmtmb::brmsfamily(f)
  sm <- fm$ord_sim_make(identity)
  pr <- tryCatch(sapply(1:nobs, function(i) vapply(seq_len(ns), function(s)
    as.numeric(sm(list(mu = mu[s, i], disc = disc[s]), list(), 1L,
                  list(tau_raw = thres[s, ]))), numeric(1))),
    error = function(e) conditionMessage(e))
  cat("sim", f, ": ", if (is.character(pr)) pr else
        paste("dim", paste(dim(pr), collapse = " x "), "range",
              paste(range(pr), collapse = "-")), "\n", sep = "")
}
