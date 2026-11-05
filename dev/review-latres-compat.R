# Attack 9: the compat rows the lane moved. Does each note state a true
# condition? Every aterm on a latent fit, against brms, sampled.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
suppressPackageStartupMessages({library(brms); library(rstan)})
set.seed(20261102)
d <- expand.grid(t = 1:5, g = factor(1:12))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(12, arima.sim(list(ar = .5), 5, sd = .5)))
d$y <- rpois(n, exp(.5 + .3 * d$x + e))
d$w <- rep(c(1, 2, 1), length.out = n)
d$cn <- sample(c("none", "right"), n, TRUE, c(.8, .2))
d$lb <- 0
d$ub <- d$y + 20
d$tr <- 6L
d$bi <- rbinom(n, 6, plogis(.3 * d$x + e))
d$bin <- rbinom(n, 1, plogis(.3 * d$x + e))
d$ex <- rep(c(1, 2), length.out = n)
d$sev <- rnorm(n, 0, 1)
d$se <- 0.2
d$pos <- rgamma(n, 3, 3 / exp(.3 * d$x + e))
d$posna <- d$pos
d$posna[sample(n, 6)] <- NA
d$ord <- factor(findInterval(.3 * d$x + e + rlogis(n), c(-1, 0, 1)) + 1L, ordered = TRUE)
d$th <- 5L
d$sub <- d$t != 3
d$ci <- pmax(d$y - 1, 0)
d$dec <- rbinom(n, 1, .5)
d$yy <- d$y
frm_try <- function(form_q, fam_q, ...) {
  ns <- asNamespace("frmtmb")
  msg_of(frm(eval(form_q, list(), ns), data = d, family = eval(fam_q, list(), ns), ...))
}
br_try <- function(form_q, fam_q) {
  ns <- asNamespace("brms")
  r <- msg_of(brms::stancode(eval(form_q, list(), ns), d, family = eval(fam_q, list(), ns)))
  if (!is.null(r$err)) return(c("REFUSES", substr(gsub("\\s+", " ", r$err), 1, 160)))
  chk <- msg_of(rstan::stanc(model_code = r$val, verbose = FALSE))
  if (!is.null(chk$err)) return(c("STANC-FAIL", substr(chk$err, 1, 120)))
  c("BUILDS", "")
}
cases <- list(
  list("weights poisson", quote(y | weights(w) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("weights cosy", quote(y | weights(w) ~ x + cosy(t, g)), quote(poisson())),
  list("cens poisson", quote(y | cens(cn) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("cens int poisson", quote(y | cens(cn, y2 = ci) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("trunc poisson", quote(y | trunc(lb = 0, ub = 100) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("trunc var", quote(y | trunc(ub = ub) ~ x + ar(t, g)), quote(poisson())),
  list("rate poisson", quote(y | rate(ex) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("trials binomial", quote(bi | trials(tr) ~ x + ar(t, g, cov = TRUE)), quote(binomial())),
  list("thres cumulative", quote(ord | thres(3) ~ x + ar(t, g, cov = TRUE)), quote(cumulative())),
  list("se lognormal", quote(pos | se(se) ~ x + ar(t, g, cov = TRUE)), quote(lognormal())),
  list("se poisson", quote(y | se(se) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("mi lognormal", quote(posna | mi() ~ x + ar(t, g, cov = TRUE)), quote(lognormal())),
  list("mi Gamma", quote(posna | mi() ~ x + ar(t, g, cov = TRUE)), quote(Gamma(link = "log"))),
  list("subset univariate", quote(y | subset(sub) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("vreal", quote(y | vreal(ex) ~ x + ar(t, g, cov = TRUE)), quote(poisson())),
  list("vint", quote(y | vint(th) ~ x + ar(t, g, cov = TRUE)), quote(poisson()))
)
for (cc in cases) {
  fr <- frm_try(cc[[2]], cc[[3]])
  br <- br_try(cc[[2]], cc[[3]])
  cat(sprintf("%-20s brms %-10s frm %s", cc[[1]], br[1], if (!is.null(fr$err)) "REFUSES" else "FITS"))
  if (is.null(fr$err)) cat(sprintf("  ll %.4f code %d latres %d autocor %d", as.numeric(logLik(fr$val)), fr$val$opt$convergence, length(fi("frame_latres")(fr$val$frame)), length(fr$val$frame$autocor)))
  cat("\n")
  if (nzchar(br[2])) cat("    brms:", br[2], "\n")
  if (!is.null(fr$err)) cat("    frm: ", substr(fr$err, 1, 220), "\n")
  if (length(fr$warn)) cat("    warn:", substr(fr$warn[1], 1, 150), "\n")
}
cat("\n== the compat table's rows\n")
ct <- frm_compat("ar()")
ct <- as.data.frame(ct)
print(table(ct$kind_b, ct$status))
print(ct[ct$kind_b == "aterm", c("feature_b", "status")], row.names = FALSE)
print(unique(substr(ct$note[ct$kind_b == "aterm"], 1, 120)))
print(ct[ct$kind_b == "family" & ct$status == "refused", c("feature_b", "status")], row.names = FALSE)
print(ct[ct$feature_b %in% c("fitted", "predict", "simulate", "residuals", "residuals_osa", "nl", "quadrature", "importance", "rescor", "loo", "conditional_effects", "emmeans"), c("feature_b", "status")], row.names = FALSE)
for (fb in c("predict", "simulate", "residuals", "dharma", "residuals_osa")) {
  r <- ct[ct$feature_b == fb, "note"]
  if (length(r)) cat(fb, ":", substr(r, 1, 340), "\n")
}
cat("\nkind:family rows on families with a joint density (structure loglik) or no mu:\n")
fams <- unique(ct$feature_b[ct$kind_b == "family"])
for (fn in fams) {
  e <- tryCatch(frmtmb:::fam_structure(frmtmb:::frm_family(fn)), error = function(z) NULL)
}
print(length(fams))
print(sort(fams))
cat("\nfcor rows\n")
ct2 <- as.data.frame(frm_compat("fcor()"))
print(table(ct2$kind_b, ct2$status))
print(head(ct2[ct2$kind_b == "family" & ct2$status != "refused", c("feature_b", "status")], 20), row.names = FALSE)
print(as.data.frame(frm_compat("cosy()", "weights()")), row.names = FALSE)
