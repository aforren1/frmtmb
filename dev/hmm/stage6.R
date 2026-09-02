# Stage 6: NA masking (probe F4), the remaining guards, frm_sample.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-hmm",
                                   quiet = TRUE))
source("C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/hmm-common.R")

## ---- probe F4's data, verbatim ---------------------------------------
K <- 2L
set.seed(64)
N4 <- 30L; T4 <- 20L; n4 <- N4 * T4
G4 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2, 2, byrow = TRUE)
mu4 <- c(0, 3); sg4 <- c(0.6, 0.6)
d4 <- do.call(rbind, lapply(seq_len(N4), function(g) {
  s <- integer(T4); s[1] <- sample.int(K, 1, prob = c(0.5, 0.5))
  for (t in seq_len(T4 - 1L)) s[t + 1L] <- sample.int(K, 1,
                                                      prob = G4[s[t], ])
  data.frame(g = g, t = seq_len(T4), y = rnorm(T4, mu4[s], sg4[s]))
}))
miss <- which(d4$t %in% c(7, 8, 14))
d4$miss <- as.integer(seq_len(n4) %in% miss)
d4$y_na <- d4$y; d4$y_na[miss] <- NA
d4$y_sent <- d4$y; d4$y_sent[miss] <- 0

fam <- hmm(K = 2, gaussian(), time = t, group = g, init = "estimated")
fm <- frm(bf(y_na ~ 1), family = fam, data = d4)
cat("masked-route logLik:", sprintf("%.9f", as.numeric(logLik(fm))), "\n")
cat("probe F4 reference : -691.400711096\n")

## the independent numeric reference at frm's own estimates
e <- unlist(fixef(fm))
lp_masked <- function(y, mu, sg, ms) {
  cbind(ifelse(ms == 1, 0, dnorm(y, mu[1], sg[1], log = TRUE)),
        ifelse(ms == 1, 0, dnorm(y, mu[2], sg[2], log = TRUE)))
}
G <- tpm_tv_num(rbind(c(0, e[["tr12.(Intercept)"]]),
                      c(0, e[["tr22.(Intercept)"]])), matrix(0, K, K), 0, K)
lpm <- lp_masked(d4$y_sent, c(e[["mu1.(Intercept)"]],
                              e[["mu2.(Intercept)"]]),
                 exp(c(e[["sigma1.(Intercept)"]],
                       e[["sigma2.(Intercept)"]])), d4$miss)
rbg <- split(seq_len(n4), d4$g)
ll <- sum(vapply(rbg, function(r)
  fwd_num_tv(lpm, function(rr) G, r,
             softmax0(fm$estimates[["hmm_ldel"]])), numeric(1)))
cat("numeric reference  :", sprintf("%.9f", ll),
    " diff", abs(ll - as.numeric(logLik(fm))), "\n")
cat("tpm masked:", round(as.vector(t(G)), 5), " (true .9 .1 .25 .75)\n")
cat("n used    :", nobs(fm), " rows in data:", n4, "\n")

## the na.omit route on the same data, for contrast: the chain SHORTENS
d4b <- d4[!is.na(d4$y_na), ]
fo <- frm(bf(y_na ~ 1), family = fam, data = d4b)
eo <- unlist(fixef(fo))
Go <- tpm_tv_num(rbind(c(0, eo[["tr12.(Intercept)"]]),
                       c(0, eo[["tr22.(Intercept)"]])),
                 matrix(0, K, K), 0, K)
cat("tpm dropped rows:", round(as.vector(t(Go)), 5),
    " (probe F4: .87158 .12842 .249 .751)\n")
cat("fitted()/residuals() at masked rows are NA:",
    all(is.na(residuals(fm)[d4$miss == 1])), "\n")
cat("hmm_probs() still defined there:",
    all(is.finite(hmm_probs(fm)[d4$miss == 1, ])), "\n\n")

## ---- the remaining refusals -------------------------------------------
d4$w <- 1; d4$cc <- 0; d4$sd <- 1
tries <- list(
  weights = quote(frm(bf(y | weights(w) ~ 1), family = fam, data = d4)),
  cens = quote(frm(bf(y | cens(cc) ~ 1), family = fam, data = d4)),
  trunc = quote(frm(bf(y | trunc(lb = -9) ~ 1), family = fam, data = d4)),
  se = quote(frm(bf(y | se(sd) ~ 1), family = fam, data = d4)),
  mi = quote(frm(bf(y_na | mi() ~ 1), family = fam, data = d4)),
  mvbf = quote(frm(mvbf(bf(y ~ 1) + hmm(2, gaussian(), time = t,
                                        group = g),
                        bf(y_sent ~ 1) + gaussian()), data = d4)),
  profile = quote(frm(bf(y ~ 1), family = fam, data = d4,
                      control = frmtmb_control(profile = TRUE))),
  dup_time = quote(frm(bf(y ~ 1),
                       family = hmm(2, gaussian(), time = miss, group = g),
                       data = d4))
)
for (nm in names(tries)) {
  r <- try(suppressMessages(eval(tries[[nm]])), silent = TRUE)
  cat(sprintf("%-9s -> %s\n", nm,
              if (inherits(r, "try-error"))
                substr(conditionMessage(attr(r, "condition")), 1, 110)
              else "NO REFUSAL"))
}

## ---- symmetric-start warning -------------------------------------------
f0 <- frm(bf(y ~ 1), family = fam, data = d4)
st <- f0$estimates
st$beta[] <- 0
w <- tryCatch(frm(bf(y ~ 1), family = fam, data = d4,
                  start = list(beta = st$beta)),
              warning = function(e) conditionMessage(e))
cat("\nsymmetric start ->", substr(w, 1, 140), "\n")

## ---- frm_sample --------------------------------------------------------
if (requireNamespace("tmbstan", quietly = TRUE)) {
  s <- try(suppressWarnings(suppressMessages(
    frm_sample(f0, chains = 1, iter = 200, refresh = 0))), silent = TRUE)
  cat("\nfrm_sample:",
      if (inherits(s, "try-error"))
        substr(conditionMessage(attr(s, "condition")), 1, 150)
      else paste("OK, class", paste(class(s), collapse = "/")), "\n")
} else {
  cat("\ntmbstan not installed\n")
}

## ---- the family= and + spellings agree ---------------------------------
fa <- frm(bf(y ~ 1), family = hmm(2, gaussian(), time = t, group = g),
          data = d4)
fp <- frm(bf(y ~ 1) + hmm(2, gaussian(), time = t, group = g), data = d4)
cat("\nfamily= vs + : logLik diff",
    abs(as.numeric(logLik(fa)) - as.numeric(logLik(fp))),
    " coef diff", max(abs(unlist(fixef(fa)) - unlist(fixef(fp)))), "\n")
