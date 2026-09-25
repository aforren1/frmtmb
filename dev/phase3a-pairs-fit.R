# Item 3.2: a 4-channel pairwise frame fitted with (1 | pair), against
# four two-channel fits. The construction is test-cross-pairs.R's last
# block, seed 45: 8192 samples, loadings 1, 0.8, 0.5, 0.25 on one white
# source plus unit white noise, segments = 16, pairs Fz-Cz, Fz-Pz,
# Cz-Pz and Pz-Oz. The true coherence of a pair (a, b) is
# a^2 b^2 / ((a^2 + 1) (b^2 + 1)), flat in frequency.
Sys.setenv(PHASE3A_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.coupling)
})
phase3a_where("frmtmb.coupling")

xp_rec <- function(n = 4096L, load = c(1, 0.8, 0.5, 0.25)) {
  src <- stats::rnorm(n)
  X <- vapply(load, function(l) l * src + stats::rnorm(n), numeric(n))
  colnames(X) <- c("Fz", "Cz", "Pz", "Oz")
  X
}
load <- c(Fz = 1, Cz = 0.8, Pz = 0.5, Oz = 0.25)
if (nzchar(Sys.getenv("PHASE3A_LOAD"))) {
  load[] <- as.numeric(strsplit(Sys.getenv("PHASE3A_LOAD"), ",")[[1L]])
}
set.seed(45)
X <- xp_rec(8192L, load = unname(load))
pr <- rbind(c("Fz", "Cz"), c("Fz", "Pz"), c("Cz", "Pz"), c("Pz", "Oz"))
lab <- paste(pr[, 1L], pr[, 2L], sep = "-")
truth <- apply(pr, 1L, function(p) {
  a <- load[[p[1L]]]^2; b <- load[[p[2L]]]^2
  a * b / ((a + 1) * (b + 1))
})
xp <- frm_cross_pairs(X, pairs = pr, segments = 16L)
cat("rows", nrow(xp), "per pair", nrow(xp) / 4, "n", unique(xp$n), "\n")
form_sep <- bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
               pow2 ~ 1, coh ~ 1, phase ~ 1)
sep <- do.call(rbind, lapply(seq_len(nrow(pr)), function(k) {
  d <- frm_cross_spectrum(X[, pr[k, 1L]], X[, pr[k, 2L]], segments = 16L)
  fit <- frm(form_sep, family = cross_wishart(), data = d)
  cbind(frm_coherence(fit)[1L, ], conv = fit$opt$convergence,
        maxgrad = diagnose(fit, quiet = TRUE)$max_grad)
}))
nd <- data.frame(pair = factor(lab, levels = lab))
w1 <- character(0)
by_pair <- withCallingHandlers(frm(
  bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 0 + pair,
     pow2 ~ 0 + pair, coh ~ 0 + pair, phase ~ 0 + pair),
  family = cross_wishart(), data = xp),
  warning = function(w) { w1 <<- c(w1, conditionMessage(w))
                          invokeRestart("muffleWarning") })
w2 <- character(0)
shrunk <- withCallingHandlers(frm(
  bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 0 + pair,
     pow2 ~ 0 + pair, coh ~ 1 + (1 | pair), phase ~ 0 + pair),
  family = cross_wishart(), data = xp),
  warning = function(w) { w2 <<- c(w2, conditionMessage(w))
                          invokeRestart("muffleWarning") })
cb <- frm_coherence(by_pair, newdata = nd)
cs <- frm_coherence(shrunk, newdata = nd)
out <- data.frame(pair = lab, truth = truth, sep = sep$.estimate,
                  sep_se = sep$.se, by_pair = cb$.estimate,
                  by_pair_se = cb$.se, shrunk = cs$.estimate,
                  shrunk_se = cs$.se)
print(out, digits = 6, row.names = FALSE)
cat(sprintf("separate fits: convergence %s, max |gradient| %s\n",
            paste(sep$conv, collapse = " "),
            paste(sprintf("%.1e", sep$maxgrad), collapse = " ")))
for (nm in c("by_pair", "shrunk")) {
  f <- get(nm)
  dg <- diagnose(f, quiet = TRUE)
  cat(sprintf("%s: convergence %d, max |gradient| %.2e, pdHess %s, logLik %.6f\n",
              nm, f$opt$convergence, dg$max_grad, isTRUE(dg$pdHess),
              as.numeric(logLik(f))))
}
cat("by_pair warnings:", length(w1), if (length(w1)) substr(w1[1], 1, 80), "\n")
cat("shrunk warnings:", length(w2), if (length(w2)) substr(w2[1], 1, 80), "\n")
cat(sprintf("by_pair vs separate on the logit scale: max |diff| / se %.3e, max |se ratio - 1| %.3e\n",
            max(abs(qlogis(cb$.estimate) - qlogis(sep$.estimate)) / sep$.se),
            max(abs(cb$.se / sep$.se - 1))))
cat(sprintf("sum of separate logLik vs by_pair logLik: see above; shrunk vs separate on the logit scale: max |diff| / se %.3f\n",
            max(abs(qlogis(cs$.estimate) - qlogis(sep$.estimate)) / sep$.se)))
print(VarCorr(shrunk))
