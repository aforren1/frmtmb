# Would an MCSE that takes the antithetic pair as the unit, or the larger of the
# two, cover where the iid one reads low? Full fits over importance_seed, with
# both numbers read from the fit's own log weights.
# Usage: Rscript dev/review-impsamp-pairmcse.R <design> <dseed> <draws> [S]
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs()
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-designs.R")
a <- commandArgs(TRUE)
nm <- a[[1]]
dseed <- as.integer(a[[2]])
N <- as.integer(a[[3]])
S <- if (length(a) >= 4) as.integer(a[[4]]) else 100L
de <- design_list[[nm]](dseed)
ns <- asNamespace("frmtmb")
dl <- split(de$data, de$data$g)
gh <- gh_rule(81L)
exact_nll <- function(par) {
  sd <- de$sdof(par)
  tot <- 0
  for (d in dl) {
    fn <- function(u) de$ll(par, d, u) + dnorm(u, 0, sd, log = TRUE)
    o <- optimize(function(u) -fn(u), c(-8 * sd - 15, 8 * sd + 15), tol = 1e-7)
    m <- o$minimum
    h <- 1e-3
    cur <- -(fn(m + h) - 2 * fn(m) + fn(m - h)) / h^2
    s <- if (is.finite(cur) && cur > 1e-8) 1 / sqrt(cur) else sd
    nodes <- m + sqrt(2) * s * gh$x
    v <- fn(nodes) + gh$x^2 + log(gh$w)
    tot <- tot - (max(v) + log(sum(exp(v - max(v)))) + log(sqrt(2) * s))
  }
  tot
}
rows <- list()
for (s in seq_len(S)) {
  fit <- suppressWarnings(frm(de$formula, family = de$family, data = de$data,
                              importance = N,
                              control = frmtmb_control(importance_seed = s)))
  prop <- ns$imp_frozen_proposal(fit)
  am <- as.matrix(prop$io$amat(ns$imp_par_list(prop$template, fit$opt$par, prop$map)))
  nd <- ncol(am)
  h <- nd %/% 2L
  # log m_g has variance var(pair mean) / (n pairs) / mean^2
  vp <- 0
  for (g in seq_len(nrow(am))) {
    w <- exp(am[g, ] - max(am[g, ]))
    pm <- (w[seq_len(h)] + w[h + seq_len(h)]) / 2
    vp <- vp + var(pm) / h / mean(w)^2
  }
  rows[[s]] <- data.frame(seed = s, err = -fit$opt$objective + exact_nll(unname(fit$opt$par)),
                          iid = fit$importance$mcse, pair = sqrt(vp),
                          ess = fit$importance$ess_min)
}
r <- do.call(rbind, rows)
mx <- pmax(r$iid, r$pair)
cat(sprintf("%s dseed %d N %d S %d: sd_err %.4g med iid %.4g med pair %.4g | cover95 iid %.2f pair %.2f max %.2f | sd/iid %.2f sd/pair %.2f | ess_min %.3f\n",
            nm, dseed, N, S, sd(r$err), median(r$iid), median(r$pair),
            mean(abs(r$err) <= 1.96 * r$iid), mean(abs(r$err) <= 1.96 * r$pair),
            mean(abs(r$err) <= 1.96 * mx), sd(r$err) / median(r$iid),
            sd(r$err) / median(r$pair), min(r$ess)))
