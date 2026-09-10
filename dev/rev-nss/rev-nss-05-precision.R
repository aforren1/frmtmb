# REVIEW of lane nss, attacks 5, 6 and 10.
#
# 5. "Truncating manufactures a precision the data do not support":
#    standard errors 1.0x to 41.6x smaller at the 264.6 h design. The
#    lane reports this for `k21` only. This script reports it for ALL
#    FIVE structural parameters, because "does it generalize beyond
#    k21" is the question the brief asks.
# 6. The record correction. `dev/reviews/2026-09-09-lincmt.md` says the
#    k21 factor of 3.16 was measured at "107 h, ii = 24"; the lane says
#    the script it names runs ke 0.1, k12 0.2, k21 0.008, which is
#    264.6 h. Checked HERE BY CONSTRUCTION, by running the same design
#    and looking for -5.54547 and -4.39225, and by running the 107 h
#    design and looking for something that is not them.
# 10. Whether the acceptance criterion discriminates.
#
# Both designs, six seeds each, `frm_lincmt()` as the vehicle exactly
# as the lane used it: n_ss = Inf is the exact limit the shipped
# default now reaches, n_ss = 20 writes out the cycles the old default
# truncated at.
#
# Script path: dev/rev-nss/rev-nss-05-precision.R
# Seeds: 101 to 106.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()

designs <- list(
  d265 = list(tag = "264.6 h", KE = 0.1,  K12 = 0.2, K21 = 0.008,
              KA = 1.0, V = 10, II = 24),
  d107 = list(tag = "107.1 h", KE = 0.15, K12 = 0.3, K21 = 0.02,
              KA = 1.0, V = 10, II = 24))
NS <- 30L

sim <- function(g, seed, sd_obs = 0.15) {
  set.seed(seed)
  tt <- g$II * c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1)
  d <- data.frame(id = factor(rep(seq_len(NS), each = length(tt))),
                  time = rep(tt, NS))
  i <- as.integer(d$id)
  lke <- log(g$KE) + rnorm(NS, 0, 0.2)
  lka <- log(g$KA) + rnorm(NS, 0, 0.3)
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = g$II,
                   addl = 0L, ss = TRUE)
  mu <- numeric(nrow(d))
  for (j in seq_len(NS)) {
    k <- which(i == j)
    mu[k] <- frm_lincmt(parms = list(ke = exp(lke[[j]]), k12 = g$K12,
                                     k21 = g$K21, ka = exp(lka[[j]]),
                                     V = g$V),
                        times = d$time[k], ncmt = 2, depot = TRUE,
                        events = ev)
  }
  d$conc <- mu + rnorm(nrow(d), 0, sd_obs)
  list(d = d, ev = ev)
}

fit <- function(g, z, nss) {
  doses <- z$ev
  main <- if (is.finite(nss))
    conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                   k21 = exp(lk21), ka = exp(lka),
                                   V = exp(lV)),
                      times = time, group = id, ncmt = 2,
                      depot = TRUE, events = doses, n_ss = 20L)
  else
    conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                   k21 = exp(lk21), ka = exp(lka),
                                   V = exp(lV)),
                      times = time, group = id, ncmt = 2,
                      depot = TRUE, events = doses)
  bd <- bf(main, lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
           lk21 ~ 1, lV ~ 1, nl = TRUE)
  frm(bd + gaussian(), data = z$d, se = TRUE,
      start = list(beta = c(log(g$KE), log(g$KA), log(g$K12),
                            log(g$K21), log(g$V))))
}

# confint() carries the estimate and the interval for every fixed
# effect; the se is the half-width over 1.96, which is how the lane's
# nss-18 reads it.
pull <- function(f) {
  ci <- confint(f)
  rn <- rownames(ci)
  keep <- grep("^(lke|lka|lk12|lk21|lV)_", rn)
  nm <- sub("_.*$", "", rn[keep])
  list(est = setNames(as.numeric(ci[keep, "est"]), nm),
       se = setNames(as.numeric(ci[keep, "upr"] - ci[keep, "lwr"]) /
                       (2 * 1.96), nm))
}

sel <- Sys.getenv("REV_DESIGN", "")
for (dn in if (nzchar(sel)) sel else names(designs)) {
  g <- designs[[dn]]
  b <- g$KE + g$K12 + g$K21
  lz <- (b - sqrt(b * b - 4 * g$KE * g$K21)) / 2
  cat(sprintf("\n\n======== design %s: ke %.3f k12 %.3f k21 %.4f ka %.1f",
              g$tag, g$KE, g$K12, g$K21, g$KA))
  cat(sprintf("\n         lambda_z %.6f -> terminal half-life %.1f h,",
              lz, log(2) / lz))
  cat(sprintf(" ii %d, lambda_z*ii %.4f ========\n", g$II, lz * g$II))
  truth <- c(lke = log(g$KE), lka = log(g$KA), lk12 = log(g$K12),
             lk21 = log(g$K21), lV = log(g$V))
  ses <- list()
  for (seed in 101:106) {
    z <- sim(g, seed)
    a <- tryCatch(pull(fit(g, z, Inf)), error = function(e) e)
    b2 <- tryCatch(pull(fit(g, z, 20L)), error = function(e) e)
    if (inherits(a, "condition") || inherits(b2, "condition")) {
      cat("seed", seed, "fit error\n"); next
    }
    nm <- names(truth)
    cat(sprintf("\nseed %d\n", seed))
    cat(sprintf("  %-6s %9s | %10s %9s %6s | %10s %9s %6s | %7s\n",
                "par", "truth", "limit est", "limit se", "cov",
                "trunc est", "trunc se", "cov", "se ratio"))
    for (p in nm) {
      ea <- a$est[[p]]; sa <- a$se[[p]]
      eb <- b2$est[[p]]; sb <- b2$se[[p]]
      ca <- abs(ea - truth[[p]]) <= 1.96 * sa
      cb <- abs(eb - truth[[p]]) <= 1.96 * sb
      cat(sprintf(paste0("  %-6s %9.4f | %10.4f %9.4f %6s |",
                         " %10.4f %9.4f %6s | %7.2f\n"),
                  p, truth[[p]], ea, sa, ca, eb, sb, cb, sa / sb))
      ses[[length(ses) + 1L]] <- data.frame(seed = seed, par = p,
                                            se_lim = sa, se_tr = sb,
                                            ratio = sa / sb)
    }
  }
  s <- do.call(rbind, ses)
  cat("\n-- se ratio (limit / truncated) by parameter, six seeds --\n")
  print(round(tapply(s$ratio, s$par, function(x)
    c(min = min(x), median = stats::median(x), max = max(x)))[
      c("lke", "lka", "lk12", "lk21", "lV")], 3))
  cat("\nse ratios above 1.5, by parameter:\n")
  print(table(s$par[s$ratio > 1.5]))
}
cat("\ndone\n")
