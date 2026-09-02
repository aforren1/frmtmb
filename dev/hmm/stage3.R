# Stage 3: emission families other than gaussian.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-hmm",
                                   quiet = TRUE))
source("C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/hmm-common.R")

sim_paths <- function(N, Tl, G, K, seed) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tl); s[1] <- sample.int(K, 1)
    for (t in 2:Tl) s[t] <- sample.int(K, 1, prob = G[s[t - 1], ])
    data.frame(g = g, t = seq_len(Tl), state = s)
  }))
}
G2 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2, 2, byrow = TRUE)

## --- poisson emissions -----------------------------------------------
dp <- sim_paths(30, 25, G2, 2L, 101)
dp$y <- rpois(nrow(dp), c(1.5, 9)[dp$state])
fp <- frm(bf(y ~ 1),
          family = hmm(K = 2, poisson(), time = t, group = g,
                       init = "estimated"), data = dp)
cat("poisson logLik:", sprintf("%.9f", as.numeric(logLik(fp))),
    " df", attr(logLik(fp), "df"), "\n")
e <- unlist(fixef(fp))
lam <- exp(c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]]))
lg <- c(e[["tr12.(Intercept)"]], e[["tr22.(Intercept)"]])
G <- rbind(c(1, exp(lg[1])) / (1 + exp(lg[1])),
           c(1, exp(lg[2])) / (1 + exp(lg[2])))
lpm <- vapply(lam, function(l) dpois(dp$y, l, log = TRUE),
              numeric(nrow(dp)))
dl <- softmax0(fp$estimates[["hmm_ldel"]])
rows <- hmm_seq_index(dp$g, dp$t)
ll <- sum(vapply(rows, function(r) fwd_num(lpm[r, , drop = FALSE], G, dl),
                 numeric(1)))
cat("poisson numeric forward:", sprintf("%.9f", ll),
    " diff", abs(ll - as.numeric(logLik(fp))), "\n")
cat("lambda:", round(lam, 5), " (true 1.5, 9)\n")
cat("decoding accuracy:", mean(hmm_viterbi(fp) == dp$state), "\n\n")

if (requireNamespace("depmixS4", quietly = TRUE)) {
  set.seed(3)
  dm <- depmixS4::depmix(y ~ 1, data = dp, nstates = 2,
                         family = poisson(),
                         ntimes = as.integer(table(dp$g)))
  best <- -Inf
  for (i in 1:6) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE, emcontrol = depmixS4::em.control(
        random.start = TRUE, tol = 1e-12, maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      best <- max(best, as.numeric(depmixS4::logLik(ff)))
    }
  }
  cat("poisson depmixS4:", sprintf("%.9f", best),
      " diff", abs(best - as.numeric(logLik(fp))), "\n\n")
}

## --- categorical (multinomial) emissions -----------------------------
dc <- sim_paths(40, 15, G2, 2L, 202)
P1 <- c(0.55, 0.25, 0.12, 0.08)
P2 <- c(0.08, 0.15, 0.32, 0.45)
cat_of <- vapply(dc$state, function(s) {
  sample.int(4, 1, prob = if (s == 1) P1 else P2)
}, integer(1))
Y <- matrix(0L, nrow(dc), 4)
Y[cbind(seq_len(nrow(dc)), cat_of)] <- 1L
dc$Y <- Y
dc$cf <- factor(cat_of)
fc <- frm(bf(Y ~ 1),
          family = hmm(K = 2, multinomial(K = 4), time = t, group = g,
                       init = "estimated"), data = dc)
cat("categorical logLik:", sprintf("%.9f", as.numeric(logLik(fc))),
    " df", attr(logLik(fc), "df"), "\n")
ec <- unlist(fixef(fc))
prob_state <- function(k) {
  eta <- c(0, vapply(2:4, function(j)
    ec[[paste0("mu", j, k, ".(Intercept)")]], numeric(1)))
  exp(eta) / sum(exp(eta))
}
cat("state 1 probs:", round(prob_state(1), 4), "\n")
cat("state 2 probs:", round(prob_state(2), 4), "\n")
cat("true         :", round(P1, 4), "|", round(P2, 4), "\n")
cat("decoding accuracy:", mean(hmm_viterbi(fc) == dc$state), "\n")

if (requireNamespace("depmixS4", quietly = TRUE)) {
  set.seed(5)
  dm <- depmixS4::depmix(cf ~ 1, data = dc, nstates = 2,
                         family = depmixS4::multinomial("identity"),
                         ntimes = as.integer(table(dc$g)))
  best <- -Inf; bp <- NULL
  for (i in 1:8) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE, emcontrol = depmixS4::em.control(
        random.start = TRUE, tol = 1e-12, maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      v <- as.numeric(depmixS4::logLik(ff))
      if (v > best) { best <- v; bp <- depmixS4::getpars(ff) }
    }
  }
  cat("categorical depmixS4:", sprintf("%.9f", best),
      " diff", abs(best - as.numeric(logLik(fc))), "\n")
  cat("depmixS4 pars:", paste(round(bp, 5), collapse = " "), "\n")
}

## --- refusals ---------------------------------------------------------
for (ex in list(quote(hmm(2, cumulative())),
                quote(hmm(2, mixture(gaussian(), gaussian()))),
                quote(hmm(1, gaussian())),
                quote(hmm(12, gaussian())))) {
  r <- try(eval(ex), silent = TRUE)
  cat("\n", deparse(ex), "->\n  ",
      conditionMessage(attr(r, "condition")), "\n")
}
