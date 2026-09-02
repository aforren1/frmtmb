# Stage 1 verification: gaussian emissions, constant transitions.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-hmm",
                                   quiet = TRUE))
source("C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/hmm-common.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/"
a1 <- readRDS(paste0(D, "probeA1.rds"))
b1 <- readRDS(paste0(D, "probeB1.rds"))
flat <- function(f) {
  e <- unlist(fixef(f)); names(e) <- sub("\\.\\(Intercept\\)$", "", names(e))
  e
}

## --- A1: one sequence of 200, free initial distribution -------------
da <- data.frame(y = a1$y, t = seq_along(a1$y), g = 1L)
fa <- frm(bf(y ~ 1),
          family = hmm(K = 2, gaussian(), time = t, group = g,
                       init = "estimated"),
          data = da)
cat("A1 frm logLik      :", sprintf("%.9f", as.numeric(logLik(fa))), "\n")
cat("A1 probe reference :", sprintf("%.9f", a1$ll), "\n")
cat("A1 diff            :", abs(as.numeric(logLik(fa)) - a1$ll), "\n")
cat("A1 par diff vs probe:",
    max(abs(unname(flat(fa))[c(1, 3, 2, 4, 5, 6)] - a1$par[1:6])), "\n\n")

## --- B1 data with CONSTANT transitions ------------------------------
dd <- b1$dat
fb <- frm(bf(y ~ 1),
          family = hmm(K = 2, gaussian(), time = t, group = g,
                       init = "estimated"),
          data = dd)
cat("B1const frm logLik:", sprintf("%.9f", as.numeric(logLik(fb))), "\n")
est <- flat(fb)
mu <- c(est[["mu1"]], est[["mu2"]])
sg <- exp(c(est[["sigma1"]], est[["sigma2"]]))
lg <- c(est[["tr12"]], est[["tr22"]])
G <- rbind(c(1, exp(lg[1])) / (1 + exp(lg[1])),
           c(1, exp(lg[2])) / (1 + exp(lg[2])))
dl <- softmax0(fb$estimates[["hmm_ldel"]])
lpm <- lpmat_gauss(dd$y, mu, sg)
rows <- hmm_seq_index(dd$g, dd$t)
llnum <- sum(vapply(rows, function(r) fwd_num(lpm[r, , drop = FALSE], G, dl),
                    numeric(1)))
cat("B1const numeric forward at the estimates:", sprintf("%.9f", llnum),
    "\n  diff:", abs(llnum - as.numeric(logLik(fb))), "\n")

if (requireNamespace("depmixS4", quietly = TRUE)) {
  set.seed(9)
  dm <- depmixS4::depmix(y ~ 1, data = dd, nstates = 2,
                         ntimes = as.integer(table(dd$g)))
  best <- -Inf; bp <- NULL
  for (i in 1:6) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE,
      emcontrol = depmixS4::em.control(random.start = TRUE,
                                       tol = 1e-12, maxit = 5000))),
      silent = TRUE)
    if (!inherits(ff, "try-error")) {
      v <- as.numeric(depmixS4::logLik(ff))
      if (v > best) { best <- v; bp <- depmixS4::getpars(ff) }
    }
  }
  cat("depmixS4 logLik:", sprintf("%.9f", best),
      " diff:", abs(best - as.numeric(logLik(fb))), "\n")
  # depmixS4 pars: delta(2), transition rows (2 x 2 probs), then
  # per-state (mu, sd)
  cat("depmixS4 pars:", paste(round(bp, 6), collapse = " "), "\n")
  cat("frm  tpm:", paste(round(as.vector(t(G)), 6), collapse = " "), "\n")
  cat("frm  mu/sd:", paste(round(c(mu[1], sg[1], mu[2], sg[2]), 6),
                           collapse = " "), "\n")
}

## --- T = 1 sequences collapse to mixture() exactly ------------------
set.seed(21)
n1 <- 400
z <- sample(1:2, n1, TRUE, prob = c(0.6, 0.4))
d1 <- data.frame(y = rnorm(n1, c(0, 3)[z], c(0.9, 0.6)[z]),
                 g = seq_len(n1), t = 1L)
fh <- try(frm(bf(y ~ 1),
              family = hmm(K = 2, gaussian(), time = t, group = g,
                           init = "estimated"), data = d1), silent = TRUE)
cat("\nall-singleton refusal:\n", conditionMessage(attr(fh, "condition")),
    "\n")
# Held at constants the transitions are identified, and the T = 1 HMM is
# then EXACTLY a two-component mixture (probe F2).
fh2 <- frm(bf(y ~ 1, tr12 = 0, tr22 = 0),
           family = hmm(K = 2, gaussian(), time = t, group = g,
                        init = "estimated"), data = d1)
fm2 <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian()), data = d1)
cat("hmm(T = 1) logLik :", sprintf("%.12f", as.numeric(logLik(fh2))), "\n")
cat("mixture()  logLik :", sprintf("%.12f", as.numeric(logLik(fm2))), "\n")
cat("diff              :",
    abs(as.numeric(logLik(fh2)) - as.numeric(logLik(fm2))), "\n")
eh <- flat(fh2); em <- flat(fm2)
cat("hmm  mu/sigma:", paste(round(unname(eh[c("mu1", "mu2", "sigma1",
                                              "sigma2")]), 8),
                            collapse = " "), "\n")
cat("mix  mu/sigma:", paste(round(unname(em[c("mu1", "mu2", "sigma1",
                                              "sigma2")]), 8),
                            collapse = " "), "\n")
cat("hmm weight1  :",
    round(1 / (1 + exp(fh2$estimates[["hmm_ldel"]])), 8), "\n")
cat("mix weight1  :", round(plogis(em[["theta1"]]), 8), "\n")
cat("df hmm/mix   :", attr(logLik(fh2), "df"), attr(logLik(fm2), "df"), "\n")
