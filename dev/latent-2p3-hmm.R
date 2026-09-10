# Lane `latent`, item 2.3 of dev/extension-gaps-plan.md.
#
# hmm() at the plan's realistic scale, 50 sequences x 500 steps, K = 3
# gaussian, scored against the simulator's truth and against a
# third-party fit of the SAME model on the SAME data.
#
#   arm A  fixed transitions, against depmixS4. depmixS4 has no random
#          effects, so this arm drops the random effect and keeps
#          everything else, and both sides estimate ONE initial
#          distribution shared across sequences.
#   arm B  `tr12 ~ (1 | id)`, against hmmTMB with the same random
#          effect. Two conventions have to be pinned or the comparison
#          is between different models, and probe
#          dev/latent-2p3-hmmtmb-probe4.R is where each was measured:
#            * hmmTMB's `initial_state = "estimated"` estimates a
#              SEPARATE initial distribution PER SEQUENCE. frmtmb's
#              `init = "estimated"` estimates one shared. Left alone
#              the two log-likelihoods sit 21.2 apart on 20 sequences
#              and nothing says why. Both sides here fix the initial
#              distribution at uniform, which is also what the
#              simulator draws.
#            * hmmTMB references the DIAGONAL of the transition matrix
#              and refuses any other reference once a formula matrix is
#              given; frmtmb references state 1. Row 1 is the same
#              under both, so `tr12` and `tr13` compare directly and
#              rows 2 and 3 are compared as a transition MATRIX.
#            * a column called `state` is read by hmmTMB as KNOWN
#              states, silently (probe D3 of dev/hmm-feasibility.md),
#              so it is dropped before hmmTMB sees the frame.
#
#   Rscript dev/latent-2p3-hmm.R <arm> <reps> [out.tsv] [ns] [tl]
#
# THE TWO SIDES ARE NOT STARTED THE SAME WAY, ON PURPOSE. frmtmb runs
# from its own shipped cold start, because that is the behaviour under
# test. depmixS4 gets four random EM starts and the best is kept;
# hmmTMB is started AT THE TRUTH. Both references are therefore closer
# to an oracle than a user would be, which is what makes the arm able to
# catch frmtmb converging to a local optimum. It also means a
# disagreement is evidence about frmtmb's starting values and not about
# either likelihood.
#
# Seeds are 20260909 + r for replicate r; depmixS4's random EM starts
# use 900000 + r.

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
source("dev/latent-hmm-sim.R")

args <- commandArgs(trailingOnly = TRUE)
ARM <- if (length(args) >= 1L) args[[1L]] else "A"
REPS <- if (length(args) >= 2L) as.integer(args[[2L]]) else 20L
OUT <- if (length(args) >= 3L) args[[3L]] else
  paste0("dev/latent-2p3-hmm-", ARM, ".tsv")
NID <- if (length(args) >= 4L) as.integer(args[[4L]]) else 8L
NS <- if (length(args) >= 5L) as.integer(args[[5L]]) else 50L
TL <- if (length(args) >= 6L) as.integer(args[[6L]]) else 500L
NDEP <- 4L                      # depmixS4 random EM starts per replicate
K <- 3L

# THE IDENTITY AND THE RECOVERY TABLE ARE DIFFERENT QUESTIONS AND COST
# DIFFERENT AMOUNTS. An identity is a property of the two codes and one
# data set establishes it; a recovery table is a property of the
# estimator and needs replicates. The third-party arm therefore runs on
# the FIRST `NID` replicates only, which are the first NID seeds and so
# are not selected on their outcome, while every replicate contributes
# to the recovery table. depmixS4's random-start EM cost 8 to 224
# seconds per replicate at this size, and one hmmTMB fit is about two
# minutes, so running either on every replicate spends the whole budget
# re-establishing something already established.

elapsed <- function(expr) {
  t0 <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# The fitted transition matrix, from the tr{i}{j} intercepts, with
# state 1 the reference cell of every row.
frm_tpm <- function(b) {
  eta <- matrix(0, K, K - 1L)
  for (i in seq_len(K)) {
    for (j in 2:K) {
      eta[i, j - 1L] <- unname(b[paste0("tr", i, j, ".(Intercept)")])
    }
  }
  hmm_tpm(eta)
}

frm_mu <- function(b) unname(b[paste0("mu", seq_len(K), ".(Intercept)")])
frm_sd <- function(b) {
  exp(unname(b[paste0("sigma", seq_len(K), ".(Intercept)")]))
}

# Wald coverage of one confint() row against a truth, as 1 or 0, plus
# the row's own half width. NA when the row is missing.
cover_of <- function(ci, row, truth) {
  i <- match(row, rownames(ci))
  if (is.na(i)) return(c(NA_real_, NA_real_, NA_real_))
  c(as.integer(truth > ci[i, "lwr"] && truth < ci[i, "upr"]),
    ci[i, "est"] - truth,
    (ci[i, "upr"] - ci[i, "lwr"]) / (2 * stats::qnorm(0.975)))
}

tr_rows <- as.vector(t(outer(seq_len(K), 2:K,
                             function(i, j) paste0("tr", i, j))))
tr_truth <- as.vector(t(hmm_truth$eta))

hdr <- c("arm", "rep", "seed", "rows", "perm", "identity_ok",
         "ll_frm", "ll_ref", "ll_rel_gap", "est_max_gap",
         "tpm_max_gap", "mu_err_max", "sd_err_max", "tpm_err_max",
         "sd_tr12", "s_frm", "s_ref", "maxgrad_rel", "pdhess",
         "nlminb_fn", "nlminb_gr",
         paste0("cov_", c(paste0("mu", 1:3), paste0("sigma", 1:3),
                          tr_rows, "theta1")),
         paste0("err_", c(paste0("mu", 1:3), paste0("sigma", 1:3),
                          tr_rows, "theta1")),
         paste0("se_", c(paste0("mu", 1:3), paste0("sigma", 1:3),
                         tr_rows, "theta1")))
if (!file.exists(OUT)) {
  cat(paste(hdr, collapse = "\t"), "\n", sep = "", file = OUT)
}

run_A <- function(seed, want_ref) {
  s <- hmm_sim_fixed(seed = seed, ns = NS, tl = TL)
  d <- s$d
  form <- bf(y ~ 1)
  fam <- hmm(K = K, gaussian(), time = t, group = id, init = "estimated")
  fit <- NULL
  s_frm <- elapsed(fit <- suppressWarnings(suppressMessages(
    frm(form, family = fam, data = d))))

  if (!want_ref) {
    return(list(fit = fit, d = d, truth = s, s_frm = s_frm,
                s_ref = NA_real_, ll_ref = NA_real_, ref = NULL))
  }
  dm <- depmixS4::depmix(y ~ 1, data = d, nstates = K,
                         ntimes = as.integer(table(d$id)))
  set.seed(900000L + seed %% 100000L)
  best <- -Inf
  bf_ <- NULL
  s_ref <- elapsed({
    for (i in seq_len(NDEP)) {
      ff <- try(suppressMessages(depmixS4::fit(
        dm, verbose = FALSE, emcontrol = depmixS4::em.control(
          random.start = TRUE, tol = 1e-12, maxit = 5000))),
        silent = TRUE)
      if (!inherits(ff, "try-error")) {
        v <- as.numeric(depmixS4::logLik(ff))
        if (v > best) { best <- v; bf_ <- ff }
      }
    }
  })
  list(fit = fit, d = d, truth = s, s_frm = s_frm, s_ref = s_ref,
       ll_ref = best, ref = bf_)
}

run_B <- function(seed, want_ref) {
  s <- hmm_sim(seed = seed, ns = NS, tl = TL)
  d <- s$d
  form <- bf(y ~ 1, tr12 ~ 1 + (1 | id))
  fam <- hmm(K = K, gaussian(), time = t, group = id, init = "uniform")
  fit <- NULL
  s_frm <- elapsed(fit <- suppressWarnings(suppressMessages(
    frm(form, family = fam, data = d))))

  if (!want_ref) {
    return(list(fit = fit, d = d, truth = s, s_frm = s_frm,
                s_ref = NA_real_, ll_ref = NA_real_, ref = NULL,
                hid = NULL))
  }
  dh <- data.frame(ID = factor(d$id), t = d$t, y = d$y)
  fmat <- matrix("~1", K, K)
  diag(fmat) <- "."
  fmat[1L, 2L] <- "~s(ID, bs = \"re\")"
  h0 <- suppressMessages(hmmTMB::MarkovChain$new(
    data = dh, n_states = K, formula = fmat,
    initial_state = "estimated"))
  ld <- h0$delta0(log = TRUE, as_matrix = FALSE)
  fp <- stats::setNames(rep(NA_real_, length(ld)), rownames(ld))
  hid <- suppressMessages(hmmTMB::MarkovChain$new(
    data = dh, n_states = K, formula = fmat,
    initial_state = "estimated", fixpar = list(delta0 = fp)))
  obs <- hmmTMB::Observation$new(
    data = dh, n_states = K, dists = list(y = "norm"),
    par = list(y = list(mean = hmm_truth$mu,
                        sd = rep(hmm_truth$sigma, K))))
  hm <- hmmTMB::HMM$new(obs = obs, hid = hid)
  s_ref <- elapsed(suppressWarnings(hm$fit(silent = TRUE)))
  list(fit = fit, d = d, truth = s, s_frm = s_frm, s_ref = s_ref,
       ll_ref = hm$llk(), ref = hm, hid = hid)
}

for (r in seq_len(REPS)) {
  seed <- 20260909L + r
  want_ref <- r <= NID
  z <- if (identical(ARM, "A")) run_A(seed, want_ref) else
    run_B(seed, want_ref)
  fit <- z$fit
  b <- unlist(fixef(fit))
  ll_frm <- as.numeric(logLik(fit))
  mu <- frm_mu(b)
  sg <- frm_sd(b)
  G <- frm_tpm(b)

  al <- hmm_align(mu, hmm_truth$mu)
  p <- al$perm
  Gp <- G[p, p, drop = FALSE]
  mu_err <- max(abs(mu[p] - hmm_truth$mu))
  sd_err <- max(abs(sg[p] - hmm_truth$sigma))
  Gt <- hmm_tpm(hmm_truth$eta)
  tpm_err <- max(abs(Gp - Gt))

  ## ---- identity against the third party ----------------------------
  est_gap <- NA_real_
  tpm_gap <- NA_real_
  ll_gap <- NA_real_
  if (!want_ref) {
    NULL
  } else if (identical(ARM, "A")) {
    pr <- depmixS4::getpars(z$ref)
    # depmixS4 lays out: K prior parameters, then K rows of K
    # transition probabilities, then per-state (mean, sd)
    off <- K
    Gr <- matrix(pr[(off + 1L):(off + K * K)], K, K, byrow = TRUE)
    off <- off + K * K
    mr <- pr[off + seq(1L, 2L * K, by = 2L)]
    sr <- pr[off + seq(2L, 2L * K, by = 2L)]
    ar <- hmm_align(mu, mr)
    est_gap <- max(abs(mu[ar$perm] - mr), abs(sg[ar$perm] - sr))
    tpm_gap <- max(abs(G[ar$perm, ar$perm, drop = FALSE] - Gr))
  } else {
    hp <- z$ref$par()
    mr <- as.numeric(hp$obspar["y.mean", , 1])
    sr <- as.numeric(hp$obspar["y.sd", , 1])
    Gr <- matrix(as.numeric(hp$tpm[, , 1]), K, K)
    ar <- hmm_align(mu, mr)
    sdr <- as.numeric(z$hid$sd_re()[1L, 1L])
    est_gap <- max(abs(mu[ar$perm] - mr), abs(sg[ar$perm] - sr),
                   abs(sqrt(VarCorr(fit)[[1L]][1L, 1L]) - sdr))
    # DO NOT READ `tpm_max_gap` AS A DISAGREEMENT ON THIS ARM. `G` is
    # the POPULATION transition matrix, with the random effect at zero;
    # hmmTMB's par()$tpm[, , 1] is the matrix at DATA ROW 1, which
    # carries sequence 1's own random intercept. The two differ by 0.04
    # to 0.18 for that reason alone. dev/latent-2p3-tpmgap.R adds the
    # conditional mode back and closes it to 1.8e-06; the coefficients
    # of row 1, which reference state 1 under both conventions, agree
    # there to 3.3e-04 of one standard error.
    tpm_gap <- max(abs(G[ar$perm, ar$perm, drop = FALSE] - Gr))
  }
  if (want_ref) ll_gap <- abs(ll_frm - z$ll_ref) / abs(ll_frm)

  ## ---- coverage of the fit's own intervals -------------------------
  ci <- suppressWarnings(stats::confint(fit))
  rows <- c(paste0("mu", 1:3, "_(Intercept)"),
            paste0("sigma", 1:3, "_(Intercept)"),
            paste0(tr_rows, "_(Intercept)"),
            if (identical(ARM, "B")) "theta_1" else NA_character_)
  truth <- c(hmm_truth$mu, rep(log(hmm_truth$sigma), 3), tr_truth,
             if (identical(ARM, "B")) log(hmm_truth$sd_tr12) else
               NA_real_)
  # coverage is only meaningful when the fitted labels ARE the true
  # ones; the permutation is recorded either way
  ident <- identical(as.integer(p), seq_len(K))
  cv <- matrix(NA_real_, length(rows), 3L)
  if (ident) {
    for (i in seq_along(rows)) {
      if (!is.na(rows[i])) cv[i, ] <- cover_of(ci, rows[i], truth[i])
    }
  }

  d1 <- frmtmb::diagnose(fit, quiet = TRUE)
  ev <- fit$opt$evaluations
  row <- c(ARM, r, seed, nrow(z$d), paste(p, collapse = ""),
           as.integer(ident),
           formatC(c(ll_frm, z$ll_ref), digits = 12, format = "f"),
           formatC(c(ll_gap, est_gap, tpm_gap, mu_err, sd_err, tpm_err,
                     if (identical(ARM, "B"))
                       sqrt(VarCorr(fit)[[1L]][1L, 1L]) else NA_real_,
                     z$s_frm, z$s_ref, d1$max_grad / abs(ll_frm)),
                   digits = 6, format = "g"),
           isTRUE(d1$pdHess),
           if (is.null(ev)) c(NA, NA) else as.integer(ev[1:2]),
           formatC(c(cv[, 1L], cv[, 2L], cv[, 3L]), digits = 6,
                   format = "g"))
  cat(paste(row, collapse = "\t"), "\n", sep = "", file = OUT,
      append = TRUE)
  cat(sprintf("%s rep %3d  ll %.4f ref %s  relgap %s  estgap %s",
              ARM, r, ll_frm,
              if (want_ref) sprintf("%.4f", z$ll_ref) else "-",
              if (want_ref) sprintf("%.2e", ll_gap) else "-",
              if (want_ref) sprintf("%.2e", est_gap) else "-"),
      sprintf("  mu_err %.4f tpm_err %.4f perm %s  %.0fs/%.0fs\n",
              mu_err, tpm_err, paste(p, collapse = ""), z$s_frm,
              if (is.na(z$s_ref)) 0 else z$s_ref))
  flush(stdout())
}
