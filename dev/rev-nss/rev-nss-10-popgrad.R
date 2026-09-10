# REVIEW of lane nss, attack 4 at population scale.
#
# The lane and the shipped NEWS both say: on a 30-subject dataset,
# truncating the run-in moves the objective by 60.66 units and turns
# `d/dlog(k21)` from +2.915 into -294.63. That number reaches users
# through NEWS.md, so it has to reproduce.
#
# It is measured through `frm_lincmt()`, and the substitution is
# licensed by the agreement between the two engines' objectives and
# gradients. Both halves are checked here in one table: all four arms
# side by side rather than two pairwise comparisons, so the licence and
# the contrast can be read off the same rows.
#
# Script path: dev/rev-nss/rev-nss-10-popgrad.R
# Seed 101. Design as dev/nss/nss-17-accept-lin.R.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
rev_env()

KE <- 0.15; K12 <- 0.3; K21 <- 0.02; KA <- 1.0; V <- 10; II <- 24
NS <- 30L
set.seed(101L)
tt <- II * c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1)
d <- data.frame(id = rep(seq_len(NS), each = length(tt)),
                time = rep(tt, NS))
lke <- log(KE) + rnorm(NS, 0, 0.2)
lka <- log(KA) + rnorm(NS, 0, 0.3)
EVL <- data.frame(time = 0, state = "depot", value = 100, ii = II,
                  addl = 0L, ss = TRUE)
EVO <- data.frame(time = 0, state = 1L, value = 100, ii = II,
                  ss = TRUE)
mu <- numeric(nrow(d))
for (j in seq_len(NS)) {
  k <- which(d$id == j)
  mu[k] <- frm_lincmt(parms = list(ke = exp(lke[[j]]), k12 = K12,
                                   k21 = K21, ka = exp(lka[[j]]),
                                   V = V),
                      times = d$time[k], ncmt = 2, depot = TRUE,
                      events = EVL)
}
y <- mu + rnorm(nrow(d), 0, 0.15)

two_oral <- function(t, y, p)
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))

mk <- function(engine, ext) function(th) {
  "c" <- RTMB::ADoverload("c")
  b <- th[1:6]
  bk <- th[6 + seq_len(NS)]
  ba <- th[6 + NS + seq_len(NS)]
  ke <- exp(b[1] + bk)[d$id]
  ka <- exp(b[2] + ba)[d$id]
  pr <- if (engine == "ode")
    frm_ode(two_oral, init = list(0, 0, 0), times = d$time,
            group = d$id,
            parms = list(ke, rep(exp(b[3]), nrow(d)),
                         rep(exp(b[4]), nrow(d)), ka),
            events = EVO, output = 2L, n_ss = 20L, ss_tol = Inf,
            ss_extrapolate = ext) / exp(b[5])
  else
    frm_lincmt(parms = list(ke = ke, k12 = rep(exp(b[3]), nrow(d)),
                            k21 = rep(exp(b[4]), nrow(d)), ka = ka,
                            V = rep(exp(b[5]), nrow(d))),
               times = d$time, group = d$id, ncmt = 2, depot = TRUE,
               events = EVL, n_ss = if (ext) Inf else 20L)
  -sum(dnorm(y, pr, exp(b[6]), log = TRUE)) +
    0.5 * sum(bk^2) / 0.04 + 0.5 * sum(ba^2) / 0.09
}

th0 <- c(log(KE), log(KA), log(K12), log(K21), log(V), log(0.15),
         lke - log(KE), lka - log(KA))

arms <- list(
  `frm_ode ss_extrapolate=TRUE`  = list("ode", TRUE),
  `frm_ode ss_extrapolate=FALSE` = list("ode", FALSE),
  `frm_lincmt n_ss=Inf`          = list("lin", TRUE),
  `frm_lincmt n_ss=20`           = list("lin", FALSE))
res <- list()
cat(sprintf("\n%-30s %18s %14s %14s\n", "arm", "objective",
            "d/dlk21", "d/dlke"))
for (an in names(arms)) {
  tp <- MakeTape(mk(arms[[an]][[1L]], arms[[an]][[2L]]), th0)
  v <- tp(th0)
  g <- as.numeric(tp$jacfun()(th0))
  res[[an]] <- list(v = v, g = g)
  cat(sprintf("%-30s %18.6f %14.4f %14.4f\n", an, v, g[4L], g[1L]))
}
lin <- function(k) res[[k]]
cat("\n-- the licence: does frm_ode match frm_lincmt arm for arm? --\n")
p <- list(c("frm_ode ss_extrapolate=TRUE", "frm_lincmt n_ss=Inf"),
          c("frm_ode ss_extrapolate=FALSE", "frm_lincmt n_ss=20"))
for (q in p) {
  a <- res[[q[[1L]]]]; b <- res[[q[[2L]]]]
  cat(sprintf("  %-30s objective rel %9.2e  gradient max rel %9.2e over %d\n",
              q[[1L]], abs(a$v - b$v) / abs(b$v),
              max(abs(a$g - b$g)) / max(abs(b$g)), length(b$g)))
}
cat("\n-- the contrast NEWS.md publishes --\n")
a <- res[["frm_lincmt n_ss=Inf"]]; b <- res[["frm_lincmt n_ss=20"]]
cat(sprintf("  objective moves %.4f units (%.6f -> %.6f)\n",
            b$v - a$v, a$v, b$v))
cat(sprintf("  d/dlk21 %+.4f -> %+.4f\n", a$g[4L], b$g[4L]))
cat("\n-- and through frm_ode itself, which NEWS does not say --\n")
a <- res[["frm_ode ss_extrapolate=TRUE"]]
b <- res[["frm_ode ss_extrapolate=FALSE"]]
cat(sprintf("  objective moves %.4f units\n", b$v - a$v))
cat(sprintf("  d/dlk21 %+.4f -> %+.4f\n", a$g[4L], b$g[4L]))
cat("\ndone\n")
