# Licensing the substitution in nss-17: does frm_lincmt() stand in for
# frm_ode() on the acceptance design?
#
# A fit is determined by its objective and its derivatives, so this
# compares those directly on the whole 30-subject dataset rather than
# fitting twice: same data, same parameters, same random-effect values,
# `frm_ode(ss_extrapolate = TRUE)` against `frm_lincmt()` at its
# default, and `frm_ode(ss_extrapolate = FALSE)` against
# `frm_lincmt(n_ss = 20)`. One tape build each instead of an hour of
# optimizer iterations.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-19-objgrad.R
# Seed 101. Design as nss-17-accept-lin.R.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode); library(RTMB)})
nss_report_env()

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

two_oral <- function(t, y, p) {
  list(c(-p[4L] * y[1L],
         p[4L] * y[1L] - (p[1L] + p[2L]) * y[2L] + p[3L] * y[3L],
         p[2L] * y[2L] - p[3L] * y[3L]))
}

# theta = (lke, lka, lk12, lk21, lV, log sigma) plus 2 * NS random
# effects, exactly the parameter vector the fit carries
mk <- function(engine, ext) {
  function(th) {
    "c" <- RTMB::ADoverload("c")
    b <- th[1:6]
    bk <- th[6 + seq_len(NS)]
    ba <- th[6 + NS + seq_len(NS)]
    ke <- exp(b[1] + bk)[d$id]
    ka <- exp(b[2] + ba)[d$id]
    pr <- if (engine == "ode") {
      frm_ode(two_oral, init = list(0, 0, 0), times = d$time,
              group = d$id,
              parms = list(ke, rep(exp(b[3]), nrow(d)),
                           rep(exp(b[4]), nrow(d)), ka),
              events = EVO, output = 2L, n_ss = 20L, ss_tol = Inf,
              ss_extrapolate = ext) / exp(b[5])
    } else {
      frm_lincmt(parms = list(ke = ke, k12 = rep(exp(b[3]), nrow(d)),
                              k21 = rep(exp(b[4]), nrow(d)), ka = ka,
                              V = rep(exp(b[5]), nrow(d))),
                 times = d$time, group = d$id, ncmt = 2, depot = TRUE,
                 events = EVL, n_ss = if (ext) Inf else 20L)
    }
    -sum(dnorm(y, pr, exp(b[6]), log = TRUE)) +
      0.5 * sum(bk^2) / 0.04 + 0.5 * sum(ba^2) / 0.09
  }
}

th0 <- c(log(KE), log(KA), log(K12), log(K21), log(V), log(0.15),
         lke - log(KE), lka - log(KA))
cmp <- function(tag, ext) {
  a <- MakeTape(mk("ode", ext), th0)
  b <- MakeTape(mk("lin", ext), th0)
  va <- a(th0); vb <- b(th0)
  ga <- as.numeric(a$jacfun()(th0))
  gb <- as.numeric(b$jacfun()(th0))
  cat(sprintf("%-28s objective %16.9f vs %16.9f  rel %9.2e\n", tag,
              va, vb, abs(va - vb) / abs(vb)))
  cat(sprintf("%-28s gradient max rel %9.2e over %d elements\n", "",
              max(abs(ga - gb)) / max(abs(gb)), length(gb)))
  cat(sprintf("%-28s d/dlk21  %14.8f vs %14.8f\n", "", ga[4L], gb[4L]))
  invisible(NULL)
}
cat("\n")
cmp("tail summed vs n_ss = Inf", TRUE)
cmp("truncated vs n_ss = 20", FALSE)

cat("\nand the two frm_lincmt arms against each other, which is the\n",
    "contrast nss-17 fits:\n", sep = "")
bi <- MakeTape(mk("lin", TRUE), th0)
b2 <- MakeTape(mk("lin", FALSE), th0)
cat(sprintf("  objective %16.9f vs %16.9f  difference %g\n",
            bi(th0), b2(th0), b2(th0) - bi(th0)))
gi <- as.numeric(bi$jacfun()(th0))
g2 <- as.numeric(b2$jacfun()(th0))
cat(sprintf("  d/dlk21   %14.8f vs %14.8f\n", gi[4L], g2[4L]))
