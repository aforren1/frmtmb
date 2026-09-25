# Lane mv, 2026-09-25: ordinal families inside a multivariate model.
#
# 1. IDENTITY: with no shared random effect the multivariate
#    log-likelihood is the sum of the univariate fits' log-likelihoods,
#    because the joint density factorizes and every parameter belongs to
#    exactly one factor. The residual is printed at full precision as a
#    numerical check of the bookkeeping (namespaced thresholds), not as
#    evidence about the density.
# 2. The ordinal part against ordinal::clm (an independent ML fit).
# 3. A shared |ID| random effect across an ordinal and a gaussian
#    response against a hand-written RTMB objective (Laplace over the
#    same group effects), at a shared parameter point and at the ML
#    optimum.
#
# Run: Rscript dev/mv-validate-ordinal.R (lane library first). Seeds are
# fixed below.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-mv")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))

set.seed(20260925)
n_g <- 60
m <- 10
n <- n_g * m
g <- factor(rep(seq_len(n_g), each = m))
x <- rnorm(n)
Su <- matrix(c(0.8^2, 0.6 * 0.8 * 0.7, 0.6 * 0.8 * 0.7, 0.7^2), 2)
U <- matrix(rnorm(n_g * 2), n_g) %*% chol(Su)
lat <- 0.9 * x + U[g, 1] + rlogis(n)
o <- cut(lat, c(-Inf, -1.2, 0, 1.1, Inf), labels = FALSE)
o2 <- cut(-0.5 * x + rlogis(n), c(-Inf, -0.5, 0.7, Inf), labels = FALSE)
y1 <- 1 + 0.5 * x + U[g, 2] + rnorm(n, 0, 0.9)
dd <- data.frame(x, g, o, o2, y1)

## 1. the sum identity, two ordinal families and a gaussian
mv <- frm(bf(o ~ x) + cumulative() + bf(o2 ~ x) + sratio() +
            bf(y1 ~ x) + gaussian(), data = dd)
u1 <- frm(bf(o ~ x) + cumulative(), data = dd)
u2 <- frm(bf(o2 ~ x) + sratio(), data = dd)
u3 <- frm(bf(y1 ~ x) + gaussian(), data = dd)
l_mv <- as.numeric(logLik(mv))
l_sum <- as.numeric(logLik(u1)) + as.numeric(logLik(u2)) +
  as.numeric(logLik(u3))
cat(sprintf("mv logLik            %.15g\n", l_mv))
cat(sprintf("sum of univariate    %.15g\n", l_sum))
cat(sprintf("residual             %.3e\n", l_mv - l_sum))
# the same identity at a shared parameter point, which does not lean on
# the three optimizers landing in the same place
set.seed(11)
pm <- mv$opt$par + rnorm(length(mv$opt$par), 0, 0.1)
lm_ <- mv$obj$env$parList(pm)
par_of <- function(fit, from) {
  tpl <- fit$obj$env$parList(fit$opt$par)
  out <- tpl
  for (nm in names(tpl)) {
    if (nm == "beta" || nm == "betad") {
      src <- names(mv$frame$par_template[[nm]])
      dst <- names(fit$frame$par_template[[nm]])
      resp <- names(fit$spec$responses)
      out[[nm]] <- lm_[[nm]][match(paste0(resp, "_", dst), src)]
    } else {
      out[[nm]] <- lm_[[paste0(names(fit$spec$responses), "_", nm)]]
    }
  }
  unlist(out)
}
f1 <- u1$obj$fn(par_of(u1))
f2 <- u2$obj$fn(par_of(u2))
f3 <- u3$obj$fn(par_of(u3))
cat(sprintf("shared point: mv %.15g  sum %.15g  residual %.3e\n",
            -mv$obj$fn(pm), -(f1 + f2 + f3),
            -mv$obj$fn(pm) + (f1 + f2 + f3)))
print(fixef(mv))

## 2. ordinal part against ordinal::clm
cl <- ordinal::clm(factor(o, ordered = TRUE) ~ x, data = dd)
cat(sprintf("clm logLik %.10f   frmtmb univariate %.10f\n",
            as.numeric(logLik(cl)), as.numeric(logLik(u1))))
cat("clm thresholds:", sprintf("%.6f", cl$alpha), " slope:",
    sprintf("%.6f", cl$beta), "\n")
fx <- fixef(mv)
cat("mv  thresholds:",
    sprintf("%.6f", fx[c("o_Intercept[1]", "o_Intercept[2]",
                         "o_Intercept[3]"), 1]),
    " slope:", sprintf("%.6f", fx["o_x", 1]), "\n")

## 3. shared |ID| random effect: frmtmb against a hand-written objective
fid <- frm(bf(o ~ x + (1 | p | g)) + cumulative() +
             bf(y1 ~ x + (1 | p | g)) + gaussian(), data = dd)
print(summary(fid))
gi <- as.integer(g)
nll_ref <- function(q) {
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  S <- frmtmb:::us_sigma(q$theta, 2L)
  Um <- RTMB::matrix(q$u, n_g, 2)
  nll <- -sum(RTMB::dmvnorm(Um, 0, S, log = TRUE))
  tau <- c(q$tr[1], q$tr[1] + exp(q$tr[2]),
           q$tr[1] + exp(q$tr[2]) + exp(q$tr[3]))
  eta <- q$bo * x + Um[gi, 1]
  Fk <- function(k) {
    if (k == 0) return(0 * eta)
    if (k == 4) return(0 * eta + 1)
    RTMB::plogis(tau[k] - eta)
  }
  # category probabilities by indicator sums: data-only selection
  pr <- 0 * eta
  for (k in 1:4) {
    pr <- pr + (o == k) * (Fk(k) - Fk(k - 1))
  }
  nll <- nll - sum(log(pr))
  mu <- q$b1[1] + q$b1[2] * x + Um[gi, 2]
  nll - sum(RTMB::dnorm(y1, mu, exp(q$ls), log = TRUE))
}
q0 <- list(bo = 0, tr = c(-1, 0, 0), b1 = c(0, 0), ls = 0,
           theta = c(0, 0, 0), u = numeric(n_g * 2))
ob <- RTMB::MakeADFun(nll_ref, q0, random = "u", silent = TRUE)
op <- nlminb(ob$par, ob$fn, ob$gr,
             control = list(iter.max = 2000, eval.max = 2000))
cat(sprintf("ML logLik frmtmb     %.10f\n", as.numeric(logLik(fid))))
cat(sprintf("ML logLik reference  %.10f\n", -op$objective))
# shared point: frmtmb's own parameter vector, mapped by name
pf <- fid$obj$env$parList(fid$opt$par)
set.seed(12)
pert <- function(v) v + rnorm(length(v), 0, 0.05)
bo <- pert(pf$beta)
bd <- pert(pf$betad)
th <- pert(pf$theta)
tr <- pert(pf[["o_tau_raw"]])
bn <- names(fid$frame$par_template$beta)
cat("frmtmb beta layout:", bn, "\n")
pfull <- fid$opt$par
pfull[names(pfull) == "beta"] <- bo
pfull[names(pfull) == "betad"] <- bd
pfull[names(pfull) == "theta"] <- th
pfull[names(pfull) == "o_tau_raw"] <- tr
qref <- c(bo = bo[match("o_x", bn)], tr = tr,
          b1 = bo[match(c("y1_(Intercept)", "y1_x"), bn)],
          ls = bd, theta = th)
cat(sprintf("shared point: frmtmb %.12f  reference %.12f  rel %.3e\n",
            -fid$obj$fn(pfull), -ob$fn(qref),
            abs(fid$obj$fn(pfull) - ob$fn(qref)) / abs(ob$fn(qref))))
cat("frmtmb rescaled RE cor:",
    sprintf("%.6f", stats::cov2cor(frmtmb:::varcorr_matrices(fid)[[1]])[1, 2]),
    " reference:",
    sprintf("%.6f", stats::cov2cor(frmtmb:::us_sigma(
      op$par[names(op$par) == "theta"], 2L))[1, 2]), "\n")
