# The stated unmeasured risk of the reviewer's `n_ss = "auto"`: how far
# does the contraction ratio r move between the starting values, where
# "auto" would have to choose, and the optimum, where the run-in is
# actually performed?
#
# r is a function of the parameters alone, so this needs no ODE solve:
# the vehicle is frm_lincmt(), which reaches the same optimum as
# frm_ode() to seven significant digits (dev/lincmt-findings.md) and
# costs seconds rather than an hour. r = exp(-lambda_z * ii) with
# lambda_z the slow disposition eigenvalue, unless ka is slower, which
# the max() below covers.
#
# Two movements are reported, because "auto" is exposed to both:
#   - the population values move from the start to the optimum;
#   - the run-in is performed PER GROUP, and at the start every random
#     effect is 0, so the worst group at the optimum is not the
#     population value at the start.
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-08-rmove.R
# Seeds: 101..106 per design, recorded in the output.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode/dev/nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
nss_report_env()

lamz <- function(ke, k12, k21) {
  b <- ke + k12 + k21
  (b - sqrt(b * b - 4 * ke * k21)) / 2
}
r_of <- function(ke, k12, k21, ka, ii)
  max(exp(-lamz(ke, k12, k21) * ii), exp(-ka * ii))
n_of <- function(r) if (is.na(r) || r <= 0 || r >= 1) NA_integer_ else
  as.integer(ceiling(log(1e-9) / log(r)))

sim_one <- function(ke, k12, k21, ka, V, ii, ns, seed, sd_obs = 0.15) {
  set.seed(seed)
  tt <- ii * c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 1)
  d <- data.frame(id = factor(rep(seq_len(ns), each = length(tt))),
                  time = rep(tt, ns))
  i <- as.integer(d$id)
  lke <- log(ke) + rnorm(ns, 0, 0.2)
  lka <- log(ka) + rnorm(ns, 0, 0.3)
  ev <- data.frame(time = 0, state = "depot", value = 100, ii = ii,
                   addl = 0L, ss = TRUE)
  mu <- numeric(nrow(d))
  for (j in seq_len(ns)) {
    k <- which(i == j)
    mu[k] <- frm_lincmt(parms = list(ke = exp(lke[[j]]), k12 = k12,
                                     k21 = k21, ka = exp(lka[[j]]),
                                     V = V),
                        times = d$time[k], ncmt = 2, events = ev,
                        depot = TRUE)
  }
  d$conc <- mu + rnorm(nrow(d), 0, sd_obs)
  list(d = d, ev = ev, lke = lke, lka = lka)
}

fit_one <- function(z, start) {
  doses <- z$ev
  bd <- bf(conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                          k21 = exp(lk21), ka = exp(lka),
                                          V = exp(lV)),
                             times = time, group = id, ncmt = 2,
                             depot = TRUE, events = doses),
           lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
           lk21 ~ 1, lV ~ 1, nl = TRUE)
  frm(bd + gaussian(), data = z$d,
      start = if (is.null(start)) NULL else list(beta = start))
}

designs <- list(
  list(tag = "2cmt t1/2 107 h, ii = 24",
       ke = 0.15, k12 = 0.3, k21 = 0.02, ka = 1.0, V = 10, ii = 24),
  list(tag = "2cmt t1/2 265 h, ii = 24",
       ke = 0.1, k12 = 0.2, k21 = 0.008, ka = 1.0, V = 10, ii = 24),
  list(tag = "2cmt t1/2 23 h, ii = 8",
       ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1, V = 10, ii = 8))

starts <- list(
  truth = function(g) c(log(g$ke), log(g$ka), log(g$k12), log(g$k21),
                        log(g$V)),
  # a start a modeller would write down without knowing the answer
  naive = function(g) c(log(0.1), log(1), log(0.1), log(0.1), log(10)),
  # frmtmb's own cold start: every nonlinear beta at 0, so every rate
  # constant is 1 (R/fit.R make_start())
  cold = function(g) c(0, 0, 0, 0, 0))

cat("\nn_start is what 'auto' would choose at the values the tape is\n",
    "built at; n_opt is what the run-in actually needed at the fitted\n",
    "values. group_opt is the worst SUBJECT at the optimum.\n\n", sep = "")
cat(sprintf("%-26s %-8s %5s %8s %8s %8s %8s %8s %9s\n", "design",
            "start", "seed", "r_start", "r_opt", "n_start", "n_opt",
            "n_grp", "short@n_s"))
for (g in designs) {
  for (sname in names(starts)) {
    for (seed in 101:106) {
      z <- sim_one(g$ke, g$k12, g$k21, g$ka, g$V, g$ii, 30L, seed)
      st <- starts[[sname]](g)
      f <- tryCatch(fit_one(z, st), error = function(e) e)
      if (inherits(f, "condition")) {
        cat(sprintf("%-26s %-8s %5d  fit error: %s\n", g$tag, sname,
                    seed, substr(conditionMessage(f), 1, 40)))
        next
      }
      fe <- fixef(f)
      e <- if (is.matrix(fe)) fe[, 1L] else unlist(fe)
      names(e) <- sub("^.*_", "", names(e))
      p <- function(nm) as.numeric(e[[grep(nm, names(e))[1L]]])
      ke <- exp(p("lke")); ka <- exp(p("lka")); k12 <- exp(p("lk12"))
      k21 <- exp(p("lk21")); vv <- exp(p("lV"))
      r_opt <- r_of(ke, k12, k21, ka, g$ii)
      # the worst group: random effects on lke and lka only
      re <- tryCatch(ranef(f), error = function(e) NULL)
      r_grp <- r_opt
      if (!is.null(re)) {
        bke <- tryCatch(as.numeric(re[["lke: 1 | id"]][, 1L]),
                        error = function(e) NULL)
        bka <- tryCatch(as.numeric(re[["lka: 1 | id"]][, 1L]),
                        error = function(e) NULL)
        if (!is.null(bke) && !is.null(bka)) {
          for (j in seq_along(bke)) {
            r_grp <- max(r_grp, r_of(exp(p("lke") + bke[[j]]), k12, k21,
                                     exp(p("lka") + bka[[j]]), g$ii))
          }
        }
      }
      r_start <- if (is.null(st)) NA_real_ else
        r_of(exp(st[1L]), exp(st[3L]), exp(st[4L]), exp(st[2L]), g$ii)
      n_s <- n_of(r_start); n_o <- n_of(r_opt); n_g <- n_of(r_grp)
      short <- if (is.na(n_s)) NA_real_ else r_grp^n_s
      cat(sprintf("%-26s %-8s %5d %8.5f %8.5f %8s %8d %8d %9.2e\n",
                  g$tag, sname, seed, r_start, r_opt,
                  if (is.na(n_s)) "-" else format(n_s), n_o, n_g,
                  short))
    }
  }
}
