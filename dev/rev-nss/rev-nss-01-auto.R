# REVIEW of lane nss, attack 1: the rejection of `n_ss = "auto"`.
#
# Three questions, in order.
#   A. Does `make_start()` really leave every nonlinear beta at 0 when
#      the user passes no `start`? Asked of par_template() directly, not
#      inferred from the source.
#   B. Does the arithmetic of dev/nss-findings.md section 2 reproduce?
#      n_start = 7 on the 23 h / ii = 8 design at the cold start, and
#      the worst subject at the optimum needing 169.
#   C. THE COMPARISON THE LANE DID NOT MAKE. Its "4.2e-01 where today's
#      blind 20 ships 3.9e-03" puts r_grp^n_start beside a MEASURED
#      trajectory error taken at the population truth. The like-for-like
#      figure for the blind 20 is r_grp^20 on the same worst subject.
#      This script prints it, and prints the floored rule
#      max(20, auto), which the lane never considered.
#
# Script path: dev/rev-nss/rev-nss-01-auto.R
# Seeds: 101..106. Vehicle frm_lincmt(), as in nss-08-rmove.R.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages({library(frmtmb); library(frmtmb.ode)})
rev_env()

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
  list(d = d, ev = ev)
}

make_bd <- function(doses)
  bf(conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                    k21 = exp(lk21), ka = exp(lka),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 2,
                       depot = TRUE, events = doses),
     lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
     lk21 ~ 1, lV ~ 1, nl = TRUE)

## ---- A. the cold start, asked of the template ---------------------
cat("\n== A. par_template() with no `start` ==\n")
z0 <- sim_one(0.2, 0.4, 0.1, 1.1, 10, 8, 30L, 101)
doses <- z0$ev
tpl <- par_template(make_bd(doses) + gaussian(), data = z0$d)
b <- tpl[["beta"]]
cat("beta length:", length(b), "\n")
print(b)
cat("all nonlinear beta exactly 0 (identical):",
    identical(as.numeric(b), rep(0, length(b))), "\n")
cat("number of exact zeros:", sum(b == 0), "of", length(b), "\n")

## ---- B/C. the decisive row, plus the two missing columns ----------
designs <- list(
  list(tag = "23 h, ii 8", ke = 0.2, k12 = 0.4, k21 = 0.1, ka = 1.1,
       V = 10, ii = 8, start = "cold"),
  list(tag = "107 h, ii 24", ke = 0.15, k12 = 0.3, k21 = 0.02,
       ka = 1.0, V = 10, ii = 24, start = "truth"))

cat("\n== B/C. n_start against the worst subject at the optimum ==\n")
cat("short@auto  = r_grp ^ n_start   (what the lane reports)\n")
cat("short@20    = r_grp ^ 20        (TODAY, same subject, same units)\n")
cat("short@floor = r_grp ^ max(20, n_start)\n\n")
cat(sprintf("%-14s %-6s %5s %8s %8s %7s %7s %10s %10s %10s\n",
            "design", "start", "seed", "r_start", "r_grp", "n_start",
            "n_grp", "short@auto", "short@20", "short@floor"))
for (g in designs) {
  for (seed in 101:106) {
    z <- sim_one(g$ke, g$k12, g$k21, g$ka, g$V, g$ii, 30L, seed)
    doses <- z$ev
    st <- if (identical(g$start, "cold")) c(0, 0, 0, 0, 0) else
      c(log(g$ke), log(g$ka), log(g$k12), log(g$k21), log(g$V))
    f <- tryCatch(frm(make_bd(doses) + gaussian(), data = z$d,
                      start = list(beta = st)),
                  error = function(e) e)
    if (inherits(f, "condition")) {
      cat(sprintf("%-14s %-6s %5d  fit error: %s\n", g$tag, g$start,
                  seed, substr(conditionMessage(f), 1, 40)))
      next
    }
    fe <- fixef(f)
    e <- if (is.matrix(fe)) fe[, 1L] else unlist(fe)
    names(e) <- sub("^.*_", "", names(e))
    p <- function(nm) as.numeric(e[[grep(nm, names(e))[1L]]])
    ke <- exp(p("lke")); ka <- exp(p("lka")); k12 <- exp(p("lk12"))
    k21 <- exp(p("lk21"))
    r_opt <- r_of(ke, k12, k21, ka, g$ii)
    re <- tryCatch(ranef(f), error = function(e) NULL)
    r_grp <- r_opt
    if (!is.null(re)) {
      bke <- tryCatch(as.numeric(re[["lke: 1 | id"]][, 1L]),
                      error = function(e) NULL)
      bka <- tryCatch(as.numeric(re[["lka: 1 | id"]][, 1L]),
                      error = function(e) NULL)
      if (!is.null(bke) && !is.null(bka))
        for (j in seq_along(bke))
          r_grp <- max(r_grp, r_of(exp(p("lke") + bke[[j]]), k12, k21,
                                   exp(p("lka") + bka[[j]]), g$ii))
    }
    r_start <- r_of(exp(st[1L]), exp(st[3L]), exp(st[4L]), exp(st[2L]),
                    g$ii)
    n_s <- n_of(r_start); n_g <- n_of(r_grp)
    nf <- if (is.na(n_s)) 20L else max(20L, n_s)
    cat(sprintf("%-14s %-6s %5d %8.5f %8.5f %7s %7s %10.2e %10.2e %10.2e\n",
                g$tag, g$start, seed, r_start, r_grp,
                if (is.na(n_s)) "-" else format(n_s),
                if (is.na(n_g)) "-" else format(n_g),
                if (is.na(n_s)) NA_real_ else r_grp^n_s,
                r_grp^20, r_grp^nf))
  }
}
cat("\ndone\n")
