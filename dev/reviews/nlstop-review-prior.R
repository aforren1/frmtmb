# Review of lane nlstop: prior_conflict_escape() false triggers, the point
# kept, the message, and the cost of the restart. Arms in one process:
# "on" is the lane's build, "off" replaces prior_conflict_escape() by the
# identity (0.81.0's path). A trace records whether the restart ran.
# Usage: Rscript dev/reviews/nlstop-review-prior.R [case ...]
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
orig <- ns$prior_conflict_escape
ran <- new.env()
cnt <- new.env(); cnt$k <- 0L
oo <- ns$optimize_obj
utils::assignInNamespace("optimize_obj", function(...) {
  cnt$k <- cnt$k + 1L
  oo(...)
}, "frmtmb")
traced <- function(obj, opt, frame, entries, control, bounds,
                   par_units = NULL, tally = NULL) {
  t0 <- proc.time()[["elapsed"]]
  n0 <- cnt$k
  o <- orig(obj, opt, frame, entries, control, bounds, par_units, tally)
  ran$t <- proc.time()[["elapsed"]] - t0
  ran$n <- cnt$k - n0
  ran$from <- opt$objective
  ran$to <- o$objective
  ran$won <- !is.null(o$prior_escape)
  o
}
set_arm <- function(arm) {
  utils::assignInNamespace("prior_conflict_escape",
    if (arm == "on") traced else function(obj, opt, ...) opt, "frmtmb")
}
args <- commandArgs(TRUE)
want <- function(k) !length(args) || k %in% args

cases <- list()
# A1: data-dominated, weak prior: slope 10, normal(0, 1), n = 5000
if (want("a1")) {
  set.seed(11); n <- 5000
  d <- data.frame(x = rnorm(n)); d$y <- 1 + 10 * d$x + rnorm(n)
  cases$a1 <- list(f = bf(y ~ x), d = d, fam = gaussian(),
                   pr = set_prior("normal(0, 1)", class = "b"))
}
# A2: brms-style student_t(3, 0, 2.5): standardized (no conflict) and an
# unstandardized predictor in units 100 times smaller (slope 300)
if (want("a2")) {
  set.seed(12); n <- 1000
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n) / 100)
  d$y <- 0.5 + 3 * d$x1 + 300 * d$x2 + rnorm(n, 0, 2)
  cases$a2s <- list(f = bf(y ~ x1), d = d, fam = gaussian(),
                    pr = set_prior("student_t(3, 0, 2.5)", class = "b"))
  cases$a2u <- list(f = bf(y ~ x1 + x2), d = d, fam = gaussian(),
                    pr = set_prior("student_t(3, 0, 2.5)", class = "b"))
  d$yb <- rbinom(n, 1, plogis(-0.3 + 1 * d$x1 + 200 * d$x2))
  cases$a2b <- list(f = bf(yb ~ x1 + x2), d = d, fam = bernoulli(),
                    pr = set_prior("student_t(3, 0, 2.5)", class = "b"))
}
# A3: many coefficients in conflict at once (20 slopes at 3, normal(0,
# 0.5)), and a mix (10 at 3, 10 at 0.2)
if (want("a3")) {
  set.seed(13); n <- 2000; p <- 20
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  d <- as.data.frame(X)
  d$y <- drop(X %*% rep(3, p)) + rnorm(n, 0, 2)
  d$ym <- drop(X %*% rep(c(3, 0.2), each = 10)) + rnorm(n, 0, 2)
  rhs <- paste(colnames(X), collapse = " + ")
  cases$a3all <- list(f = bf(as.formula(paste("y ~", rhs))), d = d,
                      fam = gaussian(),
                      pr = set_prior("normal(0, 0.5)", class = "b"))
  cases$a3mix <- list(f = bf(as.formula(paste("ym ~", rhs))), d = d,
                      fam = gaussian(),
                      pr = set_prior("normal(0, 0.5)", class = "b"))
}
# A4: a large GLMM where the restart triggers: binomial, 20000 rows, 200
# groups, 5 slopes at 2, normal(0, 0.25) on b
if (want("a4")) {
  set.seed(14); n <- 20000; ng <- 200
  d <- data.frame(g = factor(sample(ng, n, TRUE)))
  for (k in 1:5) d[[paste0("x", k)]] <- rnorm(n)
  u <- rnorm(ng, 0, 0.7)
  eta <- -0.5 + 2 * (d$x1 + d$x2 + d$x3 + d$x4 + d$x5) + u[d$g]
  d$y <- rbinom(n, 1, plogis(eta))
  cases$a4 <- list(f = bf(y ~ x1 + x2 + x3 + x4 + x5 + (1 | g)), d = d,
                   fam = bernoulli(),
                   pr = set_prior("normal(0, 0.25)", class = "b"))
  # gaussian LMM, 20000 rows, slope 5 under normal(0, 1)
  d$yg <- 1 + 5 * d$x1 + u[d$g] + rnorm(n)
  cases$a4g <- list(f = bf(yg ~ x1 + x2 + (1 | g)), d = d, fam = gaussian(),
                    pr = set_prior("normal(0, 1)", class = "b"))
}
# I1: the excluded case, a centered intercept in conflict
if (want("i1")) {
  set.seed(6201); n <- 200
  d <- data.frame(x = rnorm(n))
  d$yl <- 1e6 * (1 + 0.05 * d$x) + rnorm(n, 0, 1e3)
  cases$i1 <- list(f = bf(yl ~ 1 + x), d = d, fam = gaussian(),
                   pr = set_prior("normal(0, 1000)", class = "Intercept"))
  cases$i1nc <- list(f = bf(yl ~ 0 + Intercept + x), d = d,
                     fam = gaussian(),
                     pr = set_prior("normal(0, 1000)", class = "b",
                                    coef = "Intercept"))
}
nrep <- as.integer(Sys.getenv("REVIEW_NREP", "3"))
cat("ROWHEAD case arm rep secs llp code nwarn wclass msg ran restart_s",
    "restart_evals won\n")
for (cn in names(cases)) {
  cs <- cases[[cn]]
  for (rep in seq_len(nrep)) for (arm in if (rep %% 2) c("on", "off") else
                                    c("off", "on")) {
    set_arm(arm)
    rm(list = ls(ran), envir = ran)
    t0 <- proc.time()[["elapsed"]]
    r <- tt(frm(cs$f, family = cs$fam, data = cs$d, prior = cs$pr))
    secs <- proc.time()[["elapsed"]] - t0
    mc <- attr(r, "mc")
    cat(sprintf("ROW %s %s %d %.2f %.6f %s %d %s %s %s %s %s %s\n", cn, arm,
                rep, secs, llp(r),
                if (iserr(r)) paste("ERR", r$msg) else r$opt$convergence,
                length(attr(r, "w")), wsum(r),
                any(mc == "frmtmb_prior_multistart"),
                isTRUE(ran$n > 0),
                if (is.null(ran$t)) NA else sprintf("%.2f", ran$t),
                if (is.null(ran$n)) NA else ran$n,
                if (is.null(ran$won)) NA else ran$won))
    if (rep == 1L && arm == "on" && !iserr(r)) {
      cat("  beta:", format(r$estimates$beta, digits = 6), "\n")
      if (length(attr(r, "m"))) cat("  msg:", attr(r, "m"), sep = "\n  ")
      cat("\n")
    }
    flush.console()
  }
}
cat("DONE\n")
