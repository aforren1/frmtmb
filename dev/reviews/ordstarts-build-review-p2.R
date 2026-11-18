# Reviewer of lane ordstarts, build round. Probe 2 (questions 2, 4, 7):
# - "one kind": pairs the test calls identical that are not (priors per
#   component, fixed or predicted weights, offsets, constant disc), and the
#   weights' reference moved by a theta2 formula.
# - refit paths: a counter on mixture_swap_candidates() and on the swaps
#   run under refit(), frm_bootstrap(), influence(), frm_allfit(),
#   update(), start =, autoscale.
# - test-ordinal-mixture.R F3 case: the warnings in full on the new default.
# - test-compat.R cells: which cell warns, and starts = 10.
#   Rscript ordstarts-build-review-p2.R > ordstarts-build-review-p2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
q <- function(e) suppressWarnings(suppressMessages(e))
rec <- function(f) {
  ms <- f$opt$mixture_starts
  if (is.null(ms)) return("record NULL")
  paste(names(ms), signif(ms, 4), sep = "=", collapse = " ")
}
# the label switch of a K = 2 ordinal mixture done by hand: every mu1_/mu2_
# (and disc) coefficient, the two raw threshold vectors, and every weight
# coefficient negated
hand_swap <- function(f) {
  p <- f$opt$par
  nm <- ns$outer_par_map(f)$names
  s <- p
  for (pre in c("mu", "disc")) {
    a <- grep(paste0("^", pre, "1_"), nm)
    b <- grep(paste0("^", pre, "2_"), nm)
    ka <- sub(paste0("^", pre, "1_"), "", nm[a])
    kb <- sub(paste0("^", pre, "2_"), "", nm[b])
    if (!length(a) || !identical(sort(ka), sort(kb))) next
    s[a[order(ka)]] <- p[b[order(kb)]]
    s[b[order(kb)]] <- p[a[order(ka)]]
  }
  t1 <- which(names(p) == "tau_raw1")
  t2 <- which(names(p) == "tau_raw2")
  if (length(t1) == length(t2) && length(t1)) {
    s[t1] <- p[t2]
    s[t2] <- p[t1]
  }
  w <- grep("^theta[0-9]+_", nm)
  s[w] <- -p[w]
  s
}
one_kind <- function(label, f) {
  sw <- tryCatch(ns$mixture_ord_swap_starts(f, "y"), error = function(e) e)
  if (inherits(sw, "error")) {
    cat(label, ": swap_starts error:", conditionMessage(sw), "\n")
    return(invisible())
  }
  s <- hand_swap(f)
  v0 <- f$obj$fn(f$opt$par)
  v1 <- f$obj$fn(s)
  # a refit from the hand swap, the objective's own optimizer
  o <- tryCatch(stats::nlminb(s, f$obj$fn, f$obj$gr,
                              control = list(eval.max = 2000,
                                             iter.max = 1000)),
                error = function(e) NULL)
  cat(sprintf(paste0("%-44s | swaps the build runs %d | record: %s |",
                     " objective at the hand swap - at the fit %+.4g |",
                     " refit from it: objective - fit's %s (conv %s)\n"),
              label, length(sw), rec(f), v1 - v0,
              if (is.null(o)) "error" else sprintf("%+.4f", o$objective - v0),
              if (is.null(o)) "-" else o$convergence))
  f$obj$fn(f$opt$par)
  invisible()
}

cat("\n== 2. What the build calls one kind ==\n")
for (seed in 1:6) {
  d <- mk(seed, 300)
  d$o <- 0.8 * d$z
  cc <- mixture(cumulative(), cumulative())
  f <- q(frm(bf(y ~ x + z), family = cc, data = d))
  one_kind(paste("seed", seed, "cum_cum, no prior (control)"), f)
  for (pm in c(1.5, -0.8)) {
    pr <- set_prior(sprintf("normal(%g, 0.2)", pm), class = "b", coef = "x",
                    dpar = "mu1")
    f <- tryCatch(q(frm(bf(y ~ x + z), family = cc, data = d, prior = pr)),
                  error = function(e) e)
    if (inherits(f, "error")) {
      cat("prior fit error:", conditionMessage(f), "\n")
    } else {
      one_kind(paste("seed", seed, "cum_cum, prior N(", pm, ",0.2) on mu1 x"),
               f)
    }
  }
  if (seed > 2) next
  f <- tryCatch(q(frm(bf(y ~ x + z, theta1 = 0.8), family = cc, data = d)),
                error = function(e) e)
  if (inherits(f, "error")) cat("theta1 = 0.8:", conditionMessage(f), "\n") else
    one_kind(paste("seed", seed, "cum_cum, theta1 = 0.8 fixed"), f)
  f <- tryCatch(q(frm(bf(y ~ x + z, mu1 ~ x + z + offset(o), mu2 ~ x + z),
                      family = cc, data = d)), error = function(e) e)
  if (inherits(f, "error")) cat("offset:", conditionMessage(f), "\n") else
    one_kind(paste("seed", seed, "cum_cum, offset in mu1 only"), f)
  f <- tryCatch(q(frm(bf(y ~ x + z, disc2 = 2), family = cc, data = d)),
                error = function(e) e)
  if (inherits(f, "error")) cat("disc2 = 2:", conditionMessage(f), "\n") else
    one_kind(paste("seed", seed, "cum_cum, disc2 = 2 fixed"), f)
}

cat("\n== 2b. A theta2 formula moves the weights' reference ==\n")
d <- mk(13, 1000)
cs <- mixture(cumulative(), sratio())
for (form in list(bf(y ~ x + z), bf(y ~ x + z, theta1 ~ 1),
                  bf(y ~ x + z, theta2 ~ 1))) {
  f0 <- q(frm(form, family = cs, data = d,
              control = frmtmb_control(starts = 0)))
  nm <- ns$outer_par_map(f0)$names
  sw <- ns$mixture_ord_swap_starts(f0, "y")
  w <- grep("^theta", nm)
  fD <- q(frm(form, family = cs, data = d))
  cat(sprintf("%-28s weight %s = %.4f | swapped start's weight %.4f | D: %s | ll A0 %.4f D %.4f\n",
              paste(deparse(form$pforms), collapse = " "),
              nm[w], f0$opt$par[w], sw[[1]]$par[w], rec(fD),
              as.numeric(logLik(f0)), as.numeric(logLik(fD))))
}

cat("\n== 4. Refit paths: swap candidates built (counter) ==\n")
cnt <- new.env()
cnt$n <- 0L
orig_cand <- ns$mixture_swap_candidates
utils::assignInNamespace("mixture_swap_candidates", function(fit, resps) {
  cnt$n <- cnt$n + 1L
  orig_cand(fit, resps)
}, "frmtmb")
d <- mk(13, 300)
step <- function(label, expr) {
  n0 <- cnt$n
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(q(expr), error = function(e) e)
  cat(sprintf("%-46s swap candidate calls %d  (%.1f s)%s\n", label,
              cnt$n - n0, proc.time()[["elapsed"]] - t0,
              if (inherits(v, "error")) paste(" ERROR", conditionMessage(v))
              else ""))
  v
}
f <- step("frm() default", frm(bf(y ~ x + z), family = cs, data = d))
cat("  record:", rec(f), "\n")
step("refit(f, newresp)", refit(f, d$y))
step("frm_bootstrap(f, nsim = 3)", frm_bootstrap(f, nsim = 3, seed = 1))
step("influence(f)", influence(f))
step("frm_allfit(f) start_from_mle = TRUE", frm_allfit(
  f, optimizers = list(nlminb = "nlminb", optim = "optim")))
step("frm_allfit(f) start_from_mle = FALSE", frm_allfit(
  f, optimizers = list(nlminb = "nlminb", optim = "optim"),
  start_from_mle = FALSE))
u <- step("update(f)", update(f))
cat("  update record:", rec(u), "\n")
step("update(f, newdata = mk(14, 300))", update(f, newdata = mk(14, 300)))
st <- f$estimates[c("beta", "tau_raw1", "tau_raw2")]
step("frm(start = own)", frm(bf(y ~ x + z), family = cs, data = d,
                             start = st))
step("frm(control = starts = 0)", frm(bf(y ~ x + z), family = cs, data = d,
                                      control = frmtmb_control(starts = 0)))
d2 <- d
d2$x <- d2$x * 1000
a <- step("frm(autoscale = TRUE, x * 1000)",
          frm(bf(y ~ x + z), family = cs, data = d2,
              control = frmtmb_control(autoscale = TRUE)))
if (!inherits(a, "error")) cat("  autoscale record:", rec(a), "\n")
utils::assignInNamespace("mixture_swap_candidates", orig_cand, "frmtmb")

cat("\n== 7a. test-ordinal-mixture.R F3 case on the new default ==\n")
set.seed(1)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- stats::rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
d$y <- 1L + rowSums(outer(lat, stats::quantile(lat, 1:3 / 4), ">"))
for (ctl in list(frmtmb_control(), frmtmb_control(starts = 10))) {
  w <- character()
  fit <- withCallingHandlers(
    suppressMessages(frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()),
                         data = d, control = ctl)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  fx <- suppressWarnings(fixef(fit))
  cat("control starts =", format(ctl$starts), "| logLik",
      format(as.numeric(logLik(fit)), digits = 10), "code",
      fit$opt$convergence, "| record:", rec(fit), "\n")
  cat("  NaN SE rows:", paste(rownames(fx)[is.nan(fx[, "Est.Error"])],
                             collapse = ", "), "\n")
  cat("  warnings:\n", paste0("   - ", w, collapse = "\n"), "\n")
}

cat("\n== 7b. test-compat.R ordinal mixture cells ==\n")
set.seed(20261005)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
d$o <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
cells <- list(
  list(bf(o ~ x, disc1 ~ 0 + z), mixture(cumulative(), sratio())),
  list(bf(o ~ x), mixture(cumulative(threshold = "sum_to_zero"), sratio())),
  list(bf(o ~ x), mixture(cratio(), cumulative(threshold = "sum_to_zero"),
                          order = "mu")),
  list(bf(o ~ x), mixture(cumulative(threshold = "equidistant"), sratio())))
for (i in seq_along(cells)) {
  for (ctl in list(frmtmb_control(), frmtmb_control(starts = 10),
                   frmtmb_control(starts = 0))) {
    w <- character()
    fit <- withCallingHandlers(
      suppressMessages(frm(cells[[i]][[1]], family = cells[[i]][[2]],
                           data = d, control = ctl)),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    cat(sprintf("cell %d starts=%-4s ll %.4f code %d | %s | degen warn %s\n",
                i, format(ctl$starts), as.numeric(logLik(fit)),
                fit$opt$convergence, rec(fit),
                any(grepl("degenerate boundary", w, fixed = TRUE))))
    if (length(w)) cat("   ", substr(w, 1, 300), sep = "\n    ")
    cat("\n")
  }
}
cat("DONE p2\n")
