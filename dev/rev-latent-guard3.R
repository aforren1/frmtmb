# Reviewer, item 2.3 claim 1 and item 5: three things guard2 left open.
#
#   F. how often does hmm_starts() print "the original fit found a local
#      optimum" on a UNIMODAL fit, where the incumbent IS the answer?
#      guard2 case A hit it on the first try; this measures the rate.
#   G. a refit that converges CLEANLY to something worse, forced with a
#      large jitter around the warm d4 incumbent.
#   H. hmm_starts_scale()'s "unit" fallback: reachable or not, by
#      construction, in both directions.
#
#   Rscript dev/rev-latent-guard3.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
hr <- function(s) cat("\n======== ", s, " ========\n", sep = "")

small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}
d4_data <- function() {
  set.seed(2026)
  K <- 2L; N <- 25L; Tg <- 30L
  G <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
  mu <- c(0, 3); sg <- c(0.6, 0.6); sd_b <- c(0.7, 0.5)
  stat <- local({
    A <- rbind(t(diag(K) - G), 1)
    drop(qr.solve(A, c(rep(0, K), 1)))
  })
  b <- cbind(stats::rnorm(N, 0, sd_b[1L]), stats::rnorm(N, 0, sd_b[2L]))
  d <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tg); s[1L] <- sample.int(K, 1L, prob = stat)
    for (t in seq_len(Tg - 1L)) s[t + 1L] <- sample.int(K, 1L, prob = G[s[t], ])
    data.frame(ID = g, t = seq_len(Tg),
               y = stats::rnorm(Tg, mu[s] + b[g, s], sg[s]))
  }))
  d$gf <- factor(d$ID); d
}
d4_form <- bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))
d4_fam <- hmm(K = 2, gaussian(), time = t, group = ID, init = "stationary")

## ---- F: the false-alarm rate on a unimodal fit ---------------------
hr("F. 'the original fit found a local optimum', on a unimodal fit")
cat("Six independent data sets, one hmm_starts(n = 8, jitter = 2) each.\n")
cat("Every refit reaching the SAME optimum is the correct answer; the\n")
cat("question is what the summary then SAYS.\n\n")
res <- do.call(rbind, lapply(1:6, function(k) {
  dd <- small_data(seed = 4500L + k)
  f <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)
  ms <- suppressWarnings(hmm_starts(f, n = 8, jitter = 2, seed = 8800L + k))
  out <- utils::capture.output(print(ms))
  data.frame(data_seed = 4500L + k,
             modes = nrow(ms$modes),
             spread = ms$spread,
             gap_best_minus_orig = as.numeric(logLik(ms$best)) -
               ms$original_logLik,
             says_local_optimum =
               any(grepl("found a local optimum", out)),
             says_original_best =
               any(grepl("original fit was the best found", out)))
}))
print(res, row.names = FALSE, digits = 6)
cat("\nruns on a single-optimum surface that were told the original\n",
    "fit found a LOCAL optimum: ", sum(res$says_local_optimum), " of ",
    nrow(res), "\n", sep = "")
cat("largest gap that triggered it:",
    format(max(res$gap_best_minus_orig), digits = 3), "\n")
cat("the spread printed beside it:",
    format(max(res$spread), digits = 3), "\n")
cat("`lli > logLik(best)` and `gap > 0` carry no tolerance, so on a\n",
    "unimodal surface the incumbent loses to optimizer noise whenever\n",
    "any of the n refits lands a fraction of a ulp higher.\n", sep = "")

## ---- G: a CONVERGED refit that is worse -----------------------------
hr("G. a refit that converges cleanly to a WORSE optimum")
d4 <- d4_data()
tpl <- par_template(d4_form, family = d4_fam, data = d4)
st <- tpl
st$beta[["mu1_(Intercept)"]] <- -0.185185
st$beta[["mu2_(Intercept)"]] <- 3.121877
st$betad[["sigma1_(Intercept)"]] <- log(0.610455)
st$betad[["sigma2_(Intercept)"]] <- log(0.602054)
st$betad[["tr12_(Intercept)"]] <- log(0.173438 / 0.826562)
st$betad[["tr22_(Intercept)"]] <- log(0.828240 / 0.171760)
st$theta <- log(c(0.645297, 0.427590))
fw <- frm(d4_form, family = d4_fam, data = d4, start = st)
cat("warm incumbent logLik:", format(as.numeric(logLik(fw)), digits = 14),
    "\n")
for (jt in c(4, 8)) {
  ms <- suppressWarnings(hmm_starts(fw, n = 8, jitter = jt, seed = 4242))
  worse <- ms$table[ms$table$status == "converged" &
                      ms$table$logLik < ms$original_logLik - 1e-6, ]
  cat("\njitter", jt, ": converged", ms$n_converged, " not", ms$n_not_converged,
      " error", ms$n_error, " modes", nrow(ms$modes),
      " spread", format(ms$spread, digits = 8), "\n")
  cat("  CONVERGED refits strictly worse than the incumbent:", nrow(worse),
      "\n")
  if (nrow(worse)) {
    cat("  logLiks:", paste(format(worse$logLik, digits = 12),
                            collapse = " "), "\n")
    cat("  grad_rel:", paste(format(worse$grad_rel, digits = 3),
                             collapse = " "), "\n")
    cat("  deepest below the incumbent:",
        format(ms$original_logLik - min(worse$logLik), digits = 8), "\n")
  }
  cat("  best stayed the incumbent:",
      isTRUE(all.equal(as.numeric(logLik(ms$best)), ms$original_logLik,
                       tolerance = 1e-7)), "\n")
  print(ms$modes, row.names = FALSE, digits = 12)
}

## ---- H: the "unit" fallback ----------------------------------------
hr("H. hmm_starts_scale()'s 'unit' fallback, by construction")
cat("The branch is reached when confint() errors, returns the wrong\n")
cat("number of rows, disagrees with fit$opt$par, or reports no usable\n")
cat("standard error at all. Each is tried below.\n")

cat("\n-- H1: can the alignment test EVER fail? --\n")
cat("confint.frmtmb_fit() sets its `est` column to object$opt$par\n")
cat("(R/confint.R:383). The test compares that column with\n")
cat("fit[['opt']][['par']], so it is a vector against itself unless\n")
cat("confint()'s own length(est) != length(se) fallback fires, which\n")
cat("the core comments call defensive only.\n")
dd <- small_data()
f1 <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
          data = dd)
ci <- suppressWarnings(stats::confint(f1))
cat("identical(ci[,'est'], opt$par):",
    identical(unname(as.numeric(ci[, "est"])),
              unname(as.numeric(f1$opt$par))), "\n")
cat("nrow(ci) == length(opt$par):",
    nrow(ci) == length(f1$opt$par), "\n")

cat("\n-- H2: a fit with NO usable standard error --\n")
cat("Over-stating K on data from a single state leaves the extra\n")
cat("state's parameters unidentified. If every standard error comes\n")
cat("back non-finite, `!any(ok)` fires and scale_how is 'unit'.\n")
try_case <- function(label, expr) {
  f <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("  %-28s frm() ERROR: %s\n", label,
                substr(conditionMessage(f), 1, 70)))
    return(invisible(NULL))
  }
  ci <- tryCatch(suppressWarnings(stats::confint(f)),
                 error = function(e) e)
  if (inherits(ci, "error")) {
    cat(sprintf("  %-28s confint() ERROR -> fallback\n", label))
  } else {
    se <- (ci[, "upr"] - ci[, "lwr"]) / (2 * stats::qnorm(0.975))
    cat(sprintf("  %-28s npar %2d  se finite&pos %d of %d\n", label,
                length(f$opt$par), sum(is.finite(se) & se > 0), length(se)))
  }
  sc <- frmtmb.latent:::hmm_starts_scale(f)
  cat(sprintf("  %-28s scale_how = %s\n", label, sc$how))
  invisible(f)
}
set.seed(77)
one_state <- data.frame(id = rep(1:8, each = 25), t = rep(1:25, 8),
                        y = stats::rnorm(200, 0, 1))
try_case("K=2 on one-state data",
         frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
             data = one_state))
try_case("K=3 on one-state data",
         frm(bf(y ~ 1), family = hmm(K = 3, gaussian(), time = t, group = id),
             data = one_state))
set.seed(78)
tiny <- data.frame(id = rep(1:2, each = 4), t = rep(1:4, 2),
                   y = stats::rnorm(8))
try_case("K=3 on 8 rows",
         frm(bf(y ~ 1), family = hmm(K = 3, gaussian(), time = t, group = id),
             data = tiny))
dd2 <- small_data()
dd2$xa <- stats::rnorm(nrow(dd2))
dd2$xb <- dd2$xa
try_case("collinear covariate",
         frm(bf(y ~ xa + xb),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = dd2))
dd3 <- small_data()
dd3$const <- 1
try_case("y constant within state",
         frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = id),
             data = transform(dd3, y = round(y * 0 + rep(c(0, 3),
                                                          length.out = nrow(dd3)),
                                             10))))

cat("\n-- H3: the fallback's own arithmetic, called directly --\n")
cat("A fit whose confint() cannot be produced at all reaches the\n")
cat("branch; the branch itself is exercised here on a stub, to show\n")
cat("what it returns rather than to claim a fit reaches it.\n")
stub <- list(opt = list(par = c(a = 1, b = 2, c = 3)))
class(stub) <- "frmtmb_fit"
sc <- frmtmb.latent:::hmm_starts_scale(stub)
cat("  confint() on a stub errors -> how =", sc$how,
    " s =", paste(sc$s, collapse = " "), "\n")
