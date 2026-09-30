# Reviewer re-check, lane sampfix: attack draws_laplace_watch(), the
# probe's re-raise, and draws_as_caller(). Data seed 77, draws seed 1
# (lap_pair() as in dev/sampfix-rev-01-probe.R).
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})

lap_pair <- function(fit, n = 6L, seed = 1L, tweak = NULL) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  if (!is.null(tweak)) M <- tweak(M)
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  full <- structure(list(stanfit = NULL, draws = M, fit = fit),
                    class = "frmtmb_draws")
  lap <- full
  lap$draws <- M[, setdiff(colnames(M), inner), drop = FALSE]
  stopifnot(frmtmb.sample:::draws_is_laplace(lap))
  list(full = full, lap = lap)
}
res <- function(expr) tryCatch(suppressMessages(suppressWarnings(expr)),
                               error = function(e) e)
show <- function(lab, a, b = NULL) {
  f <- function(r) if (inherits(r, "error"))
    paste0("ERROR[", substr(gsub("\\s+", " ", conditionMessage(r)), 1, 90), "]") else
      sprintf("OK nonfinite=%d/%d", sum(!is.finite(unlist(r))), length(unlist(r)))
  cat(sprintf("%-48s laplace: %s\n%-48s full:    %s\n", lab, f(a), "",
              if (is.null(b)) "-" else f(b)))
}

set.seed(77)
G <- 8; n <- 10
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
u <- rnorm(G, 0, 0.6)[dd$g]
dd$y <- 0.5 + 0.4 * dd$x + u + rnorm(nrow(dd), 0, 0.7)
dd$yn <- 0.5 * dd$x + exp(0.3 + u)^0.8 + rnorm(nrow(dd), 0, 0.3)
q <- function(...) suppressWarnings(suppressMessages(frm(...)))

cat("== 1a. watch hole: draw 1 non-finite everywhere for its own reason\n")
fnl <- q(bf(yn ~ log(c1) * x + exp(a)^k, c1 ~ 1, a ~ 1 + (1 | g), k ~ 1,
            nl = TRUE), family = gaussian(), data = dd,
         prior = set_prior("normal(1.6, 1)", nlpar = "c1"))
cat("nl fit converged:", isTRUE(fnl$opt$convergence == 0), " c1 =",
    exp(0) * fnl$estimates$beta[1], "\n")
tw <- function(M) {
  M[1, "b_c1_Intercept"] <- -1   # log(-1) is NaN on draw 1, every row
  M[1, "b_k_Intercept"] <- 0     # NA^0 = 1 hides a on draw 1
  M
}
p <- lap_pair(fnl, tweak = tw)
a <- res(posterior_epred(p$lap)); b <- res(posterior_epred(p$full))
show("posterior_epred default", a, b)
if (!inherits(a, "error") && !inherits(b, "error")) {
  cat(sprintf("  draws 2..6: laplace non-finite %d, full non-finite %d (of %d)\n",
              sum(!is.finite(a[-1, ])), sum(!is.finite(b[-1, ])), length(a[-1, ])))
}

cat("== 1b. watch false alarm: a draw that is non-finite for its own reason\n")
p <- lap_pair(q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
sd1 <- VarCorr(p$full, summary = FALSE)$g$sd[1, 1]
sds <- VarCorr(p$full, summary = FALSE)$g$sd[, 1]
cat(sprintf("  sd draws: %s\n", paste(format(sds, digits = 5), collapse = " ")))
h <- sprintf("exp(1e5 * (sd_g__Intercept - %.10f)) > 1", sd1)
a <- res(hypothesis(p$lap, h, class = NULL)); b <- res(hypothesis(p$full, h, class = NULL))
show("hypothesis exp overflow on later draws", a, b)
h2 <- sprintf("log(sd_g__Intercept - %.10f) > -5", sd1 - 1e-9)
a <- res(hypothesis(p$lap, h2, class = NULL)); b <- res(hypothesis(p$full, h2, class = NULL))
show("hypothesis log of a negative on later draws", a, b)
# a real-model case: a lognormal expectation that overflows on one draw
fl <- q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
dd$pos <- exp(dd$y)
flog <- q(bf(pos ~ x + (1 | g)), family = lognormal(), data = dd)
pl <- lap_pair(flog, tweak = function(M) { M[3, "b_Intercept"] <- 720; M })
a <- res(posterior_epred(pl$lap, re_formula = NA))
b <- res(posterior_epred(pl$full, re_formula = NA))
show("lognormal epred re NA, draw 3 overflows", a, b)

cat("== 2. probe re-raise: calls that work on full draws\n")
p <- lap_pair(q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
for (cl in list(
  quote(posterior_epred(D, re_formula = NA)),
  quote(posterior_predict(D, re_formula = NA)),
  quote(conditional_effects(D, effects = "x", resolution = 5)),
  quote(conditional_effects(D, effects = "x", resolution = 5, re_formula = NULL, seed = 2)),
  quote(hypothesis(D, "sd_g__Intercept > 0.1", class = NULL)),
  quote(posterior_epred(D, newdata = data.frame(x = c(1, NA), g = factor(1, levels = 1:8)), re_formula = NA)),
  quote(posterior_epred(D, newdata = data.frame(x = 1, g = NA), re_formula = NA)),
  quote(fitted(D, re_formula = NA, summary = FALSE))
)) {
  a <- res(eval(cl, list(D = p$lap))); b <- res(eval(cl, list(D = p$full)))
  same <- if (!inherits(a, "error") && !inherits(b, "error")) identical(a, b) else NA
  show(paste(deparse(cl, width.cutoff = 200L), collapse = ""), a, b)
  cat("  identical:", same, "\n")
}

cat("== 3. caller names, nesting, reset\n")
st <- function() frmtmb.sample:::draws_call_state$what
first <- function(expr) {
  r <- res(expr)
  m <- if (inherits(r, "error")) gsub("\\s+", " ", conditionMessage(r)) else "NO ERROR"
  sprintf("%s | hint=%s | state after=%s", substr(m, 1, 70),
          grepl("re_formula = NA", m, fixed = TRUE), format(st()))
}
L <- p$lap
cat("loo_compare:           ", first(loo_compare(L, L)), "\n")
cat("loo:                   ", first(loo(L)), "\n")
cat("pp_check loo weights:  ", first(pp_check(L, type = "loo_pit_overlay")), "\n")
cat("pp_check dens:         ", first(pp_check(L, ndraws = 2)), "\n")
cat("hypothesis scope coef: ", first(hypothesis(L, "Intercept > 0", scope = "coef", group = "g")), "\n")
cat("residuals pearson:     ", first(residuals(L, type = "pearson")), "\n")
cat("predictive_error epred:", first(predictive_error(L, method = "posterior_epred")), "\n")
cat("fitted linear:         ", first(fitted(L, scale = "linear")), "\n")
cat("bayes_R2:              ", first(bayes_R2(L)), "\n")
cat("-- reset after an error raised by something other than the refusal\n")
cat("fitted(bad arg):       ", first(fitted(L, re_formula = NA, probs = "x")), "\n")
cat("then posterior_epred:  ", first(posterior_epred(L)), "\n")
cat("-- a user function nesting two calls\n")
uf <- function(d) { try(fitted(d), silent = TRUE); posterior_predict(d) }
cat("user fn:               ", first(uf(L)), "\n")
cat("-- interrupt-like: an error thrown from inside the wrapped expression\n")
r <- res(frmtmb.sample:::draws_as_caller("x()", stop("boom")))
cat("state after inner stop:", format(st()), "\n")
cat("DONE\n")
