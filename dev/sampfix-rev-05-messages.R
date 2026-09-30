# Reviewer, lane sampfix, script 05: the full text of every refusal of
# laplace draws, per method, with three checks on each: does it name the
# function the user called, does it give the working route (sample
# without laplace = TRUE), and does any re_formula = NA hint name an
# argument the called function accepts?
#   Rscript dev/sampfix-rev-05-messages.R   (data seed 77, draws seed 1)
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})

lap_pair <- function(fit, n = 6L, seed = 1L) {
  tpl <- fit$frame[["par_template"]]
  est <- unlist(lapply(names(tpl), function(cp) fit$estimates[[cp]]))
  lab <- frmtmb::brms_par_labels(fit)
  set.seed(seed)
  M <- matrix(rep(est, each = n) + stats::rnorm(n * length(est), 0, 0.05),
              n, dimnames = list(NULL, lab))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  inner <- setdiff(lab, frmtmb::brms_par_labels(fit, include_random = FALSE))
  lap <- structure(list(stanfit = NULL, draws = M[, setdiff(colnames(M), inner),
                                                  drop = FALSE], fit = fit),
                   class = "frmtmb_draws")
  stopifnot(frmtmb.sample:::draws_is_laplace(lap))
  lap
}
takes_re <- function(fname) {
  f <- getS3method(fname, "frmtmb_draws", optional = TRUE)
  !is.null(f) && "re_formula" %in% names(formals(f))
}
msg <- function(called, expr, model) {
  pdf(NULL)
  r <- tryCatch(suppressMessages(suppressWarnings(expr)),
                error = function(e) e)
  dev.off()
  if (!inherits(r, "error")) {
    cat(sprintf("[%s] %-22s COMPUTED (%s)\n", model, called, class(r)[1]))
    return(invisible())
  }
  m <- gsub("\\s+", " ", conditionMessage(r))
  names_called <- grepl(paste0(called, "()"), m, fixed = TRUE)
  route <- grepl("without laplace = TRUE|laplace = FALSE", m)
  hint <- grepl("re_formula = NA", m, fixed = TRUE)
  hint_ok <- !hint || takes_re(called)
  cat(sprintf("[%s] %-22s names-called=%s route=%s re_formula-hint=%s%s\n    %s\n",
              model, called, names_called, route, hint,
              if (hint && !hint_ok) " (CALLED FUNCTION HAS NO re_formula)" else "",
              m))
}

set.seed(77)
G <- 8; n <- 10
dd <- data.frame(g = factor(rep(seq_len(G), each = n)))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(G, 0, 0.6)[dd$g] + rnorm(nrow(dd), 0, 0.7)
q <- function(...) suppressWarnings(suppressMessages(frm(...)))

lap <- lap_pair(q(bf(y ~ x + (1 | g)), family = gaussian(), data = dd))
M <- "(1|g)"
msg("posterior_epred", posterior_epred(lap), M)
msg("posterior_linpred", posterior_linpred(lap), M)
msg("posterior_predict", posterior_predict(lap), M)
msg("fitted", fitted(lap), M)
msg("predict", predict(lap), M)
msg("residuals", residuals(lap), M)
msg("predictive_error", predictive_error(lap), M)
msg("predictive_interval", predictive_interval(lap), M)
msg("bayes_R2", bayes_R2(lap), M)
msg("loo_R2", loo_R2(lap), M)
msg("pp_check", pp_check(lap, ndraws = 3), M)
msg("pp_check", pp_check(lap, type = "loo_pit_overlay"), M)
msg("conditional_effects", conditional_effects(lap, effects = "x",
                                               re_formula = NULL), M)
msg("hypothesis", hypothesis(lap, "Intercept > 0", scope = "ranef",
                             group = "g"), M)
msg("ranef", ranef(lap), M)
msg("coef", coef(lap), M)
msg("log_lik", log_lik(lap), M)
msg("loo", loo(lap), M)
msg("waic", waic(lap), M)
msg("psis", psis(lap), M)
msg("loo_compare", loo_compare(lap, lap), M)
msg("loo_moment_match", loo_moment_match(lap), M)
msg("loo_subsample", loo_subsample(lap), M)

lap_s <- lap_pair(q(bf(y ~ s(x, k = 5)), family = gaussian(), data = dd))
M <- "s(x)"
msg("posterior_epred", posterior_epred(lap_s), M)
msg("posterior_epred", posterior_epred(lap_s, re_formula = NA), M)
msg("conditional_effects", conditional_effects(lap_s, effects = "x"), M)
msg("fitted", fitted(lap_s), M)

set.seed(6)
dm <- data.frame(g = factor(rep(1:8, each = 15L)))
dm$y <- c(rnorm(60, -2), rnorm(60, 3)) + rnorm(8, 0, 0.3)[dm$g]
lap_m <- lap_pair(q(bf(y ~ 1 + (1 | g)),
                    family = frmtmb::mixture(gaussian(), gaussian()), data = dm))
M <- "mixture"
msg("pp_mixture", pp_mixture(lap_m), M)
msg("pp_mixture", pp_mixture(lap_m, re_formula = NA), M)
cat("DONE\n")
