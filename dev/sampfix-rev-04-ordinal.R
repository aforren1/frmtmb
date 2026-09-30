# Reviewer, lane sampfix, script 04: the ordinal renaming on the families
# and shapes the worker did not run.
#
#   Rscript dev/sampfix-rev-04-ordinal.R ref    # saves whole ds objects
#   Rscript dev/sampfix-rev-04-ordinal.R lane   # samples again, compares
#
# Data seed 606, sampler seed 3, one chain of 150 draws. On the lane the
# reference build's SAVED draws objects (fit and all, as 0.13.0 wrote
# them) are read back as old saved draws.
arm <- commandArgs(trailingOnly = TRUE)[1L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("ARM ", arm, " frmtmb.sample from ", dirname(find.package("frmtmb.sample")),
    "\n", sep = "")
LOG <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-rev-log"

set.seed(606L)
n <- 300L
do <- data.frame(x = rnorm(n), z = rnorm(n),
                 fc = factor(sample(c("a", "b", "c"), n, TRUE)),
                 grp = factor(sample(c("p", "q"), n, TRUE)),
                 g = factor(sample(1:10, n, TRUE)))
ct <- function(eta, cuts = c(-Inf, -0.5, 0.7, 1.5, Inf)) {
  factor(cut(eta + rlogis(n), cuts, labels = FALSE), ordered = TRUE)
}
ug <- rnorm(10, 0, 0.5)
do$yo <- ct(0.8 * do$x + ug[do$g])
do$yo2 <- ct(-0.5 * do$z)
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
fits <- list(
  acat = q(bf(yo ~ x), family = acat(), data = do),
  cratio = q(bf(yo ~ x), family = cratio(), data = do),
  sratio = q(bf(yo ~ x), family = sratio(), data = do),
  cratio_csfac = q(bf(yo ~ x + cs(fc)), family = cratio(), data = do),
  acat_thresgr = q(bf(yo | thres(gr = grp) ~ x), family = acat(), data = do),
  cum_thresgr_probit = q(bf(yo | thres(gr = grp) ~ x),
                         family = cumulative("probit"), data = do),
  mv_cum_sratio = q(bf(yo ~ x, family = cumulative()) +
                      bf(yo2 ~ z, family = sratio()), data = do),
  cum_re = q(bf(yo ~ x + (1 | g)), family = cumulative(), data = do)
)
mix <- tryCatch(q(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
                  data = do), error = function(e) e)
if (inherits(mix, "error")) {
  cat("mixture(cumulative, cumulative) refused by frm(): ",
      substr(conditionMessage(mix), 1, 150), "\n", sep = "")
  do$yc <- c(rnorm(n / 2, -2), rnorm(n / 2, 2)) + 0.3 * do$x
  fits$mix_gauss <- q(bf(yc ~ x), family = mixture(gaussian(), gaussian()),
                      data = do)
} else fits$mix_cum <- mix

smp <- function(fit, ...) suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3, ...)))
out <- list()
for (nm in names(fits)) {
  ds <- tryCatch(smp(fits[[nm]]), error = function(e) e)
  if (inherits(ds, "error")) {
    cat(nm, " SAMPLE ERROR ", conditionMessage(ds), "\n"); next
  }
  out[[nm]] <- ds
}
# a prior on the thresholds, brms spelling
out$acat_prior <- smp(fits$acat,
                      prior = set_prior("normal(0, 1)", class = "Intercept"))
out$cum_re_lap <- tryCatch(smp(fits$cum_re, laplace = TRUE),
                           error = function(e) e)

if (identical(arm, "ref")) {
  saveRDS(out, file.path(LOG, "04-ordinal-ref.rds"))
  cat("saved\nDONE\n")
  quit(save = "no")
}

ref <- readRDS(file.path(LOG, "04-ordinal-ref.rds"))
ep <- function(ds) {
  rs <- names(ds$fit$spec$responses)
  re <- if (frmtmb.sample:::draws_is_laplace(ds)) NA else NULL
  lapply(rs, function(r) posterior_epred(ds, ndraws = 15, resp = r,
                                        re_formula = re))
}
for (nm in names(out)) {
  a <- out[[nm]]; b <- ref[[nm]]
  cat("\n-- ", nm, "\n", sep = "")
  if (inherits(a, "error")) { cat("  lane error: ", conditionMessage(a), "\n"); next }
  fit <- a$fit
  vf <- variables(fit)
  cat("  variables(ds):  ", paste(setdiff(variables(a), "lp__"), collapse = " "), "\n", sep = "")
  cat("  variables(fit) minus sd_/cor_ all in draws: ",
      all(vf[!grepl("^(sd|cor)_", vf)] %in% variables(a)), "\n", sep = "")
  nc <- frmtmb.sample:::draws_natural_cols(fit)
  bm <- b$draws
  for (o in nc$ordinal) {
    j <- match(o$internal, colnames(bm))
    bm[, j] <- frmtmb.sample:::draws_rowmap(bm[, j, drop = FALSE], o$map)
    colnames(bm)[j] <- o$names
  }
  cat(sprintf("  ref draws mapped: names same %s, max|diff| %.3g\n",
              identical(colnames(bm), colnames(a$draws)),
              max(abs(bm - a$draws))))
  ea <- ep(a); eb <- ep(b)
  # the old object read by the new code
  eo <- tryCatch(ep(b), error = function(e) e)
  cat(sprintf("  epred lane vs ref max|diff| %.3g\n",
              max(mapply(function(x, y) max(abs(x - y)), ea, eb))))
  fa <- fixef(a); fb <- fixef(b)
  cat(sprintf("  fixef rows same %s, max|diff| %.3g\n",
              identical(rownames(fa), rownames(fb)), max(abs(fa - fb))))
  ra <- tryCatch(posterior_predict(a, ndraws = 5), error = function(e) e)
  cat("  posterior_predict: ", if (inherits(ra, "error")) conditionMessage(ra) else "OK", "\n", sep = "")
  # old saved object (fit and all from the reference build) under lane
  cat("  OLD saved object: variables ", paste(setdiff(variables(b), "lp__"), collapse = " "), "\n", sep = "")
  fo <- tryCatch(fixef(b), error = function(e) e)
  so <- tryCatch(summary(b), error = function(e) e)
  po <- tryCatch(posterior_summary(b), error = function(e) e)
  cat("  OLD: fixef ", if (inherits(fo, "error")) conditionMessage(fo) else
    sprintf("rows %s", paste(rownames(fo), collapse = ",")),
    "; summary ", if (inherits(so, "error")) conditionMessage(so) else "OK",
    "; posterior_summary ", if (inherits(po, "error")) conditionMessage(po) else "OK",
    "\n", sep = "")
  # readers by name on the new object
  ex <- function(lab, e) {
    r <- tryCatch(suppressWarnings(e), error = function(e) e)
    cat(sprintf("  %-34s %s\n", lab, if (inherits(r, "error"))
      paste("ERROR:", substr(gsub("\n", " ", conditionMessage(r)), 1, 140)) else "OK"))
    invisible(r)
  }
  thr <- grep("Intercept\\[", colnames(a$draws), value = TRUE)
  cat("  threshold columns: ", paste(thr, collapse = " "), "\n", sep = "")
  if (length(thr) >= 2) {
    h1 <- sub("^b_", "", thr[1]); h2 <- sub("^b_", "", thr[2])
    r <- ex(sprintf("hypothesis('%s < %s')", h1, h2),
            hypothesis(a, paste(h1, "<", h2)))
    if (!inherits(r, "error")) {
      direct <- mean(a$draws[, thr[1]] < a$draws[, thr[2]])
      cat(sprintf("    Post.Prob %.4f, direct from the columns %.4f\n",
                  r$hypothesis$Post.Prob, direct))
    }
  }
  ex("as_draws_df", posterior::as_draws_df(a))
  ex("posterior_summary", posterior_summary(a))
  if (length(thr)) {
    ex("posterior_summary(variable = thr[1])", posterior_summary(a, variable = thr[1]))
    pdf(NULL); ex("mcmc_plot(variable = thr)", print(mcmc_plot(a, variable = thr))); dev.off()
    ex("as.matrix(variable = thr)", as.matrix(a, variable = thr))
  }
  ex("summary", summary(a))
  ex("conditional_effects", conditional_effects(a, effects = "x"))
}
cat("\n-- named-list prior spelled as the new draws column\n")
r <- tryCatch(smp(fits$acat, prior = list("b_Intercept[1]" = prior_normal(0, 1))),
              error = function(e) e)
cat("  ", if (inherits(r, "error")) paste("ERROR:", substr(conditionMessage(r), 1, 200)) else
  "OK (accepted)", "\n", sep = "")
r <- tryCatch(smp(fits$acat, prior = list(tau_raw = prior_normal(0, 1))),
              error = function(e) e)
cat("  internal name tau_raw: ", if (inherits(r, "error")) paste("ERROR:", substr(conditionMessage(r), 1, 200)) else
  "OK (accepted)", "\n", sep = "")
cat("\nDONE\n")
