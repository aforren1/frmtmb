# Reviewer, claim 1: the four existing ordinal families must not move.
# Fits a battery of cumulative/sratio/cratio/acat models and captures
# every post-fit output in one arm; dev/fams2-rev-ord-cmp.R compares the
# two arms with identical().
#   Rscript dev/fams2-rev-ord.R base|lane
arm <- commandArgs(TRUE)[1]
base_lib <- c("C:/Users/adf44/source/r/rellib-r3",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-fams2-lib",
                               base_lib) else base_lib)
suppressPackageStartupMessages({
  library(frmtmb)
  library(emmeans)
})
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
out_dir <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-out"
dir.create(out_dir, showWarnings = FALSE)

# drop environments and functions so identical() compares values only
strip <- function(x, depth = 0) {
  if (depth > 12) return("<deep>")
  if (is.function(x)) return("<fn>")
  if (is.environment(x)) return("<env>")
  if (inherits(x, "formula")) return(deparse1(x))
  if (inherits(x, "ggplot")) {
    return(list(data = strip(x$data, depth + 1),
                layers = lapply(x$layers, function(l) {
                  tryCatch(strip(l$data, depth + 1), error = function(e) NULL)
                })))
  }
  if (is.list(x)) {
    at <- attributes(x)
    y <- lapply(x, strip, depth = depth + 1)
    at$names <- NULL
    at <- lapply(at, strip, depth = depth + 1)
    attributes(y) <- c(list(names = names(x)), at)
    return(y)
  }
  at <- attributes(x)
  if (length(at)) {
    for (nm in setdiff(names(at), c("dim", "dimnames", "names", "levels",
                                     "class", "row.names"))) {
      attr(x, nm) <- strip(at[[nm]], depth + 1)
    }
  }
  x
}

cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(list(msg = conditionMessage(e)), class = "cap_error")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw))
      invokeRestart("muffleWarning")
    })
  list(value = strip(v), warnings = w)
}

post <- function(fit, d, nd = NULL, emm = NULL, cs_var = NULL) {
  r <- list()
  r$logLik <- cap(logLik(fit))
  r$fixef <- cap(fixef(fit))
  r$coef <- cap(coef(fit))
  r$ranef <- cap(ranef(fit))
  r$vcov <- cap(vcov(fit))
  r$summary <- cap(summary(fit))
  r$variables <- cap(variables(fit))
  r$estimates <- cap(fit$estimates)
  r$fitted <- cap(fitted(fit))
  r$fitted_lin <- cap(fitted(fit, scale = "linear"))
  r$fitted_draws <- cap({ set.seed(11); fitted(fit, summary = FALSE,
                                                ndraws = 20) })
  r$fitted_re0 <- cap(fitted(fit, re_formula = NA))
  r$linpred <- cap(frm_linpred(fit, type = "response"))
  r$predict <- cap({ set.seed(12); predict(fit) })
  r$predict_draws <- cap({ set.seed(13); predict(fit, summary = FALSE,
                                                  ndraws = 20) })
  if (!is.null(nd)) {
    r$fitted_nd <- cap(fitted(fit, newdata = nd))
    r$predict_nd <- cap({ set.seed(14); predict(fit, newdata = nd) })
    r$sim_nd <- cap(simulate(fit, nsim = 3, seed = 15, newdata = nd))
  }
  r$simulate <- cap(simulate(fit, nsim = 5, seed = 16))
  for (ty in c("response", "ordinary", "pearson", "deviance", "osa")) {
    r[[paste0("resid_", ty)]] <- cap(residuals(fit, type = ty))
  }
  r$ce <- cap(conditional_effects(fit))
  r$ce_cat <- cap(conditional_effects(fit, categorical = TRUE))
  r$ce_predict <- cap({ set.seed(17); conditional_effects(fit,
                                                           method = "predict") })
  r$default_prior <- cap(default_prior(fit))
  r$prior_summary <- cap(prior_summary(fit))
  r$dharma <- cap({
    dh <- dharma_residuals(fit, nsim = 50, seed = 18)
    list(res = dh$scaledResiduals, fit = dh$fittedPredictedResponse)
  })
  r$pp_check <- cap({ set.seed(19); pp_check(fit, type = "bars", ndraws = 5) })
  r$log_lik <- cap(log_lik(fit))
  r$confint <- cap(confint(fit))
  r$plot_fitted <- cap({
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    plot(fit, which = 1)
    "plotted"
  })
  if (!is.null(emm)) {
    r$emm <- cap(as.data.frame(emmeans::emmeans(fit, emm)))
    r$emm_epred <- cap(as.data.frame(emmeans::emmeans(fit, emm, epred = TRUE)))
  }
  r
}

set.seed(20260929)
n <- 400
ng <- 20
g <- factor(rep(seq_len(ng), length.out = n))
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
grp <- factor(sample(c("u", "v"), n, TRUE))
ue <- rnorm(ng, 0, 0.7)[g]
lat <- 0.9 * x + ue + 0.5 * (fc == "b") - 0.4 * (fc == "c") + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1, 0.2, 1.3, Inf)))   # 1..4
lat2 <- -0.6 * x + rlogis(n)
y2 <- as.integer(cut(lat2, c(-Inf, -0.5, 0.8, Inf)))     # 1..3
d <- data.frame(y = y, y2 = y2, x = x, g = g, fc = fc, grp = grp)
d$y0 <- d$y - 1L                                         # 0..3
d$yf <- factor(c("lo", "mid", "hi", "top")[d$y],
               levels = c("lo", "mid", "hi", "top"), ordered = TRUE)
# an unused level between two used ones and one at the top
d$yu <- factor(c("lo", "mid", "hi", "top")[d$y],
               levels = c("lo", "never", "mid", "hi", "top", "nobody"),
               ordered = TRUE)
nd <- data.frame(x = c(-1, 0, 1.5), g = factor(c(1, 2, 3), levels = levels(g)),
                 fc = factor(c("a", "b", "c"), levels = levels(fc)),
                 grp = factor(c("u", "v", "u"), levels = levels(grp)))
pr <- set_prior("normal(0, 3)", class = "Intercept")

specs <- list(
  cum_logit_re = list(f = y ~ x + fc + (1 | g), fam = quote(cumulative()),
                      emm = "fc"),
  cum_probit_re = list(f = y ~ x + fc + (1 | g),
                       fam = quote(cumulative(link = "probit")), emm = "fc"),
  cum_cloglog = list(f = y ~ x, fam = quote(cumulative(link = "cloglog"))),
  cum_thresK = list(f = y | thres(5) ~ x, fam = quote(cumulative()),
                    prior = pr),
  sratio_thresK = list(f = y | thres(5) ~ x, fam = quote(sratio()),
                       prior = pr),
  cum_thres_gr = list(f = y | thres(gr = grp) ~ x, fam = quote(cumulative())),
  acat_thres_gr = list(f = y | thres(gr = grp) ~ x, fam = quote(acat())),
  sratio_cs = list(f = y ~ x + cs(fc), fam = quote(sratio()), emm = "fc"),
  cratio_cs = list(f = y ~ x + cs(fc), fam = quote(cratio()), emm = "fc"),
  acat_cs = list(f = y ~ x + cs(fc), fam = quote(acat()), emm = "fc"),
  sratio_cs_num = list(f = y ~ cs(x) + (1 | g), fam = quote(sratio())),
  cum_ofactor = list(f = yf ~ x + (1 | g), fam = quote(cumulative())),
  acat_ofactor = list(f = yf ~ x, fam = quote(acat())),
  cum_unused = list(f = yu ~ x, fam = quote(cumulative())),
  cratio_unused = list(f = yu ~ x, fam = quote(cratio())),
  cum_zero = list(f = y0 ~ x, fam = quote(cumulative())),
  sratio_zero = list(f = y0 ~ x, fam = quote(sratio())),
  cum_yplus1 = list(f = y ~ x, fam = quote(cumulative()))
)

res <- list()
for (nm in names(specs)) {
  s <- specs[[nm]]
  cat("fitting", nm, "\n")
  fit <- cap(frm(s$f, data = d, family = eval(s$fam), prior = s$prior))
  fitv <- tryCatch(frm(s$f, data = d, family = eval(s$fam), prior = s$prior),
                   error = function(e) NULL, warning = function(w) {
                     suppressWarnings(frm(s$f, data = d, family = eval(s$fam),
                                          prior = s$prior))
                   })
  res[[nm]] <- list(fit = if (inherits(fit$value, "cap_error")) fit else
                      list(warnings = fit$warnings))
  if (!is.null(fitv)) {
    res[[nm]]$post <- post(fitv, d, nd = if (!grepl("thres_gr|unused|zero|ofactor",
                                                     nm)) nd, emm = s$emm)
  }
}

# a multivariate ordinal response
cat("fitting mv\n")
mv <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + cumulative() +
                             bf(y2 ~ x) + sratio(), data = d))
res$mv <- list(post = list(
  logLik = cap(logLik(mv)), fixef = cap(fixef(mv)), coef = cap(coef(mv)),
  variables = cap(variables(mv)), summary = cap(summary(mv)),
  fitted = cap(fitted(mv)), fitted_y2 = cap(fitted(mv, resp = "y2")),
  predict = cap({ set.seed(21); predict(mv) }),
  predict_draws = cap({ set.seed(22); predict(mv, summary = FALSE,
                                              ndraws = 10) }),
  simulate = cap(simulate(mv, nsim = 3, seed = 23)),
  resid = cap(residuals(mv, type = "response")),
  resid_y2 = cap(residuals(mv, resp = "y2", type = "pearson")),
  ce = cap(conditional_effects(mv)),
  ce_cat = cap(conditional_effects(mv, categorical = TRUE)),
  default_prior = cap(default_prior(mv)),
  fitted_nd = cap(fitted(mv, newdata = nd))
))
# the default-prior tables before any fit
res$dp_formula <- list(
  cum = cap(default_prior(y ~ x + (1 | g), data = d, family = cumulative())),
  hc_cs = cap(default_prior(y ~ x + cs(fc), data = d, family = sratio())),
  thres = cap(default_prior(y | thres(gr = grp) ~ x, data = d,
                            family = cumulative()))
)
saveRDS(res, file.path(out_dir, paste0("ord-", arm, ".rds")))
cat("saved", length(res), "entries\n")
