# Review round 2: frm_sample() through the tmp_xi chart against brms on
# the round-1 data (n = 40, y ~ 1; seeds 1401, 1403, 1402), default
# priors in both, 4 chains x 2000, seed 1; and the warnings for cens(),
# trunc() and two GEV responses.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb from", find.package("frmtmb"), "; sample from",
    find.package("frmtmb.sample"), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
divs <- function(sf) {
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0))
}
qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
qrow <- function(v) paste(sprintf("%8.4f", quantile(v, qs)), collapse = " ")
for (cfg in list(c(-0.45, 1401), c(0.2, 1403), c(-0.85, 1402))) {
  set.seed(cfg[2])
  d <- data.frame(y = rg(40, 0, 1, cfg[1]))
  fit <- suppressWarnings(frm(y ~ 1, family = gen_extreme_value(), data = d))
  msgs <- character(0); ws <- character(0)
  dr <- withCallingHandlers(frm_sample(fit, chains = 4, iter = 2000, seed = 1,
                                       refresh = 0),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  dm <- as.data.frame(dr)
  bf <- suppressWarnings(suppressMessages(brms::brm(y ~ 1, data = d,
          family = brms::gen_extreme_value(), chains = 4, iter = 2000,
          seed = 1, backend = "rstan", refresh = 0, silent = 2)))
  bd <- brms::as_draws_df(bf)
  cat(sprintf("\n=== xi true %.2f, seed %d: divergences frm %d, brms %d of 4000\n",
              cfg[1], cfg[2], divs(dr$stanfit), divs(bf$fit)))
  cat("frm_sample warnings:", if (length(ws)) unique(substr(ws, 1, 90)) else "none", "\n")
  cat("prior lines:", grep("xi|sigma|Intercept", unlist(strsplit(msgs, "\n")), value = TRUE), sep = "\n  ")
  cat("quantiles            2.5%      25%      50%      75%    97.5%\n")
  cat("xi        frm  ", qrow(dm$xi), "\n          brms ", qrow(bd$xi), "\n")
  cat("Intercept frm  ", qrow(dm$b_Intercept), "\n          brms ", qrow(bd$b_Intercept), "\n")
  cat("sigma     frm  ", qrow(dm$sigma), "\n          brms ", qrow(bd$sigma), "\n")
  cat(sprintf("brms P(xi < -1) = %.4f; frm min xi %.6f; brms min xi %.4f\n",
              mean(bd$xi < -1), min(dm$xi), min(bd$xi)))
  # brms's posterior conditioned on xi > -1, the region frmtmb holds
  if (mean(bd$xi < -1) > 0) {
    k <- bd$xi > -1
    cat("brms | xi > -1:  xi", qrow(bd$xi[k]), "\n                 Intercept",
        qrow(bd$b_Intercept[k]), "\n                 sigma", qrow(bd$sigma[k]), "\n")
  }
}

cat("\n=== warnings where the chart is not used\n")
set.seed(1401)
d <- data.frame(y = rg(40, 0, 1, -0.45), y2 = rg(40, 0, 1, -0.3))
d$cc <- as.numeric(d$y > quantile(d$y, 0.9))
d$lb <- min(d$y) - 0.5
show <- function(lab, fit) {
  ws <- character(0)
  withCallingHandlers(suppressMessages(frm_sample(fit, chains = 1, iter = 300,
                                                  seed = 1, refresh = 0)),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  hit <- grep("tmp_xi", ws, value = TRUE)
  cat(lab, ":", if (length(hit)) substr(hit, 1, 200) else "NO tmp_xi WARNING", "\n")
}
show("cens()", frm(y | cens(cc) ~ 1, family = gen_extreme_value(), data = d))
show("trunc(lb)", frm(y | trunc(lb = lb) ~ 1, family = gen_extreme_value(), data = d))
show("two GEV responses",
     suppressWarnings(frm(bf(y ~ 1) + gen_extreme_value() + bf(y2 ~ 1) +
                            gen_extreme_value(), data = d)))
show("xi ~ 1 (log1p by design, no chart)",
     frm(bf(y ~ 1, xi ~ 1), family = gen_extreme_value(), data = d))
show("mixture of two GEVs",
     suppressWarnings(frm(y ~ 1, family = mixture(gen_extreme_value(),
                                                  gen_extreme_value()), data = d)))
