# Reviewer, lane sampfix, script 03: pp_check()'s loo_* types beyond the
# worker's 2-chain brms comparison.
#   (a) 4 chains of 250 (sampler seed 9, data seed 41): the lw that
#       bayesplot receives equals brms's rule computed by hand (chain ids
#       when every draw is used, one chain otherwise), and differs from
#       the other rule, so the check can tell them apart.
#   (b) every loo_* type bayesplot::available_ppc() lists draws.
#   (c) refusals: laplace draws, cov = TRUE ARMA; dens_overlay still
#       draws on both.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(frmtmb.sample)})

set.seed(41)
dd <- data.frame(g = factor(rep(1:6, each = 8)), t = rep(1:8, 6))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(6, 0, 0.5)[dd$g] + rnorm(nrow(dd), 0, 0.7)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 4, iter = 500, refresh = 0, seed = 9)))
cat("nchains", nchains(ds), "ndraws", ndraws(ds), "\n")

cap <- new.env()
fns <- c("ppc_loo_pit_overlay", "ppc_loo_pit_qq", "ppc_loo_pit_ecdf",
         "ppc_loo_calibration", "ppc_loo_calibration_grouped",
         "ppc_loo_intervals", "ppc_loo_ribbon")
for (fn in fns) {
  fm <- names(formals(get(fn, asNamespace("bayesplot"))))
  lwx <- if ("lw" %in% fm) quote(if (!missing(lw)) lw) else NULL
  psx <- if ("psis_object" %in% fm) quote(if (!missing(psis_object))
    psis_object) else NULL
  trace(fn, where = asNamespace("bayesplot"), print = FALSE,
        tracer = bquote({
          assign(.(fn), list(lw = .(lwx), ps = .(psx)), envir = cap)
        }))
}
lw_of <- function(c) if (!is.null(c$lw)) c$lw else
  stats::weights(c$ps, log = TRUE)
hand <- function(ll, cid) {
  r <- loo::relative_eff(exp(ll), chain_id = cid)
  stats::weights(loo::psis(-ll, r_eff = r), log = TRUE)
}
pdf(NULL)
ll <- log_lik(ds)
n <- nrow(ll)
lw_chains <- hand(ll, rep(1:4, each = n / 4))
lw_one <- hand(ll, rep(1L, n))
cat(sprintf("(a) all draws: chain-id rule vs one-chain rule differ by %.3g\n",
            max(abs(lw_chains - lw_one))))
for (ty in c("loo_pit_overlay", "loo_intervals")) {
  rm(list = ls(cap), envir = cap)
  suppressMessages(suppressWarnings(pp_check(ds, type = ty)))
  got <- lw_of(get(paste0("ppc_", ty), cap))
  cat(sprintf("    %-16s all draws: max|lw - chain rule| %.3g, max|lw - one-chain rule| %.3g\n",
              ty, max(abs(got - lw_chains)), max(abs(got - lw_one))))
}
ids <- seq(1, n, by = 3)
lw_sub <- hand(ll[ids, ], rep(1L, length(ids)))
rm(list = ls(cap), envir = cap)
suppressMessages(suppressWarnings(pp_check(ds, type = "loo_pit_qq",
                                           draw_ids = ids)))
got <- lw_of(cap$ppc_loo_pit_qq)
cat(sprintf("    loo_pit_qq draw_ids by 3: max|lw - one-chain rule on subset| %.3g\n",
            max(abs(got - lw_sub))))
rm(list = ls(cap), envir = cap)
suppressMessages(suppressWarnings(pp_check(ds, type = "loo_pit_qq",
                                           ndraws = n)))
got <- lw_of(cap$ppc_loo_pit_qq)
cat(sprintf("    loo_pit_qq ndraws = all: max|lw - chain rule| %.3g\n",
            max(abs(got - lw_chains))))

cat("(b) every loo_* type bayesplot lists\n")
lt <- grep("^loo_", sub("^ppc_", "", bayesplot::available_ppc("")),
           value = TRUE)
for (ty in lt) {
  args <- list(ds, type = ty)
  if (grepl("grouped", ty)) args$group <- "g"
  r <- tryCatch(suppressMessages(suppressWarnings(do.call(pp_check, args))),
                error = function(e) e)
  if (!inherits(r, "error")) r <- tryCatch({print(r); r},
                                           error = function(e) e)
  cat(sprintf("    %-26s %s\n", ty, if (inherits(r, "error"))
    paste("ERROR:", substr(conditionMessage(r), 1, 110)) else
      paste("OK", class(r)[1])))
}

cat("(c) refusals\n")
lap <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 9,
             laplace = TRUE)))
cat("    laplace draws detected:", frmtmb.sample:::draws_is_laplace(lap), "\n")
show <- function(lab, expr) {
  r <- tryCatch(suppressMessages(suppressWarnings(expr)),
                error = function(e) e)
  cat(sprintf("    %-40s %s\n", lab, if (inherits(r, "error"))
    paste("ERROR:", substr(gsub("\n", " ", conditionMessage(r)), 1, 230)) else
      paste("OK", class(r)[1])))
}
show("laplace loo_pit_overlay", pp_check(lap, type = "loo_pit_overlay"))
show("laplace dens_overlay re NA", pp_check(lap, ndraws = 5, re_formula = NA))
fa <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = gaussian(), data = dd)
dsa <- suppressWarnings(suppressMessages(
  frm_sample(fa, chains = 2, iter = 300, refresh = 0, seed = 9)))
show("ARMA cov=TRUE log_lik", log_lik(dsa))
show("ARMA cov=TRUE loo_pit_overlay", pp_check(dsa, type = "loo_pit_overlay"))
show("ARMA cov=TRUE loo_intervals", pp_check(dsa, type = "loo_intervals"))
show("ARMA cov=TRUE dens_overlay", pp_check(dsa, ndraws = 5))
show("loo_pit (deprecated) on full draws", pp_check(ds, type = "loo_pit"))
dev.off()
cat("DONE\n")
