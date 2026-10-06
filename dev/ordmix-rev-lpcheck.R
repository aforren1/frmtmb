# The log density of frmtmb's ordinal mixtures against brms 2.23.0's
# compiled Stan program, at matched parameter values.
#
# Usage: Rscript dev/ordmix-lpcheck.R <case> [lane|base]
# Output: one block of LP lines per case on stdout (the driver
# dev/ordmix-lpcheck.sh writes dev/ordmix-lpcheck-log/<case>.txt).
#
# For each case: simulate data (seed 20261005 + case number), fit frmtmb
# by maximum likelihood, build brms's program with every prior flat (so
# its log density with adjust_transform = FALSE is the log likelihood),
# and evaluate both at the optimum and at three perturbed parameter
# vectors (seed 99 + case). Each point's frmtmb value is -obj$fn(p); the
# brms value is rstan::log_prob() at the translated parameters. Also
# reports brms's gradient at frmtmb's optimum, which vanishes only when
# the translation is right and the optimum is shared.
#
# `patch_delta = TRUE` rewrites brms's equidistant mixture program the
# way its transformed parameters read it (brms defect 4: the parameters
# block declares a bare `delta` per component, which stanc refuses, and
# the body reads `delta_mu<k>`), so the intended model can be compiled.
args <- commandArgs(TRUE)
case <- if (length(args)) args[1] else "cum2"
arm <- if (length(args) >= 2) args[2] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cache_dir <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-rev-stan-cache"
stan_mod <- function(code) {
  f <- tempfile(fileext = ".stan")
  writeLines(c(code, paste0("// rstan ", packageVersion("rstan"))), f)
  key <- unname(tools::md5sum(f))
  path <- file.path(cache_dir, paste0(key, ".rds"))
  if (file.exists(path)) {
    m <- try(readRDS(path), silent = TRUE)
    if (!inherits(m, "try-error")) return(m)
  }
  m <- rstan::stan_model(model_code = code, save_dso = TRUE)
  saveRDS(m, path)
  m
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## data ---------------------------------------------------------------
cases <- c("cum2", "probit_sratio", "disc_theta", "mu_cum_sratio",
           "mu_stz", "stz_flex", "cs_sratio_acat", "gr_cratio_acat",
           "gr_mu", "thres5", "cum3", "hurdle2", "mu2z", "acat_probit",
           "equi_flex", "mu_disc", "hurdle_mu", "mu_flex_stz",
           "hu_gr_logit", "hu_gr_probit", "hu_gr_disc", "hu_gr_equi",
           "hu_gr_stz", "hu_cs_probit", "hu_cs_disc", "hu_cs_logit",
           "mix_hu_gr", "mix_hu_cs", "acat_probit_1", "cum_cloglog_1",
           "r_three_mixlink", "r_disc2x", "r_theta3", "r_gr_uneq",
           "r_gr_uneq_mu", "r_weights", "r_mu2cs", "r_hu_mix_disc",
           "r_offset", "r_mu_three", "r_gr_hu_uneq", "r_theta3b")
ci <- match(case, cases)
if (is.na(ci)) stop("unknown case ", case)
set.seed(20261005 + ci)
n <- 400
x <- rnorm(n)
z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
if (case == "cum3") {
  # three classes for three components: on two-class data the third
  # component collapses a category, and the cumulative density then
  # meets its own cancellation (dev/ordmix-findings.md)
  cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
  lat <- c(1.5, -0.8, 0)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
}
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
yh <- ifelse(runif(n) < 0.2, 0L, y)
w <- runif(n, 0.5, 2)
if (case %in% c("r_three_mixlink", "r_theta3", "r_mu_three", "r_theta3b")) {
  cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
  lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
  y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
}
# group b observes only categories 1..3: per-group threshold counts
if (case %in% c("r_gr_uneq", "r_gr_uneq_mu", "r_gr_hu_uneq")) {
  y <- ifelse(g == "b", pmin(y, 3L), y)
  yh <- ifelse(runif(n) < 0.2, 0L, y)
}
d <- data.frame(y, yh, x, z, g, w)

spec <- switch(case,
  cum2 = list(f = y ~ x, fam = function(p) p$mixture(p$cumulative(),
                                                     p$cumulative())),
  probit_sratio = list(f = y ~ x, fam = function(p) {
    p$mixture(p$cumulative("probit"), p$sratio())
  }),
  disc_theta = list(f = quote(bf(y ~ x, disc1 ~ 0 + z, theta1 ~ z)),
                    fam = function(p) p$mixture(p$cumulative(),
                                                p$cumulative())),
  mu_cum_sratio = list(f = y ~ x, fam = function(p) {
    p$mixture(p$cumulative(), p$sratio(), order = "mu")
  }),
  mu_stz = list(f = y ~ x, fam = function(p) {
    p$mixture(p$cumulative(threshold = "sum_to_zero"),
              p$cumulative(threshold = "sum_to_zero"), order = "mu")
  }),
  stz_flex = list(f = y ~ x, fam = function(p) {
    p$mixture(p$sratio(threshold = "sum_to_zero"), p$cumulative())
  }),
  cs_sratio_acat = list(f = y ~ cs(x), fam = function(p) {
    p$mixture(p$sratio(), p$acat())
  }),
  gr_cratio_acat = list(f = y | thres(gr = g) ~ x, fam = function(p) {
    p$mixture(p$cratio(), p$acat())
  }),
  gr_mu = list(f = y | thres(gr = g) ~ x, fam = function(p) {
    p$mixture(p$cumulative(), p$sratio(), order = "mu")
  }),
  thres5 = list(f = y | thres(4) ~ x, fam = function(p) {
    p$mixture(p$cumulative(), p$cumulative())
  }),
  cum3 = list(f = y ~ x, fam = function(p) {
    p$mixture(p$cumulative(), p$cumulative(), p$cumulative())
  }),
  hurdle2 = list(f = quote(bf(yh ~ x, hu1 ~ z)), fam = function(p) {
    p$mixture(p$hurdle_cumulative("probit"), p$hurdle_cumulative("probit"))
  }),
  mu2z = list(f = quote(bf(y ~ x, mu2 ~ z)), fam = function(p) {
    p$mixture(p$cumulative(), p$cumulative())
  }),
  acat_probit = list(f = y ~ x, fam = function(p) {
    p$mixture(p$acat("probit"), p$cumulative("cloglog"))
  }),
  equi_flex = list(f = y ~ x, patch_delta = TRUE, fam = function(p) {
    p$mixture(p$cumulative(threshold = "equidistant"),
              p$sratio(threshold = "equidistant"))
  }),
  mu_disc = list(f = quote(bf(y ~ x, disc2 ~ 0 + z)), fam = function(p) {
    p$mixture(p$cumulative(), p$cumulative(), order = "mu")
  }),
  hurdle_mu = list(f = yh ~ x, fam = function(p) {
    p$mixture(p$hurdle_cumulative(), p$hurdle_cumulative(), order = "mu")
  }),
  mu_flex_stz = list(f = y ~ x, fam = function(p) {
    p$mixture(p$acat(), p$cratio(threshold = "sum_to_zero"),
              order = "mu")
  }),
  # hurdle_cumulative() with thres(gr = ) and cs(). With cs(), or with
  # disc modeled, brms's logit program reads the top category at
  # thres[nthres + 1] (upstream brms-1), so those rows use the probit,
  # and hu_cs_logit records the failure
  hu_gr_logit = list(f = quote(bf(yh | thres(gr = g) ~ x, hu ~ z)),
                     fam = function(p) p$hurdle_cumulative()),
  hu_gr_probit = list(f = yh | thres(gr = g) ~ x,
                      fam = function(p) p$hurdle_cumulative("probit")),
  hu_gr_disc = list(f = quote(bf(yh | thres(gr = g) ~ x, disc ~ 0 + z)),
                    fam = function(p) p$hurdle_cumulative("probit")),
  hu_gr_equi = list(f = yh | thres(gr = g) ~ x, fam = function(p) {
    p$hurdle_cumulative(threshold = "equidistant")
  }),
  hu_gr_stz = list(f = yh | thres(gr = g) ~ x, fam = function(p) {
    p$hurdle_cumulative("probit", threshold = "sum_to_zero")
  }),
  hu_cs_probit = list(f = yh ~ cs(x),
                      fam = function(p) p$hurdle_cumulative("probit")),
  hu_cs_disc = list(f = quote(bf(yh ~ cs(x), disc ~ 0 + z, hu ~ z)),
                    fam = function(p) p$hurdle_cumulative("probit")),
  hu_cs_logit = list(f = yh ~ cs(x),
                     fam = function(p) p$hurdle_cumulative()),
  mix_hu_gr = list(f = yh | thres(gr = g) ~ x, fam = function(p) {
    p$mixture(p$hurdle_cumulative("probit"), p$hurdle_cumulative("probit"))
  }),
  mix_hu_cs = list(f = yh ~ cs(x), fam = function(p) {
    p$mixture(p$hurdle_cumulative("probit"), p$hurdle_cumulative("probit"),
              order = "mu")
  }),
  # the two components of acat_probit alone, to see which one brms's
  # NaN gradient at the optimum comes from
  acat_probit_1 = list(f = y ~ x, fam = function(p) p$acat("probit")),
  cum_cloglog_1 = list(f = y ~ x, fam = function(p) p$cumulative("cloglog")),
  r_three_mixlink = list(f = y ~ x, fam = function(p) {
    p$mixture(p$cumulative("probit"), p$sratio("cloglog"), p$acat())
  }),
  r_disc2x = list(f = quote(bf(y ~ x, disc2 ~ z)), fam = function(p) {
    p$mixture(p$cumulative(), p$cratio())
  }),
  r_theta3 = list(f = quote(bf(y ~ x, theta1 ~ z, theta2 ~ 1)),
                  fam = function(p) {
    p$mixture(p$cumulative(), p$cumulative(), p$sratio())
  }),
  r_gr_uneq = list(f = y | thres(gr = g) ~ x, fam = function(p) {
    p$mixture(p$cumulative("probit"), p$cumulative())
  }),
  r_gr_uneq_mu = list(f = y | thres(gr = g) ~ x, fam = function(p) {
    p$mixture(p$cumulative(), p$cratio(threshold = "sum_to_zero"),
              order = "mu")
  }),
  r_weights = list(f = y | weights(w) ~ x, fam = function(p) {
    p$mixture(p$cumulative(), p$sratio())
  }),
  r_mu2cs = list(f = quote(bf(y ~ x, mu2 ~ x + cs(z))), fam = function(p) {
    p$mixture(p$cumulative(), p$sratio())
  }),
  r_hu_mix_disc = list(f = quote(bf(yh ~ x, hu1 ~ z, hu2 ~ x,
                                    disc1 ~ 0 + z)),
                       fam = function(p) {
    p$mixture(p$hurdle_cumulative("probit"), p$hurdle_cumulative("probit"))
  }),
  r_offset = list(f = y ~ x + offset(z), fam = function(p) {
    p$mixture(p$cumulative(), p$acat())
  }),
  r_theta3b = list(f = quote(bf(y ~ x, theta1 ~ z, theta2 ~ x)),
                   fam = function(p) {
    p$mixture(p$cumulative(), p$sratio(), p$acat())
  }),
  r_mu_three = list(f = quote(bf(y ~ x, theta1 ~ z)), fam = function(p) {
    p$mixture(p$cumulative(), p$sratio(), p$acat(), order = "mu")
  }),
  r_gr_hu_uneq = list(f = quote(bf(yh | thres(gr = g) ~ x, hu ~ z)),
                      fam = function(p) p$hurdle_cumulative())
)
pk_frm <- list(mixture = frmtmb::mixture, cumulative = frmtmb::cumulative,
               sratio = frmtmb::sratio, cratio = frmtmb::cratio,
               acat = frmtmb::acat,
               hurdle_cumulative = frmtmb::hurdle_cumulative)
pk_brm <- list(mixture = brms::mixture, cumulative = brms::cumulative,
               sratio = brms::sratio, cratio = brms::cratio,
               acat = brms::acat,
               hurdle_cumulative = brms::hurdle_cumulative)
ff <- if (!inherits(spec$f, "formula")) eval(spec$f, list(bf = frmtmb::bf)) else
  frmtmb::bf(spec$f)
fb <- if (!inherits(spec$f, "formula")) eval(spec$f, list(bf = brms::bf)) else
  brms::bf(spec$f)
fam_b <- suppressMessages(spec$fam(pk_brm))
cat("case", case, "arm", arm, "frmtmb", format(packageVersion("frmtmb")),
    "lib", find.package("frmtmb"), "\n")
fit <- frm(ff, family = spec$fam(pk_frm), data = d,
           control = frmtmb_control(grad_tol = 1e-8))
cat("frmtmb logLik", format(as.numeric(logLik(fit)), digits = 15), "\n")

## brms program -------------------------------------------------------
pr <- brms::get_prior(fb, data = d, family = fam_b)
pr$prior <- ""
code <- as.character(brms::stancode(fb, data = d, family = fam_b,
                                    prior = pr))
if (isTRUE(spec$patch_delta)) {
  # brms defect 4: one `delta` per component in the parameters block,
  # read as delta_mu<k> in transformed parameters. Rename the k-th
  # declaration to what the body reads.
  lines <- strsplit(code, "\n")[[1]]
  di <- grep("^  real<lower=0> delta;|^  real delta;", lines)
  for (j in seq_along(di)) {
    lines[di[j]] <- sub("delta;", paste0("delta_mu", j, ";"), lines[di[j]])
  }
  code <- paste(lines, collapse = "\n")
}
sdat <- brms::standata(fb, data = d, family = fam_b, prior = pr)
mod <- stan_mod(code)
sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))

## translation ---------------------------------------------------------
# the names declared in the Stan parameters block, in order
stan_par_names <- function(code) {
  lines <- sub("//.*$", "", strsplit(code, "\n", fixed = TRUE)[[1]])
  start <- grep("^parameters[[:space:]]*[{]", lines)
  out <- character()
  depth <- 0L
  for (i in seq(start[[1]], length(lines))) {
    ln <- lines[[i]]
    depth <- depth + lengths(regmatches(ln, gregexpr("[{]", ln))) -
      lengths(regmatches(ln, gregexpr("[}]", ln)))
    m <- regmatches(ln, regexpr("([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*;[[:space:]]*$", ln))
    if (length(m)) out <- c(out, trimws(sub(";.*$", "", m)))
    if (i > start[[1]] && depth <= 0L) break
  }
  out
}
fam <- fit$spec$responses[[1]]$family
mx <- fam$mix$ord
is_mix <- !is.null(mx)
K <- mx$K
rn <- fit$spec$responses[[1]]$resp_name
lpk <- function(dp) fit$frame$linpreds[[frmtmb:::linpred_key(rn, dp)]]
coef_of <- function(est, dp) {
  lp <- lpk(dp)
  if (!is.null(lp$constant)) return(NULL)
  v <- est[[lp$par]][lp$idx]
  names(v) <- colnames(lp$X)[seq_along(v)]
  v
}
# brms's suffix of one dpar: none for a plain family's mu
dsx <- function(dp) if (dp == "mu" && !is_mix) "" else paste0("_", dp)
# brms's X for one dpar, with its column names
bX <- function(dp) {
  X <- sdat[[paste0("X", dsx(dp))]]
  if (is.null(X)) return(NULL)
  cn <- colnames(X)
  X <- as.matrix(X)
  colnames(X) <- cn
  X
}
shift_of <- function(est, dp) {
  if (is.null(sdat[[paste0("Kc", dsx(dp))]])) return(0)
  X <- bX(dp)
  b <- coef_of(est, dp)
  cn <- setdiff(colnames(X), "Intercept")
  sum(colMeans(X)[cn] * unname(b[cn]))
}
comp_tau <- function(est, k) {
  if (!is_mix) {
    return(frmtmb:::ord_threshold_values(fam, est[["tau_raw"]]))
  }
  raw <- est[[mx$tau_names[k]]]
  frmtmb:::ord_threshold_values(mx$views[[k]], raw)
}
nth <- fam[["thres"]][["nthres"]] %||% length(comp_tau(fit$estimates, 1L))
lay <- frmtmb:::thres_layout(nth)
slice <- function(tau, gi) {
  if (length(nth) == 1L) tau else tau[lay$start[gi]:lay$end[gi]]
}
# component number and threshold group of an Intercept-like name
kg_of <- function(nm, pre) {
  rest <- sub(paste0("^", pre), "", nm)
  k <- if (grepl("^_mu[0-9]+", rest)) {
    as.integer(sub("^_mu([0-9]+).*$", "\\1", rest))
  } else 1L
  rest <- sub("^_mu[0-9]+", "", rest)
  gi <- if (grepl("^_[0-9]+$", rest)) as.integer(sub("^_", "", rest)) else 1L
  list(k = k, gi = gi, dp = if (is_mix) paste0("mu", k) else "mu")
}
translate <- function(est) {
  need <- stan_par_names(code)
  out <- list()
  for (nm in need) {
    if (grepl("^b(_(mu|disc|hu|theta)[0-9]*)?$", nm)) {
      dp <- if (nm == "b") "mu" else sub("^b_", "", nm)
      b <- coef_of(est, dp)
      cn <- colnames(bX(dp))
      cn <- setdiff(cn, if (!is.null(sdat[[paste0("Kc", dsx(dp))]])) {
        "Intercept"
      })
      cn[cn == "Intercept"] <- "(Intercept)"
      out[[nm]] <- array(unname(b[cn]), length(cn))
    } else if (grepl("^Intercept_(disc|hu|theta)[0-9]*$", nm)) {
      dp <- sub("^Intercept_", "", nm)
      b <- coef_of(est, dp)
      out[[nm]] <- unname(b[["(Intercept)"]]) + {
        X <- bX(dp)
        cn <- setdiff(colnames(X), "Intercept")
        sum(colMeans(X)[cn] * unname(b[cn]))
      }
    } else if (grepl("^Intercept(_mu[0-9]+)?(_[0-9]+)?$", nm)) {
      a <- kg_of(nm, "Intercept")
      out[[nm]] <- as.numeric(slice(comp_tau(est, a$k), a$gi) -
                                shift_of(est, a$dp))
    } else if (grepl("^first_Intercept(_mu[0-9]+)?(_[0-9]+)?$", nm)) {
      a <- kg_of(nm, "first_Intercept")
      out[[nm]] <- slice(comp_tau(est, a$k), a$gi)[1] - shift_of(est, a$dp)
    } else if (grepl("^delta(_mu[0-9]+)?(_[0-9]+)?$", nm)) {
      a <- kg_of(nm, "delta")
      tk <- slice(comp_tau(est, a$k), a$gi)
      out[[nm]] <- tk[2] - tk[1]
    } else if (grepl("^fixed_Intercept(_[0-9]+)?$", nm)) {
      gi <- if (grepl("_[0-9]+$", nm)) as.integer(sub("^.*_", "", nm)) else 1L
      # the shared vector itself: a flexible component reads it as is
      kf <- which(!vapply(mx$comps, function(cp) {
        identical(cp$threshold, "sum_to_zero")
      }, NA))[1]
      out[[nm]] <- as.numeric(slice(comp_tau(est, if (is.na(kf)) 1L else kf),
                                    gi))
    } else if (grepl("^bcs(_mu[0-9]+)?$", nm)) {
      lp <- lpk(if (nm == "bcs") "mu" else sub("^bcs_", "", nm))
      v <- unlist(lapply(lp$cs, function(ct) est[[ct$par]]))
      out[[nm]] <- matrix(v, nrow = length(lp$cs), byrow = TRUE)
    } else if (identical(nm, "theta")) {
      eta <- vapply(seq_len(K), function(k) {
        if (k == fam$mix$ref) return(0)
        unname(coef_of(est, paste0("theta", k))[["(Intercept)"]])
      }, 0)
      out[[nm]] <- exp(eta - max(eta)) / sum(exp(eta - max(eta)))
    } else if (grepl("^hu[0-9]*$", nm)) {
      out[[nm]] <- plogis(unname(coef_of(est, nm)[["(Intercept)"]]))
    } else {
      stop("no translation rule for Stan parameter ", nm)
    }
  }
  out
}
# brms keeps its dirichlet(1) on a simplex theta under a flat prior, and
# its density on K shares is the constant Gamma(K): 1 for two components
lconst <- if ("theta" %in% stan_par_names(code)) lgamma(K) else 0
cat("brms constant added to frmtmb:", lconst, "\n")
lp_at <- function(p) {
  est <- fit$obj$env$parList(p)
  up <- rstan::unconstrain_pars(sf, translate(est))
  list(brms = rstan::log_prob(sf, up, adjust_transform = FALSE,
                              gradient = FALSE),
       frm = -fit$obj$fn(p) + lconst, up = up)
}
r0 <- lp_at(fit$opt$par)
g0 <- rstan::grad_log_prob(sf, r0$up, adjust_transform = FALSE)
ulp <- function(a, b) abs(a - b) / (.Machine$double.eps * max(abs(a), abs(b)))
cat(sprintf("LP case=%s point=opt brms=%.15g frm=%.15g reldiff=%.3g ulp=%.1f\n",
            case, r0$brms, r0$frm, abs(r0$brms - r0$frm) / abs(r0$frm),
            ulp(r0$brms, r0$frm)))
cat(sprintf("GRAD case=%s max|grad brms at frmtmb opt|=%.3g npar=%d\n",
            case, max(abs(g0)), length(g0)))
set.seed(99 + ci)
for (j in 1:3) {
  p <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.3)
  r <- lp_at(p)
  cat(sprintf("LP case=%s point=perturb%d brms=%.15g frm=%.15g reldiff=%.3g ulp=%.1f\n",
              case, j, r$brms, r$frm, abs(r$brms - r$frm) / abs(r$frm),
              ulp(r$brms, r$frm)))
}
gf <- fit$obj$gr(fit$opt$par)
cat(sprintf("GRADFRM case=%s max|grad frmtmb at its opt|=%.3g\n", case,
            max(abs(gf))))
