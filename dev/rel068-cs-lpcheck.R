# The log density of cs() on cumulative() (lifted at 0.68.0) against
# brms 2.23.0's compiled Stan program, at matched parameter values. The
# translation and the comparison below the data block are
# dev/ordmix-lpcheck.R's, copied unchanged (that script's lines 169
# to the end).
#
# Usage: Rscript dev/rel068-cs-lpcheck.R <case>
#
# For each case: simulate data (seed 20261006 + case number), fit
# frmtmb by maximum likelihood, build brms's program with every prior
# flat (so its log density with adjust_transform = FALSE is the log
# likelihood), and evaluate both at the optimum and at three perturbed
# parameter vectors (seed 99 + case). Each point's frmtmb value is
# -obj$fn(p); the brms value is rstan::log_prob() at the translated
# parameters. Also reports brms's gradient at frmtmb's optimum, which
# vanishes only when the translation is right and the optimum is
# shared. A perturbed point whose offsets cross two thresholds on some
# row has a NaN density on both sides; the line says so.
args <- commandArgs(TRUE)
case <- if (length(args)) args[1] else "cum_cs_logit"
arm <- "release"
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
  library(rstan)
})
cache_dir <- "C:/Users/adf44/source/r/frmtmb-wt-release/dev/stan-cache"
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
cases <- c("cum_cs_logit", "cum_cs_probit_disc", "cum_cs_factor",
           "cum_cs_equi", "mix_cum_cs", "cum_cs_probit2")
ci <- match(case, cases)
if (is.na(ci)) stop("unknown case ", case)
set.seed(20261006 + ci)
n <- 400
x <- rnorm(n)
z <- rnorm(n)
g <- factor(sample(c("a", "b", "c"), n, TRUE))
# threshold-specific slopes small against the gaps between thresholds,
# so that no row's thresholds cross at the optimum
u <- rlogis(n) + 0.5 * x + 0.4 * z
y <- 1L + (u > -1 + 0.15 * x) + (u > 0.3) + (u > 1.5 - 0.15 * x +
                                                0.2 * (g == "b"))
yh <- y
d <- data.frame(y, yh, x, z, g)

spec <- switch(case,
  cum_cs_logit = list(f = y ~ cs(x), fam = function(p) p$cumulative()),
  cum_cs_probit_disc = list(f = quote(bf(y ~ z + cs(x), disc ~ 0 + z)),
                            fam = function(p) p$cumulative("probit")),
  cum_cs_factor = list(f = y ~ x + z + cs(g),
                       fam = function(p) p$cumulative()),
  cum_cs_equi = list(f = y ~ z + cs(x), fam = function(p) {
    p$cumulative(threshold = "equidistant")
  }),
  mix_cum_cs = list(f = y ~ cs(x), fam = function(p) {
    p$mixture(p$cumulative(), p$sratio())
  }),
  cum_cs_probit2 = list(f = y ~ cs(x) + cs(z),
                        fam = function(p) p$cumulative("probit"))
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
