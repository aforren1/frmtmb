# Item 3.6: how often the widened rp_floored() fires on ORDINARY
# right-censored designs, at the df values flexsurv users pick.
#
# Usage: Rscript phase3a-falsealarm36.R <part> <out tsv>
#   part "sim"  simulated designs, hazard scale
#   part "real" nine real datasets, all three scales
#
# Every row is one fit. `old` is the 0.7.0 check (event rows only),
# `new_only` counts fits the widening refuses that 0.7.0 did not. Both
# counts come from ONE report on the lane build: n_nonmonotone is the
# 0.7.0 quantity unchanged, and n_nonmonotone_censored is the addition.
# Seeds 20260923 upward, one per replicate, recorded in the row.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
Sys.setenv(PHASE3A_ARM = "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})
phase3a_where("frmtmb.spline")
out <- args[2L]

emit <- function(row) {
  utils::write.table(as.data.frame(row), out, sep = "\t", quote = FALSE,
                     row.names = FALSE, append = file.exists(out),
                     col.names = !file.exists(out))
}

# What rp_floored(action = "error") does on a fit, measured by calling
# it: "refuse", "warn" or "silent". Punch round 2 split the two.
gate_of <- function(f) {
  w <- FALSE
  e <- tryCatch(withCallingHandlers(rp_floored(f), warning = function(c) {
    w <<- TRUE
    invokeRestart("muffleWarning")
  }), error = function(err) "refuse")
  if (identical(e, "refuse")) "refuse" else if (w) "warn" else "silent"
}

one_fit <- function(tag, seed, d, form, df, scale, flab) {
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(form, family = royston_parmar(df = df, scale = scale), data = d),
    warning = function(cnd) {
      warns <<- c(warns, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  base <- c(list(design = tag, seed = seed, n = nrow(d),
                 cens_frac = round(mean(d$censored != 0), 3), df = df,
                 scale = scale, form = flab))
  if (inherits(f, "error")) {
    emit(c(base, list(status = "fit_error", conv = NA, n_nonmono = NA,
                      n_nonmono_cens = NA, old = NA, new = NA,
                      new_only = NA, fitend_warn_cens = NA, gate = NA, max_rise = NA,
                      min_detadx_cens = NA, elapsed = round(el, 2),
                      msg = gsub("[\t\n]", " ",
                                 substr(conditionMessage(f), 1, 100)))))
    return(invisible(NULL))
  }
  r <- rp_floored(f, action = "report")
  # the lowest derivative on any censored row, so a near miss is visible
  fam <- stats::family(f)
  fx <- frmtmb.spline:::sp_rp_fitted(f, fam)
  mdc <- if (any(fx$cens != 0)) min(fx$detadx[fx$cens != 0]) else NA
  emit(c(base, list(
    status = "ok", conv = f$opt$convergence,
    n_nonmono = r[["n_nonmonotone"]],
    n_nonmono_cens = r[["n_nonmonotone_censored"]],
    old = as.integer(r[["n_nonmonotone"]] > 0),
    new = as.integer(r[["n_nonmonotone"]] + r[["n_nonmonotone_censored"]] > 0),
    new_only = as.integer(r[["n_nonmonotone"]] == 0 &&
                            r[["n_nonmonotone_censored"]] > 0),
    fitend_warn_cens = as.integer(any(grepl("fitted survival rises",
                                            warns))),
    min_detadx_cens = signif(mdc, 4), elapsed = round(el, 2),
    gate = gate_of(f), max_rise = signif(r[["max_survival_rise"]], 3),
    msg = "")))
}

# ------------------------------------------------------------ simulated
# Four baseline shapes that cover the hazards a survival paper meets:
# decreasing (Weibull 0.8), increasing (Weibull 1.5), rise then fall
# (log-logistic 2.5) and lognormal. A treatment arm shifts the time
# scale by exp(-0.5), which is proportional hazards only for the
# Weibull, so the PH fits on the other two are misspecified the way
# real data usually are.
sim_times <- function(shape, n, trt) {
  sc <- exp(0.5 * trt)
  switch(shape,
    weib08 = stats::rweibull(n, 0.8, sc),
    weib15 = stats::rweibull(n, 1.5, sc),
    llogis = { u <- stats::runif(n); sc * (u / (1 - u))^(1 / 2.5) },
    lnorm = stats::rlnorm(n, log(sc), 1))
}
# light: administrative at the 90th percentile plus uniform dropout;
# heavy: administrative at the median plus uniform dropout
sim_design <- function(shape, n, cens, seed) {
  set.seed(seed)
  trt <- stats::rbinom(n, 1L, 0.5)
  t <- sim_times(shape, n, trt)
  tau <- unname(stats::quantile(t, if (cens == "light") 0.9 else 0.5))
  cc <- pmin(tau, stats::runif(n, 0, if (cens == "light") 3 * tau else
    1.5 * tau))
  data.frame(t = pmin(t, cc), censored = as.integer(t > cc), trt = trt)
}

if (args[1L] == "sim") {
  ph <- bf(t | cens(censored) ~ trt)
  nph <- bf(t | cens(censored) ~ trt, gamma1 ~ trt)
  rep_n <- 10L
  for (shape in c("weib08", "weib15", "llogis", "lnorm")) {
    for (n in c(200L, 1000L)) {
      for (cens in c("light", "heavy")) {
        for (i in seq_len(rep_n)) {
          seed <- 20260923L + i
          d <- sim_design(shape, n, cens, seed)
          for (df in 1:5) {
            one_fit(paste(shape, cens, sep = "_"), seed, d, ph, df,
                    "hazard", "ph")
            one_fit(paste(shape, cens, sep = "_"), seed, d, nph, df,
                    "hazard", "gamma1~x")
          }
        }
      }
    }
  }
}

# ----------------------------------------------------------------- real
if (args[1L] == "real") {
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$t <- bc$recyrs; bc$censored <- 1L - bc$censrec; bc$x <- bc$group
  sv <- function(nm) getExportedValue("survival", nm)
  lung <- stats::na.omit(sv("lung")[, c("time", "status", "sex")])
  lung <- data.frame(t = lung$time, censored = as.integer(lung$status == 1),
                     x = factor(lung$sex))
  vet <- sv("veteran")
  vet <- data.frame(t = vet$time, censored = 1L - vet$status,
                    x = factor(vet$trt))
  col <- sv("colon"); col <- col[col$etype == 2, ]
  col <- data.frame(t = col$time, censored = 1L - col$status, x = col$rx)
  rot <- sv("rotterdam")
  rot <- data.frame(t = pmax(rot$dtime, 1), censored = 1L - rot$death,
                    x = factor(rot$hormon))
  gb <- sv("gbsg")
  gb <- data.frame(t = gb$rfstime, censored = 1L - gb$status,
                   x = factor(gb$hormon))
  pb <- sv("pbc"); pb <- pb[!is.na(pb$trt), ]
  pb <- data.frame(t = pb$time, censored = as.integer(pb$status != 2),
                   x = factor(pb$trt))
  ov <- sv("ovarian")
  ov <- data.frame(t = ov$futime, censored = 1L - ov$fustat,
                   x = factor(ov$rx))
  mg <- sv("mgus2")
  mg <- data.frame(t = pmax(mg$futime, 1), censored = 1L - mg$death,
                   x = factor(mg$sex))
  sets <- list(bc = bc[, c("t", "censored", "x")], lung = lung,
               veteran = vet, colon_death = col, rotterdam_death = rot,
               gbsg_rfs = gb, pbc_death = pb, ovarian = ov,
               mgus2_death = mg)
  ph <- bf(t | cens(censored) ~ x)
  nph <- bf(t | cens(censored) ~ x, gamma1 ~ x)
  for (nm in names(sets)) {
    d <- sets[[nm]]
    for (scale in c("hazard", "odds", "normal")) {
      for (df in 1:5) {
        one_fit(nm, NA_integer_, d, ph, df, scale, "ph")
        one_fit(nm, NA_integer_, d, nph, df, scale, "gamma1~x")
      }
    }
  }
}
cat("done", args[1L], "\n")
