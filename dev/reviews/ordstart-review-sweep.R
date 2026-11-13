# Review of lane ordstart: fit a set of ordinal-mixture cases on one arm
# (base = rellib-r20, 0.81.0; lane = wt-ordstart-lib) and record, per
# fit, the logLik, the warnings, the degeneracy row per component and the
# outer parameters. On the lane arm the candidates of mixture_starts()
# are recorded too (objective and tier from mixture_starts_pick(), and
# the degeneracy row and max |par| at each from mixture_ord_start_class()).
#   Rscript ordstart-review-sweep.R arm set out_prefix [seeds]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
set <- args[2]
out <- args[3]
seeds_arg <- if (length(args) > 3) eval(parse(text = args[4])) else NULL
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r20",
  lane = c("C:/Users/adf44/source/r/wt-ordstart-lib",
           "C:/Users/adf44/source/r/rellib-r20"),
  lane0 = c("C:/Users/adf44/source/r/wt-ordstart-lib",
            "C:/Users/adf44/source/r/rellib-r20"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "set", set, "frmtmb", format(packageVersion("frmtmb")),
    "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
lane <- arm %in% c("lane", "lane0")
rec <- new.env()
if (lane) {
  rec$cands <- list()
  rec$pick <- list()
  suppressMessages(trace("mixture_ord_start_class", where = ns, print = FALSE,
    exit = quote({
      assign("cands", c(get("cands", envir = .GlobalEnv$rec),
        list(list(obj = fit$opt$objective,
                  maxpar = max(abs(fit$opt$par)),
                  sharpen = dg$sharpen, near = dg$near,
                  collapsed = dg$collapsed))), envir = .GlobalEnv$rec)
    })))
  suppressMessages(trace("mixture_starts_pick", where = ns, print = FALSE,
    exit = quote({
      assign("pick", list(obj_v = obj_v, tier = tier,
                          keep = returnValue()), envir = .GlobalEnv$rec)
    })))
}
rec <- rec
assign("rec", rec, envir = .GlobalEnv)

# two latent classes; class k: lat = b_k * x + c_k + noise; K categories
# cut at the quantiles of the pooled latent values (each category used)
gen <- function(seed, n = 300, K = 4, b1 = 1.5, b2 = -0.8, c1 = 1, c2 = -1,
                p = 0.4, noise = "logit", nx = 0, single = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, p)
  if (single) cls[] <- 1L
  e <- switch(noise, logit = rlogis(n), probit = rnorm(n),
              cloglog = log(rexp(n)))
  lat <- ifelse(cls == 1, b1 * d$x + c1, b2 * d$x + c2) + e
  if (nx > 0) {
    X <- matrix(rnorm(n * nx), n, nx)
    colnames(X) <- paste0("w", seq_len(nx))
    lat <- lat + drop(X %*% rep(0.2, nx))
    d <- cbind(d, X)
  }
  cuts <- stats::quantile(lat, seq_len(K - 1) / K)
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}

cases <- list()
add <- function(name, data, formula, family, control = NULL, prior = NULL) {
  cases[[name]] <<- list(data = data, formula = formula, family = family,
                         control = control, prior = prior)
}
if (set == "steep") {
  # well-separated classes with steep slopes: a genuinely steep proper
  # maximum, the false-positive risk of the step classifier
  ss <- seeds_arg %||% 1:2
  for (n in c(300, 3000, 20000)) for (s in c(3, 5, 8)) for (K in c(3, 6)) {
    for (sd in if (n == 20000) ss[1] else ss) {
      add(sprintf("steep_n%d_s%g_K%d_seed%d", n, s, K, sd),
          quote(gen(SEED, n = N, K = KK, b1 = S, b2 = -S, c1 = 2, c2 = -2)),
          bf(y ~ x + z), quote(mixture(cumulative(), cumulative())))
      cases[[length(cases)]]$data <- eval(substitute(
        quote(gen(SEED, n = N, K = KK, b1 = S, b2 = -S, c1 = 2, c2 = -2)),
        list(SEED = sd, N = n, KK = K, S = s)))
    }
  }
  # one class, a strong predictor, fitted with two components
  for (n in c(300, 3000)) for (s in c(3, 8)) for (sd in ss) {
    add(sprintf("single_n%d_s%g_seed%d", n, s, sd), NULL, bf(y ~ x + z),
        quote(mixture(cumulative(), cumulative())))
    cases[[length(cases)]]$data <- eval(substitute(
      quote(gen(SEED, n = N, K = 4, b1 = S, c1 = 0, single = TRUE)),
      list(SEED = sd, N = n, S = s)))
  }
} else if (set == "sweep") {
  ss <- seeds_arg %||% 1:8
  defs <- list(
    cc_probit = list(K = 5, noise = "probit", f = bf(y ~ x + z),
                     fam = quote(mixture(cumulative("probit"),
                                         cumulative("probit")))),
    sr_cr_cloglog = list(K = 5, noise = "cloglog", f = bf(y ~ x + z),
                         fam = quote(mixture(sratio("cloglog"), cratio()))),
    acat_acat_K3 = list(K = 3, f = bf(y ~ x + z),
                        fam = quote(mixture(acat(), acat()))),
    cum3_K5 = list(K = 5, f = bf(y ~ x + z),
                   fam = quote(mixture(cumulative(), cumulative(),
                                       cumulative()))),
    theta_formula = list(K = 4, f = bf(y ~ x + z, theta2 ~ z),
                         fam = quote(mixture(cumulative(), cumulative()))),
    dpar_formula = list(K = 4, f = bf(y ~ 1, mu1 ~ x, mu2 ~ x + z),
                        fam = quote(mixture(cumulative(), sratio()))),
    cs_sr_acat = list(K = 4, f = bf(y ~ cs(x) + z),
                      fam = quote(mixture(sratio(), acat()))),
    disc_free = list(K = 5, f = bf(y ~ x + z, disc2 ~ 1),
                     fam = quote(mixture(cumulative(), cumulative()))),
    order_mu = list(K = 4, f = bf(y ~ x + z),
                    fam = quote(mixture(cumulative(), sratio(),
                                        order = "mu"))),
    prior_b = list(K = 4, f = bf(y ~ x + z),
                   fam = quote(mixture(cumulative(), cumulative())),
                   prior = quote(c(set_prior("normal(0, 5)", class = "b",
                                             dpar = "mu1"),
                                   set_prior("normal(0, 5)", class = "b",
                                             dpar = "mu2")))),
    cr_cr_probit_K6 = list(K = 6, noise = "probit", f = bf(y ~ x + z),
                           fam = quote(mixture(cratio("probit"),
                                               cratio("probit")))),
    acat_cumcll_K7 = list(K = 7, f = bf(y ~ x + z),
                          fam = quote(mixture(acat(),
                                              cumulative("cloglog")))),
    cc_logit_K7 = list(K = 7, f = bf(y ~ x + z),
                       fam = quote(mixture(cumulative(), cumulative()))))
  for (nm in names(defs)) for (sd in ss) {
    df <- defs[[nm]]
    add(sprintf("%s_seed%d", nm, sd), NULL, df$f, df$fam,
        prior = df$prior)
    cases[[length(cases)]]$data <- eval(substitute(
      quote(gen(SEED, n = 300, K = KK, noise = NZ)),
      list(SEED = sd, KK = df$K, NZ = df$noise %||% "logit")))
  }
} else if (set == "wide") {
  ss <- seeds_arg %||% 1:3
  for (sd in ss) {
    add(sprintf("wide20_K7_n1000_seed%d", sd), NULL,
        bf(as.formula(paste("y ~ x + z +", paste0("w", 1:20,
                                                  collapse = " + ")))),
        quote(mixture(cumulative(), cumulative())))
    cases[[length(cases)]]$data <- eval(substitute(
      quote(gen(SEED, n = 1000, K = 7, nx = 20)), list(SEED = sd)))
  }
}

ctl_for <- function(cs) {
  if (arm == "lane0") return(frmtmb_control(starts = 0))
  cs$control
}
rows <- list()
pars <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- eval(cs$data)
  w <- character()
  rec$cands <- list()
  rec$pick <- NULL
  t0 <- proc.time()[["elapsed"]]
  a <- list(cs$formula, family = eval(cs$family), data = d)
  if (!is.null(ctl_for(cs))) a$control <- ctl_for(cs)
  if (!is.null(cs$prior)) a$prior <- eval(cs$prior)
  f <- tryCatch(withCallingHandlers(do.call(frm, a), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("CASE %s ERROR %s\n", nm, conditionMessage(f)))
    rows[[nm]] <- data.frame(case = nm, ll = NA, code = NA, elapsed = el,
                             degen = NA, stepw = NA, ms = NA, nwarn = NA,
                             maxsharp = NA, warn = "",
                             err = conditionMessage(f))
    next
  }
  dg <- tryCatch(frmtmb:::mixture_ord_degeneracy(f, names(f$spec$responses)[1]),
                 error = function(e) NULL)
  ms <- f$opt$mixture_starts
  pk <- rec$pick
  cand_s <- if (length(rec$cands)) {
    paste(vapply(seq_along(rec$cands), function(i) {
      cc <- rec$cands[[i]]
      sprintf("%.4f/%s/%.3g/%.3g", -cc$obj,
              if (!is.null(pk)) pk$tier[i] else "?",
              max(cc$sharpen, na.rm = TRUE), cc$maxpar)
    }, ""), collapse = " ")
  } else ""
  rows[[nm]] <- data.frame(
    case = nm, ll = as.numeric(logLik(f)), code = f$opt$convergence,
    elapsed = el,
    degen = any(grepl("ended at a degenerate boundary", w)),
    stepw = any(grepl("jittered starts reached a higher", w)),
    ms = if (is.null(ms)) "" else paste(names(ms), signif(ms, 4), sep = "=",
                                         collapse = ","),
    nwarn = length(w),
    maxsharp = if (is.null(dg)) NA else max(dg$sharpen, na.rm = TRUE),
    warn = paste(gsub("[\t\n]", " ", w), collapse = " || "),
    err = "")
  pars[[nm]] <- f$opt$par
  cat(sprintf("CASE %s ll %.4f code %d el %.2f degen %s stepw %s ms %s keep %s cands %s\n",
              nm, as.numeric(logLik(f)), f$opt$convergence, el,
              rows[[nm]]$degen, rows[[nm]]$stepw, rows[[nm]]$ms,
              if (is.null(pk)) "-" else pk$keep, cand_s))
}
R <- do.call(rbind, rows)
utils::write.table(R, paste0(out, ".tsv"), sep = "\t", row.names = FALSE,
                   quote = TRUE)
saveRDS(pars, paste0(out, "-pars.rds"))
cat("DONE", nrow(R), "\n")
