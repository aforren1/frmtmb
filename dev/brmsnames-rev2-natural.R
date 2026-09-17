## Reviewer recheck, MAJOR 3: the natural-scale dpar columns. Every reader
## that hands a draw back to the model must see the SAME link value it saw
## at base, so each reader's output must be identical (or equal to the
## round-trip ulp) between arms on the same stanfit. Run once per arm, then
## compare, pairing by POSITION, never by name.
##   Rscript dev/brmsnames-rev2-natural.R base
##   Rscript dev/brmsnames-rev2-natural.R lane
##   Rscript dev/brmsnames-rev2-natural.R compare
## Data seed 41, sampler seed 9, chains 2, iter 200.
arm <- commandArgs(trailingOnly = TRUE)[1L]
only <- commandArgs(trailingOnly = TRUE)[-1L]
libs <- list(
  base = c("C:/Users/adf44/source/r/rellib-r3"),
  lane = c("C:/Users/adf44/source/r/brmsnames-lib",
           "C:/Users/adf44/source/r/rellib-r3"),
  compare = c("C:/Users/adf44/source/r/rellib-r3"))
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
outf <- function(a, nm) {
  sprintf("dev/stan-cache/brmsnames-rev2-natural-%s-%s.rds", a, nm)
}
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}

if (arm != "compare") {
  q(library(frmtmb)); q(library(frmtmb.sample))
  cat("arm", arm, "frmtmb from", dirname(system.file(package = "frmtmb")),
      " sample from", dirname(system.file(package = "frmtmb.sample")), "\n")
  set.seed(41)
  n <- 160; G <- 8
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(G), length.out = n)))
  u <- rnorm(G, 0, 0.7)
  d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(n, 0, 4)          # sigma 4
  d$yc <- rnbinom(n, mu = exp(1.5 + 0.3 * d$x), size = 3)  # shape 3
  lam <- exp(1 + 0.4 * d$x)
  d$yz <- ifelse(runif(n) < 0.3, 0L, rpois(n, lam))         # zi 0.3
  d$yt <- 1 + 0.5 * d$x + 2 * rt(n, df = 4)                 # nu 4
  mb <- plogis(-0.3 + 0.5 * d$x)
  d$yb <- rbeta(n, mb * 12, (1 - mb) * 12)                  # phi 12
  d$yh <- ifelse(runif(n) < 0.35, 0,
                 rgamma(n, shape = 2.5, rate = 2.5 / exp(0.5 + 0.3 * d$x)))
  d$ya <- 2 + 0.7 * d$x + rnorm(n, 0, 3)
  d$y2 <- -1 + 0.2 * d$x + rnorm(n, 0, 0.3)
  pri <- set_prior("normal(0, 1)", class = "b")
  ms <- list(
    gauss = list(f = bf(y ~ x + (1 | g)), fam = gaussian(), dp = "sigma",
                 nat = "sigma > 3"),
    negbin = list(f = bf(yc ~ x), fam = negbinomial(), dp = "shape",
                  nat = "shape > 2"),
    zip = list(f = bf(yz ~ x), fam = zero_inflated_poisson(), dp = "zi",
               nat = "zi > 0.2"),
    student = list(f = bf(yt ~ x), fam = student(), dp = "nu",
                   nat = "nu > 3"),
    beta = list(f = bf(yb ~ x), fam = Beta(), dp = "phi", nat = "phi > 8"),
    hurdle = list(f = bf(yh ~ x), fam = hurdle_gamma(), dp = "hu",
                  nat = "hu > 0.3"),
    mv = list(f = mvbf(bf(ya ~ x), bf(y2 ~ x), rescor = FALSE),
              fam = gaussian(), dp = "sigma", nat = "sigma_ya > 2",
              resp = "ya"),
    gauss_prior = list(f = bf(y ~ x + (1 | g)), fam = gaussian(),
                       dp = "sigma", nat = "sigma > 3", prior = pri),
    gauss_laplace = list(f = bf(y ~ x + (1 | g)), fam = gaussian(),
                         dp = "sigma", nat = "sigma > 3", laplace = TRUE),
    gauss_ncp = list(f = bf(y ~ x + (1 | g)), fam = gaussian(),
                     dp = "sigma", nat = "sigma > 3", ncp = TRUE)
  )
  if (length(only)) ms <- ms[only]
  for (nm in names(ms)) {
    M <- ms[[nm]]
    fit <- q(frm(M$f, family = M$fam, data = d))
    args <- list(fit, chains = 2, iter = 200, refresh = 0, seed = 9)
    if (!is.null(M$prior)) args$prior <- M$prior
    if (isTRUE(M$laplace)) args$laplace <- TRUE
    if (isTRUE(M$ncp)) args$reparameterize <- TRUE
    # the mode-anchored inits jitter with R's RNG, so the sampler seed
    # alone does not fix the draws; reset per model, independent of what
    # the readers of the previous model consumed
    set.seed(sum(utf8ToInt(nm)))
    ds <- q(do.call(frm_sample, args))
    rs <- M$resp
    ss <- function(expr) { set.seed(1); q(expr) }
    lane <- arm == "lane"
    r <- list(
      draws = unname(ds$draws), names = colnames(ds$draws),
      epred = try1(q(posterior_epred(ds, resp = rs))),
      linpred = try1(q(posterior_linpred(ds, resp = rs))),
      linpred_dpar = try1(q(posterior_linpred(ds, dpar = M$dp, resp = rs))),
      epred_dpar = try1(q(posterior_epred(ds, dpar = M$dp, resp = rs))),
      predict = try1(ss(posterior_predict(ds, resp = rs))),
      loglik = try1(q(log_lik(ds, resp = rs))),
      loo = try1(q(loo(ds, resp = rs))$estimates),
      waic = try1(q(waic(ds, resp = rs))$estimates),
      r2 = try1(q(bayes_R2(ds, resp = rs))),
      perr = try1(ss(predictive_error(ds, resp = rs))),
      pint = try1(ss(predictive_interval(ds, resp = rs))),
      ce = try1(q(conditional_effects(ds, "x", resp = rs))[[1]][
        , c("estimate__", "se__", "lower__", "upper__")]),
      ce_pred = try1(ss(conditional_effects(ds, "x", resp = rs,
                                            method = "posterior_predict"))[[1]][
        , c("estimate__", "se__", "lower__", "upper__")]),
      ce_nonrobust = try1(q(do.call(conditional_effects, c(
        list(ds, "x", resp = rs), if (lane) list(robust = FALSE))))[[1]][
        , c("estimate__", "se__", "lower__", "upper__")]),
      ce_dpar_nonrobust = try1(q(do.call(conditional_effects, c(
        list(ds, "x", dpar = M$dp, resp = rs),
        if (lane) list(robust = FALSE))))[[1]][
        , c("estimate__", "se__", "lower__", "upper__")]),
      ce_dpar = try1(q(conditional_effects(ds, "x", dpar = M$dp,
                                           resp = rs))[[1]][
        , c("estimate__", "se__", "lower__", "upper__")]),
      ppc = try1(ss(pp_check(ds, type = "dens_overlay", ndraws = 5,
                             resp = rs))$data),
      ppc_stat = try1(ss(pp_check(ds, type = "stat", ndraws = NULL,
                                  resp = rs))$data),
      vc = try1(q(VarCorr(ds))),
      hyp_x = try1(q(hypothesis(ds, if (!is.null(rs)) "ya_x > 0" else
        "x > 0")$hypothesis)),
      hyp_x_point = try1(q(hypothesis(ds, if (!is.null(rs)) "ya_x = 0" else
        "x = 0")$hypothesis)),
      hyp_nat = try1(q(if (lane) hypothesis(ds, M$nat, class = NULL)$hypothesis
                       else hypothesis(ds, M$nat))),
      summary = try1(q(summary(ds))),
      print = try1(capture.output(print(ds)))
    )
    if (nm %in% c("gauss", "gauss_prior")) {
      r$hyp_mix <- try1(q(if (lane)
        hypothesis(ds, "sigma - sd_g__Intercept > 0", class = NULL)$hypothesis
        else hypothesis(ds, "sigma - sd_g__Intercept > 0")))
      r$ranef <- try1(q(ranef(ds, summary = FALSE)))
      r$cl <- try1(q(check_laplace(fit, chains = 2, iter = 200,
                                   refresh = 0, seed = 11)))
    }
    saveRDS(r, outf(arm, nm))
    cat(nm, "done; natural-looking columns:",
        grep("sigma|shape|zi|nu|phi|hu", colnames(ds$draws), value = TRUE),
        "\n")
  }
} else {
  num <- function(v) {
    if (is.character(v) && length(v) == 1L && grepl("^ERROR", v)) return(v)
    if (is.list(v) && !is.data.frame(v)) v <- unlist(lapply(v, num))
    if (is.data.frame(v)) v <- unlist(Filter(is.numeric, v))
    as.numeric(unlist(v))
  }
  rel <- function(a, b) {
    if (length(a) != length(b)) return(sprintf("LENGTH %d vs %d",
                                               length(a), length(b)))
    ok <- is.finite(a) & is.finite(b)
    if (any(is.finite(a) != is.finite(b))) return("FINITENESS DIFFERS")
    if (!any(ok)) return("no finite")
    m <- max(abs(a[ok] - b[ok]) / pmax(abs(a[ok]), abs(b[ok]), 1e-300))
    sprintf("%.3g", m)
  }
  fl <- list.files("dev/stan-cache", "^brmsnames-rev2-natural-base-")
  for (f in fl) {
    nm <- sub("^brmsnames-rev2-natural-base-(.*)[.]rds$", "\\1", f)
    B <- readRDS(file.path("dev/stan-cache", f))
    L <- readRDS(outf("lane", nm))
    cat("\n== model", nm, "==\n")
    dif <- which(B$names != L$names)
    cat("  column renames (base -> lane):",
        paste(B$names[dif], L$names[dif], sep = "->")[
          !grepl("^r_|^b\\[", B$names[dif])], "\n")
    same <- vapply(seq_len(ncol(B$draws)), function(j)
      identical(B$draws[, j], L$draws[, j]), TRUE)
    cat("  draws columns identical:", sum(same), "of", length(same),
        "; differing:", L$names[!same], "\n")
    for (j in which(!same)) {
      cat(sprintf("    %s: base mean %.4g, lane mean %.4g, lane vs exp(base) %s, vs plogis(base) %s\n",
                  L$names[j], mean(B$draws[, j]), mean(L$draws[, j]),
                  rel(L$draws[, j], exp(B$draws[, j])),
                  rel(L$draws[, j], plogis(B$draws[, j]))))
    }
    for (k in setdiff(names(L), c("draws", "names"))) {
      b <- num(B[[k]]); l <- num(L[[k]])
      if (is.character(b) || is.character(l)) {
        cat(sprintf("  %-13s base: %s | lane: %s\n", k,
                    substr(paste(b, collapse = ""), 1, 90),
                    substr(paste(l, collapse = ""), 1, 90)))
        next
      }
      cat(sprintf("  %-13s identical %-5s max rel diff %s\n", k,
                  identical(b, l), rel(b, l)))
    }
  }
}
