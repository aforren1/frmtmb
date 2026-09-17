## Output match against brms, on the SAME draws: every draws-side call
## lane brmsnames changed, answered by frmtmb.sample and by brms's own
## installed method, compared with identical().
##
##   Rscript dev/brmsnames-match.R base > dev/brmsnames-log/match-base.txt
##   Rscript dev/brmsnames-match.R lane > dev/brmsnames-log/match-lane.txt
##
## How brms is made to answer on these draws. brms's methods read a
## brmsfit's `$fit`, a stanfit, through `x$fit@sim$samples`. The shim is
## brms's own `brm(..., empty = TRUE)` object for the same formula and
## data, given a copy of the sampler's stanfit whose `sim` holds the
## post-warmup draws under brms's names. brms then runs its own compiled
## bodies (fixef.brmsfit, ranef.brmsfit, VarCorr.brmsfit, ...) on the
## same numbers.
##
## What the shim cannot supply: brms's sd_ draws. The sampler stores log
## standard deviations (theta_1); the shim's `sd_g__Intercept` is the
## standard deviation computed by frmtmb's varcorr_values() from those
## draws. So the VarCorr() rows test brms's layout and summary path on
## the same derived numbers, and not the log-to-sd transform, which is
## exp() and is checked separately below against exp(theta_1).
##
## Draws: dev/stan-cache/brmsnames-draws-<arm>.rds (dev/brmsnames-draws.R,
## data seed 9, sampler seed 20260915), and a correlated-slope model
## sampled here: data seed 11, frm_sample(chains = 2, iter = 600,
## seed = 20260916), cached as dev/stan-cache/brmsnames-slope-<arm>.rds.
arm <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/brmsnames-libs.R")
brmsnames_libs(arm)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(posterior))
q(requireNamespace("brms"))
cat("arm", arm, " frmtmb", format(packageVersion("frmtmb")),
    " frmtmb.sample", format(packageVersion("frmtmb.sample")), " from",
    dirname(system.file(package = "frmtmb.sample")), "\n")
cat("brms", format(packageVersion("brms")), " posterior",
    format(packageVersion("posterior")), "\n\n")

o <- readRDS(sprintf("dev/stan-cache/brmsnames-draws-%s.rds", arm))
slope_path <- sprintf("dev/stan-cache/brmsnames-slope-%s.rds", arm)
if (file.exists(slope_path)) {
  o2 <- readRDS(slope_path)
} else {
  set.seed(11)
  d2 <- data.frame(x = stats::rnorm(200), g = factor(rep(1:10, 20)))
  u <- cbind(stats::rnorm(10, 0, 0.7), stats::rnorm(10, 0, 0.4))
  d2$y <- stats::rnorm(200, 1 + 0.5 * d2$x + u[d2$g, 1] +
                         u[d2$g, 2] * d2$x, 1)
  f2 <- frm(frmtmb::bf(y ~ x + (1 + x | g)), family = gaussian(), data = d2)
  ds2 <- q(frm_sample(f2, chains = 2, iter = 600, refresh = 0,
                      seed = 20260916))
  o2 <- list(ds = ds2, fit = f2, data = d2)
  saveRDS(o2, slope_path)
}

## brms's names for a draws column, from the lane's labeler run on the
## fit. The base arm's draws carry other names; the shim needs brms's
## in both arms, and the VALUES are bitwise the same in both arms
## (checked below), so the base arm's shim is the lane arm's shim.
lane_labels <- function(fit) {
  lp <- c("C:/Users/adf44/source/r/brmsnames-lib", .libPaths())
  f <- gsub("\\\\", "/", tempfile(fileext = ".rds"))
  saveRDS(fit, f)
  out <- gsub("\\\\", "/", tempfile(fileext = ".rds"))
  code <- sprintf(paste0(
    ".libPaths(c(%s)); suppressMessages(library(frmtmb));",
    "fit <- readRDS('%s');",
    "lay <- varcorr_layout(fit);",
    "saveRDS(list(lab = c(brms_par_labels(fit), 'lp__'), lay = lay), '%s')"),
    paste0("'", lp, "'", collapse = ","), f, out)
  system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)))
  readRDS(out)
}

## derived sd, cor and sigma draws, from the lane's core, same reason
lane_derived <- function(ds) {
  lp <- c("C:/Users/adf44/source/r/brmsnames-lib", .libPaths())
  f <- gsub("\\\\", "/", tempfile(fileext = ".rds"))
  saveRDS(ds, f)
  out <- gsub("\\\\", "/", tempfile(fileext = ".rds"))
  code <- sprintf(paste0(
    ".libPaths(c(%s)); suppressMessages({library(frmtmb);",
    "library(frmtmb.sample)});",
    "ds <- readRDS('%s'); fit <- ds$fit; ",
    "lab <- c(brms_par_labels(fit), 'lp__'); colnames(ds$draws) <- lab;",
    "lay <- varcorr_layout(fit);",
    "per <- frmtmb.sample:::draws_varcorr_values(ds, lay);",
    "saveRDS(list(per = per, lay = lay), '%s')"),
    paste0("'", lp, "'", collapse = ","), f, out)
  system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)))
  readRDS(out)
}

## A stanfit carrying `M` (post-warmup draws, chain-major, named columns)
make_shim <- function(bfit, sf, M) {
  nc <- sf@sim$chains
  n <- nrow(M) %/% nc
  sim <- sf@sim
  sim$samples <- lapply(seq_len(nc), function(ch) {
    rows <- (ch - 1L) * n + seq_len(n)
    stats::setNames(lapply(seq_len(ncol(M)), function(j) {
      unname(M[rows, j])
    }), colnames(M))
  })
  sim$iter <- n
  sim$warmup <- 0
  sim$warmup2 <- rep(0, nc)
  sim$n_save <- rep(n, nc)
  sim$permutation <- lapply(seq_len(nc), function(i) seq_len(n))
  sim$pars_oi <- colnames(M)
  sim$dims_oi <- stats::setNames(rep(list(integer(0)), ncol(M)),
                                 colnames(M))
  sim$fnames_oi <- colnames(M)
  sim$n_flatnames <- ncol(M)
  sf@sim <- sim
  bfit$fit <- sf
  bfit
}

res <- list()
cmp <- function(label, ours, theirs) {
  a <- tryCatch(q(ours), error = function(e) e)
  b <- tryCatch(q(theirs), error = function(e) e)
  st <- if (inherits(b, "error")) {
    "brms-error"
  } else if (inherits(a, "error")) {
    "ERROR"
  } else if (identical(a, b)) {
    "identical"
  } else if (isTRUE(all.equal(a, b, check.attributes = TRUE))) {
    "all.equal"
  } else {
    "DIFFERS"
  }
  why <- if (inherits(a, "error")) conditionMessage(a) else
    if (inherits(b, "error")) conditionMessage(b) else if (st == "DIFFERS")
      paste(utils::head(all.equal(a, b), 2L), collapse = " | ") else ""
  res[[length(res) + 1L]] <<- data.frame(model = model, call = label,
                                         result = st,
                                         note = substr(gsub("\n", " ", why),
                                                       1, 90))
  cat(sprintf("%-9s %-46s %-10s %s\n", model, label, st,
              substr(gsub("\n", " ", why), 1, 70)))
}

for (model in c("intercept", "slope")) {
  obj <- if (model == "intercept") o else o2
  ds <- obj$ds
  info <- lane_labels(obj$fit)
  lab <- info$lab
  # the lane stores sigma with no formula on its natural scale, as brms
  # does; the base build stored its log, so the base arm's column is put
  # on the lane's scale before either arm is compared or derived from
  if (arm == "base") {
    j <- which(lab == "sigma")
    ds$draws[, j] <- exp(ds$draws[, j])
  }
  der <- lane_derived(ds)
  other <- if (arm == "base") "lane" else "base"
  op <- sprintf(if (model == "intercept") "dev/stan-cache/brmsnames-draws-%s.rds"
                else "dev/stan-cache/brmsnames-slope-%s.rds", other)
  if (file.exists(op)) {
    od <- readRDS(op)$ds$draws
    if (other == "base") {
      j <- which(lab == "sigma")
      od[, j] <- exp(od[, j])
    }
    cat("draw values identical() to the ", other, " arm's, sigma on one ",
        "scale: ", identical(unname(ds$draws), unname(od)), "\n", sep = "")
  }
  M <- ds$draws
  colnames(M) <- lab
  keyg <- names(der$lay$groups)[1L]
  sdm <- do.call(rbind, lapply(der$per, function(v) v[[keyg]]$sd))
  rn <- der$lay$groups[[keyg]]$rnames
  colnames(sdm) <- paste0("sd_g__", rn)
  corm <- NULL
  if (!is.null(der$per[[1L]][[keyg]]$cor)) {
    corm <- do.call(rbind, lapply(der$per, function(v) {
      C <- v[[keyg]]$cor
      C[lower.tri(C)]
    }))
    colnames(corm) <- "cor_g__Intercept__x"
  }
  sigm <- vapply(der$per, function(v) v[["residual__"]]$sd, 1)
  cat("== ", model, ": ", nrow(M), " x ", ncol(M), " draws ==\n", sep = "")
  cat("columns of this arm's draws: ",
      paste(utils::head(colnames(ds$draws), 6L), collapse = " "), "\n")
  cat("exp(theta_1) vs derived sd, max |diff|: ",
      format(max(abs(exp(M[, "theta_1"]) - sdm[, 1L]))), "\n", sep = "")

  d <- obj$data
  # brms's side is the model frmtmb fitted, formula for formula (punch
  # round 1, MAJOR 3): sigma has no formula in either, so both have the
  # natural-scale `sigma` column
  fml_simple <- if (model == "intercept") y ~ x + (1 | g) else
    y ~ x + (1 + x | g)
  b_simple <- q(brms::brm(fml_simple, data = d, empty = TRUE, backend = "rstan"))
  b_sig <- b_simple
  shimA <- make_shim(b_sig, ds$stanfit, M)
  # brms's own ranef() reshapes its r_ draws by position, relying on
  # brms's storage order (every level of one coefficient, then the next
  # coefficient). frmtmb stores level by level. A shim in frmtmb's order
  # would scramble brms's answer, so the group-level comparisons use a
  # shim in brms's order; the variable ORDER itself is a divergence
  # recorded in dev/brmsnames-findings.md
  rcols <- grep("^r_", colnames(M))
  rn_r <- colnames(M)[rcols]
  lev <- sub("^r_[^[]*[[]([^,]*),.*$", "\\1", rn_r)
  cf <- sub("^r_[^[]*[[][^,]*,(.*)[]]$", "\\1", rn_r)
  ord <- order(match(cf, unique(cf)), as.integer(lev))
  Mb <- M
  Mb[, rcols] <- M[, rcols[ord]]
  colnames(Mb)[rcols] <- rn_r[ord]
  shimB <- make_shim(b_sig, ds$stanfit, Mb)
  shimV <- make_shim(b_sig, ds$stanfit, cbind(M, sdm, corm))
  shimR <- shimV

  cat("\n-- controls: the instrument must be able to say no --\n")
  # rvars carry a cache ENVIRONMENT, so two constructions of the same
  # draws are never identical(); brms against itself shows that ceiling
  cmp("CONTROL brms as_draws_rvars twice",
      brms:::as_draws_rvars.brmsfit(shimA),
      brms:::as_draws_rvars.brmsfit(shimA))
  # one draw of one coefficient moved by 1e-9, brms against brms, so the
  # control does not depend on which arm is running: must read DIFFERS
  Mp <- M
  Mp[1L, "b_x"] <- Mp[1L, "b_x"] + 1e-9
  cmp("CONTROL brms fixef(summary = FALSE), perturbed 1e-9",
      brms:::fixef.brmsfit(shimA, summary = FALSE),
      brms:::fixef.brmsfit(make_shim(b_sig, ds$stanfit, Mp),
                           summary = FALSE))
  # the same at 1e-3, past all.equal()'s tolerance: must read DIFFERS
  Mp[1L, "b_x"] <- M[1L, "b_x"] + 1e-3
  cmp("CONTROL brms fixef(summary = FALSE), perturbed 1e-3",
      brms:::fixef.brmsfit(shimA, summary = FALSE),
      brms:::fixef.brmsfit(make_shim(b_sig, ds$stanfit, Mp),
                           summary = FALSE))

  cat("\n-- names --\n")
  cmp("variables()", variables(ds), brms:::variables.brmsfit(shimA))

  cat("\n-- accessors --\n")
  cmp("as.matrix(x)", as.matrix(ds), brms:::as.matrix.brmsfit(shimA))
  cmp("as.matrix(x, variable = 'b_x')", as.matrix(ds, variable = "b_x"),
      brms:::as.matrix.brmsfit(shimA, variable = "b_x"))
  cmp("as.matrix(x, 'b_x') [pars]", as.matrix(ds, "b_x"),
      brms:::as.matrix.brmsfit(shimA, "b_x"))
  cmp("as.matrix(x, draw = 1:5)", as.matrix(ds, draw = 1:5),
      brms:::as.matrix.brmsfit(shimA, draw = 1:5))
  cmp("as.array(x)", as.array(ds), brms:::as.array.brmsfit(shimA))
  cmp("as.array(x, variable = '^b_', regex)",
      as.array(ds, variable = "^b_", regex = TRUE),
      brms:::as.array.brmsfit(shimA, variable = "^b_", regex = TRUE))
  cmp("as.data.frame(x)", as.data.frame(ds),
      brms:::as.data.frame.brmsfit(shimA))
  cmp("as_draws_array(x)", as_draws_array(ds),
      brms:::as_draws_array.brmsfit(shimA))
  cmp("as_draws_df(x, 'b_x')", as_draws_df(ds, "b_x"),
      brms:::as_draws_df.brmsfit(shimA, "b_x"))
  cmp("as_draws_matrix(x)", as_draws_matrix(ds),
      brms:::as_draws_matrix.brmsfit(shimA))
  cmp("as_draws_list(x)", as_draws_list(ds),
      brms:::as_draws_list.brmsfit(shimA))
  cmp("as_draws_rvars(x)", as_draws_rvars(ds),
      brms:::as_draws_rvars.brmsfit(shimA))
  cmp("as_draws(x)", as_draws(ds), brms:::as_draws.brmsfit(shimA))

  cat("\n-- summaries --\n")
  cmp("posterior_summary(x)", posterior_summary(ds),
      brms:::posterior_summary.brmsfit(shimA))
  cmp("posterior_summary(x, '^b_')", posterior_summary(ds, "^b_"),
      brms:::posterior_summary.brmsfit(shimA, "^b_"))
  cmp("posterior_summary(x, robust, probs)",
      posterior_summary(ds, probs = c(0.1, 0.9), robust = TRUE),
      brms:::posterior_summary.brmsfit(shimA, probs = c(0.1, 0.9),
                                       robust = TRUE))
  cmp("posterior_interval(x)", posterior_interval(ds),
      brms:::posterior_interval.brmsfit(shimA))
  cmp("rhat(x)", rhat(ds), brms:::rhat.brmsfit(shimA))
  cmp("neff_ratio(x)", neff_ratio(ds), brms:::neff_ratio.brmsfit(shimA))
  cmp("bayes_R2(x, summary = FALSE) columns",
      colnames(bayes_R2(ds, summary = FALSE)), "R2")
  r2 <- tryCatch(q(bayes_R2(ds, summary = FALSE)), error = function(e) NULL)
  cmp("bayes_R2(x, NULL, TRUE, TRUE) on its draws",
      bayes_R2(ds, NULL, TRUE, TRUE),
      brms:::posterior_summary.default(r2, robust = TRUE))

  cat("\n-- fixef / ranef / coef --\n")
  cmp("fixef(x)", fixef(ds), brms:::fixef.brmsfit(shimA))
  cmp("fixef(x, FALSE)", fixef(ds, FALSE), brms:::fixef.brmsfit(shimA, FALSE))
  cmp("fixef(x, robust, probs)",
      fixef(ds, robust = TRUE, probs = c(0.1, 0.9)),
      brms:::fixef.brmsfit(shimA, robust = TRUE, probs = c(0.1, 0.9)))
  cmp("fixef(x, pars = 'x')", fixef(ds, pars = "x"),
      brms:::fixef.brmsfit(shimA, pars = "x"))
  cmp("ranef(x)", ranef(ds), brms:::ranef.brmsfit(shimB))
  cmp("ranef(x, FALSE)", ranef(ds, FALSE), brms:::ranef.brmsfit(shimB, FALSE))
  cmp("ranef(x, groups = 'g', pars = 'Intercept')",
      ranef(ds, groups = "g", pars = "Intercept"),
      brms:::ranef.brmsfit(shimB, groups = "g", pars = "Intercept"))
  cmp("coef(x)", coef(ds), brms:::coef.brmsfit(shimB))
  cmp("coef(x, FALSE)", coef(ds, FALSE), brms:::coef.brmsfit(shimB, FALSE))
  cmp("coef(x, robust = TRUE)", coef(ds, robust = TRUE),
      brms:::coef.brmsfit(shimB, robust = TRUE))

  cat("\n-- VarCorr --\n")
  vs <- tryCatch(q(VarCorr(ds)), error = function(e) e)
  cmp("VarCorr(x)$g", VarCorr(ds)[["g"]], brms:::VarCorr.brmsfit(shimV)$g)
  cmp("VarCorr(x, NULL, FALSE)$g", VarCorr(ds, NULL, FALSE)[["g"]],
      brms:::VarCorr.brmsfit(shimV, NULL, FALSE)$g)
  cmp("VarCorr(x, robust, probs)$g",
      VarCorr(ds, robust = TRUE, probs = c(0.1, 0.9))[["g"]],
      brms:::VarCorr.brmsfit(shimV, robust = TRUE, probs = c(0.1, 0.9))$g)
  cmp("names(VarCorr(x))", names(VarCorr(ds)),
      names(brms:::VarCorr.brmsfit(shimR)))
  cmp("VarCorr(x)$residual__", VarCorr(ds)[["residual__"]],
      brms:::VarCorr.brmsfit(shimR)$residual__)

  cat("\n-- hypothesis --\n")
  cmp("hypothesis(x, 'x > 0')$hypothesis",
      hypothesis(ds, "x > 0")$hypothesis,
      brms:::hypothesis.brmsfit(shimA, "x > 0")$hypothesis)
  cmp("hypothesis(x, 'x = 0.5')$hypothesis[1:5]",
      hypothesis(ds, "x = 0.5")$hypothesis[1:5],
      brms:::hypothesis.brmsfit(shimA, "x = 0.5")$hypothesis[1:5])
  cmp("hypothesis(x, 'x < Intercept', robust)$hypothesis",
      hypothesis(ds, "x < Intercept", robust = TRUE)$hypothesis,
      brms:::hypothesis.brmsfit(shimA, "x < Intercept",
                                robust = TRUE)$hypothesis)
  cmp("names(hypothesis(x, 'x > 0'))", names(hypothesis(ds, "x > 0")),
      names(brms:::hypothesis.brmsfit(shimA, "x > 0")))
  cmp("hypothesis(x, 'x > 0')$samples",
      hypothesis(ds, "x > 0")$samples,
      brms:::hypothesis.brmsfit(shimA, "x > 0")$samples)
  cmp("hypothesis(x, 'x > 0')$class", hypothesis(ds, "x > 0")$class,
      brms:::hypothesis.brmsfit(shimA, "x > 0")$class)
  cmp("hypothesis(x, ..., scope = 'coef')$hypothesis",
      hypothesis(ds, "Intercept > 0", scope = "coef",
                 group = "g")$hypothesis,
      brms:::hypothesis.brmsfit(shimB, "Intercept > 0", scope = "coef",
                                group = "g")$hypothesis)
  cmp("hypothesis(x, ..., scope = 'ranef')$hypothesis",
      hypothesis(ds, "Intercept < 0", scope = "ranef",
                 group = "g")$hypothesis,
      brms:::hypothesis.brmsfit(shimB, "Intercept < 0", scope = "ranef",
                                group = "g")$hypothesis)
  cat("\n")
}

tab <- do.call(rbind, res)
utils::write.table(tab, sprintf("dev/brmsnames-log/match-%s.tsv", arm),
                   sep = "\t", row.names = FALSE, quote = FALSE)
ctl <- startsWith(tab$call, "CONTROL")
cat("== controls, arm ", arm, " ==\n", sep = "")
print(tab[ctl, c("model", "call", "result")], row.names = FALSE)
cat("\n== counts, arm ", arm, " (controls excluded) ==\n", sep = "")
print(table(tab$result[!ctl]))
cat("comparisons:", sum(!ctl), "\n")
cat("DONE\n")
