# Lane wt-mvprior: nothing but the refused priors may move. Fits every
# design below, none of which carries a prefix-less prior on a
# multivariate model or a location prior without nlpar on a nonlinear
# one, under ML, REML = TRUE and control(profile = TRUE), and saves
# vcov(), fixef(), logLik(), summary(), and predict() at a fixed seed.
# The univariate designs DO carry prefix-less priors, which must stay
# untouched. Also frm_sample()'s default priors on the multivariate
# designs: the resolved entries and the negative log prior at the
# estimate, which is what said nothing there relied on the broadcast.
#   MVPRIOR_ARM=base Rscript dev/mvprior-bitwise.R -> bitwise-base.rds
#   MVPRIOR_ARM=lane Rscript dev/mvprior-bitwise.R -> bitwise-lane.rds
#   Rscript dev/mvprior-bitwise.R compare          -> prints the verdict
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-log"
`%||%` <- function(a, b) if (is.null(a)) b else a
if (identical(commandArgs(trailingOnly = TRUE)[1], "compare")) {
  a <- readRDS(file.path(OUT, "bitwise-base.rds"))
  b <- readRDS(file.path(OUT, "bitwise-lane.rds"))
  stopifnot(identical(names(a), names(b)))
  n_ok <- 0L
  n_all <- 0L
  for (k in names(a)) {
    if (!identical(names(a[[k]]), names(b[[k]]))) {
      # one arm fitted and the other refused: print both outcomes
      cat(sprintf("%-36s %-10s DIFFERS base: %s | lane: %s
", k, "outcome",
                  substr(paste(a[[k]][["error"]] %||% "fitted"), 1, 90),
                  substr(paste(b[[k]][["error"]] %||% "fitted"), 1, 90)))
      n_all <- n_all + 1L
      next
    }
    for (q in names(a[[k]])) {
      # environments are pointers that two processes never share, so a
      # closure is compared by its code and values
      same <- identical(a[[k]][[q]], b[[k]][[q]], ignore.environment = TRUE,
                        ignore.bytecode = TRUE, ignore.srcref = TRUE)
      cat(sprintf("%-36s %-10s %s\n", k, q,
                  if (same) "identical" else "DIFFERS"))
      n_ok <- n_ok + same
      n_all <- n_all + 1L
    }
  }
  cat("identical:", n_ok, "of", n_all, "\n")
  quit(save = "no")
}
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-prelude.R")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(frmtmb.latent)
  library(frmtmb.eam)
})
cat("arm", mvprior_arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), " frmtmb.sample from",
    dirname(find.package("frmtmb.sample")), "\n")

set.seed(20260923)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, 20)),
                t = rep(1:20, each = 10))
u <- rnorm(10, 0, 0.6)
d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(n, 0, exp(-0.2 + 0.3 * d$z))
d$y1 <- d$y
d$y2 <- -1 + 0.3 * d$x + 0.5 * u[d$g] + rnorm(n)
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + u[d$g]))
d$ax <- abs(d$x)
d$ynl <- (2 + rnorm(10, 0, 0.3)[d$g]) * exp(-0.5 * d$ax) + rnorm(n, 0, 0.1)
dmi <- d
dmi$xm <- 0.5 * dmi$z + rnorm(n)
dmi$ym <- 1 + 0.7 * dmi$xm + rnorm(n)
dmi$xm[sample.int(n, 25)] <- NA
# appended after every draw the earlier designs use, so their data are
# the same as in the first battery
d$cat <- factor(c("a", "b", "c")[1 + (d$x + rlogis(n) > 0) +
                                   (d$z + rlogis(n) > 0.5)])
d$ym <- ifelse(rbinom(n, 1, 0.5) == 1, 3 + d$x, -1 + d$x) + rnorm(n, 0, 0.5)

# punch round 2: the extension families with several locations and the
# cs() ordinal models, drawn after everything above so the earlier
# designs keep their data
d$ord <- factor(cut(d$x + 0.5 * d$z + rlogis(n), c(-Inf, -1, 0, 1.2, Inf),
                    labels = FALSE), ordered = TRUE)
dlca <- local({
  cl <- rbinom(n, 1, 0.4) + 1
  pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
  Y <- matrix(0L, n, 4)
  for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
  dl <- data.frame(x = rnorm(n))
  dl$Y <- Y
  dl
})
dhmm <- do.call(rbind, lapply(1:10, function(id) {
  st <- integer(20)
  st[1] <- 1L
  for (t in 2:20) {
    st[t] <- sample.int(2, 1, prob = if (st[t - 1] == 1) c(.9, .1) else
      c(.2, .8))
  }
  data.frame(id = id, t = 1:20, x = rnorm(20), y = rnorm(20, c(0, 3)[st]))
}))
dlba <- lba_simulate(300, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.4, ndt = 0.2)
dlba$x <- rnorm(nrow(dlba))
drdm <- rdm_simulate(300, v = c(2.4, 1.6, 1.0), A = 0.5, k = 0.5, ndt = 0.2)
drdm$x <- rnorm(nrow(drdm))
dmvn <- local({
  clm <- rbinom(n, 1, 0.4)
  dm <- data.frame(x = rnorm(n))
  dm$Y <- cbind(rnorm(n, ifelse(clm == 1, 0, 3)),
                rnorm(n, ifelse(clm == 1, 0, 4)))
  dm
})

p <- function(...) set_prior(...)
designs <- list(
  # univariate, WITH prefix-less priors: must not move
  uni_prefixless = list(
    bf(y ~ x + (1 | g)) + gaussian(), d,
    p("normal(0, 1)", class = "b") + p("normal(0, 5)", class = "Intercept") +
      p("exponential(1)", class = "sd") +
      p("student_t(3, 0, 2.5)", class = "sigma")),
  uni_dpar = list(
    bf(y ~ x, sigma ~ z) + gaussian(), d,
    p("normal(0, 1)", class = "b", dpar = "sigma") +
      p("normal(0, 2)", class = "Intercept", dpar = "sigma") +
      p("normal(0, 1)", class = "b", coef = "x")),
  uni_poisson = list(
    bf(cnt ~ x + (1 | g)) + poisson(), d,
    p("normal(0, 1)", class = "b", lb = -2)),
  uni_noprior = list(bf(y ~ x + (1 | g)) + gaussian(), d, NULL),
  # nonlinear, with nlpar and a prefix-less sigma: must not move
  nl = list(
    bf(ynl ~ a * exp(-b * ax), a ~ 1 + (1 | g), b ~ 1, nl = TRUE) +
      gaussian(), d,
    p("normal(2, 1)", nlpar = "a") + p("normal(0.5, 1)", nlpar = "b") +
      p("student_t(3, 0, 2.5)", class = "sigma")),
  # multivariate, every prior with resp: must not move
  mv_resp = list(
    bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)) + set_rescor(FALSE) +
      gaussian(), d,
    p("normal(0, 1)", class = "b", resp = "y1") +
      p("normal(0, 2)", class = "b", resp = "y2") +
      p("normal(0, 5)", class = "Intercept", resp = "y1") +
      p("student_t(3, 0, 2.5)", class = "sigma", resp = "y2") +
      p("exponential(1)", class = "sd", resp = "y1")),
  mv_rescor = list(
    bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE) + gaussian(), d,
    p("lkj(2)", class = "rescor") + p("normal(0, 1)", class = "b",
                                      resp = "y2")),
  mv_mixed = list(
    bf(y1 ~ x, family = gaussian()) + bf(cnt ~ x, family = poisson()) +
      set_rescor(FALSE), d,
    p("normal(0, 1)", class = "b", resp = "cnt") +
      p("student_t(3, 0, 2.5)", class = "sigma", resp = "y1")),
  mv_noprior = list(
    bf(y1 ~ x + (1 | g)) + bf(y2 ~ x) + set_rescor(FALSE) + gaussian(), d,
    NULL),
  mv_mi = list(
    bf(ym ~ mi(xm) + z) + bf(xm | mi() ~ z) + set_rescor(FALSE) +
      gaussian(), dmi,
    p("normal(0, 1)", class = "b", resp = "ym")),
  # several locations, every prior naming its dpar: must not move
  cat_dpar = list(
    bf(cat ~ x + (1 | g)) + categorical(), d,
    p("normal(0, 1)", class = "b", dpar = "mub") +
      p("normal(0, 2)", class = "Intercept", dpar = "muc") +
      p("exponential(1)", class = "sd", dpar = "mub")),
  cat_noprior = list(bf(cat ~ x) + categorical(), d, NULL),
  cat_re_noprior = list(bf(cat ~ x + (1 | g)) + categorical(), d, NULL),
  mix_noprior = list(bf(ym ~ x) + mixture(gaussian(), gaussian()), d, NULL),
  mix_dpar = list(
    bf(ym ~ x) + mixture(gaussian(), gaussian()), d,
    p("normal(0, 1)", class = "b", dpar = "mu1") +
      p("normal(-1, 2)", class = "Intercept", dpar = "mu1") +
      p("normal(3, 2)", class = "Intercept", dpar = "mu2") +
      p("student_t(3, 0, 2.5)", class = "sigma1")),
  mvcat_dpar = list(
    bf(y1 ~ x, family = gaussian()) +
      bf(cat ~ x, family = categorical()) + set_rescor(FALSE), d,
    p("normal(0, 1)", class = "b", dpar = "mub", resp = "cat") +
      p("normal(0, 1)", class = "b", resp = "y1")),
  lca_dpar = list(
    bf(Y ~ x) + lca(K = 2), dlca,
    p("normal(0, 1)", class = "b", dpar = "theta1")),
  hmm_dpar = list(
    bf(y ~ x) + hmm(K = 2, gaussian(), time = t, group = id), dhmm,
    p("normal(0, 1)", class = "b", dpar = "mu1") +
      p("normal(3, 2)", class = "Intercept", dpar = "mu2")),
  lba_dpar = list(
    bf(rt | vint(choice) ~ x) + lba(3), dlba,
    p("normal(0, 1)", class = "b", dpar = "v1")),
  rdm_dpar = list(
    bf(rt | vint(choice) ~ x) + rdm(3), drdm,
    p("normal(0, 1)", class = "b", dpar = "v2")),
  mvn_dpar = list(
    bf(Y ~ x) + mixture_mvn(K = 2, D = 2), dmvn,
    p("normal(0, 1)", class = "b", dpar = "mu1d1")),
  # cs() with no class b prior: must not move
  cs_noprior = list(bf(ord ~ z + cs(x)) + acat(), d, NULL),
  cs_intercept = list(bf(ord ~ z + cs(x)) + sratio(), d,
                      p("normal(0, 2)", class = "Intercept")),
  mv_ar = list(
    bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)) + bf(y2 ~ x) +
      set_rescor(FALSE) + gaussian(), d,
    p("normal(0, 0.5)", class = "ar", resp = "y1"))
)
settings <- list(ML = list(), REML = list(REML = TRUE),
                 profile = list(control = frmtmb_control(profile = TRUE)))
safe <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}
grab <- function(fit) {
  out <- list(vcov = safe(vcov(fit)), fixef = safe(fixef(fit)),
              logLik = safe(logLik(fit)), summary = safe(summary(fit)))
  set.seed(7)
  out$predict <- safe(predict(fit))
  out
}
res <- list()
for (nm in names(designs)) {
  ds <- designs[[nm]]
  for (st in names(settings)) {
    fit <- tryCatch(suppressWarnings(suppressMessages(do.call(
      frm, c(list(ds[[1]], data = ds[[2]], prior = ds[[3]]),
             settings[[st]])))),
      error = function(e) e)
    key <- paste(nm, st)
    res[[key]] <- if (inherits(fit, "error")) {
      list(error = conditionMessage(fit))
    } else {
      grab(fit)
    }
    cat(key, if (inherits(fit, "error")) conditionMessage(fit) else "ok",
        "\n")
  }
}
# frm_sample()'s defaults, resolved against the multivariate designs
for (nm in c("mv_resp", "mv_rescor", "mv_mixed", "mv_noprior", "mv_mi",
             "mv_ar", "cat_dpar", "cat_noprior", "cat_re_noprior",
             "mix_dpar", "mix_noprior", "mvcat_dpar")) {
  ds <- designs[[nm]]
  uf <- suppressMessages(frm(ds[[1]], data = ds[[2]],
                             dry_run = "objective"))
  fit <- suppressWarnings(suppressMessages(frm(ds[[1]], data = ds[[2]])))
  ri <- safe(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(uf, ds[[3]]))$ri)
  res[[paste("sample defaults", nm)]] <- list(
    entries = if (is.list(ri)) ri$entries else ri,
    nlp = if (is.list(ri)) neg_log_prior_fn(ri$entries)(fit$estimates) else
      ri)
  cat("sample defaults", nm, "ok\n")
}
saveRDS(res, file.path(OUT, paste0("bitwise-", mvprior_arm, ".rds")))
cat("saved", length(res), "designs\n")
