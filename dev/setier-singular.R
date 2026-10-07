# Lane setier: how often an ordinary lme4-style fit at a variance
# boundary gets frmtmb's boundary message, its standard-error warning, or
# nothing, against lme4's isSingular() on the same data.
#   Rscript dev/setier-singular.R <lib or "base"> <out .tsv> [seeds]
# Designs (seeds 1..S each, default 100):
#   ri6   y ~ x + (1 | g), 6 groups of 5, true group sd 0
#   ri20  y ~ x + (1 | g), 20 groups of 5, true group sd 0
#   ri20s y ~ x + (1 | g), 20 groups of 5, true group sd 0.3 (mostly not
#         singular: the false-alarm arm)
#   rs20  y ~ x + (1 + x | g), 20 groups of 6, true slope sd 0
#   bin15 cbind(k, 10 - k) ~ x + (1 | g), binomial, 15 groups of 4, sd 0
#   bin15s the same with group sd 0.5 (false-alarm arm)
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages({
  library(frmtmb)
  library(lme4)
})
out <- args[2]
S <- if (length(args) >= 3) as.integer(args[3]) else 100L
cat("lib:", find.package("frmtmb"), " BLAS probe (s):",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
ns <- asNamespace("frmtmb")
conds <- function(expr) {
  w <- character(0); m <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) {
    m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage")
  })
  list(value = val, warnings = w, messages = m)
}
gen <- function(design, seed) {
  set.seed(seed)
  if (design %in% c("ri6", "ri20", "ri20s")) {
    G <- if (design == "ri6") 6 else 20
    sdg <- if (design == "ri20s") 0.3 else 0
    d <- data.frame(g = factor(rep(seq_len(G), each = 5)),
                    x = rnorm(G * 5))
    d$y <- 1 + 0.5 * d$x + rnorm(G, 0, sdg)[d$g] + rnorm(G * 5)
    list(d = d, ff = y ~ x + (1 | g), fam = gaussian())
  } else if (design == "rs20") {
    d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
    d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
    list(d = d, ff = y ~ x + (1 + x | g), fam = gaussian())
  } else {
    sdg <- if (design == "bin15s") 0.5 else 0
    d <- data.frame(g = factor(rep(1:15, each = 4)), x = rnorm(60))
    p <- plogis(-0.5 + 0.5 * d$x + rnorm(15, 0, sdg)[d$g])
    d$k <- rbinom(60, 10, p)
    d$nk <- 10 - d$k
    list(d = d, ff = cbind(k, nk) ~ x + (1 | g), fam = binomial())
  }
}
rows <- list()
for (design in c("ri6", "ri20", "ri20s", "rs20", "bin15", "bin15s")) {
  for (seed in seq_len(S)) {
    g <- gen(design, seed)
    lm4 <- suppressMessages(suppressWarnings(
      if (identical(g$fam$family, "gaussian")) {
        lmer(g$ff, data = g$d, REML = FALSE)
      } else glmer(g$ff, data = g$d, family = g$fam)))
    t0 <- proc.time()[["elapsed"]]
    r <- conds(frm(g$ff, family = g$fam, data = g$d))
    t1 <- proc.time()[["elapsed"]]
    fit <- r$value
    s <- conds(summary(fit))
    v <- conds(VarCorr(fit))
    sdr <- ns$sdr_of(fit)
    lost <- sdr$se_lost
    th <- fit$estimates[["theta"]]
    rows[[length(rows) + 1L]] <- data.frame(
      design = design, seed = seed,
      lme4_singular = isSingular(lm4),
      code = fit$opt$convergence,
      min_theta = min(th),
      n_lost = length(lost),
      lost = paste(names(lost), lost, sep = ":", collapse = ";"),
      boundary_msg = sum(grepl("^Boundary [(]singular[)] fit",
                               c(r$messages, s$messages))),
      se_warn = sum(grepl("Standard errors are not available",
                          c(r$warnings, s$warnings))),
      other_warn = sum(!grepl("Standard errors are not available",
                              c(r$warnings, s$warnings))),
      varcorr_warn = length(v$warnings),
      se_theta_max = suppressWarnings(max(sqrt(diag(sdr$cov.fixed))[
        grepl("^theta", rownames(sdr$cov.fixed))])),
      fit_s = t1 - t0)
  }
}
res <- do.call(rbind, rows)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", nrow(res), "rows to", out, "\n")
