# Item 3. brms is the tiebreaker, and the question is whether rhat()
# on frmtmb.sample draws returns the SAME R-HAT DEFINITION brms
# returns on a fit. Three things are measured:
#   1. which generic brms's own rhat() reaches, and where rhat.brmsfit
#      is registered;
#   2. what each implementation computes, read off the bodies;
#   3. the NUMBERS, on one set of draws, three ways.
.libPaths(c("C:/Users/adf44/source/r/sgrev-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
q(library(posterior)); q(library(bayesplot)); q(library(rstan))

cat("== 1. where brms's rhat comes from ==\n")
q(requireNamespace("brms", quietly = TRUE))
bg <- get("rhat", envir = asNamespace("brms"))
cat("brms's rhat generic environment: ",
    environmentName(environment(bg)), "\n")
for (p in c("posterior", "bayesplot")) {
  t <- get(".__S3MethodsTable__.", envir = asNamespace(p),
           inherits = FALSE)
  cat("  rhat.brmsfit in ", p, "'s table: ",
      exists("rhat.brmsfit", envir = t, inherits = FALSE), "\n", sep = "")
}
m <- get("rhat.brmsfit", envir = get(".__S3MethodsTable__.",
         envir = asNamespace("posterior"), inherits = FALSE))
cat("brms's rhat.brmsfit body:\n")
cat(paste(deparse(body(m)), collapse = "\n"), "\n\n")

cat("== 2. what each implementation computes ==\n")
cat("-- frmtmb.sample's rhat.frmtmb_draws --\n")
t <- get(".__S3MethodsTable__.", envir = asNamespace("frmtmb.sample"),
         inherits = FALSE)
fm <- get("rhat.frmtmb_draws", envir = t, inherits = FALSE)
cat(paste(deparse(body(fm)), collapse = "\n"), "\n")
cat("-- bayesplot's rhat.stanfit --\n")
bt <- get(".__S3MethodsTable__.", envir = asNamespace("bayesplot"),
          inherits = FALSE)
cat(paste(grep("^rhat", ls(bt, all.names = TRUE), value = TRUE),
          collapse = " "), "\n")
bs <- tryCatch(get("rhat.stanfit", envir = bt, inherits = FALSE),
               error = function(e) NULL)
if (!is.null(bs)) cat(paste(deparse(body(bs)), collapse = "\n"), "\n")
cat("-- posterior's rhat.default / rhat.draws --\n")
pt <- get(".__S3MethodsTable__.", envir = asNamespace("posterior"),
          inherits = FALSE)
cat(paste(grep("^rhat", ls(pt, all.names = TRUE), value = TRUE),
          collapse = " "), "\n")
pd <- get("rhat.default", envir = pt, inherits = FALSE)
cat(paste(deparse(body(pd)), collapse = "\n"), "\n\n")

cat("== 3. the numbers ==\n")
cache <- "dev/stan-cache/sgrev-draws.rds"
if (file.exists(cache)) {
  ds <- readRDS(cache)
} else {
  set.seed(9)
  dd <- data.frame(x = stats::rnorm(120), g = factor(rep(1:6, 20)))
  dd$y <- stats::rnorm(120, 1 + 0.5 * dd$x +
                         stats::rnorm(6, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ds <- q(frm_sample(fit, chains = 4, iter = 1000, refresh = 0,
                     seed = 20260915))
  saveRDS(ds, cache)
}
a <- rhat(ds)
cat("rhat(ds) length ", length(a), " names ",
    paste(head(names(a), 6), collapse = ","), "\n")

arr <- as_draws_array(ds)
post <- vapply(posterior::variables(arr), function(v)
  posterior::rhat(posterior::extract_variable_matrix(arr, v)), 0)
sd <- posterior::summarise_draws(arr, "rhat")
pb <- setNames(sd$rhat, sd$variable)

common <- intersect(names(a), names(pb))
cat("variables in common: ", length(common), "\n")
cat(sprintf("%-22s %12s %12s %12s %10s\n", "variable",
            "rhat(ds)", "posterior", "post(extract)", "rel.diff"))
d <- numeric(0)
for (v in common) {
  r1 <- a[[v]]; r2 <- pb[[v]]
  rel <- abs(r1 - r2) / abs(r2)
  d <- c(d, rel)
  cat(sprintf("%-22s %12.8f %12.8f %12.8f %10.3e\n", v, r1, r2,
              post[[v]], rel))
}
cat("\nmax relative difference rhat(ds) vs posterior::rhat: ",
    format(max(d), digits = 6), "\n")
cat("identical():   ", identical(unname(a[common]), unname(pb[common])),
    "\n")
cat("how far from 1 is the DIAGNOSTIC itself: max |posterior - 1| = ",
    format(max(abs(pb[common] - 1)), digits = 6), "\n")
cat("difference as a fraction of (posterior rhat - 1): ",
    format(max(abs(a[common] - pb[common]) / abs(pb[common] - 1)),
           digits = 6), "\n")

# what a rank-normalized vs classic split-Rhat difference looks like
# on the SAME draws, so the size above has a yardstick
cat("\n-- yardstick: rstan's classic split-Rhat on the same stanfit --\n")
sf <- ds$stanfit
if (!is.null(sf)) {
  s <- rstan::summary(sf)$summary
  rr <- s[, "Rhat"]
  cm <- intersect(names(a), names(rr))
  cat("max |rhat(ds) - rstan Rhat| over ", length(cm), " vars: ",
      format(max(abs(a[cm] - rr[cm])), digits = 6), "\n")
  cat("max |posterior - rstan Rhat|: ",
      format(max(abs(pb[cm] - rr[cm])), digits = 6), "\n")
}
cat("DONE\n")
