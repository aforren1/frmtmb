# Reviewer: newdata, frm_sample(), saved fit without frmtmb.survival, and
# a third-party option term written from the help page alone.
#   Rscript dev/review-bhaz-compose2.R > dev/review-bhaz-log/compose2.txt 2>&1
WT <- "C:/Users/adf44/source/r/frmtmb-wt-bhaz"
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(frmtmb.survival))
cat("lib:", find.package("frmtmb"), find.package("frmtmb.survival"), "\n")
set.seed(2); n <- 300; x <- rnorm(n); z <- rnorm(n)
tt <- stats::rexp(n, exp(-0.5 + 0.7 * x)); ct <- stats::rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 z = z)
f6b <- frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = cox(),
           data = dd)
f6c <- frm(bf(time | cens(cens) ~ x), family = cox(df = 6), data = dd)
nd <- dd[c(3, 50, 120), ]
cat("frm_linpred(newdata) identical:",
    identical(frm_linpred(f6b, newdata = nd),
              frm_linpred(f6c, newdata = nd)), "\n")
nd2 <- nd; nd2$time <- NULL; nd2$cens <- NULL
r <- tryCatch(frm_linpred(f6b, newdata = nd2), error = function(e)
  conditionMessage(e))
cat("frm_linpred(newdata without response):", if (is.character(r)) r else
  paste("ok, identical to cox(df=6):",
        identical(r, tryCatch(frm_linpred(f6c, newdata = nd2),
                              error = function(e) "err"))), "\n")
cat("frm_curve on bhaz fit:", tryCatch({ g <- data.frame(x = c(-1, 0, 1));
  a <- frm_curve(f6b, g); b <- frm_curve(f6c, g)
  dif <- names(a)[!vapply(names(a), function(k) identical(a[[k]], b[[k]]),
                          NA)]
  paste("identical", identical(a, b), "| columns that differ:",
        paste(dif, collapse = ","), "| attrs differ:",
        paste(setdiff(union(names(attributes(a)), names(attributes(b))),
                      names(attributes(a))[vapply(names(attributes(a)),
          function(k) identical(attr(a, k), attr(b, k)), NA)]),
          collapse = ","))
}, error = function(e) conditionMessage(e)), "\n")
cat("cox_baseline identical:",
    identical(cox_baseline(f6b), cox_baseline(f6c)), "\n")

cat("\n== frm_sample ==\n")
suppressPackageStartupMessages(library(frmtmb.sample))
Sys.setenv(FRMTMB_STAN_CACHE = file.path(WT, "dev/stan-cache"))
sb <- suppressWarnings(suppressMessages(frm_sample(
  f6b, chains = 1, iter = 300, seed = 1, refresh = 0)))
sc <- suppressWarnings(suppressMessages(frm_sample(
  f6c, chains = 1, iter = 300, seed = 1, refresh = 0)))
mb <- as.matrix(sb); mc <- as.matrix(sc)
cat("variables:", paste(colnames(mb)[grepl("sbhaz|b_", colnames(mb))],
                        collapse = " "), "\n")
cat("draws identical:", identical(mb, mc), " max abs diff",
    max(abs(mb - mc)), "\n")
cat("draws log_lik identical:", identical(log_lik(sb), log_lik(sc)), "\n")
eb <- tryCatch(posterior_epred(sb, newdata = nd), error = function(e)
  conditionMessage(e))
cat("posterior_epred(newdata):", if (is.character(eb)) substr(eb, 1, 120)
    else identical(eb, posterior_epred(sc, newdata = nd)), "\n")
gp <- tryCatch(get_prior(bf(time | cens(cens) + bhaz(df = 6) ~ x),
                         family = cox(), data = dd),
               error = function(e) conditionMessage(e))
print(as.data.frame(gp)[as.data.frame(gp)$class == "sbhaz", ])

cat("\n== saved fit, frmtmb.survival absent ==\n")
rds <- file.path(WT, "dev/review-bhaz-fit.rds")
saveRDS(f6b, rds)
corelib <- file.path(WT, "dev/review-bhaz-corelib")
child <- function(libs, lines) {
  f <- tempfile(fileext = ".R")
  writeLines(c(sprintf(".libPaths(%s)", paste(deparse(libs), collapse = "")),
               lines), f)
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  c("--vanilla", f), stdout = TRUE,
                                  stderr = TRUE))
  cat(paste0("   ", out), sep = "\n")
}
child(c(corelib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"),
      c("suppressMessages(library(frmtmb))",
        "cat('survival installed:', nzchar(system.file(package = 'frmtmb.survival')), '\\n')",
        sprintf("f <- tryCatch(readRDS('%s'), error = function(e) e)", rds),
        "if (inherits(f, 'error')) cat('readRDS:', conditionMessage(f), '\\n')",
        "for (nm in c('frm_linpred', 'summary', 'logLik')) {",
        "  r <- tryCatch({get(nm)(f); 'ok'}, error = function(e) conditionMessage(e))",
        "  cat(nm, ':', substr(r, 1, 300), '\\n')",
        "}",
        "r <- tryCatch(frm(bf(time | cens(cens) + bhaz(df = 6) ~ x), family = 'cox', data = data.frame(time = 1:5, cens = 0, x = 1:5), dry_run = 'spec'), error = function(e) conditionMessage(e))",
        "cat('bhaz() without survival:', substr(r, 1, 300), '\\n')"))
child(.libPaths(),
      c("suppressMessages(library(frmtmb.survival))",
        sprintf("f <- readRDS('%s')", rds),
        "nd <- f$data[c(3, 50, 120), ]",
        "cat('reloaded frm_linpred(newdata):', format(frm_linpred(f, newdata = nd), digits = 10), '\n')",
        "d2 <- f$data; d2$z <- rnorm(nrow(d2))",
        "u <- suppressWarnings(update(f, . ~ . + z, newdata = d2))",
        "cat('reloaded update(+z): K', length(cox_baseline(u)), '\\n')"))
unlink(rds)

cat("\n== third-party option term from the help page ==\n")
# A normal family with a known sd that the term sets: y | sdk(2) ~ x.
# The help says option(family, call, env) returns the family the fit
# uses; this one rebuilds it.
knownsd <- function(sd = 1) {
  fam <- frmtmb_family(
    "knownsd", dpars = "mu", links = list(mu = "identity"),
    accepts_aterms = c("sdk"),
    lpdf = function(y, dpars, aterms, extra) {
      stats::dnorm(y, dpars$mu, sd, log = TRUE)
    })
  fam$knownsd_sd <- sd
  fam
}
calls <- 0L
frmtmb_local_registries()
r <- tryCatch(frmtmb_register_aterm("sdk", option = function(family, call,
                                                               env) {
  calls <<- calls + 1L
  knownsd(eval(call[[2L]], env))
}), error = function(e) conditionMessage(e))
cat("register:", if (is.null(r)) "ok" else r, "\n")
set.seed(4); dk <- data.frame(x = rnorm(100)); dk$y <- 1 + dk$x +
  rnorm(100, 0, 2)
ft <- tryCatch(frm(bf(y | sdk(2) ~ x), family = knownsd(), data = dk),
               error = function(e) e)
if (inherits(ft, "error")) {
  cat("fit with option term: ERROR", conditionMessage(ft), "\n")
} else {
  ref <- sum(stats::dnorm(dk$y, fitted(lm(y ~ x, dk)), 2, log = TRUE))
  cat(sprintf("fit with option term: logLik %.10f, reference %.10f\n",
              as.numeric(logLik(ft)), ref))
  cat("option calls during frm():", calls, "\n")
  calls <- 0L
  u <- update(ft, . ~ . + I(x^2))
  cat("update calls:", calls, " stored family sd after fit:",
      family(ft)$knownsd_sd, "\n")
}
# a non-idempotent option: multiplies the width each time it runs
frmtmb_register_aterm("sdm", option = function(family, call, env) {
  family$knownsd_sd2 <- (family$knownsd_sd2 %||% 1) * eval(call[[2L]], env)
  family
})
kf <- knownsd(); kf$accepts_aterms <- NULL
r <- tryCatch({
  fm <- frm(bf(y | sdm(2) ~ x), family = knownsd(), data = dk)
  u <- update(fm, . ~ . + I(x^2))
  c(fit = family(fm)$knownsd_sd2 %||% NA, upd = family(u)$knownsd_sd2 %||% NA)
}, error = function(e) conditionMessage(e))
cat("non-idempotent option, sd2 after fit / after update:", r, "\n")
cat("compat rows for sdk:\n")
print(tryCatch(frm_compat("sdk"), error = function(e) conditionMessage(e)))
