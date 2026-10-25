## mo() of integers that do not start at 0: brms 2.23.0 subtracts the
## minimum (get_mo_values()); frmtmb codes D = max(x). Seed 5.
##   Rscript dev/review-brmsnames-mo.R lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5)
n <- 300
d <- data.frame(k = sample(1:4, n, replace = TRUE))
d$k2 <- d$k + 1L
d$k0 <- d$k - 1L
zeta <- c(0.6, 0.3, 0.1)
d$y <- 2 + 1.5 * c(0, cumsum(zeta))[d$k] * 3 + rnorm(n)
cat("== brms 2.23.0 standata\n")
for (v in c("k0", "k", "k2")) {
  sd <- q(brms::standata(as.formula(paste0("y ~ mo(", v, ")")), data = d))
  cat(sprintf("mo(%s): range %s, Jmo = %d, Xmo_1 range %s\n", v,
              paste(range(d[[v]]), collapse = ".."), sd$Jmo[1],
              paste(range(sd$Xmo_1), collapse = "..")))
}
q(library(frmtmb))
cat("== frmtmb", arm, dirname(system.file(package = "frmtmb")), "\n")
fits <- list()
for (v in c("k0", "k", "k2")) {
  f <- as.formula(paste0("y ~ mo(", v, ")"))
  w <- character(0)
  fit <- withCallingHandlers(frm(f, data = d), warning = function(w0) {
    w <<- c(w, conditionMessage(w0))
    invokeRestart("muffleWarning")
  })
  fits[[v]] <- fit
  mi <- fit$frame$mo_info %||% fit$frame$dpars[[1]]$mo_info
  fe <- fixef(fit)
  cat(sprintf("mo(%s): logLik %.6f df %d; warnings %d\n", v,
              as.numeric(logLik(fit)), attr(logLik(fit), "df"), length(w)))
  print(round(fe[, 1:2], 4))
  cat("simo:", paste(round(fit$estimates[grep("^zeta",
      names(fit$estimates))][[1]], 4), collapse = " "), "\n")
}
nd <- data.frame(k = 1:4, k0 = 0:3, k2 = 2:5)
cat("fitted at k = 1..4, from mo(k0), mo(k), mo(k2):\n")
print(round(sapply(fits, function(f) fitted(f, newdata = nd)[, 1]), 6))
cat("newdata with only k = 2..4 (brms would re-read the minimum unless",
    "stored):\n")
print(round(fitted(fits$k, newdata = data.frame(k = 2:4))[, 1], 6))
if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
  q(library(frmtmb.sample))
  Sys.setenv(FRMTMB_STAN_CACHE = file.path(getwd(), "dev/stan-cache"))
  for (v in c("k0", "k")) {
    w <- character(0)
    ds <- withCallingHandlers(
      frm_sample(fits[[v]], chains = 2, iter = 1000, refresh = 0, seed = 11),
      warning = function(w0) {
        w <<- c(w, conditionMessage(w0))
        invokeRestart("muffleWarning")
      })
    s <- fixef(ds)
    cat(sprintf("frm_sample mo(%s): warnings %d\n", v, length(w)))
    if (length(w)) cat("  ", substr(w, 1, 160), sep = "\n  ")
    print(round(s[, 1:4], 3)); sv <- as.matrix(ds)[, grep("^simo", variables(ds)), drop = FALSE]; print(round(colMeans(sv), 3))
    cat("simo vars:", grep("^simo", variables(ds), value = TRUE), "\n")
  }
}
