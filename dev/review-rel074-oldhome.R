# Release review of 0.74.0, short check: old_core_family_home() on real
# 0.73.0 fits, on 0.74.0 fits, and on user families built at top level.
#   Rscript dev/review-rel074-oldhome.R save-wiener <dir>   (rellib-r11)
#   Rscript dev/review-rel074-oldhome.R read <dir> <libs;..> <attach;..|none>
#   Rscript dev/review-rel074-oldhome.R user
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "save-wiener") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r11", user))
  suppressMessages({library(frmtmb); library(frmtmb.eam)})
  set.seed(3)
  dw <- data.frame(x = rnorm(120))
  dw$rt <- 0.3 + rgamma(120, 2, 6)
  dw$resp <- rbinom(120, 1, 0.6)
  fw <- frm(bf(rt | dec(resp) ~ x), family = wiener(), data = dw)
  cat("fwien logLik", format(as.numeric(logLik(fw)), digits = 12),
      "| package attr", attr(fw$spec$responses[[1]]$family, "frmtmb_package"),
      "\n")
  saveRDS(list(fit = fw, data = dw, ref = list(vcov = vcov(fw))),
          file.path(a[2], "fwien.rds"))
  quit(save = "no")
}
if (a[1] == "user") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r12", user))
  suppressMessages(library(frmtmb))
  och <- getFromNamespace("old_core_family_home", "frmtmb")
  rmp <- getFromNamespace("require_model_packages", "frmtmb")
  lp <- function(y, dpars, aterms) RTMB::dnorm(y, dpars$mu, exp(dpars$sigma), log = TRUE)
  set.seed(4)
  d <- data.frame(x = runif(80)); d$y <- 1 + 2 * d$x + rnorm(80, 0, 0.3)
  for (nm in c("cox", "wiener", "mytfam")) {
    for (ctor in c("frmtmb_family", "custom_family")) {
      fam <- tryCatch(do.call(ctor, list(family = nm, dpars = c("mu", "sigma"),
                                         links = c(mu = "identity", sigma = "identity"),
                                         lpdf = lp)),
                      error = function(e) e)
      if (inherits(fam, "error")) {
        cat(sprintf("%-14s %-8s build ERROR %s\n", ctor, nm,
                    substr(conditionMessage(fam), 1, 120)))
        next
      }
      env <- environmentName(topenv(environment(fam$lpdf)))
      home <- och(fam)
      fit <- tryCatch({
        f <- frm(bf(y ~ x), family = fam, data = d)
        sprintf("fits, logLik %.6f; summary %s", as.numeric(logLik(f)),
                tryCatch({utils::capture.output(summary(f)); "ok"},
                         error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 90))))
      }, error = function(e) paste("frm ERROR", substr(conditionMessage(e), 1, 120)))
      cat(sprintf("%-14s %-8s pkg=%s lpdf topenv=%s old_home=%s | %s\n",
                  ctor, nm, frmtmb::frm_family_package(fam), env,
                  if (is.null(home)) "NULL" else home, fit))
    }
  }
  quit(save = "no")
}
.libPaths(c(strsplit(a[3], ";", fixed = TRUE)[[1]], user))
suppressMessages(library(frmtmb))
att <- strsplit(a[4], ";", fixed = TRUE)[[1]]
for (p in setdiff(att, "none")) suppressMessages(library(p, character.only = TRUE))
och <- getFromNamespace("old_core_family_home", "frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), "| attached:", a[4], "\n")
for (nm in c("fcar", "fcox", "fmvn", "fwien", "fmvn-074")) {
  f <- file.path(a[2], paste0(nm, ".rds"))
  if (!file.exists(f)) next
  s <- readRDS(f); fit <- s$fit
  fam <- fit$spec$responses[[1]]$family
  h <- och(fam)
  cat("==", nm, "| family", substr(fam$family, 1, 30), "| pkg",
      frmtmb::frm_family_package(fam), "| old_core_family_home:",
      if (is.null(h)) "NULL" else h, "\n")
  calls <- list(
    print = quote(utils::capture.output(print(fit))),
    summary = quote(utils::capture.output(summary(fit))),
    logLik = quote(logLik(fit)),
    vcov = quote(vcov(fit)),
    fitted = quote(fitted(fit)),
    simulate = quote(simulate(fit, nsim = 1, seed = 1)),
    update = quote(update(fit, data = s$data)))
  for (k in names(calls)) {
    r <- tryCatch({ eval(calls[[k]]); "ok" }, error = function(e) {
      kk <- if (inherits(e, "frmtmb_extension_missing")) "REFUSED[ext]" else "ERR"
      paste(kk, substr(gsub("\\s+", " ", conditionMessage(e)), 1, 130))
    })
    cat(sprintf("  %-10s %s\n", k, r))
  }
}
