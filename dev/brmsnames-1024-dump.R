## Raw dump of every case: brms's names beside each frmtmb arm's.
##   Rscript dev/brmsnames-1024-dump.R > dev/brmsnames-1024-log/dump.txt
out <- "dev/brmsnames-1024-out"
P <- "^(b|bs|bsp|bcs|sds|sd|simo|cor)_"
for (i in seq_along(list.files(out, "^brms-"))) {
  b <- readRDS(sprintf("%s/brms-%02d.rds", out, i))
  cat("=====", i, b$id, "\n")
  cat("brms X       :", unlist(b$X), "\n")
  cat("brms vars    :", grep(P, b$variables, value = TRUE), "\n")
  cat("brms fixef   :", b$fixef, "\n")
  cat("brms coef    :", unlist(b$coef), "\n")
  cat("brms ce      :", b$ce, "\n")
  cat("brms ranef   :", unlist(b$ranef), "\n")
  cat("brms varcorr :", unlist(b$varcorr), "\n")
  cat("brms prior   :", paste0(b$prior$class, "/", b$prior$coef)[
    nzchar(b$prior$coef)], "\n")
  for (arm in c("base", "lane")) {
    f <- readRDS(sprintf("%s/frm-%s-%02d.rds", out, arm, i))
    if (!is.null(f$error)) {
      cat(arm, "ERROR:", f$error, "\n")
      next
    }
    cat(arm, "internal :", f$internal, "\n")
    cat(arm, "vars     :", grep(P, f$variables, value = TRUE), "\n")
    cat(arm, "fixef    :", f$fixef, "\n")
    cat(arm, "vcov     :", f$vcov, "\n")
    cat(arm, "summary  :", f$summary, "\n")
    cat(arm, "confint  :", f$confint, "\n")
    cat(arm, "coef     :", unlist(f$coef), "\n")
    cat(arm, "ce       :", f$ce, "\n")
    cat(arm, "ranef    :", unlist(f$ranef), "\n")
    cat(arm, "varcorr  :", unlist(f$varcorr), "\n")
    cat(arm, "prior    :", paste0(f$prior$class, "/", f$prior$coef)[
      nzchar(f$prior$coef)], "\n")
    cat(arm, "prior_ok :", f$prior_ok, "\n")
    cat(arm, "hyp      :", paste0(names(f$hyp), "=", f$hyp), "\n")
    cat(arm, "draws    :", grep(P, f$draws, value = TRUE), "\n")
  }
}
