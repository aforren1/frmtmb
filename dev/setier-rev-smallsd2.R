# Reviewer of lane setier: detail on lme4-singular fits of
# dev/setier-rev-smallsd.R that get no boundary message: warnings, theta
# SEs, the theta rows, the at-edge probe, logLik against lme4.
#   Rscript dev/setier-rev-smallsd2.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
src <- readLines("dev/setier-rev-smallsd.R")
eval(parse(text = src[grep("^gen <- function", src):
                        (grep("^rows <- list", src) - 1L)]))
for (cs in list(c("s1", 10), c("s1", 39), c("p1", 16), c("c95", 8),
                c("c99", 17), c("c99", 6))) {
  g <- gen(cs[1], as.integer(cs[2]))
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(g$ff, family = g$fam, data = g$d),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  lm4 <- suppressMessages(suppressWarnings(
    if (g$fam$family == "gaussian") lmer(g$ff, data = g$d, REML = FALSE)
    else glmer(g$ff, data = g$d, family = g$fam)))
  nm <- ns$outer_par_names(f)
  th <- grep("^theta", nm)
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  hc <- f$cache$hessian_fixed
  cat(sprintf("%s seed %s: code %d (%s) | logLik frmtmb %.6f lme4 %.6f\n",
              cs[1], cs[2], f$opt$convergence, f$opt$message,
              as.numeric(logLik(f)), as.numeric(logLik(lm4))))
  cat("   theta", signif(f$opt$par[th], 5), "SE", signif(se[th], 4), "\n")
  if (!is.null(hc)) {
    for (j in th) cat("   row", nm[j], "max", signif(max(abs(hc$H[j, ])), 3),
                      "diag", signif(hc$H[j, j], 3), "noise",
                      signif(max(hc$E[j, ]), 3),
                      " at_edge", if (exists("se_at_edge", ns))
                        ns$se_at_edge(f, nm[j], if (j == max(th) &&
                          length(th) == 3) sign(f$opt$par[j]) else -1)
                      else NA, "\n")
  }
  if (length(w)) cat("   W:", substr(w, 1, 160), sep = "\n   W: ")
  if (length(m)) cat("   M:", substr(m, 1, 100), sep = "\n   M: ")
  cat("\n")
}
