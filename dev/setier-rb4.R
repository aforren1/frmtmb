# Lane setier, punch round 1 (item 8): the RB4 designs of
# dev/setier-rev-a-rev3-falseloss-sweep.R (y ~ x + f + (1 + x | g2),
# 40-level f, no g2 variation, and y ~ x + (1 + x | g2)), seeds 1 to 20:
# per fit, the optimizer code, the largest gradient, what the user is
# told (convergence warning, boundary message, SE warning) and the lost
# reasons, against lme4's isSingular().
#   Rscript dev/setier-rb4.R <lib or "base">
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
tab <- NULL
for (design in c("f40", "plain")) {
  for (s in 1:20) {
    set.seed(s)
    n <- if (design == "f40") 480 else 120
    d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                    g2 = factor(rep(1:6, length.out = n)))
    d$y <- 1 + 0.5 * d$x + (if (design == "f40") rnorm(40, 0, 0.5)[d$f]
                            else 0) + rnorm(n)
    fo <- if (design == "f40") y ~ x + f + (1 + x | g2) else
      y ~ x + (1 + x | g2)
    w <- character()
    m <- character()
    fit <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    })
    s1 <- withCallingHandlers(summary(fit), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    })
    lost <- ns$sdr_of(fit)$se_lost
    lm4 <- suppressMessages(suppressWarnings(
      lme4::lmer(fo, data = d, REML = FALSE)))
    g <- max(abs(fit$obj$gr(fit$opt$par)))
    tab <- rbind(tab, data.frame(
      design, seed = s, code = fit$opt$convergence, grad = signif(g, 3),
      conv_warn = any(grepl("did not report convergence", w)),
      boundary_msg = any(grepl("^Boundary", m)),
      se_warn = any(grepl("Standard errors are not", w)),
      lost = paste(names(lost), lost, sep = ":", collapse = ","),
      singular = lme4::isSingular(lm4),
      dll = signif(as.numeric(logLik(fit)) - as.numeric(logLik(lm4)), 3)))
  }
}
print(tab, row.names = FALSE)
cat("\ncode != 0:", sum(tab$code != 0), "| of them boundary message:",
    sum(tab$code != 0 & tab$boundary_msg), "| convergence warning:",
    sum(tab$code != 0 & tab$conv_warn), "| silenced with a non-boundary",
    "loss:", sum(tab$code != 0 & !tab$conv_warn &
                   grepl("flat|concave|bound:|nonfinite", tab$lost)), "\n")
cat("lme4 singular:", sum(tab$singular), "| boundary message on them:",
    sum(tab$singular & tab$boundary_msg), "| on the others:",
    sum(!tab$singular & tab$boundary_msg), "\n")
