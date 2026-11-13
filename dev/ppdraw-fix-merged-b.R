# Fix round, merged build only (lane + boundfam): the two kinds of
# "limit" through predict(). Kind (a), a family parameter at its end
# (3269, zi ~ x, se(sigma = TRUE)): the side law draws it. Kind (b), a
# mean left undetermined by a fit with no maximum (boundfam's all-zero
# poisson and all-right-censored gaussian, seed 3): the fallback warning,
# and no draw from the raw row. Usage: Rscript dev/ppdraw-fix-merged-b.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), " se_limit_draw_held present:",
    exists("se_limit_draw_held", asNamespace("frmtmb")), "\n")
show <- function(label, fit, nd = NULL) {
  sdr <- frmtmb:::sdr_of(fit)
  lost <- sdr$se_lost
  cat("\n==", label, "\nse_lost:", if (length(lost)) {
    paste(names(lost), lost, collapse = "; ")
  } else "(none)", "\n")
  if (length(lost)) {
    cat("quiet:", frmtmb:::draw_limit_quiet(fit, lost), " held:",
        paste(names(lost), frmtmb:::draw_limit_held(fit, lost),
              collapse = "; "), "\n")
  }
  set.seed(1)
  said <- character(0)
  p <- withCallingHandlers(
    predict(fit, newdata = nd, ndraws = 400, summary = FALSE),
    message = function(m) {
      said <<- c(said, paste("M:", substr(conditionMessage(m), 1, 110)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      said <<- c(said, paste("W:", substr(conditionMessage(w), 1, 110)))
      invokeRestart("muffleWarning")
    })
  cat(unique(said), sep = "\n")
  cat(sprintf("predict mean %.4g sd of draw means %.4g\n", mean(p),
              stats::sd(rowMeans(p))))
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
set.seed(7)
d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
show("(a) 3269 zib", q(frm(bf(y | trials(N) ~ x),
                           family = zero_inflated_binomial(), data = d)))
show("(a) zib zi ~ x", q(frm(bf(y | trials(N) ~ x, zi ~ x),
                             family = zero_inflated_binomial(), data = d)))
set.seed(2)
de <- data.frame(x = stats::rnorm(50), y = 0L)
show("(b) poisson, all-zero response",
     q(frm(y ~ x, family = poisson(), data = de)))
set.seed(3)
dc <- data.frame(x = stats::rnorm(60))
dc$y <- 1 + dc$x + stats::rnorm(60)
dc$cen <- "right"
show("(b) gaussian, every row right-censored",
     q(frm(y | cens(cen) ~ x, data = dc)))
