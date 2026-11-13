# Review C: what lane boundfam's SE check makes of the lane's three fits,
# and what predict_par_drawer() then does (the merge question).
# Usage: Rscript ppdraw-review-c.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-ppdraw/dev/reviews/ppdraw-review-common.R")
cat("se_edge_limit present:",
    exists("se_edge_limit", envir = asNamespace("frmtmb")), "\n")
loud <- function(expr) withCallingHandlers(expr,
  message = function(m) {
    cat("  M[", class(m)[1], "]:", substr(conditionMessage(m), 1, 200), "\n")
    invokeRestart("muffleMessage")
  },
  warning = function(w) {
    cat("  W:", substr(conditionMessage(w), 1, 200), "\n")
    invokeRestart("muffleWarning")
  })
rep_fit <- function(label, fit, newdata = NULL) {
  cat("\n==", label, "\n")
  sdr <- frmtmb:::sdr_of(fit)
  cat("se_lost:", if (length(sdr$se_lost)) {
    paste(names(sdr$se_lost), sdr$se_lost, collapse = "; ")
  } else "(none)", "\n")
  ds <- frmtmb:::fit_draw_space(fit)
  cat("jc non-null:", !is.null(ds$jc), " boundary_only:",
      isTRUE(ds$boundary_only), "\n")
  set.seed(1)
  p <- loud(predict(fit, newdata = newdata, ndraws = 400, summary = FALSE))
  m <- rowMeans(p)
  fm <- mean(fitted(fit, newdata = newdata)[, "Estimate"])
  cat(sprintf("predict mean %.4f fitted %.4f sd of draw means %.4g NA %.4f\n",
              mean(m), fm, stats::sd(m), mean(is.na(p))))
}

set.seed(7)
d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
f1 <- loud(frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
               data = d))
rep_fit("3269 zib", f1)

set.seed(7)
dg <- data.frame(x = stats::rnorm(200), s = stats::runif(200, 0.2, 1))
dg$y <- 1 + dg$x + stats::rnorm(200, 0, dg$s)
f2 <- loud(frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(),
               data = dg))
rep_fit("se(sigma = TRUE)", f2,
        newdata = data.frame(x = c(0, 0), s = c(0.05, 5)))

set.seed(7)
do <- data.frame(x = stats::rnorm(400))
k <- stats::rbinom(400, 1, 0.3)
do$yo <- cut(do$x + ifelse(k == 1, 2, 0) + stats::rlogis(400),
             c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
f3 <- loud(frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
               data = do))
cat("\n== 3277 ordinal mixture\n")
sdr <- frmtmb:::sdr_of(f3)
cat("se_lost:", if (length(sdr$se_lost)) {
  paste(names(sdr$se_lost), sdr$se_lost, collapse = "; ")
} else "(none)", "\n")
set.seed(13)
p <- loud(predict(f3, newdata = data.frame(x = c(-1, 0, 1)), ndraws = 4000,
                  summary = FALSE))
cat("NA share", mean(is.na(p)), "\n")
