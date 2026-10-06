# Defect 9: frm_allfit() on nonlinear fits that need start =, and on a
# linear control where every optimizer should agree. Same constructions
# as dev/vigport-rev2-allfit.R in the vigport worktree (fit_loss) plus
# brms_nonlinear's first model (b1 * exp(b2 * x)) and a flat-start
# control. Prints the frm_allfit() table of the build under test, and
# for the lane build both start_from_mle settings.
#   Rscript dev/nanse-allfit.R [lib]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
lane <- "start_from_mle" %in% names(formals(frm_allfit))
show <- function(tag, f, ...) {
  cat("\n==", tag, "\n")
  a <- tryCatch(frm_allfit(f, ...), error = function(e) {
    cat("frm_allfit ERROR:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(a)) return(invisible())
  print(a)
  if (!lane) {
    for (nm in names(a$fits)) {
      g <- a$fits[[nm]]
      cat(sprintf("  %-13s %s\n", nm, if (is.null(g)) "NULL (refit failed)"
                  else sprintf("logLik %.4f  code %s",
                               as.numeric(logLik(g)), g$opt$convergence)))
    }
  }
  invisible(a)
}
loss <- brms::loss
fit_loss <- frm(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
                   ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
                data = loss, family = gaussian(),
                start = list(beta = c(5000, 1, 45)))
cat("fit_loss: code", fit_loss$opt$convergence, " logLik",
    format(as.numeric(logLik(fit_loss)), digits = 10), "\n")
show("nonlinear fit_loss (default)", fit_loss)
if (lane) show("nonlinear fit_loss, start_from_mle = FALSE", fit_loss,
               start_from_mle = FALSE)
# brms_nonlinear's first model
set.seed(1234)
b <- c(2, 0.75)
x <- rnorm(100)
y <- rnorm(100, mean = b[1] * exp(b[2] * x))
dat1 <- data.frame(x, y)
fit1 <- frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE), data = dat1,
            start = list(beta = c(1, 0.5)))
cat("\nfit1: code", fit1$opt$convergence, " logLik",
    format(as.numeric(logLik(fit1)), digits = 10), "\n")
show("nonlinear b1 * exp(b2 * x) (default)", fit1)
if (lane) show("nonlinear b1 * exp(b2 * x), start_from_mle = FALSE", fit1,
               start_from_mle = FALSE)
# control: a linear mixed model, where the default start is fine
d <- brms::epilepsy
fl <- frm(count ~ zAge + zBase * Trt + (1 | patient), data = d,
          family = poisson())
show("linear control (epilepsy poisson)", fl)
if (lane) show("linear control, start_from_mle = FALSE", fl,
               start_from_mle = FALSE)
