# Reviewer of lane setier: dev/setier-singular.R's rs20 design
# (y ~ x + (1 + x | g), slope sd 0), seeds 1..100. For every fit nlminb
# ends at a code other than 0: logLik against lme4, the largest gradient,
# and what the standard-error check would say if it ran (its analysis,
# asked directly), so item 8's recommendation rests on numbers.
#   Rscript dev/setier-rev-rs20.R <lib> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:100
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
X <- list()
for (s in seeds) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(120))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g] + rnorm(120)
  f <- suppressMessages(suppressWarnings(frm(y ~ x + (1 + x | g), data = d)))
  m <- suppressMessages(suppressWarnings(lmer(y ~ x + (1 + x | g), data = d,
                                              REML = FALSE)))
  if (f$opt$convergence == 0) next
  g <- f$obj$gr(f$opt$par)
  an <- tryCatch({
    h <- ns$fit_outer_hessian(f)
    a <- ns$cov_from_hessian(f, h$H, h$E)
    ns$se_relabel(f, a$lost)
  }, error = function(e) NULL)
  X[[length(X) + 1L]] <- data.frame(
    seed = s, code = f$opt$convergence, msg = f$opt$message,
    dll = as.numeric(logLik(f)) - as.numeric(logLik(m)),
    lme4_singular = isSingular(m), maxgrad = max(abs(g)),
    theta = paste(signif(f$opt$par[grepl("^theta", names(f$opt$par))], 4),
                  collapse = ","),
    se_check = if (is.null(an)) "error" else if (!length(an)) "none" else
      paste(names(an), an, sep = ":", collapse = ","))
}
X <- do.call(rbind, X)
print(table(X$msg))
cat("code != 0 fits:", nrow(X), "| lme4 singular among them:",
    sum(X$lme4_singular), "\n")
cat("logLik frmtmb - lme4: min", signif(min(X$dll), 3), "max",
    signif(max(X$dll), 3), "| below lme4 by more than 1e-4:",
    sum(X$dll < -1e-4), "\n")
cat("max |gradient|: median", signif(median(X$maxgrad), 3), "max",
    signif(max(X$maxgrad), 3), "\n")
print(table(X$se_check))
