# Reviewer, fix round 1: the df check, the knownsd() example, and the
# compat statuses against rellib-r12.
#   Rscript dev/review-bhaz-fix1.R lane > dev/review-bhaz-log/fix1-lane.txt
#   Rscript dev/review-bhaz-fix1.R base > dev/review-bhaz-log/fix1-base.txt
arm <- commandArgs(trailingOnly = TRUE)[[1L]]
libs <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-bhaz-lib", libs)
.libPaths(libs)
exts <- c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
          "frmtmb.ode", "frmtmb.sample", "frmtmb.spatial", "frmtmb.survival")
suppressPackageStartupMessages(for (p in exts) library(p,
                                                       character.only = TRUE))
cat("lib:", find.package("frmtmb"), find.package("frmtmb.survival"), "\n")
cc <- frm_compat_rules()
key <- paste(cc$feature_a, cc$feature_b, sep = " x ")
saveRDS(data.frame(key = key, status = cc$status),
        sprintf("C:/Users/adf44/source/r/frmtmb-wt-bhaz/dev/review-bhaz-log/compat-%s.rds",
                arm))
cat("compat rules:", nrow(cc), "\n")
if (arm == "base") quit(save = "no")

cat("\n== finding 5 ==\n")
set.seed(2); n <- 300; x <- rnorm(n)
tt <- stats::rexp(n, exp(-0.5 + 0.7 * x)); ct <- stats::rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
q <- function(e) tryCatch(e, error = function(e) structure(conditionMessage(e),
                                                       class = "err"))
a <- q(frm(bf(time | cens(cens) + bhaz(intercept = FALSE) ~ x),
           family = cox(df = 3), data = dd))
b <- q(frm(bf(time | cens(cens) ~ x), family = cox(df = 3, intercept = FALSE),
           data = dd))
cat("cox(df = 3) + bhaz(intercept = FALSE):",
    if (inherits(a, "err")) a else sprintf(
      "K %d, identical logLik %s estimates %s vcov %s", length(cox_baseline(a)),
      identical(logLik(a), logLik(b)), identical(a$estimates, b$estimates),
      identical(vcov(a), vcov(b))), "\n")
a2 <- q(frm(bf(time | cens(cens) + bhaz(degree = 2) ~ x), family = cox(df = 3),
            data = dd))
b2 <- q(frm(bf(time | cens(cens) ~ x), family = cox(df = 3, degree = 2),
            data = dd))
cat("cox(df = 3) + bhaz(degree = 2): identical logLik",
    identical(logLik(a2), logLik(b2)), "\n")
f3 <- q(cox(df = 3))
cat("cox(df = 3) constructs:", !inherits(f3, "err"), "\n")
calls <- 0L
tr <- function(...) calls <<- calls + 1L
r <- q(frm(bf(time | cens(cens) ~ x), family = cox(df = 3), data = dd,
           dry_run = "frame"))
cat("cox(df = 3), dry_run = 'frame':", if (inherits(r, "err")) r else "NO ERROR",
    "\n")
r <- q(frm(bf(time | cens(cens) ~ x), family = cox(df = 3), data = NULL,
           dry_run = "spec"))
cat("cox(df = 3), dry_run = 'spec':", if (inherits(r, "err")) r else "no error",
    "\n")
t0 <- proc.time()[[3]]
r <- q(frm(bf(time | cens(cens) ~ x), family = cox(df = 3), data = dd))
cat("cox(df = 3) fit:", if (inherits(r, "err")) r else "FIT", " in",
    round(proc.time()[[3]] - t0, 2), "s\n")
r <- q(frm(bf(time | cens(cens) + bhaz(df = 3) ~ x), family = cox(),
           data = NULL, dry_run = "spec"))
cat("bhaz(df = 3), dry_run = 'spec':", if (inherits(r, "err")) r else
  "no error", "\n")
r <- q(frm(bf(time | cens(cens) + bhaz(degree = 4) ~ x), family = cox(df = 4),
           data = NULL, dry_run = "spec"))
cat("cox(df = 4) + bhaz(degree = 4):", if (inherits(r, "err")) r else
  "no error", "\n")

cat("\n== ?frmtmb_register_aterm example ==\n")
frmtmb_local_registries()
ex <- utils::example("frmtmb_register_aterm", package = "frmtmb",
                     give.lines = TRUE)
cat(ex, sep = "\n")
r <- q(utils::example("frmtmb_register_aterm", package = "frmtmb",
                      ask = FALSE, echo = FALSE))
cat("example() ran:", !inherits(r, "err"), if (inherits(r, "err")) r, "\n")
if (exists("fit")) {
  ref <- sum(stats::dnorm(d$y, fitted(lm(y ~ x, d)), 2, log = TRUE))
  cat(sprintf("example fit logLik %.10f, independent reference %.10f\n",
              as.numeric(logLik(fit)), ref))
}
