# Reviewer: the BREAKING trunc() checks on the cases brms accepts, and
# what brms does on each (standata(), no compile). ARM=lane or base.
arm <- Sys.getenv("ARM", "lane")
lp <- c("C:/Users/adf44/source/r/rellib-r10",
        "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") lp <- c("C:/Users/adf44/source/r/wt-cdfs-rev-lib", lp)
.libPaths(lp)
suppressPackageStartupMessages({
  library(frmtmb)
  library(brms)
})
options(frmtmb.notices = FALSE)
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(5)
n <- 60
x <- rnorm(n)
yc <- round(1 + 0.5 * x + rnorm(n), 3)
yp <- rpois(n, exp(1 + 0.3 * x))
d <- data.frame(x = x, yc = yc, yp = yp)
d$lbc <- min(yc)              # a continuous response ON its lower bound
d$ubc <- max(yc)              # and on its upper bound
d$yp1 <- pmax(yp, 1)          # counts >= 1, some exactly 1
d$lbr <- ifelse(seq_len(n) %% 3 == 0, -Inf, d$yc - 0.5)  # per-row bounds
d$ubr <- ifelse(seq_len(n) %% 4 == 0, Inf, d$yc + 0.5)
d$lbp <- pmin(d$yp1, 2)       # per-row discrete lb, y == lb on many rows
d$ubp <- d$yp1 + ifelse(seq_len(n) %% 2 == 0, 0, 3)   # y == ub on half
d$lbna <- ifelse(seq_len(n) == 7, NA, d$yc - 1)
d$tr <- 12
d$yb <- rbinom(n, 12, 0.4)
d$yb1 <- pmax(d$yb, 1)
d$cc <- ifelse(seq_len(n) %% 5 == 0, 1, 0)
run <- function(lab, f, fam, data = d) {
  fr <- tryCatch({
    fit <- suppressWarnings(frm(f, family = fam, data = data))
    sprintf("fits, logLik %.4f", as.numeric(logLik(fit)))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 110)))
  br <- tryCatch({
    brms::standata(f, data = data, family = fam)
    "accepts"
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 70)))
  cat(sprintf("%-44s frmtmb: %s\n%-44s brms:   %s\n", lab, fr, "", br))
}
run("gaussian y on lb (lb = min y)", yc | trunc(lb = lbc) ~ x, gaussian())
run("gaussian y on ub (ub = max y)", yc | trunc(ub = ubc) ~ x, gaussian())
run("gaussian lb = -Inf, ub = Inf", yc | trunc(lb = -Inf, ub = Inf) ~ x, gaussian())
run("gaussian lb = -Inf, ub = max", yc | trunc(lb = -Inf, ub = ubc) ~ x, gaussian())
run("gaussian per-row lb/ub with +-Inf rows", yc | trunc(lb = lbr, ub = ubr) ~ x,
    gaussian())
run("gaussian per-row lb with an NA", yc | trunc(lb = lbna) ~ x, gaussian())
run("student y on lb", yc | trunc(lb = lbc) ~ x, student())
run("poisson y >= 1, trunc(lb = 1), y == lb", yp1 | trunc(lb = 1) ~ x, poisson())
run("poisson trunc(lb = 0) (y has 0)", yp | trunc(lb = 0) ~ x, poisson())
run("poisson trunc(lb = 0, ub = max)", yp | trunc(lb = 0, ub = 30) ~ x, poisson())
run("poisson trunc(ub = max y), y == ub", yp | trunc(ub = max(yp)) ~ x, poisson())
run("poisson lb = -Inf, ub = 30", yp | trunc(lb = -Inf, ub = 30) ~ x, poisson())
run("negbinomial per-row lb/ub, y == lb, y == ub",
    yp1 | trunc(lb = lbp, ub = ubp) ~ x, negbinomial())
run("binomial trunc(lb = 1, ub = 12), y == 1", yb1 | trials(tr) + trunc(lb = 1, ub = 12) ~ x,
    binomial())
run("binomial trunc(ub = 12) (y has 0)", yb | trials(tr) + trunc(ub = 12) ~ x,
    binomial())
run("zero_inflated_poisson trunc(ub = 30)", yp | trunc(ub = 30) ~ x,
    zero_inflated_poisson())
run("hurdle_poisson trunc(lb = 1) (zeros present)", yp | trunc(lb = 1) ~ x,
    hurdle_poisson())
run("gaussian cens + trunc, y == ub censored", yc | cens(cc) + trunc(ub = ubc) ~ x,
    gaussian())
run("gaussian y outside (ub below max)", yc | trunc(ub = 1) ~ x, gaussian())
run("gaussian lb >= ub", yc | trunc(lb = 2, ub = 2) ~ x, gaussian())
d2 <- d
d2$yc[3] <- NA
run("gaussian NA response row, y on lb", yc | trunc(lb = lbc) ~ x, gaussian(), d2)
cat("== multivariate: one response truncated\n")
r <- tryCatch({
  fit <- suppressWarnings(frm(bf(yc | trunc(lb = lbc) ~ x) + bf(yp ~ x, family = poisson()),
                              data = d, family = gaussian()))
  sprintf("fits, logLik %.4f", as.numeric(logLik(fit)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("mvbf:", r, "\n")
r <- tryCatch({
  fit <- suppressWarnings(frm(bf(yc | trunc(ub = 1) ~ x) + bf(yp ~ x, family = poisson()),
                              data = d, family = gaussian()))
  sprintf("fits, logLik %.4f", as.numeric(logLik(fit)))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("mvbf, y outside ub:", r, "\n")
