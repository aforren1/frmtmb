# Reviewer of lane ndaterm, fix round 1: predict() dpars on newdata and
# under re_formula, against fitted() and simulate(), for every family
# whose reporting scale can differ from the natural dpar scale.
# Usage: Rscript dev/ndaterm-rev-r1dpar.R <arm>   (mine | r9). Seed 7.
# Log: dev/ndaterm-rev-out/r1dpar-<arm>.txt (redirected by the caller)
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste("ERROR:", conditionMessage(e)), class = "err")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [", label, "] WARN:", substr(unique(w), 1, 120),
                     "\n")
  if (inherits(r, "err")) cat("  [", label, "]", substr(r, 1, 300), "\n")
  invisible(r)
}
ok <- function(r) !inherits(r, "err") && !is.null(r)
show <- function(label, x) {
  cat(sprintf("%-56s %s\n", label,
              paste(format(x, digits = 4), collapse = "  ")))
}
# predictive means, propagate_error TRUE and FALSE, with their MC se,
# against fitted() and simulate()
cmp <- function(label, fit, nd, re_formula = NULL, nsim = 4000) {
  f <- tr(paste(label, "fitted"),
          fitted(fit, newdata = nd, re_formula = re_formula))
  for (pe in c(TRUE, FALSE)) {
    set.seed(11)
    p <- tr(paste(label, "predict", pe),
            predict(fit, newdata = nd, re_formula = re_formula,
                    summary = FALSE, ndraws = nsim, propagate_error = pe))
    if (ok(p)) {
      show(paste0(label, ": predict mean (pe = ", pe, ")"), colMeans(p))
      show(paste0(label, ":   6 MC se"), 6 * apply(p, 2, sd) / sqrt(nsim))
    }
  }
  if (ok(f)) show(paste0(label, ": fitted"), f[, "Estimate"])
  if (!is.null(nd) && is.null(re_formula)) {
    s <- tr(paste(label, "simulate"),
            as.matrix(simulate(fit, nsim = nsim, newdata = nd, seed = 11)))
    if (ok(s)) show(paste0(label, ": simulate(newdata) mean"), rowMeans(s))
  }
}

set.seed(7)
n <- 300
d <- data.frame(x = rnorm(n), g = gl(15, 20))
k <- rbinom(n, 1, plogis(-1.4 + 0.8 * d$x))
d$y <- ifelse(k == 1, 6, 0) + d$x + rnorm(15, 0, 0.5)[d$g] +
  rnorm(n, 0, 0.5)
nd <- data.frame(x = c(-1, 0, 1), g = "1")

cat("\n== 1. mixture(gaussian, gaussian), theta1 the predicted weight\n")
f1 <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = d)
cmp("mix", f1, nd)
cat("\n== 2. the same, theta2 ~ x (component 1 the reference)\n")
f2 <- tr("theta2 fit", frm(bf(y ~ x, theta2 ~ x),
                          family = mixture(gaussian(), gaussian()),
                          data = d))
if (ok(f2)) {
  cat("dpars:", names(f2$spec$responses[[1]]$dpars), "\n")
  cmp("mix theta2~x", f2, nd)
}
cat("\n== 3. mixture with (1 | g): newdata, and re_formula = NA in sample\n")
f3 <- tr("re fit", frm(bf(y ~ x + (1 | g)),
                       family = mixture(gaussian(), gaussian()), data = d))
if (ok(f3)) {
  cmp("mix RE newdata", f3, nd)
  cmp("mix RE newdata re_formula NA", f3, nd, re_formula = NA)
  f <- fitted(f3, re_formula = NA)[, "Estimate"]
  set.seed(12)
  p <- predict(f3, re_formula = NA, summary = FALSE, ndraws = 2000)
  show("mix RE in sample re_formula NA: mean |predict - fitted|",
       mean(abs(colMeans(p) - f)))
  show("  typical 6 MC se", 6 * mean(apply(p, 2, sd)) / sqrt(2000))
}
cat("\n== 4. ordinal mixture(cumulative, cumulative)\n")
set.seed(7)
do <- data.frame(x = rnorm(400))
ko <- rbinom(400, 1, 0.3)
do$yo <- cut(d$x[1] * 0 + do$x + ifelse(ko == 1, 2, 0) + rlogis(400),
             c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
f4 <- tr("ordinal mixture fit", frm(bf(yo ~ x), family =
                                     mixture(cumulative(), cumulative()),
                                   data = do))
if (ok(f4)) {
  ndo <- data.frame(x = c(-1, 0, 1))
  pr <- fitted(f4, newdata = ndo)
  cat("fitted dims:", dim(pr), "\n")
  pm <- if (length(dim(pr)) == 3) pr[, "Estimate", ] else pr
  ec <- as.numeric(pm %*% seq_len(ncol(pm)))
  for (pe in c(TRUE, FALSE)) {
    set.seed(13)
    p <- tr("ord mix predict", predict(f4, newdata = ndo, summary = FALSE,
                                       ndraws = 4000, propagate_error = pe))
    if (ok(p)) {
      pn <- matrix(as.numeric(p), nrow(p))
      show(paste0("ord mix: predict mean category (pe = ", pe, ")"),
           colMeans(pn))
      show("ord mix:   6 MC se", 6 * apply(pn, 2, sd) / sqrt(4000))
    }
  }
  show("ord mix: fitted expected category", ec)
}
cat("\n== 5. categorical\n")
set.seed(7)
dc <- data.frame(x = rnorm(400))
eta <- cbind(0, 0.5 + 0.8 * dc$x, -0.3 - 0.6 * dc$x)
pc <- exp(eta) / rowSums(exp(eta))
dc$yc <- apply(pc, 1, function(p) sample(3, 1, prob = p))
f5 <- frm(bf(yc ~ x), family = categorical(), data = dc)
ndc <- data.frame(x = c(-1, 0, 1))
fp <- fitted(f5, newdata = ndc)
pm <- if (length(dim(fp)) == 3) fp[, "Estimate", ] else fp
for (pe in c(TRUE, FALSE)) {
  set.seed(14)
  p <- predict(f5, newdata = ndc, summary = FALSE, ndraws = 4000,
               propagate_error = pe)
  pn <- matrix(as.numeric(p), nrow(p))
  show(paste0("categorical: share of category 2 (pe = ", pe, ")"),
       colMeans(pn == 2))
}
show("categorical: fitted P(category 2)", pm[, 2])
cat("\n== 6. gaussian se(s) without sigma = TRUE (sigma reports 0)\n")
set.seed(7)
dg <- data.frame(x = rnorm(200), s = runif(200, 0.2, 1))
dg$yg <- 1 + dg$x + rnorm(200, 0, dg$s)
f6 <- frm(bf(yg | se(s) ~ x), family = gaussian(), data = dg)
set.seed(15)
p <- predict(f6, newdata = data.frame(x = 0, s = c(0.05, 5)),
             summary = FALSE, ndraws = 4000)
show("se(s): predict sd at s = 0.05, 5", apply(p, 2, sd))
cat("\n== 7. nonlinear dpar and a constant dpar\n")
set.seed(7)
dn <- data.frame(x = runif(200, 0, 4))
dn$yn <- 3 * exp(-0.5 * dn$x) + rnorm(200, 0, 0.2)
f7 <- tr("nl fit", frm(bf(yn ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE),
                       family = gaussian(), data = dn))
if (ok(f7)) cmp("nl", f7, data.frame(x = c(0, 2)))
f8 <- tr("constant sigma fit", frm(bf(yn ~ x, sigma = 0.3),
                                   family = gaussian(), data = dn))
if (ok(f8)) {
  set.seed(16)
  p <- predict(f8, newdata = data.frame(x = c(0, 2)), summary = FALSE,
               ndraws = 4000)
  show("sigma = 0.3 constant: predict sd (want ~0.3)", apply(p, 2, sd))
}
cat("done\n")
