# Reviewer of lane defects: second behavioral probe (follow-ups of
# dev/defects-rev-probe.R). Run per arm:
#   Rscript dev/defects-rev-probe2.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
.libPaths(libs)
suppressMessages(library(frmtmb))
hdr <- function(...) cat("\n=====", ..., "\n")
try_msg <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [warnings:", paste(unique(substr(w, 1, 120)), collapse = " | "), "]\n")
  v
}
say <- function(label, v) {
  if (is.character(v) && length(v) == 1) cat(" ", label, ":", substr(v, 1, 300), "\n")
  else { cat(" ", label, ":\n"); print(v) }
}

hdr("A. residuals of an mi() response at its missing rows")
set.seed(21)
n <- 100
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 1 + d$x + rnorm(n)
d$xmi <- ifelse(runif(n) < 0.2, NA, d$x)
miss <- is.na(d$xmi)
fmi <- try_msg(frm(bf(y ~ mi(xmi)) + bf(xmi | mi() ~ z) + set_rescor(FALSE), data = d))
if (!is.character(fmi)) {
  say("frame$y[['xmi']] at missing rows (first 5)", head(fmi$frame$y[["xmi"]][miss], 5))
  r <- try_msg(residuals(fmi, resp = "xmi"))
  say("residual Estimate at missing rows (first 5)", if (is.character(r)) r else head(r[miss, "Estimate"], 5))
  say("fitted at missing rows (first 5)", head(fitted(fmi, resp = "xmi")[miss, "Estimate"], 5))
  say("estimates names", names(fmi$estimates))
}
d$ymi <- ifelse(runif(n) < 0.2, NA, d$y)
fu <- try_msg(frm(ymi | mi() ~ x, data = d))
say("univariate y | mi() fit", if (is.character(fu)) fu else "fitted")
if (!is.character(fu)) {
  say("nobs", nobs(fu))
  ru <- try_msg(residuals(fu))
  say("univariate residuals dim / NA count", if (is.character(ru)) ru else c(dim(ru), sum(is.na(ru[, 1]))))
}

hdr("B. cs() layers: identity with thresholds, and a delta-method SE")
set.seed(11)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
lat <- 0.6 * d$x + rnorm(n)
d$y <- cut(lat, c(-Inf, -0.7, 0.2, 1, Inf), labels = FALSE)
fit <- frm(y ~ x + cs(z), family = sratio(), data = d)
if (identical(arm, "lane")) {
  L <- fitted(fit, scale = "linear")
  fe <- fixef(fit)
  say("fixef", fe[, 1:2])
  tau <- fe[c("Intercept[1]", "Intercept[2]", "Intercept[3]"), "Estimate"]
  cs <- fe[c("z[1]", "z[2]", "z[3]"), "Estimate"]
  bx <- fe["x", "Estimate"]
  eta_hand <- sapply(1:3, function(k) bx * d$x + cs[k] * d$z)
  say("max |layer - (bx x + cs_k z)|", max(abs(L[, "Estimate", ] - eta_hand)))
  P <- fitted(fit)
  p1 <- plogis(tau[1] - eta_hand[, 1])
  p2 <- (1 - p1) * plogis(tau[2] - eta_hand[, 2])
  say("max |P(Y=1) - F(tau1 - eta1)|, |P(Y=2) - ...|", c(max(abs(P[, "Estimate", 1] - p1)), max(abs(P[, "Estimate", 2] - p2))))
  # delta method from vcov of fixef: eta_k = bx x + cs_k z
  V <- try_msg(vcov(fit))
  say("vcov names", if (is.character(V)) V else rownames(V))
  if (!is.character(V)) {
    nm <- rownames(V)
    ix <- which(nm == "x"); iz <- match(c("z[1]", "z[2]", "z[3]"), nm)
    if (length(ix) == 1 && !anyNA(iz)) {
      se_hand <- sapply(1:3, function(k) {
        g <- cbind(d$x, d$z)
        idx <- c(ix, iz[k])
        sqrt(rowSums((g %*% V[idx, idx]) * g))
      })
      say("max relative |Est.Error - delta SE|", max(abs(L[, "Est.Error", ] - se_hand) / se_hand))
    }
  }
}

hdr("C. residuals(newdata =) on log(y) ~ x")
set.seed(22)
d <- data.frame(x = rnorm(60)); d$pos <- exp(0.2 + 0.3 * d$x + rnorm(60, 0, 0.4))
fl <- frm(log(pos) ~ x, data = d)
if (identical(arm, "lane")) {
  r0 <- residuals(fl); r1 <- try_msg(residuals(fl, newdata = d))
  say("log(pos): newdata = d equals in-sample", if (is.character(r1)) r1 else max(abs(r0 - r1)))
  say("names(fl$data)", names(fl$data))
}

hdr("D. old_levels with correlated slopes, dpar and by")
set.seed(31)
n <- 300
d <- data.frame(g = factor(rep(1:10, 30)), x = rnorm(n))
u <- c(-6, -4, -2, 0, 1, 2, 3, 4, 5, 6)
d$y <- u[d$g] + (u[d$g] / 3) * d$x + rnorm(n, 0, 0.3)
f2 <- try_msg(frm(y ~ x + (1 + x | g), data = d))
if (!is.character(f2) && identical(arm, "lane")) {
  re <- ranef(f2)$g
  fe <- fixef(f2)[, "Estimate"]
  nd2 <- data.frame(g = factor(c("A", "A")), x = c(0, 1))
  set.seed(5)
  pp <- predict(f2, newdata = nd2, allow_new_levels = TRUE,
                sample_new_levels = "old_levels", ndraws = 3000,
                summary = FALSE)
  icp <- fe[1] + re[, 1]; slp <- fe[2] + re[, 2]
  lev <- which.min(abs(icp - mean(pp[, 1])))
  say("picked level, its intercept and slope", c(lev, icp[lev], slp[lev]))
  say("mean at x=0 and slope implied", c(mean(pp[, 1]), mean(pp[, 2]) - mean(pp[, 1])))
}
d$grp <- factor(ifelse(as.integer(d$g) <= 5, "lo", "hi"))
d$y <- u[d$g] + rnorm(n, 0, 0.3)
f4 <- try_msg(frm(y ~ 1 + (1 | gr(g, by = grp)), data = d))
if (!is.character(f4) && identical(arm, "lane")) {
  nd4 <- data.frame(g = factor("A"), grp = factor("lo", levels = c("hi", "lo")))
  set.seed(9)
  mm <- replicate(30, mean(try_msg(predict(f4, newdata = nd4, allow_new_levels = TRUE,
    sample_new_levels = "old_levels", ndraws = 30, summary = FALSE))))
  say("by = grp, new level in 'lo' (lo levels have means -6..1): picked means", round(sort(unique(round(mm))), 1))
}
d$y <- u[d$g] + rnorm(n, 0, exp(0.3 * (as.integer(d$g) - 5) / 5))
f3 <- try_msg(frm(bf(y ~ 1 + (1 | g), sigma ~ 1 + (1 | g)), data = d))
if (!is.character(f3) && identical(arm, "lane")) {
  nd <- data.frame(g = factor(c("A", "B")))
  say("dpar re old_levels dim", dim(try_msg(predict(f3, newdata = nd, allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 50, summary = FALSE))))
}
# multi-membership and a nested term
d$g2 <- factor(rep(1:5, 60))
f5 <- try_msg(frm(y ~ 1 + (1 | g) + (1 | g2), data = d))
if (!is.character(f5) && identical(arm, "lane")) {
  nd5 <- data.frame(g = factor(c("A", "1")), g2 = factor(c("1", "Z")))
  say("two factors, new level in each", try_msg(predict(f5, newdata = nd5, allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 200)))
}
f6 <- try_msg(frm(y ~ 1 + (1 | mm(g, g2)), data = d))
if (!is.character(f6) && identical(arm, "lane")) {
  nd6 <- data.frame(g = factor("A"), g2 = factor("B"))
  say("mm() new levels old_levels", try_msg(predict(f6, newdata = nd6, allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 50)))
}
fg <- try_msg(frm(y ~ 1 + gp(x) + (1 | g), data = d))
if (!is.character(fg) && identical(arm, "lane")) {
  say("gp + new level old_levels", try_msg(predict(fg, newdata = data.frame(x = 0, g = factor("A")), allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 50)))
}
fc <- try_msg(frm(y ~ 1 + cs(1 | g), data = d))
if (!is.character(fc) && identical(arm, "lane")) {
  say("cs(1|g) covariance + old_levels", try_msg(predict(fc, newdata = data.frame(g = factor("A")), allow_new_levels = TRUE, sample_new_levels = "old_levels", ndraws = 50)))
}

hdr("E. ranef, data_name, gp summary, categorical, spellings")
set.seed(91)
n <- 200
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, 20)), h = factor(rep(1:5, 40)))
d$y <- rnorm(10)[d$g] + rnorm(5)[d$h] + d$x + rnorm(n)
f <- frm(bf(y ~ x + (1 + x | g) + (1 | h), sigma ~ (1 | g)), data = d)
if (identical(arm, "lane")) {
  r <- ranef(f)
  say("names / colnames g", list(names(r), colnames(r$g)))
  say("pars = 'x'", lapply(ranef(f, pars = "x"), colnames))
  say("pars = 'Intercept'", lapply(ranef(f, pars = "Intercept"), colnames))
  say("pars = 'sigma_Intercept'", lapply(ranef(f, pars = "sigma_Intercept"), colnames))
  say("groups = 'h'", names(ranef(f, groups = "h")))
  say("groups = 'zz'", length(ranef(f, groups = "zz")))
  r5 <- ranef(f, pars = "x", condVar = TRUE)
  say("condSD dim", dim(attr(r5$g, "condSD")))
  say("pars = 1", try_msg(ranef(f, pars = 1)))
}
dd <- data.frame(y = rnorm(40), x = rnorm(40))
f <- frm(y ~ x, data = dd)
if (identical(arm, "lane")) {
  f2 <- frm(y ~ x, data = dd[1:30, ])
  new_data <- dd[5:40, ]
  say("names: plain, call, update(newdata), update(data)",
      c(attr(f$data, "data_name"), attr(f2$data, "data_name"),
        attr(update(f, newdata = new_data)$data, "data_name"),
        attr(update(f, data = new_data)$data, "data_name")))
  f3 <- do.call(frm, list(y ~ x, data = dd))
  say("do.call by value", is.null(attr(f3$data, "data_name")))
  say("summary Data lines", c(grep("Data", capture.output(summary(f3)), value = TRUE), grep("Data", capture.output(summary(f2)), value = TRUE)))
  g <- function(dat) frm(y ~ x, data = dat)
  say("inside a function", attr(g(dd)$data, "data_name"))
}
