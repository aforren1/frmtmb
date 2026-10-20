# Reviewer: multinomial fitted() against brms over every cell with
# varying trials; drop1 and confint aliases; multinomial refcat
# consumers. Lane library plus the user library for brms.
# Log dev/catmulti-rev-out/consumers3.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
source("dev/catmulti-rev-data.R")
TRY <- function(nm, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err")),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") },
    message = function(cm) { cat("  MESSAGE:", conditionMessage(cm)); invokeRestart("muffleMessage") })
  cat(nm, ":\n"); if (inherits(r, "err")) cat("  ERROR:", r, "\n") else print(r)
  if (length(w)) cat("  WARNINGS:", substr(unique(w), 1, 200), sep = "\n    ")
  invisible(if (inherits(r, "err")) NULL else r)
}
cat("\n==== multinomial fitted vs brms, varying trials (seed 31) ====\n")
dm <- rev_mult_data(31)
fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dm)
bm <- readRDS("dev/catmulti-rev-out/brms-mult-vt.rds")
ff <- fitted(fm); fb <- suppressWarnings(brms:::fitted.brmsfit(bm))
z <- (ff[, "Estimate", ] - fb[, "Estimate", ]) / fb[, "Est.Error", ]
sr <- ff[, "Est.Error", ] / fb[, "Est.Error", ]
cat("dimnames equal:", identical(dimnames(ff), dimnames(fb)), "\n")
cat(sprintf("|z| median %.3f max %.3f; se ratio median %.3f range %.3f to %.3f over %d cells\n",
            median(abs(z)), max(abs(z)), median(sr), min(sr), max(sr), length(sr)))
r1 <- which(dm$N == 1)
cat(sprintf("rows with N = 1 (%d): se ratio median %.3f\n", length(r1), median(sr[r1, ])))
cat("estimates sum to trials (max |rowSums - N|):", max(abs(rowSums(ff[, "Estimate", ]) - dm$N)), "\n")
# brms predict vs frmtmb predict names
pb <- suppressWarnings(brms:::predict.brmsfit(bm)); pf <- predict(fm)
cat("predict dims brms", dim(pb), "frmtmb", dim(pf), "names equal:", identical(dimnames(pb), dimnames(pf)), "\n")

cat("\n==== multinomial refcat consumers ====\n")
fr <- frm(bf(Yc | trials(N) ~ x), family = multinomial(refcat = "y3"), data = dm)
cat("refcat y3 dpars:", names(fixef_by_dpar(fr)), "logLik", as.numeric(logLik(fr)),
    "default", as.numeric(logLik(fm)), "\n")
cat("max |fitted y3 - fitted default|:", max(abs(fitted(fr)[, "Estimate", ] - ff[, "Estimate", ])), "\n")
s <- simulate(fr, nsim = 50, seed = 3)
TRY("simulate class/dim", c(class(s), dim(as.matrix(s[[1]]))))
TRY("simulate first", head(s[[1]], 3))
fu <- TRY("update refcat NA", update(fr, family = multinomial(refcat = NA)))
if (!is.null(fu)) cat("update dpars:", names(fixef_by_dpar(fu)), "\n")
TRY("family(fit)$refcat", family(fr)$refcat)
TRY("print(multinomial(refcat = 'y2'))", print(multinomial(refcat = "y2")))
TRY("categorical(K = 3, refcat = '2') dpars", categorical(K = 3, refcat = "2")$dpars)
TRY("categorical(levels = c('x','y','z'), refcat = 'z') dpars",
    categorical(levels = c("x", "y", "z"), refcat = "z")$dpars)
TRY("categorical(refcat = c('a','b'))", categorical(refcat = c("a", "b")))
TRY("categorical(refcat = 'q') on a,b,c", frm(bf(y ~ x), family = categorical(refcat = "q"),
                                            data = rev_cat_data(21, 3)))
TRY("drop1 multinomial", drop1(fm, test = "Chisq"))

cat("\n==== drop1 ====\n")
d3 <- rev_cat_data(21, 3)
d3$z <- rnorm(nrow(d3))
f1 <- frm(bf(y ~ x + z, muc ~ x), family = categorical(refcat = "b"), data = d3)
tb <- TRY("drop1 with muc own formula", drop1(f1, test = "Chisq"))
f1z <- frm(bf(y ~ x, muc ~ x), family = categorical(refcat = "b"), data = d3)
cat("refit LRT for z:", 2 * (as.numeric(logLik(f1)) - as.numeric(logLik(f1z))), "\n")
fna <- suppressWarnings(frm(bf(y ~ x + z), family = categorical(refcat = NA), data = d3))
TRY("drop1 refcat NA", suppressWarnings(drop1(fna, test = "Chisq")))
fa <- frm(bf(y ~ x + z), family = categorical(), data = d3)
TRY("drop1 refcat a", drop1(fa, test = "Chisq"))
set.seed(8)
dmx <- data.frame(x = rnorm(300), z = rnorm(300))
dmx$y <- ifelse(runif(300) < 0.4, rnorm(300, -2 + dmx$x), rnorm(300, 2 + 0.5 * dmx$x))
fmx <- frm(bf(y ~ x + z), family = mixture(gaussian(), gaussian()), data = dmx)
TRY("drop1 mixture", drop1(fmx, test = "Chisq"))
fmx2 <- frm(bf(y ~ x), family = mixture(gaussian(), gaussian()), data = dmx)
cat("mixture refit LRT for z:", 2 * (as.numeric(logLik(fmx)) - as.numeric(logLik(fmx2))), "\n")

cat("\n==== confint aliases ====\n")
set.seed(2)
d2 <- data.frame(x = stats::rnorm(200))
e <- matrix(stats::rnorm(400), 200) %*% chol(matrix(c(1, .5, .5, 1), 2))
d2$y1 <- 1 + d2$x + e[, 1]; d2$y2 <- -1 + 0.5 * d2$x + e[, 2]
fr2 <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), data = d2)
ci <- TRY("confint rescor alias", confint(fr2, parm = "rescor__y1__y2"))
cip <- TRY("confint rescor alias profile", confint(fr2, parm = "rescor__y1__y2", method = "profile"))
V <- vcov(fr2, full = TRUE)
TRY("vcov full names", head(rownames(V), 12))
th <- fr2$estimates$thetar
cat("thetar", th, "r", th / sqrt(1 + th^2), "\n")
k <- grep("thetar", rownames(V))
if (length(k)) cat("Wald from vcov:", th + c(-1, 1) * qnorm(0.975) * sqrt(V[k, k]), "\n")
TRY("summary rescor row", summary(fr2)$rescor_pars)
TRY("cor of residuals", cor(d2$y1 - fitted(fr2)[, "Estimate", "y1"], d2$y2 - fitted(fr2)[, "Estimate", "y2"]))
cat("DONE\n")
