# Reviewer probe 2: targeted false-accept and false-refuse cases for
# draws_newdata_is_fitted()
.libPaths(c("C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb.sample"), "\n")
isf <- frmtmb.sample:::draws_newdata_is_fitted
mk <- function(fit) structure(list(stanfit = NULL, draws = matrix(0, 1, 1),
                                   fit = fit), class = "frmtmb_draws")
q <- function(expr) suppressWarnings(suppressMessages(expr))
show <- function(lbl, val) cat(sprintf("%-55s %s\n", lbl, format(val)))

set.seed(2)
n <- 40
dd <- data.frame(x = rnorm(n), z = rnorm(n),
                 g = factor(rep(letters[1:5], each = 8)),
                 cc = rep(c(0, 0, 0, 2), 10), u = runif(n))
dd$y <- 1 + 0.5 * dd$x + rnorm(n)
dd$y2 <- dd$y + 1
dd$lo <- dd$y - 2
dd$hi <- dd$y + 2
dd$xo <- factor(rep(1:4, 10), ordered = TRUE)

# 1. interval-censored upper bound, row 4 is interval censored (cc = 2)
fc <- q(frm(y | cens(cc, y2) ~ x, data = dd))
d1 <- dd; d1$y2[4] <- d1$y2[4] + 3
show("cens: upper bound of an interval row changed, accepted", isf(mk(fc), d1))
show("cens: names(fit$data)", paste(names(fc$data), collapse = ","))
# the fitted log-likelihood does read y2: a refit at the same data but
# the changed bound gives another logLik at the ML optimum
fc1 <- q(frm(y | cens(cc, y2) ~ x, data = d1))
show("cens: logLik fitted vs with changed bound", sprintf("%.6f vs %.6f",
     as.numeric(logLik(fc)), as.numeric(logLik(fc1))))

# 2. truncation bounds from variables
ft <- tryCatch(q(frm(y | trunc(lb = lo, ub = hi) ~ x, data = dd)),
               error = function(e) e)
if (!inherits(ft, "error")) {
  d2 <- dd; d2$hi[3] <- d2$hi[3] + 1
  show("trunc: ub changed, accepted", isf(mk(ft), d2))
  show("trunc: names(fit$data)", paste(names(ft$data), collapse = ","))
} else show("trunc fit error", conditionMessage(ft))

# 3. monotonic: change one ordered level
fm <- q(frm(y ~ mo(xo), data = dd))
d3 <- dd; d3$xo[1] <- factor(3, levels = 1:4, ordered = TRUE)
show("mo: xo[1] 1 -> 3, accepted", isf(mk(fm), d3))
# same labels, reversed level order (brms re-levels newdata by label)
d3b <- dd; d3b$xo <- factor(as.character(dd$xo), levels = 4:1, ordered = TRUE)
show("mo: same labels, reversed level order, accepted", isf(mk(fm), d3b))

# 4. benign variants of the true fitted data
fr <- q(frm(y ~ x + (1 | g), data = dd))
ds <- mk(fr)
show("re: data itself", isf(ds, dd))
show("re: fit$data", isf(ds, fr$data))
show("re: tibble", if (requireNamespace("tibble", quietly = TRUE))
  isf(ds, tibble::as_tibble(dd)) else NA)
dch <- dd; dch$g <- as.character(dch$g)
show("re: g as character", isf(ds, dch))
show("re: extra column", isf(ds, cbind(dd, extra = 1)))
dlev <- dd; dlev$g <- factor(dlev$g, levels = rev(levels(dd$g)))
show("re: g levels reordered, same labels", isf(ds, dlev))
show("re: rows reordered", isf(ds, dd[c(2, 1, 3:n), ]))
show("re: data.table", if (requireNamespace("data.table", quietly = TRUE))
  isf(ds, data.table::as.data.table(dd)) else NA)
# integer vs double
di <- data.frame(y = dd$y, k = rep(1:4, 10), g = dd$g)
fi <- q(frm(y ~ k, data = di))
dk <- di; dk$k <- as.double(dk$k)
show("int fitted, double newdata", isf(mk(fi), dk))
dkf <- di; dkf$k <- factor(dkf$k)
show("numeric fitted, factor newdata with same labels", isf(mk(fi), dkf))
dkc <- di; dkc$k <- as.character(dkc$k)
show("numeric fitted, character newdata", isf(mk(fi), dkc))

# 5. NA rows the fit dropped
dna <- dd; dna$x[c(5, 9)] <- NA
fna <- q(frm(y ~ x, data = dna))
show("NA: nobs", nrow(fna$data))
show("NA: raw data with NA rows", isf(mk(fna), dna))
show("NA: na.omit(data)", isf(mk(fna), na.omit(dna)))
dna2 <- dna[-c(5, 9), ]; dna2$x[1] <- dna2$x[1] + 1
show("NA: dropped rows, x changed", isf(mk(fna), dna2))
# a variable with an NA the fit kept? (in an unused column) and NA in y
dny <- dd; dny$u[3] <- NA
show("NA in an unused column", isf(ds, dny))

# 6. formula-environment variable
zz <- rnorm(n)
fe <- q(frm(y ~ x + zz, data = dd))
show("env var: data without zz", isf(mk(fe), dd))
show("env var: data with zz", isf(mk(fe), cbind(dd, zz = zz)))
show("env var: fit$data", isf(mk(fe), fe$data))

# 7. Date and logical predictors
dd7 <- dd; dd7$dt <- as.Date("2020-01-01") + 0:(n - 1); dd7$lg <- dd$x > 0
f7 <- q(frm(y ~ lg, data = dd7))
d7 <- dd7; d7$lg <- as.numeric(d7$lg)
show("logical fitted, 0/1 newdata", isf(mk(f7), d7))

# 8. fast path with a raw variable outside the frame: fit$data plus a
# changed x for a poly() model
fp <- q(frm(y ~ poly(x, 2), data = dd))
np <- fp$data; np$x <- dd$x + 1
show("poly: fit$data plus a different x column", isf(mk(fp), np))
