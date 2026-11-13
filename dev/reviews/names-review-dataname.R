# Review of lane names, 1550 and 1543: the Data: and Formula: lines.
# Usage: Rscript dev/reviews/names-review-dataname.R <base|lane>
a <- commandArgs(TRUE)
arm <- a[1]
lib <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
         "C:/Users/adf44/source/r/rellib-r20",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
dline <- function(fit) {
  s <- capture.output(summary(fit))
  p <- capture.output(print(fit))
  c(summary = grep("Data:", s, value = TRUE)[1],
    print = grep("Data:", p, value = TRUE)[1])
}
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(1)
d <- data.frame(x = rnorm(300), g = factor(rep(1:10, 30)))
d$y <- 1 + d$x + rnorm(300)
out <- list()
out$name <- dline(q(frm(y ~ x, data = d)))
out$docall <- dline(q(do.call(frm, list(y ~ x, data = d))))
out$docall_quoted <- dline(q(do.call(frm, list(y ~ x, data = quote(d)))))
dd <- d
out$variable <- dline(q(frm(y ~ x, data = dd)))
out$pipe <- dline(q(d |> frm(formula = y ~ x, data = _)))
out$subset_call <- dline(q(frm(y ~ x, data = d[d$x > -10, ])))
out$long_call <- dline(q(frm(y ~ x, data = transform(d, zzzzzzzzzzzzzzzz = x * 2, wwwwwwwwwwwwwwwwwwwwww = x))))
if (requireNamespace("tibble", quietly = TRUE)) {
  tb <- tibble::as_tibble(d)
  out$tibble_docall <- dline(q(do.call(frm, list(y ~ x, data = tb))))
}
f0 <- q(frm(y ~ x, data = d))
out$update_newdata_name <- dline(q(update(f0, newdata = dd)))
out$update_newdata_value <- dline(q(do.call(update, list(f0, newdata = d))))
out$lapply <- dline(q(lapply(list(d), function(z) frm(y ~ x, data = z))[[1]]))
out$list_data <- dline(q(do.call(frm, list(y ~ x, data = as.list(d)))))
for (nm in names(out)) cat(sprintf("%-22s %s | %s\n", nm, out[[nm]][1],
                                   out[[nm]][2]))
# cost: independent of rows?
big <- data.frame(x = rnorm(1e6))
big$y <- big$x + rnorm(1e6)
t_name <- system.time(fb <- q(frm(y ~ x, data = big)))[["elapsed"]]
t_dc <- system.time(fdc <- q(do.call(frm, list(y ~ x, data = big))))[["elapsed"]]
t_sn <- system.time(for (i in 1:5) s1 <- summary(fb))[["elapsed"]] / 5
t_sd <- system.time(for (i in 1:5) s2 <- summary(fdc))[["elapsed"]] / 5
cat(sprintf("1e6 rows: frm by name %.2f s, by do.call %.2f s; summary by name %.3f s, by do.call %.3f s\n",
            t_name, t_dc, t_sn, t_sd))
cat("do.call 1e6 Data line:", dline(fdc), "\n")
if (arm == "lane") {
  dno <- getFromNamespace("data_name_of", "frmtmb")
  tm <- vapply(c(1e3, 1e5, 1e6), function(nr) {
    dfr <- big[seq_len(nr), ]
    system.time(for (i in 1:200) dno(dfr, dfr))[["elapsed"]] / 200
  }, 0)
  cat("data_name_of seconds per call at 1e3, 1e5, 1e6 rows:",
      format(tm, digits = 3), "\n")
}
# summary of a fit whose call holds the frame (saved fit, no attribute)
fdc2 <- fdc
attr(fdc2$data, "data_name") <- NULL
cat("no attribute, call by value:", dline(fdc2)[1], "\n")
# Formula lines
tq <- function(e) tryCatch(q(e), error = function(err) err)
fl <- list(
  mix = tq(frm(bf(y ~ x, sigma1 = "sigma2"),
               family = mixture(gaussian(), gaussian()), data = d)),
  mixth = tq(frm(bf(y ~ x, theta1 = 0.3),
                 family = mixture(gaussian(), gaussian()), data = d)),
  mv_pf = tq(frm(bf(y ~ x, sigma ~ x) + bf(x ~ 1) + set_rescor(FALSE),
                 data = d)),
  nlf = tq(frm(bf(y ~ a + b * x, nlf(a ~ 1 + 0 * x), b ~ 1, nl = TRUE),
               data = d)),
  lf = tq(frm(bf(y ~ x) + lf(sigma ~ x), data = d)),
  ar = tq(frm(bf(y ~ x + ar(p = 1)), data = transform(d, t = seq_len(300)))))
for (nm in names(fl)) {
  if (inherits(fl[[nm]], "error")) {
    cat(nm, "ERROR", conditionMessage(fl[[nm]]), "\n")
    next
  }
  s <- capture.output(summary(fl[[nm]]))
  i <- grep("^Formula:", s)
  j <- grep("Data:", s)
  cat("--", nm, "\n")
  writeLines(s[i:(j - 1)])
}
cat("\nDONE\n")
