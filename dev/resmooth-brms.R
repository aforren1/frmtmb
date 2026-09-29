# Lane wt-resmooth. What brms does with a smooth indexed by a grouping
# factor under re_formula = NA, and with a NEW level of that factor in
# newdata. Stan compiles fresh on this machine, so the four fits are
# cached in dev/resmooth-brms-fits.rds and the analysis below reruns
# without them.
#   Rscript dev/resmooth-brms.R > dev/resmooth-brms.txt
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")),
    "rstan", as.character(packageVersion("rstan")),
    "mgcv", as.character(packageVersion("mgcv")), "\n")

set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)

forms <- list(
  re = y ~ s(x) + s(g, bs = "re"),
  fs = y ~ s(x, g, bs = "fs", k = 5),
  t2 = y ~ t2(x, g, bs = c("cr", "re")),
  sg = y ~ s(x) + (1 | g)
)
cache <- "dev/resmooth-brms-fits.rds"
fits <- if (file.exists(cache)) readRDS(cache) else list()
for (nm in names(forms)) {
  if (!is.null(fits[[nm]])) next
  cat("fitting", nm, format(Sys.time()), "\n"); flush.console()
  fits[[nm]] <- brm(forms[[nm]], data = d, chains = 2, iter = 1000,
                    seed = 7, refresh = 0, silent = 2)
  saveRDS(fits, cache)
  cat("done", nm, format(Sys.time()), "\n"); flush.console()
}

msg <- function(e) substr(conditionMessage(e), 1, 200)
tryv <- function(expr) tryCatch(expr, error = function(e) e)

cat("\n== does re_formula = NA change brms's expectation? ==\n")
for (nm in names(fits)) {
  a <- posterior_epred(fits[[nm]], re_formula = NA)
  b <- posterior_epred(fits[[nm]], re_formula = NULL)
  cat(sprintf("%-4s identical(NA, NULL) %-5s  max|NA-NULL| %.4g  ",
              nm, identical(a, b), max(abs(a - b))))
  # which parameters the two draws read: brms's own bookkeeping
  pp <- posterior_summary(fits[[nm]])
  cat("pars:", paste(grep("^(bs_|sds_|sd_|s_|r_)", rownames(pp),
                          value = TRUE)[1:4], collapse = " "), "\n")
}

cat("\n== brms with a NEW level of g in newdata ==\n")
nd2 <- d[1:5, ]
nd2$g <- factor("99", levels = c(levels(d$g), "99"))
nd3 <- nd2[, setdiff(names(nd2), "g"), drop = FALSE]
for (nm in names(fits)) {
  for (anl in c(FALSE, TRUE)) {
    for (rf in list(NA, NULL)) {
      p <- tryv(posterior_epred(fits[[nm]], newdata = nd2, re_formula = rf,
                                allow_new_levels = anl))
      cat(sprintf("%-4s re_formula=%-4s allow_new_levels=%-5s %s\n", nm,
                  if (is.null(rf)) "NULL" else "NA", anl,
                  if (inherits(p, "error")) paste("ERROR:", msg(p)) else
                    paste("OK mean", sprintf("%.4f", mean(colMeans(p))))))
    }
  }
  p <- tryv(posterior_epred(fits[[nm]], newdata = nd3, re_formula = NA,
                            allow_new_levels = TRUE))
  cat(sprintf("%-4s no g column, re_formula=NA   %s\n", nm,
              if (inherits(p, "error")) paste("ERROR:", msg(p)) else "OK"))
}

cat("\n== mgcv PredictMat at an unseen factor level ==\n")
for (spec in c("s(g, bs = \"re\")", "s(x, g, bs = \"fs\", k = 5)",
               "t2(x, g, bs = c(\"cr\", \"re\"))")) {
  sm <- mgcv::smoothCon(eval(str2lang(spec)), data = d,
                        absorb.cons = TRUE)[[1L]]
  M <- tryv(mgcv::PredictMat(sm, nd2))
  cat(sprintf("%-32s %s\n", spec,
              if (inherits(M, "error")) paste("ERROR:", msg(M)) else
                paste("dim", nrow(M), "x", ncol(M), "rowsum",
                      paste(sprintf("%.3g", rowSums(abs(M))),
                            collapse = " "))))
}

cat("\n== design pieces for the parity check ==\n")
out <- list(data = d)
for (nm in c("re", "fs", "t2", "sg")) {
  sd_ <- brms::standata(fits[[nm]])
  ps <- posterior_summary(fits[[nm]])
  pieces <- sd_[grep("^(Xs|Zs_)", names(sd_))]
  cat(nm, ":", paste(names(pieces),
                     vapply(pieces, function(m) paste(dim(as.matrix(m)),
                                                      collapse = "x"), ""),
                     collapse = " "), "\n")
  out[[nm]] <- list(pieces = pieces, post = ps[, "Estimate"],
                    epred_na = colMeans(posterior_epred(fits[[nm]],
                                                        re_formula = NA)),
                    epred_null = colMeans(posterior_epred(fits[[nm]],
                                                          re_formula = NULL)))
}
saveRDS(out, "dev/resmooth-brms-design.rds")
cat("saved dev/resmooth-brms-design.rds\n")
