# Reviewer, claim 2, part 2: the public newdata route on an fs smooth
# with a by = factor (the base arm answered NA without an error and the
# lane arm refuses), and a frmtmb.sample draws object at re_formula = NA.
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "| frmtmb:", find.package("frmtmb"), "\n")
cat("frmtmb.sample:", find.package("frmtmb.sample"), "\n")
msg <- function(e) substr(conditionMessage(e), 1, 160)
tryv <- function(x) tryCatch(x, error = function(e) e)
say <- function(lab, r) {
  cat(sprintf("  %-46s %s\n", lab,
              if (inherits(r, "error")) paste("ERROR:", msg(r)) else
                paste("OK", paste(sprintf("%.4f", utils::head(as.vector(r), 3)),
                                  collapse = " "))))
}

set.seed(5)
n <- 300
d <- data.frame(x = runif(n), z = runif(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(1:10, each = 30)))
d$y <- sin(2 * pi * d$x) + d$z^2 + c(0, 1, -1)[d$f] * d$x +
  rnorm(10, 0, 0.5)[d$g] + rnorm(n, 0, 0.3)

cat("== s(x, g, bs = 'fs', by = f), the PUBLIC newdata route ==\n")
fby <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5, by = f)),
                           family = gaussian(), data = d))
nd <- d[1:5, ]
say("frm_linpred(newdata, NA)", tryv(frm_linpred(fby, newdata = nd,
                                                 re_formula = NA)))
say("frm_linpred(newdata, NULL)", tryv(frm_linpred(fby, newdata = nd,
                                                   re_formula = NULL)))
say("fitted(newdata, NA)", tryv(fitted(fby, newdata = nd,
                                       re_formula = NA)[, "Estimate"]))
say("frm_linpred(in sample, NA)", tryv(frm_linpred(fby, re_formula = NA)))
say("conditional_effects(effects = x)",
    tryv(range(conditional_effects(fby, effects = "x",
                                   resolution = 8)[[1]]$estimate__)))
say("simulate(nsim = 2, newdata, NA)",
    tryv(as.matrix(simulate(fby, nsim = 2, seed = 1, newdata = nd,
                            re_formula = NA))))

cat("\n== frmtmb.sample draws at re_formula = NA ==\n")
suppressMessages(library(frmtmb.sample))
cat("frmtmb.sample version:",
    as.character(packageVersion("frmtmb.sample")), "\n")
for (lab in c("fs", "re", "t2")) {
  form <- switch(lab,
    fs = bf(y ~ s(x, g, bs = "fs", k = 5)),
    re = bf(y ~ s(x) + s(g, bs = "re")),
    t2 = bf(y ~ t2(x, g, bs = c("cr", "re"))))
  fit <- suppressWarnings(frm(form, family = gaussian(), data = d))
  dr <- tryv(suppressMessages(frm_sample(fit, chains = 1, iter = 400,
                                         refresh = 0, seed = 3)))
  if (inherits(dr, "error")) {
    cat(sprintf("%-4s frm_sample ERROR: %s\n", lab, msg(dr))); next
  }
  ea <- tryv(posterior_epred(dr, re_formula = NA))
  en <- tryv(posterior_epred(dr, re_formula = NULL))
  pa <- tryv(posterior_predict(dr, re_formula = NA, seed = 4))
  pn <- tryv(posterior_predict(dr, re_formula = NULL, seed = 4))
  cat(sprintf("%-4s epred: identical(NA,NULL) %-5s max|NA-NULL| %9.3g\n",
              lab, identical(ea, en),
              if (inherits(ea, "error")) NA_real_ else max(abs(ea - en))))
  # the population epred the draws object should give: the smooth kept at
  # each draw. Reference: the per-row sd of predict draws over sigma.
  sg <- as.vector(frm_linpred(fit, dpar = "sigma", type = "response",
                              re_formula = NA))
  if (length(sg) == 1L) sg <- rep(sg, n)
  cat(sprintf("     ppred rowsd/sigma at NA: %.3f | at NULL: %.3f\n",
              if (inherits(pa, "error")) NA_real_ else
                stats::median(apply(pa, 2, stats::sd) / sg),
              if (inherits(pn, "error")) NA_real_ else
                stats::median(apply(pn, 2, stats::sd) / sg)))
  # an unseen level through the draws route
  nd2 <- d[1:5, ]
  nd2$g <- factor("99", levels = c(levels(d$g), "99"))
  say(paste(lab, "epred newdata new level, NA"),
      tryv(colMeans(posterior_epred(dr, newdata = nd2, re_formula = NA))))
  nd3 <- nd2[, setdiff(names(nd2), "g"), drop = FALSE]
  say(paste(lab, "epred newdata no g column, NA"),
      tryv(colMeans(posterior_epred(dr, newdata = nd3, re_formula = NA))))
}
cat("DONE\n")
