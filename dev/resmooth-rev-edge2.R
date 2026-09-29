# Reviewer: edge cases the change could have broken without a test
# noticing.
#  (a) an ordinal fit with s(x) + (1 | g): newdata with NO g column at
#      re_formula = NA. The smooth needs no level, so this must still
#      answer, and its Est.Error must be finite.
#  (b) the same with a group-indexed smooth beside (1 | g): predicting
#      at an unseen level with allow_new_levels = TRUE.
#  (c) a gp() fit on newdata OFF the fitted positions.
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", ARM, "\n")
tryv <- function(x) tryCatch(x, error = function(e) e)
msg <- function(e) substr(conditionMessage(e), 1, 130)
say <- function(lab, r) cat(sprintf("  %-44s %s\n", lab,
  if (inherits(r, "error")) paste("ERROR:", msg(r)) else
    paste(sprintf("%.5f", utils::head(as.vector(r), 4)), collapse = " ")))

set.seed(31)
ng <- 8; per <- 30
d <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                x = stats::runif(ng * per, -2, 2))
lat <- sin(1.5 * d$x) + stats::rnorm(ng, 0, 0.8)[d$g] +
  stats::rlogis(nrow(d))
d$y <- factor(cut(lat, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
              ordered = TRUE)
f1 <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                          family = cumulative(), data = d))
cat("(a) ordinal s(x) + (1 | g), newdata with NO g column\n")
nd <- data.frame(x = c(-1, 0, 1))
say("fitted(nd, NA)[, Estimate, 1]",
    tryv(fitted(f1, newdata = nd, re_formula = NA)[, "Estimate", 1]))
say("fitted(nd, NA)[, Est.Error, 1]",
    tryv(fitted(f1, newdata = nd, re_formula = NA)[, "Est.Error", 1]))
say("fitted(nd, NULL)", tryv(fitted(f1, newdata = nd,
                                    re_formula = NULL)[, "Estimate", 1]))
say("frm_linpred(nd, NA, se.fit)$se.fit",
    tryv(frm_linpred(f1, newdata = nd, re_formula = NA,
                     se.fit = TRUE)$se.fit))

cat("(b) gaussian fs smooth beside (1 | h), unseen fs level\n")
set.seed(41)
d2 <- data.frame(g = factor(rep(1:8, each = 30)),
                 h = factor(rep(1:6, length.out = 240)),
                 x = stats::runif(240))
d2$y <- sin(2 * pi * d2$x) * stats::rnorm(8, 0, 0.6)[d2$g] +
  stats::rnorm(6, 0, 0.5)[d2$h] + stats::rnorm(240, 0, 0.3)
f2 <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | h)),
                          family = gaussian(), data = d2))
nd2 <- data.frame(x = c(0.2, 0.6),
                  g = factor(rep("new", 2), levels = c(levels(d2$g), "new")),
                  h = d2$h[1:2])
say("frm_linpred(nd, allow_new_levels)",
    tryv(frm_linpred(f2, newdata = nd2, allow_new_levels = TRUE)))
say("frm_linpred(nd, NA, allow_new_levels)",
    tryv(frm_linpred(f2, newdata = nd2, re_formula = NA,
                     allow_new_levels = TRUE)))
say("frm_linpred(nd, NA) no opt-in",
    tryv(frm_linpred(f2, newdata = nd2, re_formula = NA)))

cat("(c) exact gp() off the fitted positions\n")
set.seed(47)
d3 <- data.frame(x = stats::runif(120, -2, 2))
d3$y <- sin(2 * d3$x) + stats::rnorm(120, 0, 0.3)
f3 <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d3))
say("frm_linpred(off-grid, NA, se.fit)$se.fit",
    tryv(frm_linpred(f3, newdata = data.frame(x = c(-1.234, 0.777)),
                     re_formula = NA, se.fit = TRUE)$se.fit))
cat("DONE\n")
