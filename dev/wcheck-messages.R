# Lane wcheck: the full refusal text on representative cases, and the
# block labels the message reads in a multivariate model.
#   Rscript dev/wcheck-messages.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
msg <- function(label, expr) {
  m <- tryCatch(suppressWarnings(suppressMessages({
    expr
    "no refusal"
  })), error = function(e) conditionMessage(e))
  cat("== ", label, "\n", paste(strwrap(m, 78), collapse = "\n"), "\n\n",
      sep = "")
}
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, 0.8)
d$y2 <- d$y + rnorm(n)
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
G <- gaussian()
fr <- suppressWarnings(frm(mvbf(bf(y | weights(w) ~ x + (1 | p | obs)),
                                bf(y2 | weights(w) ~ x + (1 | p | obs))),
                           family = G, data = d, dry_run = "frame"))
for (bk in fr$re_blocks) {
  cat("block written label:", frmtmb:::block_written_label(bk), "| term:",
      frmtmb:::block_term_as_written(bk), "| components:",
      paste(vapply(bk$components, function(cp) {
        paste0(cp$lp_key, " = ", cp$label)
      }, ""), collapse = "; "), "\n")
}
cat("\n")
msg("mv, (1 | obs) in both, weights on y2 only",
    frm(mvbf(bf(y ~ x + (1 | obs)), bf(y2 | weights(w) ~ x + (1 | obs))),
        family = G, data = d))
msg("mv, (1 | p | obs) in both",
    frm(mvbf(bf(y | weights(w) ~ x + (1 | p | obs)),
             bf(y2 | weights(w) ~ x + (1 | p | obs))), family = G, data = d))
msg("univariate, sigma ~ 0 + fac, normal on faca only",
    frm(bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac), family = G,
        data = d, prior = set_prior("normal(0, 1)", class = "b",
                                    dpar = "sigma", coef = "faca")))
msg("Gamma, student_t(3, 0, 2.5) on shape",
    frm(bf(yg | weights(w) ~ x + (1 | obs)), family = Gamma("log"), data = d,
        prior = set_prior("student_t(3, 0, 2.5)", class = "shape")))
q <- quantile(d$y, 0.8)
d$cc <- ifelse(d$y > q, "right", "none")
d$yc <- pmin(d$y, q)
d$wc <- ifelse(d$cc == "none", 1 + 2 / sum(d$cc == "none"), 1)
msg("gaussian, 16 right censored, uncensored excess 2",
    frm(bf(yc | weights(wc) + cens(cc) ~ x + (1 | obs)), family = G,
        data = d))
ql <- quantile(d$y, 0.1)
d$c2 <- ifelse(d$y > q, "right", ifelse(d$y < ql, "left", "none"))
d$y2c <- pmin(pmax(d$y, ql), q)
d$w2 <- ifelse(d$c2 == "none", 1 + 8 / sum(d$c2 == "none"), 1)
msg("gaussian, 16 right and 8 left censored, uncensored excess 8",
    frm(bf(y2c | weights(w2) + cens(c2) ~ x + (1 | obs)), family = G,
        data = d))
cat("DONE\n")
