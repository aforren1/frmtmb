# Reviewer, lane ordinal: every new refusal, with the case where its
# condition is absent, plus nearby shapes. Data seed 20261005.
# Output: dev/ordinal-rev-log-refusals.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261005)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                h = factor(sample(c("p", "q"), n, TRUE)))
u <- stats::rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$y2 <- 1L + (u > 0) + (stats::rlogis(n) > 0.5)
d$yb <- 1L + (u > 0)                      # two categories
d$y3 <- pmin(d$y, 3L)                     # three categories
d$ym <- d$y; d$ym[d$ym == 3L] <- 2L       # unused middle category
d$yf <- factor(d$ym, levels = 1:5, ordered = TRUE)
d$ygr <- ifelse(d$h == "p", pmin(d$y, 3L), d$y)   # 2 vs 4 thresholds
d$y12h <- d$yb                             # every level one threshold

probe <- function(lab, expr, expect) {
  r <- tryCatch({
    w <- character()
    v <- withCallingHandlers(expr, warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
    paste0("FITS", if (length(w)) paste0(" [warn: ", substr(w[1], 1, 90), "]"))
  }, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 200)))
  ok <- (expect == "fit") == startsWith(r, "FITS")
  cat(sprintf("%-4s expect %-6s | %s\n   -> %s\n", if (ok) "ok" else "BAD",
              expect, lab, r))
}

cat("## sum_to_zero and class Intercept\n")
probe("stz + Intercept prior on dpar disc (disc ~ z)",
      frm(bf(y ~ x, disc ~ z), family = cumulative(threshold = "sum_to_zero"),
          data = d, prior = set_prior("normal(0, 1)", class = "Intercept",
                                      dpar = "disc")), "fit")
probe("stz + class b prior", frm(y ~ x, family = acat(threshold = "sum_to_zero"),
      data = d, prior = set_prior("normal(0, 1)", class = "b")), "fit")
probe("stz + class Intercept prior (brms fits this)",
      frm(y ~ x, family = acat(threshold = "sum_to_zero"), data = d,
          prior = set_prior("normal(0, 5)", class = "Intercept")), "refuse")
probe("mv: flexible y + stz y2, Intercept prior resp = y",
      frm(bf(y ~ x) + bf(y2 ~ x), data = d,
          family = list(cumulative(), sratio(threshold = "sum_to_zero")),
          prior = set_prior("normal(0, 5)", class = "Intercept", resp = "y")),
      "fit")
probe("mv: flexible y + stz y2, Intercept prior resp = y2",
      frm(bf(y ~ x) + bf(y2 ~ x), data = d,
          family = list(cumulative(), sratio(threshold = "sum_to_zero")),
          prior = set_prior("normal(0, 5)", class = "Intercept", resp = "y2")),
      "refuse")
probe("stz + prior = list(tau_raw = ) escape hatch",
      frm(y ~ x, family = sratio(threshold = "sum_to_zero"), data = d,
          prior = list(tau_raw = "normal(0, 5)")), "fit")

cat("\n## too few thresholds\n")
probe("2 categories, flexible", frm(yb ~ x, family = cumulative(), data = d), "fit")
probe("2 categories, equidistant (brms: int<lower=2> nthres)",
      frm(yb ~ x, family = sratio(threshold = "equidistant"), data = d), "refuse")
probe("2 categories, sum_to_zero",
      frm(yb ~ x, family = cratio(threshold = "sum_to_zero"), data = d), "refuse")
probe("3 categories, equidistant",
      frm(y3 ~ x, family = acat(threshold = "equidistant"), data = d), "fit")
probe("3 categories, sum_to_zero, hurdle",
      frm(y3 ~ x, family = hurdle_cumulative(threshold = "sum_to_zero"),
          data = d), "fit")
probe("thres(gr): levels of 2 and 4 thresholds, equidistant",
      frm(ygr | thres(gr = h) ~ x, family = cumulative(threshold = "equidistant"),
          data = d), "fit")
probe("thres(gr): every level ONE threshold, sum_to_zero (brms: int<lower=1>)",
      frm(y12h | thres(gr = h) ~ x, family = cumulative(threshold = "sum_to_zero"),
          data = d), "fit")
probe("thres(x = 5) on a 3-category response, equidistant",
      frm(y3 | thres(4) ~ x, family = cumulative(threshold = "equidistant"),
          data = d), "fit")
probe("thres(x = 5) on a 3-category response, sum_to_zero",
      frm(y3 | thres(4) ~ x, family = sratio(threshold = "sum_to_zero"),
          data = d), "fit")

cat("\n## class delta\n")
probe("delta on flexible", frm(y ~ x, family = cumulative(), data = d,
      prior = set_prior("normal(1, 1)", class = "delta")), "refuse")
probe("delta on sum_to_zero", frm(y ~ x, family = cumulative(threshold = "sum_to_zero"),
      data = d, prior = set_prior("normal(1, 1)", class = "delta")), "refuse")
probe("delta group = p on grouped equidistant",
      frm(y | thres(gr = h) ~ x, family = sratio(threshold = "equidistant"),
          data = d, prior = set_prior("normal(1, 1)", class = "delta",
                                      group = "p")), "fit")
probe("delta no group on grouped equidistant",
      frm(y | thres(gr = h) ~ x, family = sratio(threshold = "equidistant"),
          data = d, prior = set_prior("normal(1, 1)", class = "delta")), "fit")
probe("delta group = zz (no such level)",
      frm(y | thres(gr = h) ~ x, family = sratio(threshold = "equidistant"),
          data = d, prior = set_prior("normal(1, 1)", class = "delta",
                                      group = "zz")), "refuse")
probe("mv: delta resp = y2 (equidistant) beside flexible y",
      frm(bf(y ~ x) + bf(y2 ~ x), data = d,
          family = list(cumulative(), cumulative(threshold = "equidistant")),
          prior = set_prior("normal(1, 1)", class = "delta", resp = "y2")),
      "fit")
probe("delta with a bound only, lb = 0.9 (ML delta near 1)",
      { f <- frm(y ~ x, family = cumulative(threshold = "equidistant"), data = d,
                 prior = set_prior("", class = "delta", lb = 1.2))
        tau <- frmtmb:::ord_threshold_values(family(f), f$estimates$tau_raw)
        cat("      delta with lb 1.2:", diff(tau)[1], "\n"); f }, "fit")

cat("\n## ordered factor with an unused level, and brmsfamily()\n")
probe("ordered factor, empty level 3, equidistant",
      { f <- frm(yf ~ x, family = cumulative(threshold = "equidistant"), data = d)
        cat("      nthres:", family(f)$thres$nthres, " brms nthres:",
            brms::standata(brms::bf(yf ~ x), data = d,
                           family = brms::cumulative(threshold = "equidistant"))$nthres,
            "\n"); f }, "fit")
probe("brmsfamily('cumulative', threshold = 'sum_to_zero')",
      frm(y ~ x, family = brmsfamily("cumulative", threshold = "sum_to_zero"),
          data = d), "fit")
probe("brmsfamily('acat', 'probit', threshold = 'equidistant')",
      frm(y ~ x, family = brmsfamily("acat", "probit", threshold = "equidistant"),
          data = d), "fit")

cat("\n## the disc-intercept warning\n")
probe("disc ~ 0 + h (cell means span the intercept; no warning expected by the guard)",
      frm(bf(y ~ x, disc ~ 0 + h), family = sratio(), data = d), "fit")
probe("mv: disc ~ z on y, Intercept prior on its thresholds (pinned wording)",
      frm(bf(y ~ x, disc ~ z) + bf(y2 ~ x), data = d,
          family = list(cumulative(), sratio()),
          prior = set_prior("normal(0, 3)", class = "Intercept", resp = "y")),
      "fit")

cat("\n## brms's Stan data declarations for thres(gr = ) counts\n")
sc <- brms::stancode(brms::bf(y12h | thres(gr = h) ~ x), data = d,
                     family = brms::cumulative(threshold = "sum_to_zero"))
cat(grep("nthres", strsplit(sc, "\n")[[1]], value = TRUE)[1:3], sep = "\n")
sc <- brms::stancode(brms::bf(y ~ x), data = d,
                     family = brms::cumulative(threshold = "sum_to_zero"))
cat(grep("int.*nthres", strsplit(sc, "\n")[[1]], value = TRUE)[1], sep = "\n")
