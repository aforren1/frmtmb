# Reviewer re-check (punch round 1), lane ordinal: m6, the disc-intercept
# warning for a design that spans the intercept without an intercept
# column. Where it must fire, where it must not, and what silences it.
# Data seed 20261016. Output: dev/ordinal-rev2-log-m6.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordinal-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261016)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                h = factor(sample(c("p", "q", "r"), n, TRUE)),
                g = factor(sample(letters[1:6], n, TRUE)))
d$zc <- 1 - d$w          # w + zc spans the constant exactly
d$zs <- 100 * d$z        # a large-scale column, no constant in its span
u <- stats::rlogis(n) / exp(0.2 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
probe <- function(lab, expr, expect) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) paste("ERROR:", conditionMessage(e)))
  hit <- any(grepl("disc has an intercept", w))
  ok <- (expect == "warn") == hit
  cat(sprintf("%-4s expect %-6s | %s\n", if (ok) "ok" else "BAD", expect, lab))
  if (length(w)) cat("     ", substr(w[grepl("disc", w)][1], 1, 260), "\n")
  if (is.character(r)) cat("     ", r, "\n")
}
probe("disc ~ 0 + h (cell means)",
      frm(bf(y ~ x, disc ~ 0 + h), family = sratio(), data = d), "warn")
probe("disc ~ 0 + h, prior class b dpar disc (all columns)",
      frm(bf(y ~ x, disc ~ 0 + h), family = sratio(), data = d,
          prior = set_prior("normal(0, 1)", class = "b", dpar = "disc")), "quiet")
probe("disc ~ 0 + h, prior on one column (coef hp)",
      frm(bf(y ~ x, disc ~ 0 + h), family = sratio(), data = d,
          prior = set_prior("normal(0, 1)", class = "b", coef = "hp",
                            dpar = "disc")), "quiet")
probe("disc ~ 0 + h, prior class b on mu only (pinned wording)",
      frm(bf(y ~ x, disc ~ 0 + h), family = sratio(), data = d,
          prior = set_prior("normal(0, 1)", class = "b")), "warn")
probe("disc ~ 0 + w + zc (continuous columns that add to 1)",
      frm(bf(y ~ x, disc ~ 0 + w + zc), family = cumulative(), data = d), "warn")
probe("disc ~ 0 + z (continuous)",
      frm(bf(y ~ x, disc ~ 0 + z), family = cumulative(), data = d), "quiet")
probe("disc ~ 0 + zs (scale 100)",
      frm(bf(y ~ x, disc ~ 0 + zs), family = cumulative(), data = d), "quiet")
probe("disc ~ 0 + h:z (interaction, no constant)",
      frm(bf(y ~ x, disc ~ 0 + h:z), family = acat(), data = d), "quiet")
probe("disc ~ 0 + z + (1 | g)",
      frm(bf(y ~ x, disc ~ 0 + z + (1 | g)), family = cumulative(), data = d),
      "quiet")
probe("disc ~ z (intercept column)",
      frm(bf(y ~ x, disc ~ z), family = cratio(), data = d), "warn")
probe("disc ~ z, Intercept prior dpar disc",
      frm(bf(y ~ x, disc ~ z), family = cratio(), data = d,
          prior = set_prior("normal(0, 1)", class = "Intercept", dpar = "disc")),
      "quiet")
probe("default disc (held at 1)", frm(y ~ x, family = cratio(), data = d), "quiet")
probe("user-fixed disc = 2", frm(bf(y ~ x, disc = 2), family = cratio(), data = d),
      "quiet")
probe("hurdle, disc ~ 0 + h",
      frm(bf(y ~ x, disc ~ 0 + h), family = hurdle_cumulative(),
          data = transform(d, y = ifelse(runif(n) < 0.2, 0L, y))), "warn")
