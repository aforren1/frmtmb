# Reviewer, claim 4: emmeans on a group-smooth fit at NA and NULL, and
# identical() between the arms on fits with NO group-indexed smooth.
ARM <- if (identical(Sys.getenv("REVLIB"), "base")) "base" else "lane"
LIB <- if (ARM == "base") "C:/Users/adf44/source/r/rellib-r3" else
  "C:/Users/adf44/source/r/wt-resmooth-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
suppressMessages(library(emmeans))
cat("ARM:", ARM, "| frmtmb:", find.package("frmtmb"),
    "| emmeans:", as.character(packageVersion("emmeans")), "\n")
tryv <- function(x) tryCatch(x, error = function(e) e)
msg <- function(e) substr(conditionMessage(e), 1, 140)

set.seed(5)
n <- 300L
d <- data.frame(x = stats::runif(n),
                f = factor(rep(c("a", "b"), length.out = n)),
                g = factor(rep(1:10, each = 30L)),
                z = stats::runif(n, 0.1, 3))
d$y <- sin(2 * pi * d$x) + 0.7 * (d$f == "b") +
  stats::rnorm(10, 0, 0.5)[d$g] * d$x + stats::rnorm(n, 0, 0.3)
d$yo <- ordered(cut(d$y + stats::rlogis(n), 3, labels = FALSE))
d$y2 <- 0.6 * d$x + stats::rnorm(n, 0, 0.4)
d$ynl <- 2 * exp(-0.8 * d$z) + 0.4 * (d$f == "b") + stats::rnorm(n, 0, 0.15)

cat("\n== emmeans(~ x) on an fs fit, at NA and at NULL ==\n")
for (lab in c("fs", "re", "t2")) {
  fm <- switch(lab,
    fs = bf(y ~ f + s(x, g, bs = "fs", k = 5)),
    re = bf(y ~ f + s(x) + s(g, bs = "re")),
    t2 = bf(y ~ f + t2(x, g, bs = c("cr", "re"))))
  fit <- suppressWarnings(frm(fm, family = gaussian(), data = d))
  rg <- tryv(emmeans::ref_grid(fit))
  cat(sprintf("%-3s ref_grid vars: %s\n", lab,
              if (inherits(rg, "error")) paste("ERROR:", msg(rg)) else
                paste(names(rg@levels), collapse = ",")))
  for (rf in list(NULL, NA)) {
    lb <- if (is.null(rf)) "NULL" else "NA"
    e <- tryv(if (is.null(rf)) emmeans::emmeans(fit, ~ x) else
      emmeans::emmeans(fit, ~ x, re_formula = rf))
    cat(sprintf("  ~x re_formula=%-4s %s\n", lb,
                if (inherits(e, "error")) paste("ERROR:", msg(e)) else
                  paste(sprintf("%.5f",
                                as.data.frame(e)$emmean), collapse = " ")))
    e2 <- tryv(if (is.null(rf)) emmeans::emmeans(fit, ~ f) else
      emmeans::emmeans(fit, ~ f, re_formula = rf))
    cat(sprintf("  ~f re_formula=%-4s %s\n", lb,
                if (inherits(e2, "error")) paste("ERROR:", msg(e2)) else
                  paste(sprintf("%.5f",
                                as.data.frame(e2)$emmean), collapse = " ")))
  }
}

cat("\n== fits with NO group-indexed smooth: saved for identical() ==\n")
out <- list()
grab <- function(nm, e) {
  out[[nm]] <<- if (inherits(e, "error")) paste("ERROR:", msg(e)) else
    as.data.frame(e)
  cat(sprintf("  %-22s %s\n", nm,
              if (inherits(e, "error")) "ERROR" else "ok"))
}
fg <- suppressWarnings(frm(bf(y ~ f + s(x) + (1 | g)), family = gaussian(),
                          data = d))
grab("gauss_f", tryv(emmeans::emmeans(fg, ~ f)))
grab("gauss_x", tryv(emmeans::emmeans(fg, ~ x)))
grab("gauss_f_null", tryv(emmeans::emmeans(fg, ~ f, re_formula = NULL)))
fo <- suppressWarnings(frm(bf(yo ~ f + s(x)), family = cumulative(),
                          data = d))
grab("ord_f", tryv(emmeans::emmeans(fo, ~ f)))
fn <- suppressWarnings(frm(bf(ynl ~ a * exp(-b * z), a ~ f, b ~ 1,
                             nl = TRUE), family = gaussian(), data = d))
grab("nl_f", tryv(emmeans::emmeans(fn, ~ f, nlpar = "a")))
grab("nl_f_plain", tryv(emmeans::emmeans(fn, ~ f)))
fmv <- suppressWarnings(frm(bf(y ~ f + s(x)) + bf(y2 ~ f + x),
                           family = gaussian(), data = d))
grab("mv_f", tryv(emmeans::emmeans(fmv, ~ f, resp = "y")))
grab("mv_f2", tryv(emmeans::emmeans(fmv, ~ f, resp = "y2")))
saveRDS(out, file.path("dev", paste0("resmooth-rev-emm-", ARM, ".rds")))
cat("DONE\n")
