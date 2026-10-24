# Reviewer: wrappers around a bar term, refused and not. Arg 1: libs.
.libPaths(c(strsplit(commandArgs(TRUE)[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
set.seed(1)
d <- data.frame(g = gl(10, 6), x = rnorm(60), t = factor(rep(1:3, 20)),
                g1 = gl(10, 6), g2 = factor(sample(1:10, 60, TRUE)),
                h = gl(2, 30))
d$y <- rnorm(60) + rnorm(10)[d$g]
cases <- c(
  "y ~ x + foo(1 | g)",
  "y ~ x + foo((1 | g))",
  "y ~ x + foo(bar(1 | g))",
  "y ~ x + I(1 | g)",
  "y ~ x + us(1 | g)",
  "y ~ x + us(1 + x | g)",
  "y ~ x + diag(1 + x | g)",
  "y ~ x + cs(0 + t | g)",
  "y ~ x + homdiag(1 + x | g)",
  "y ~ x + ar1(0 + t | g)",
  "y ~ x + exp(0 + t | g)",
  "y ~ x + (1 | gr(g))",
  "y ~ x + (1 + x | gr(g, cor = FALSE))",
  "y ~ x + us(1 + x | gr(g, cor = FALSE))",
  "y ~ x + (1 | mm(g1, g2))",
  "y ~ x + diag(1 | mm(g1, g2))",
  "y ~ x + (1 + x || g)",
  "y ~ x + (1 | g1:h)",
  "y ~ x + (1 | h/g)",
  "y ~ x + (1 | gr(g, by = h))",
  "y ~ x + (1 | g) + diag(0 + x | g)",
  "y ~ x + (1 | gr(g, dist = 'student'))",
  "y ~ x + us_t(1 | g)")
for (fo in cases) {
  r <- tryCatch({
    f <- suppressWarnings(frm(stats::as.formula(fo), data = d))
    paste("FIT", paste(vapply(f$frame$re_blocks, `[[`, "", "covstruct"),
                       collapse = ","),
          "logLik", format(as.numeric(logLik(f)), digits = 10))
  }, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 160)))
  cat(sprintf("%-42s -> %s\n", fo, r))
}
dpar <- tryCatch({
  f <- suppressWarnings(frm(bf(y ~ x + (1 | g), sigma ~ foo(1 | g)), data = d))
  paste("FIT", paste(vapply(f$frame$re_blocks, `[[`, "", "covstruct"),
                     collapse = ","))
}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 160)))
cat("bf(y ~ x + (1|g), sigma ~ foo(1 | g)) ->", dpar, "\n")
nl <- tryCatch({
  f <- suppressWarnings(frm(bf(y ~ a + b * x, a ~ 1 + foo(1 | g), b ~ 1,
                               nl = TRUE), data = d,
                            prior = c(set_prior("normal(0, 5)", nlpar = "a"),
                                      set_prior("normal(0, 5)", nlpar = "b"))))
  paste("FIT", paste(vapply(f$frame$re_blocks, `[[`, "", "covstruct"),
                     collapse = ","))
}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 160)))
cat("nl a ~ 1 + foo(1 | g) ->", nl, "\n")
f <- frm(y ~ x + (1 | g), data = d)
rf <- tryCatch({
  p <- predict(f, re_formula = ~ foo(1 | g))
  paste("OK, first row", format(p[1, 1], digits = 8), "vs full",
        format(predict(f)[1, 1], digits = 8))
}, error = function(e) paste("ERR", substr(conditionMessage(e), 1, 160)))
cat("predict(re_formula = ~ foo(1 | g)) ->", rf, "\n")
