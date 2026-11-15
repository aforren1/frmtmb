# Lane wcheck2, backlog 4993: how check_re_structure() names a term of
# a multivariate model in its observation-level and single-level
# warnings, against the formula as written.
#   Rscript dev/wcheck2-olre-label.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
set.seed(20261111)
n <- 40
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), one = factor(1),
                g = factor(rep(1:8, each = 5)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$y2 <- 0.5 * d$x + rnorm(n)
d$y_a <- d$y
warns <- function(expr) {
  ws <- character()
  withCallingHandlers(tryCatch(expr, error = function(e) {
    ws <<- c(ws, paste("ERROR", conditionMessage(e)))
  }), warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  ws
}
show <- function(lab, fo, ...) {
  fr <- q(frm(fo, data = d, dry_run = "frame", ...))
  labs <- vapply(fr$re_blocks, function(bk) bk$term_label %||% "", "")
  cat(lab, "\n  term labels:", paste0("[", labs, "]", collapse = " "), "\n")
  ws <- warns(frm(fo, data = d, ...))
  ws <- ws[grepl("own random effect|single level", ws)]
  for (w in ws) cat("  warning:", substr(gsub("\\s+", " ", w), 1, 90), "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a
show("univariate (1 | obs)", bf(y ~ x + (1 | obs)), family = gaussian())
show("mvbf y, y2 each (1 | obs)",
     mvbf(bf(y ~ x + (1 | obs)), bf(y2 ~ x + (1 | g))),
     family = gaussian())
show("mvbf y_a, y2, (1 | obs) on y2",
     mvbf(bf(y_a ~ x + (1 | g)), bf(y2 ~ x + (1 | obs))),
     family = gaussian())
show("mvbf y, y2, (1 | one) on y (single level)",
     mvbf(bf(y ~ x + (1 | one)), bf(y2 ~ x + (1 | g))),
     family = gaussian())
show("mvbf y, y2, sigma ~ (1 | g) on y",
     mvbf(bf(y ~ x + (1 | obs), sigma ~ (1 | one)), bf(y2 ~ x)),
     family = gaussian())
cat("DONE\n")
