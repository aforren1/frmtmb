## Punch round 1, MAJOR 4, the refusal half: what ranef() and coef() do
## on draws with no group-level draws, with and without the refusal in
## draws_ranef_fill(). Real frm_sample(laplace = TRUE) draws, and the
## column-removal construction test-draws-methods.R uses.
##   Rscript dev/brmsnames-probe-laplace.R      data seed 5, sampler seed 8
source("dev/brmsnames-libs.R")
brmsnames_libs("lane")
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
q(library(frmtmb.sample))
set.seed(5)
d <- data.frame(x = stats::rnorm(200), g = factor(rep(1:10, 20)))
d$y <- 1 + d$x + stats::rnorm(10)[d$g] + stats::rnorm(200)
fit <- q(frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d))
lap <- q(frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 8,
                    laplace = TRUE))
full <- q(frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 8))
cut <- full
cut$draws <- cut$draws[, !startsWith(colnames(cut$draws), "r_"),
                       drop = FALSE]
show <- function(lbl, e) {
  r <- tryCatch({
    v <- q(e)
    sprintf("returned, NA cells %d of %d", sum(is.na(unlist(v))),
            length(unlist(v)))
  }, error = function(err) paste("ERROR:", conditionMessage(err)))
  cat(sprintf("%-34s %s\n", lbl, substr(gsub("\n", " ", r), 1, 160)))
}
run <- function(tag) {
  cat("--", tag, "--\n")
  cat("laplace real:", frmtmb.sample:::draws_is_laplace(lap),
      " column cut:", frmtmb.sample:::draws_is_laplace(cut), "\n")
  show("ranef(real laplace)", ranef(lap))
  show("coef(real laplace)", coef(lap))
  show("ranef(column cut)", ranef(cut))
  show("coef(column cut)", coef(cut))
}
run("as built")
ns <- asNamespace("frmtmb.sample")
src <- deparse(get("draws_ranef_fill", ns), width.cutoff = 500L)
hit <- grep("if (draws_is_laplace(x)) {", src, fixed = TRUE)
stopifnot(length(hit) == 1L)
src[hit] <- sub("if (draws_is_laplace(x)) {", "if (FALSE) {", src[hit],
                fixed = TRUE)
g <- eval(parse(text = src))
environment(g) <- ns
assignInNamespace("draws_ranef_fill", g, "frmtmb.sample")
run("refusal removed")
cat("DONE\n")
