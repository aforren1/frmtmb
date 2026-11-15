# Review of lane gpfield: does lane predraw2's draw check (its library,
# wt-predraw2-lib, read-only) fire on gp() fits? And, on the gpfield lane,
# does predict() leave the fit's objective state as it found it?
# Usage: Rscript gpfield-review-predraw2.R <predraw2|lane>
a <- commandArgs(TRUE)[1]
lib <- if (a == "predraw2") "C:/Users/adf44/source/r/wt-predraw2-lib" else
  "C:/Users/adf44/source/r/wt-gpfield-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
ns <- asNamespace("frmtmb")
cat("has draw_check:", exists("draw_check", ns), " has gp_field_draw:",
    exists("gp_field_draw", ns), "\n")
mk <- function(seed, n) {
  set.seed(seed)
  x <- sort(stats::runif(n, 0, 10))
  data.frame(x = x, yg = 2 * sin(x) + stats::rnorm(n, 0, 0.2),
             yp = stats::rpois(n, exp(1 + sin(x))),
             yb = stats::rbinom(n, 1, stats::plogis(1.5 * sin(x))))
}
d <- mk(9200, 30)
d2 <- mk(9201, 80)
rows <- list(
  g30 = list(bf(yg ~ gp(x)), gaussian(), d),
  p80 = list(bf(yp ~ gp(x)), poisson(), d2),
  b80 = list(bf(yb ~ gp(x)), bernoulli(), d2),
  h80 = list(bf(yg ~ gp(x, k = 15)), gaussian(), d2),
  m32 = list(bf(yg ~ gp(x, cov = "matern32")), gaussian(), d))
for (nm in names(rows)) {
  if (a == "predraw2" && nm == "m32") next
  r <- rows[[nm]]
  f <- tryCatch(suppressWarnings(frm(r[[1]], data = r[[3]], family = r[[2]])),
                error = function(e) e)
  if (inherits(f, "error")) { cat(nm, "ERROR", conditionMessage(f), "\n"); next }
  env <- f$obj$env
  st0 <- mget(intersect(c("last.par", "last.par.best", "value.best",
                          "last.par.ok"), ls(env, all.names = TRUE)),
              envir = env)
  msgs <- character(0); wns <- character(0)
  set.seed(1)
  p <- withCallingHandlers(
    predict(f, newdata = data.frame(x = c(2.5, 12)), summary = FALSE,
            ndraws = 500),
    message = function(m) {
      msgs <<- c(msgs, paste(class(m)[1], substr(conditionMessage(m), 1, 90)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      wns <<- c(wns, paste(class(w)[1], substr(conditionMessage(w), 1, 90)))
      invokeRestart("muffleWarning")
    })
  st1 <- mget(names(st0), envir = env)
  same <- vapply(names(st0), function(k) identical(st0[[k]], st1[[k]]), NA)
  chk <- if (is.environment(f$cache)) f$cache[["draw_check"]]
  cat(sprintf("PD %s %s | sd %s | draw_check fire %s excess %s rel %s | msgs [%s] | warns [%s] | state kept: %s\n",
              a, nm, paste(round(apply(p, 2, sd), 3), collapse = " "),
              format(chk$fire), format(chk$excess, digits = 3),
              format(chk$rel, digits = 3), paste(msgs, collapse = " || "),
              paste(wns, collapse = " || "),
              paste(names(same), same, sep = "=", collapse = ",")))
}
