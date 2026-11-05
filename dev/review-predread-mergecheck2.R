# Item 8: localize simulate(newdata) on a latent term on the trial merge: is the error
# me()'s path, latres's, or the two together? Data seed 8101.
#   REV_LIB=<merge lib> Rscript dev/review-predread-mergecheck2.R
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(8101)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(sprintf("g%02d", 1:30), each = 10)),
                lat = rnorm(n), sx = 0.4)
d$t <- rep(1:10, 30)
d$xo <- d$lat + rnorm(n, 0, 0.4)
d$cnt <- rpois(n, exp(0.5 + 0.7 * d$lat))
nd <- data.frame(x = c(-1, 0, 1), xo = c(-1, 0.5, 1), sx = c(0.3, 0.6, 0.9), t = 1:3,
                 g = factor("g01", levels = levels(d$g)))
for (f in list(bf(cnt ~ x + ar(time = t, gr = g, cov = TRUE)),
               bf(cnt ~ me(xo, sx) + ar(time = t, gr = g, cov = TRUE)),
               bf(cnt ~ me(xo, sx) + (1 | g)),
               bf(cnt ~ me(xo, sx)))) {
  fit <- quiet(frm(f, family = poisson(), data = d))
  cat("---", deparse(f$formula), "\n")
  r <- tryCatch(as.matrix(simulate(fit, nsim = 2, seed = 1, newdata = nd)), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    cat(paste(utils::head(sapply(sys.calls(), function(cl) substr(deparse(cl)[1], 1, 70)), 0), collapse = "\n"))
    NULL
  })
  if (!is.null(r)) print(r)
  tb <- tryCatch(withCallingHandlers(simulate(fit, nsim = 1, seed = 1, newdata = nd), error = function(e) {
    cs <- sys.calls()
    cat("CALLS:\n", paste(utils::tail(vapply(cs, function(cl) substr(deparse(cl)[1], 1, 60), ""), 8), collapse = "\n "), "\n")
  }), error = function(e) NULL)
}
