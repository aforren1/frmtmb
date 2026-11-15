# Review of lane quadseam2, attack 3: false positives of the limit rule
# on the one-block quadrature route. Identified fits with a level of
# mostly zeros, a binomial level with all or nearly all successes, small
# groups, rare counts; each by Laplace and by quadrature, seeds 1 to 8.
# A verdict is the warning classes plus which SEs are not finite.
#   Rscript quadseam2-review-limit.R <libpath|base>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
verdict <- function(fo, fam, d, quad) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(fo, family = fam, data = d, quadrature = quad),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("m:", class(x)[1])); invokeRestart("muffleMessage")
    }), error = function(e) e)
  if (inherits(f, "error")) return(list(txt = "ERROR", secs = NA))
  t1 <- proc.time()[["elapsed"]]
  se <- tryCatch(withCallingHandlers(sqrt(diag(vcov(f))),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }, message = function(x) {
      w <<- c(w, paste0("m:", class(x)[1])); invokeRestart("muffleMessage")
    }), error = function(e) NA)
  t2 <- proc.time()[["elapsed"]]
  bad <- names(se)[!is.finite(se)]
  lost <- tryCatch(paste(unique(frmtmb:::sdr_of(f)$se_lost), collapse = "+"),
                   error = function(e) "?")
  list(txt = sprintf("nonfinite[%s] lost[%s] w[%s] maxse %.3g",
                     paste(bad, collapse = ","), lost,
                     paste(sort(unique(w)), collapse = ","),
                     max(se[is.finite(se)], -Inf)),
       secs = c(fit = t1 - t0, vcov = t2 - t1))
}
cases <- list(
  # poisson level c with exactly one nonzero count
  pois1 = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
    ic <- which(d$f == "c")
    d$y[ic] <- 0L
    d$y[ic[1]] <- 1L
    list(fo = y ~ f + x + (1 | g), fam = poisson(), d = d)
  },
  # poisson level c with a low but nonzero rate
  poislow = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rpois(180, exp(0.5 + 0.3 * d$x - 3.5 * (d$f == "c") +
                            rnorm(20, 0, 0.5)[d$g]))
    list(fo = y ~ f + x + (1 | g), fam = poisson(), d = d)
  },
  # bernoulli level c all ones: separation, both routes should say so
  bernsep = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rbinom(180, 1, plogis(0.2 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g]))
    d$y[d$f == "c"] <- 1L
    list(fo = y ~ f + x + (1 | g), fam = bernoulli(), d = d)
  },
  # bernoulli level c with one zero among 60
  bern59 = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:20, each = 9)),
                    f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
    d$y <- rbinom(180, 1, plogis(0.2 + 0.5 * d$x + rnorm(20, 0, 0.7)[d$g]))
    ic <- which(d$f == "c")
    d$y[ic] <- 1L
    d$y[ic[3]] <- 0L
    list(fo = y ~ f + x + (1 | g), fam = bernoulli(), d = d)
  },
  # small groups: 40 groups of 2, bernoulli with a factor
  bernsmall = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:40, each = 2)),
                    f = factor(rep(c("a", "b"), 40)), x = rnorm(80))
    d$y <- rbinom(80, 1, plogis(-0.5 + 0.8 * (d$f == "b") + 0.5 * d$x +
                                  rnorm(40, 0, 1)[d$g]))
    list(fo = y ~ f + x + (1 | g), fam = bernoulli(), d = d)
  },
  # rare poisson counts, mean about 0.1
  poisrare = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:30, each = 5)), x = rnorm(150))
    d$y <- rpois(150, exp(-2.3 + 0.4 * d$x + rnorm(30, 0, 0.5)[d$g]))
    list(fo = y ~ x + (1 | g), fam = poisson(), d = d)
  },
  # binomial with trials, level c 10 of 10 in every row
  binsep = function(s) {
    set.seed(s)
    d <- data.frame(g = factor(rep(1:15, each = 6)),
                    f = factor(rep(c("a", "b", "c"), 30)), x = rnorm(90))
    d$n <- 10L
    d$y <- rbinom(90, 10, plogis(0.3 + 0.4 * d$x + rnorm(15, 0, 0.5)[d$g]))
    d$y[d$f == "c"] <- 10L
    list(fo = y | trials(n) ~ f + x + (1 | g), fam = binomial(), d = d)
  }
)
for (nm in names(cases)) for (s in 1:8) {
  cs <- cases[[nm]](s)
  for (quad in c(FALSE, TRUE)) {
    v <- verdict(cs$fo, cs$fam, cs$d, quad)
    cat(sprintf("ROW\t%s\t%d\t%s\t%.2f\t%.2f\t%s\n", nm, s,
                if (quad) "QUAD" else "LAP", v$secs[1], v$secs[2], v$txt))
  }
}
cat("END\n")
