# Reviewer of lane ordmix: default_prior() rows of ordinal mixtures, the
# lane build against brms 2.23.0, as sets of (class, coef, group, dpar).
# Seed 20261094.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
})
set.seed(20261094)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sample(1:4, n, TRUE)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
key <- function(p) {
  p <- as.data.frame(p)
  sort(unique(paste(p$class, p$coef, p$group, p$dpar, sep = "|")))
}
cases <- list(
  none = list(f = "y ~ x", fam = "mixture(cumulative(), sratio())"),
  mu = list(f = "y ~ x", fam = "mixture(cumulative(), sratio(), order = 'mu')"),
  disc = list(f = "bf(y ~ x, disc1 ~ z)",
              fam = "mixture(cumulative(), cumulative())"),
  gr = list(f = "y | thres(gr = g) ~ x", fam = "mixture(cumulative(), acat())"),
  equi = list(f = "y ~ x", fam = paste0("mixture(cumulative(threshold = ",
                                        "'equidistant'), sratio())")),
  hurdle = list(f = "bf(yh ~ x, hu1 ~ z)",
                fam = "mixture(hurdle_cumulative(), hurdle_cumulative())"),
  theta = list(f = "bf(y ~ x, theta1 ~ z)",
               fam = "mixture(cumulative(), sratio())")
)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  ff <- eval(parse(text = cs$f), list(bf = frmtmb::bf))
  if (inherits(ff, "formula")) ff <- frmtmb::bf(ff)
  fb <- eval(parse(text = cs$f), list(bf = brms::bf))
  if (inherits(fb, "formula")) fb <- brms::bf(fb)
  ours <- tryCatch(key(default_prior(ff, data = d,
                                     family = eval(parse(text = cs$fam),
                                                   asNamespace("frmtmb")))),
                   error = function(e) paste("ERROR", conditionMessage(e)))
  theirs <- key(suppressMessages(brms::default_prior(
    fb, data = d, family = eval(parse(text = cs$fam), asNamespace("brms")))))
  cat(sprintf("== %s: ours %d rows, brms %d rows\n", nm, length(ours),
              length(theirs)))
  cat("   only brms:", paste(setdiff(theirs, ours), collapse = "  "), "\n")
  cat("   only ours:", paste(setdiff(ours, theirs), collapse = "  "), "\n")
}
