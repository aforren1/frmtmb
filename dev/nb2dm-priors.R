# frm_sample()'s default priors (default_prior(route = "sample")) against
# brms 2.23.0's default_prior() for negbinomial2 and
# dirichlet_multinomial, on the same data.
#   Rscript dev/nb2dm-priors.R [lib] > dev/nb2dm-out/priors.txt 2>&1
a <- commandArgs(trailingOnly = TRUE)
lib <- if (length(a)) a[1] else "C:/Users/adf44/source/r/wt-nb2dm-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r10",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(system.file(package = "frmtmb")),
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "\n")
rows <- function(tab) {
  tab <- as.data.frame(tab)
  stats::setNames(tab$prior, paste(tab$class, tab$coef, tab$group, tab$resp,
                                   tab$dpar, sep = "|"))
}
set.seed(5)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(6, 10))
d$y <- rnbinom(n, mu = exp(1 + 0.3 * d$x), size = 2)
Y <- t(sapply(seq_len(n), function(i) rmultinom(1, 12, c(0.2, 0.3, 0.5))))
colnames(Y) <- c("a", "b", "c")
d$Y <- Y
d$tr <- rowSums(Y)
cases <- list(
  nb2 = list(f = y ~ x, fb = "negbinomial2"),
  nb2_sigma = list(f = frmtmb::bf(y ~ x, sigma ~ z), fb = "negbinomial2",
                   bb = brms::bf(y ~ x, sigma ~ z)),
  nb2_re = list(f = y ~ x + (1 | g), fb = "negbinomial2"),
  dm = list(f = Y | trials(tr) ~ x, fb = "dirichlet_multinomial"),
  dm_phi = list(f = frmtmb::bf(Y | trials(tr) ~ x, phi ~ z),
                bb = brms::bf(Y | trials(tr) ~ x, phi ~ z),
                fb = "dirichlet_multinomial"),
  dm_re = list(f = Y | trials(tr) ~ x + (1 | g), fb = "dirichlet_multinomial")
)
bad <- 0L
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fr <- tryCatch(rows(frmtmb::default_prior(cs$f, d,
                                            family = frmtmb::brmsfamily(cs$fb),
                                            route = "sample")),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  bf_ <- cs$bb %||% cs$f
  br <- rows(brms::default_prior(bf_, d, family = brms::brmsfamily(cs$fb)))
  cat("\n==", nm, "\n")
  if (length(fr) == 1L && startsWith(fr, "ERROR")) {
    cat(fr, "\n")
    bad <- bad + 1L
    next
  }
  keys <- union(names(br), names(fr))
  tab <- data.frame(key = keys, brms = unname(br[keys]),
                    frmtmb = unname(fr[keys]))
  # brms stores "" for a flat slot and for a row that inherits from the
  # class row above it; the rows this lane can move are the ones where
  # brms states a density
  dens <- !is.na(tab$brms) & nzchar(tab$brms)
  tab$same <- ifelse(dens, tab$brms == tab$frmtmb, NA)
  print(tab, row.names = FALSE)
  cat("  rows where brms states a density:", sum(dens), " equal:",
      sum(tab$same, na.rm = TRUE), "\n")
  bad <- bad + sum(!tab$same, na.rm = TRUE)
}
cat("\nrows with a brms density that differ, all cases:", bad, "\n")
