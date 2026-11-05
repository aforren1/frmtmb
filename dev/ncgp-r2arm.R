# Bernoulli and low-count poisson stop short on rows 0.78.0 reaches: does taking the lower of
# the working start and the scan's start, by the objective, close them without losing the
# rows the working start won? Patched in process, against the best of 18 starts in
# dev/ncgp-out/r2/ng.txt. Seeds as dev/review-ncgp-r2-ng.R.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
patch <- function(nm, f) { environment(f) <- ns; unlockBinding(nm, ns); assign(nm, f, envir = ns); lockBinding(nm, ns) }
scan0 <- ns$gp_start_scan
arm <- function(obj, template, frame, start, to_beta = identity) {
  bks <- Filter(function(bk) identical(bk[["covstruct"]], "gp"), frame[["re_blocks"]] %||% list())
  if (!length(bks) || !is.null(start[["theta"]])) return(obj)
  val <- function(tpl) {
    p <- outer_from_template(tpl, obj, frame)
    if (is.null(p)) return(NA_real_)
    v <- tryCatch(obj$fn(p), error = function(e) NA_real_)
    if (isTRUE(is.finite(v))) v else NA_real_
  }
  m <- gp_marginal_max(frame)
  if (!is.null(m)) return(scan0(obj, template, frame, start, to_beta))
  w <- tryCatch(gp_start_working(frame), error = function(e) NULL)
  covered <- integer(0)
  tw <- NULL
  if (!is.null(w)) {
    tpl <- template
    tpl[["theta"]][w$idx] <- w$theta
    if (!is.na(val(tpl))) { tw <- tpl; covered <- w$idx }
  }
  obj <- gp_start_scan_units(obj, template, frame, val, integer(0))
  ps <- obj$par
  vs <- tryCatch(obj$fn(ps), error = function(e) NA_real_)
  if (!is.null(tw)) {
    vw <- val(tw)
    if (is.na(vs) || (!is.na(vw) && vw <= vs)) {
      obj <- gp_start_scan_units(obj, tw, frame, val, covered)
    }
  }
  obj
}
patch("gp_start_scan", arm)
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x, two = function(x) sin(x) + 0.4 * sin(10 * x))
mk <- function(case, x, f, n) {
  eta <- f(x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))), fml = y ~ gp(x), fam = poisson()),
    pois_high = list(d = data.frame(x = x, y = rpois(n, exp(7 + eta))), fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fml = y ~ gp(x), fam = bernoulli()),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20, mu = exp(4 + eta))), fml = y ~ gp(x), fam = negbinomial()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fml = y ~ gp(x), fam = Gamma(link = "log")),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()))
}
ref <- readLines("dev/ncgp-out/r2/ng.txt")
for (cs in c("pois_low", "pois_high", "bern", "negbin", "gamma", "student")) for (cn in names(curves)) for (n in c(60L, 120L)) {
  set.seed(56000 + 13 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  m <- mk(cs, x, curves[[cn]], n)
  r <- rev_fit(m$fml, m$d, family = m$fam)
  ll <- as.numeric(logLik(r$fit))
  line <- grep(sprintf("^RNG2 lane %s %s n %d seed 1 ", cs, cn, n), ref, value = TRUE)
  bst <- as.numeric(sub(".*\\| best ([-0-9.]+) theta.*", "\\1", line))
  cat(sprintf("RAR %s %s n %d ll %.4f best %.4f gap %.2e conv %d w %d\n", cs, cn, n, ll, bst, ll - max(bst, ll), r$fit$opt$convergence, length(r$warn)))
  flush.console()
}
