# Review of lane arstep: does the save/restore around reading the second
# start (gp_start_scan()) matter, and does obj_settle() in
# gp_second_descent()? Same process, interleaved, on the 40 rows of the
# ncgp review's 150 where the lane's trace shows a second descent
# (dev/arstep-out/trace-ng-rows.txt), data as dev/arstep-ng.R, seed 1.
#   lane        as shipped
#   norestore   gp_start_scan() without its on.exit() restore
#   nosettle    gp_second_descent() without obj_settle()
# Prints ll and the max relative difference of estimates$b and fitted()
# from the lane arm.
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
rv_libs("lane")
ns <- asNamespace("frmtmb")
edit_fn <- function(f, from, to) {
  txt <- deparse(f)
  i <- grep(from, txt, fixed = TRUE)
  stopifnot(length(i) == 1L)
  txt[i] <- gsub(from, to, txt[i], fixed = TRUE)
  eval(parse(text = txt))
}
o_scan <- ns$gp_start_scan
o_sd <- ns$gp_second_descent
m_scan <- edit_fn(o_scan, "on.exit(obj_state_restore(obj, saved), add = TRUE)", "invisible(NULL)")
m_sd <- edit_fn(o_sd, "obj_settle(obj, opt, best, saved)", "invisible(NULL)")
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, wig = function(x) sin(7 * x) + 0.2 * x,
               two = function(x) sin(x) + 0.4 * sin(10 * x), smooth = function(x) sin(x),
               step = function(x) 2 * (x > 5) - 1)
mk <- function(case, x, f, n) {
  eta <- f(x)
  switch(case,
    pois_low = list(d = data.frame(x = x, y = rpois(n, exp(eta / 2))), fml = y ~ gp(x), fam = poisson()),
    pois_mid = list(d = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))), fml = y ~ gp(x), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fml = y ~ gp(x), fam = bernoulli()),
    binom5 = list(d = data.frame(x = x, tr = 5, y = rbinom(n, 5, plogis(eta - mean(eta)))), fml = y | trials(tr) ~ gp(x), fam = binomial()),
    negbin = list(d = data.frame(x = x, y = rnbinom(n, size = 20, mu = exp(4 + eta))), fml = y ~ gp(x), fam = negbinomial()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fml = y ~ gp(x), fam = Gamma(link = "log")),
    lognorm = list(d = data.frame(x = x, y = exp(eta + rnorm(n, sd = 0.02))), fml = y ~ gp(x), fam = lognormal()),
    student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
    gsigma = list(d = data.frame(x = x, z = rnorm(n), y = eta + 0.02 * rnorm(n)), fml = bf(y ~ gp(x), sigma ~ z), fam = gaussian()))
}
tr <- readLines(file.path(RV_WT, "dev/arstep-out/trace-ng-rows.txt"))
tr <- tr[grepl(" gap ", tr)]
fe <- function(f) { v <- fitted(f); as.numeric(if (is.matrix(v) || is.data.frame(v)) v[, 1] else v) }
for (ln in tr) {
  s <- strsplit(trimws(ln), " +")[[1]]
  case <- s[1]; cn <- s[2]; n <- as.integer(s[3])
  if (n > 80L) next  # keep the run short: the n 40 and 80 rows
  set.seed(74000 + 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  m <- mk(case, x, curves[[cn]], n)
  res <- list()
  for (arm in c("lane", "norestore", "nosettle")) {
    rv_put("gp_start_scan", if (arm == "norestore") m_scan else o_scan)
    rv_put("gp_second_descent", if (arm == "nosettle") m_sd else o_sd)
    q <- rv_quiet(frm(m$fml, data = m$d, family = m$fam))
    res[[arm]] <- q$value
  }
  rv_put("gp_start_scan", o_scan); rv_put("gp_second_descent", o_sd)
  L <- res$lane
  dd <- function(f) {
    if (inherits(f, "error")) return("ERROR")
    sprintf("ll %.8f db %.1e dfit %.1e", as.numeric(logLik(f)),
            max(abs(f$estimates$b - L$estimates$b)) / max(abs(L$estimates$b)),
            max(abs(fe(f) - fe(L))) / max(abs(fe(L))))
  }
  cat(sprintf("RS %-8s %-6s %3d lane ll %.8f second %s | norestore %s | nosettle %s\n",
              case, cn, n, as.numeric(logLik(L)),
              paste(signif(L$opt$gp_second, 4), collapse = ","),
              dd(res$norestore), dd(res$nosettle)))
  flush.console()
}
cat("DONE\n")
