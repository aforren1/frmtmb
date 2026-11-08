# Review of lane arstep: the fit object around the second descent.
# Rows from the lane's trace (dev/arstep-findings.md section 9) where the
# second descent ran: some lost, some won. For each row:
#   default    frm() as shipped on this arm
#   again      the same call a second time (determinism in one process)
#   nosecond   (lane only) gp_second_descent() replaced by "drop the alt,
#              return opt": what the fit would be with the alt read but no
#              second descent, i.e. the base's path plus the alt's restore
#   nosettle   (lane only) gp_second_descent() without obj_settle(): the
#              last.par contamination mutant
#   atopt      (lane only) a fit started at the kept optimum (theta, beta,
#              betad given), which no second start reaches
# Snapshots go to RDS for review-arstep-statecmp.R.
# Usage: Rscript review-arstep-state.R lane|base
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
arm <- commandArgs(TRUE)[1]
rv_libs(arm)
ns <- asNamespace("frmtmb")
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x,
               wig = function(x) sin(7 * x) + 0.2 * x,
               two = function(x) sin(x) + 0.4 * sin(10 * x),
               step = function(x) 2 * (x > 5) - 1)
mk <- function(case, cn, n, seed = 1) {
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  eta <- curves[[cn]](x)
  switch(case,
    pois_mid = list(d = data.frame(x = x, y = rpois(n, exp(1.5 + eta / 2))), fam = poisson()),
    bern = list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta - mean(eta)))), fam = bernoulli()),
    gamma = list(d = data.frame(x = x, y = rgamma(n, 400, 400 / exp(eta))), fam = Gamma(link = "log")),
    binom5 = list(d = data.frame(x = x, tr = 5, y = rbinom(n, 5, plogis(eta - mean(eta)))), fam = binomial()))
}
fml <- function(case) if (case == "binom5") y | trials(tr) ~ gp(x) else y ~ gp(x)
rows <- list(c("bern", "mix", 40), c("bern", "step", 40), c("binom5", "mix", 40),
             c("pois_mid", "wig", 40), c("gamma", "two", 40),
             c("bern", "two", 160))
snap <- function(fit, warns) {
  env <- fit$obj$env
  v <- tryCatch(vcov(fit), error = function(e) conditionMessage(e))
  list(opt = fit$opt, est = fit$estimates,
       last.par = env$last.par, last.par.best = env$last.par.best,
       value.best = env$value.best,
       vcov = v, summ = capture.output(print(summary(fit))),
       pred = tryCatch(predict(fit), error = function(e) conditionMessage(e)),
       fitted = tryCatch(fitted(fit), error = function(e) conditionMessage(e)),
       ll = as.numeric(logLik(fit)), start = fit$start, warns = warns,
       sdr_par = tryCatch(fit$sdr$par.fixed %||% ns$sdr_of(fit)$par.fixed,
                          error = function(e) NULL),
       alt_left = exists(".frmtmb_gp_alt", envir = env, inherits = FALSE))
}
run <- function(m, case, start = NULL) {
  a <- list(fml(case), data = m$d, family = m$fam)
  if (!is.null(start)) a$start <- start
  r <- rv_quiet(do.call(frm, a))
  if (inherits(r$value, "error")) stop(conditionMessage(r$value))
  snap(r$value, r$warnings)
}
orig <- ns$gp_second_descent
out <- list()
for (rw in rows) {
  case <- rw[1]; cn <- rw[2]; n <- as.integer(rw[3])
  key <- paste(case, cn, n)
  m <- mk(case, cn, n)
  res <- list()
  res$default <- run(m, case)
  res$again <- run(m, case)
  if (arm == "lane") {
    rv_put("gp_second_descent", function(obj, opt, ...) {
      if (exists(".frmtmb_gp_alt", envir = obj$env, inherits = FALSE)) {
        rm(".frmtmb_gp_alt", envir = obj$env)
      }
      opt
    })
    res$nosecond <- run(m, case)
    txt <- deparse(orig)
    i <- grep("obj_settle(obj, opt, best, saved)", txt, fixed = TRUE)
    stopifnot(length(i) == 1L)
    txt[i] <- "invisible(NULL)"
    rv_put("gp_second_descent", eval(parse(text = txt)))
    res$nosettle <- run(m, case)
    rv_put("gp_second_descent", orig)
    est <- res$default$est
    st <- list(theta = est$theta, beta = est$beta)
    if (length(est$betad)) st$betad <- est$betad
    res$atopt <- run(m, case, start = st)
  }
  cat(sprintf("ST %s %-18s ll default %.8f again %.8f%s gp_second %s\n", arm, key,
              res$default$ll, res$again$ll,
              if (arm == "lane") sprintf(" nosecond %.8f nosettle %.8f atopt %.8f",
                                         res$nosecond$ll, res$nosettle$ll,
                                         res$atopt$ll) else "",
              paste(format(res$default$opt$gp_second, digits = 4), collapse = ",")))
  flush.console()
  out[[key]] <- res
}
saveRDS(out, file.path(RV_OUT, paste0("state-", arm, ".rds")))
cat("DONE\n")
