# Short check of fix round 2, items 3 and 5, as data for the merge rules. Not a review
# of lane ncgp. On wt-ncgp-lib: (a) frm_sample(fit) from the mode where frm() has hit
# the 1e-6 floor: does the stuck chain of F2 come back? (b) the warning's rule read on
# that tape's chart, r = nu * diag(L^-T (H - I) L^-1), against the closed value
# 1e-12 sdgp^2 / sigma^2. Seed 8901, n 60 for (b), seed 6601, n 15 for (a).
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-ncgpinit.R ncgpinit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb:", dirname(find.package("frmtmb")), " sample:", dirname(find.package("frmtmb.sample")), "\n")
set.seed(8901)
n <- 60
x <- round(runif(n, 0, 6), 3)
e <- rnorm(n)
for (s in c(1e-3, 1e-4, 3e-5, 1e-5, 5e-6, 3e-6, 1e-6, 1e-7, 0)) {
  d <- data.frame(x = x, y = sin(2 * x) + s * e)
  fr <- suppressMessages(suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d)))
  bk <- Filter(function(b) identical(b$covstruct, "gp"), fr$frame$re_blocks)[[1]]
  th <- fr$estimates$theta[bk$theta_idx]
  env <- fr$obj$env
  par <- env$last.par.best
  ri <- which(names(par)[env$random] == "b")[bk$b_idx]
  r <- tryCatch({
    H <- as.matrix(env$spHess(par, random = TRUE))[ri, ri]
    L <- frmtmb:::gp_factor(th, bk)
    Li <- solve(L)
    W <- t(Li) %*% (H - diag(nrow(H))) %*% Li
    stats::quantile(1e-12 * pmax(diag(W), 0), 0.9, names = FALSE)
  }, error = function(er) NA_real_)
  sg <- as.numeric(sigma(fr)); sd <- exp(th[[1]])
  cat(sprintf("RCHART s %.0e | frm sigma %.3g sdgp %.3g conv %d | r90 on the chart %.3g, 1e-12 sd^2/sigma^2 %.3g, rule %.3f, fires %s\n",
              s, sg, sd, fr$opt$convergence, r, 1e-12 * sd^2 / sg^2, sqrt(2 / n), isTRUE(r >= sqrt(2 / n))))
}
set.seed(6601)
n <- 15
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
for (s in c(1e-3, 1e-6, 0)) {
  d <- data.frame(x = x, y = sin(2 * x) + s * e0)
  fr <- suppressMessages(suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d)))
  for (ini in list(NULL, "random")) {
    t0 <- proc.time()[["elapsed"]]
    msg <- character()
    args <- list(fr, chains = 4, iter = 1000, warmup = 500, seed = 11, refresh = 0, cores = 1,
                 control = list(max_treedepth = 10))
    if (!is.null(ini)) args$init <- ini
    fs <- withCallingHandlers(do.call(frm_sample, args),
      message = function(m) {
        msg <<- c(msg, conditionMessage(m))
        invokeRestart("muffleMessage")
      }, warning = function(w) invokeRestart("muffleWarning"))
    dv <- posterior::as_draws_df(fs)
    m <- posterior::extract_variable_matrix(dv, "sigma")
    sp <- rstan::get_sampler_params(fs$stanfit, inc_warmup = FALSE)
    cat(sprintf("NINIT s %.0e frm sigma %.3g conv %d init %-14s says random %-5s | sigma q50 %.3g rhat %.3f stepsize %.2g | divergent %d | %.0f s\n",
                s, as.numeric(sigma(fr)), fr$opt$convergence, if (is.null(ini)) "default" else ini,
                any(grepl("start at random", msg, fixed = TRUE)), median(as.numeric(m)),
                posterior::rhat(m), mean(vapply(sp, function(q) mean(q[, "stepsize__"]), 0)),
                sum(vapply(sp, function(q) sum(q[, "divergent__"]), 0)),
                proc.time()[["elapsed"]] - t0))
  }
}
