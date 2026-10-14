# Probe frmtmb.eam's wiener (worktree source via load_all) for the brms
# parity audit: which mu links it accepts, the dpar links it reports,
# and the default priors frmtmb::default_prior() lists on both routes.
# Output: dev/famuni-out/wiener-eamprobe.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
suppressMessages(library(frmtmb.sample))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "\n")
set.seed(20261013)
n <- 80
d <- data.frame(x = rnorm(n))
sim <- RWiener::rwiener(n, alpha = 1.4, tau = 0.25, beta = 0.5, delta = 0.8)
d$rt <- sim$q
d$d <- as.integer(sim$resp == "upper")

cat("\n== formals(frmtmb.eam::wiener)\n"); print(names(formals(wiener)))
cat("\n== mu links accepted by wiener(link = )\n")
lks <- c("identity", "log", "softplus", "squareplus", "logit", "inverse",
         "probit", "sqrt")
for (lk in lks) {
  r <- tryCatch({
    f <- wiener(link = lk)
    fit <- suppressWarnings(frm(bf(rt | dec(d) ~ x), family = f, data = d))
    sprintf("ok  link=%s  logLik=%.6f  conv=%s", f$links$mu$name %||%
              as.character(f$links$mu), as.numeric(logLik(fit)),
            fit$opt$convergence %||% NA)
  }, error = function(e) paste("ERR", conditionMessage(e)))
  cat(sprintf("%-10s %s\n", lk, substr(r, 1, 120)))
}
cat("\n== link_* arguments\n")
for (arg in c("link_bs", "link_ndt", "link_bias")) {
  r <- tryCatch({do.call(wiener, setNames(list("log"), arg)); "ok"},
                error = function(e) paste("ERR", conditionMessage(e)))
  cat(arg, ":", substr(r, 1, 120), "\n")
}

for (route in c("fit", "sample")) {
  for (fo in list(bf(rt | dec(d) ~ x, bs ~ x, ndt ~ 1, bias ~ 1),
                  bf(rt | dec(d) ~ 1, bs ~ 1, ndt ~ 1, bias ~ 1),
                  bf(rt | dec(d) ~ x))) {
    cat("\n== default_prior route =", route, " ", deparse(fo$formula), "+",
        paste(names(fo$pforms), collapse = ","), "\n")
    r <- tryCatch(print(frmtmb::default_prior(fo, data = d,
                                              family = wiener(),
                                              route = route)),
                  error = function(e) cat("ERR", conditionMessage(e), "\n"))
  }
}
