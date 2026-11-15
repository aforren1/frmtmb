# Lane gpfield: test-gp-field.R's assertions on y * (1 + 2e-15 k), k = 1 to
# 20 (and k = 0), lane build. Prints one line per row and k: the
# assertion's numbers and its verdict. Args: rows (comma) of
# pois_two, approx_b, pow_a, st80, st40, logis_pois, exp_bern
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
ns <- asNamespace("frmtmb")
rows <- strsplit(a[1], ",")[[1]]
ks <- 0:20
spec <- list(
  pois_two = list(tag = "ng:pois_mid:two:80:1",
                  start = c(log(0.4), log(0.1))),
  approx_b = list(tag = "gate:approx_b", start = c(log(0.05), log(0.056))),
  pow_a = list(tag = "sweep:pow_a:gaussian:120:2:12",
               start = c(log(0.8), log(0.9))),
  st80 = list(tag = "ng:student:step:80:1"),
  st40 = list(tag = "ng:student:step:40:1"),
  logis_pois = list(tag = "sweep:logis_a:poisson:40:1"),
  exp_bern = list(tag = "sweep:exp_a:bernoulli:120:2"))
for (rn in rows) {
  s <- spec[[rn]]
  r0 <- gf_row(s$tag)
  ref <- NULL
  if (!is.null(s$start)) {
    ref <- gf_fit(r0$fml, data = r0$data, family = r0$family,
                  start = list(theta = s$start))$fit
  }
  fails <- 0L
  for (k in ks) {
    d <- r0$data
    if (rn %in% c("pois_two", "logis_pois", "exp_bern")) {
      # a count or a binary response has no ulps to move: perturb x
      d$x <- d$x * (1 + 2e-15 * k)
    } else {
      d$y <- d$y * (1 + 2e-15 * k)
    }
    q <- gf_fit(r0$fml, data = d, family = r0$family)
    f <- q$fit
    cls <- paste(unique(sub(":.*", "", q$warn)), collapse = ";")
    ok <- switch(rn,
      pois_two = , approx_b = , pow_a = {
        l1 <- as.numeric(logLik(ref))
        as.numeric(logLik(f)) >= l1 - 1e-6 * abs(l1)
      },
      st80 = , exp_bern = grepl("frmtmb_gp_laplace", cls),
      st40 = any(grepl("inner Hessian is near singular", q$warn)),
      logis_pois = {
        tk <- ns$gp_laplace_tk(f)
        abs(tk) < 1 && !grepl("frmtmb_gp_laplace", cls)
      })
    if (!isTRUE(ok)) fails <- fails + 1L
    cat(sprintf("ULPS %s k %2d ll %.6f conv %d [%s] ok %s\n", rn, k,
                as.numeric(logLik(f)), f$opt$convergence, cls, ok))
    flush.console()
  }
  cat(sprintf("ULPSUM %s draws %d fail %d\n", rn, length(ks), fails))
}
cat("DONE\n")
