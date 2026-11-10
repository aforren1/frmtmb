# Lane support, review finding 1: the random-effect fits of the
# reviewer's dev/review-support-is.R and is2.R, with opt$par and logLik()
# in hexadecimal, so that two builds compare bit for bit. Run once per
# build, pinned to one core (dev/support-pin.bat), and diff.
#
#   Rscript dev/support-fitsame-re.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
rgev <- function(n, mu, sigma, xi) {
  mu + sigma * ((-log(runif(n)))^(-xi) - 1) / xi
}
hex <- function(v) paste(sprintf("%a", v), collapse = " ")
run <- function(label, expr) {
  r <- tryCatch(suppressWarnings(expr), error = function(e) e)
  if (inherits(r, "error")) {
    cat(label, "ERROR", substr(conditionMessage(r), 1, 160), "\n")
  } else {
    cat(label, "par", hex(r$opt$par), "\n")
    cat(label, "logLik", hex(as.numeric(logLik(r))),
        sprintf("(%.12g)", as.numeric(logLik(r))), "\n")
  }
}
sim <- function(seed, xi) {
  set.seed(seed)
  g <- rep(seq_len(15), each = 8)
  u <- rnorm(15, 0, 0.7)
  d <- data.frame(g = factor(g), y = rgev(120, 1 + u[g], 1, xi))
  d$h <- factor(paste(d$g, rep(1:2, length.out = nrow(d))))
  d
}
gev <- gen_extreme_value()
for (seed in 1:3) for (xi in c(0.5, 0.8)) {
  d <- sim(seed, xi)
  lab <- sprintf("seed %d xi %.1f", seed, xi)
  run(paste(lab, "laplace"), frm(y ~ 1 + (1 | g), family = gev, data = d))
  run(paste(lab, "importance"), frm(y ~ 1 + (1 | g), family = gev, data = d,
                                    importance = 500L))
  run(paste(lab, "quad1"), frm(y ~ 1 + (1 | g), family = gev, data = d,
                               quadrature = TRUE))
  run(paste(lab, "quadnested"), frm(y ~ 1 + (1 | g) + (1 | h), family = gev,
                                    data = d, quadrature = TRUE))
}
for (seed in c(2, 4)) {
  d <- sim(seed, 1.5)
  run(sprintf("seed %d xi 1.5 importance", seed),
      frm(y ~ 1 + (1 | g), family = gev, data = d, importance = 500L))
}
