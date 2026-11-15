# Items 4846 and 561, route (a'): the error of the one-block Gauss-Kronrod
# rule at its anchor, TMB's raw rule against a "whitened" one. The
# whitened tape integrates u with b = bhat + s * u, bhat and s the
# Laplace conditional mode and scale at the anchor (constants on the
# tape), so the integrand TMB calibrates (by differences of step dx = 1,
# mode search stopped at a gain of ytol = 1e-2) has unit width and mode 0.
# Both rules are taped at each fit's own estimate and read there, against
# stats::integrate() per group (ref_ll() of
# dev/reviews/quadseam-review-acc.R). The fit itself is the build's own.
#   Rscript dev/quadseam2-white.R <libpath|base> <design,...> <seeds>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
designs <- strsplit(a[2], ",")[[1]]
seeds <- eval(parse(text = a[3]))
ns <- asNamespace("frmtmb")
src <- readLines("dev/reviews/quadseam-review-acc.R")
# gen(), fo, frm_fam and ref_ll() from the review's script
eval(parse(text = src[grep("^gen <- function", src):
                        (grep("^for \\(des in designs\\)", src) - 1L)]))

rule_at <- function(f, white) {
  frame <- f$frame
  nll <- ns$build_objective(frame)
  est <- f$estimates
  map <- frame[["map"]]
  ictl <- ns$predictor_inner_control(ns$predictor_units(frame))
  q <- f$opt$par
  tpl <- ns$imp_par_list(est, q, map)
  if (white) {
    lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = map, silent = TRUE,
                           inner.control = ictl)
    lap$fn(q)
    full <- lap$env$last.par
    bhat <- as.numeric(full[lap$env$random])
    h <- Matrix::diag(lap$env$spHess(full, random = TRUE))
    s <- ifelse(is.finite(h) & h > 0, 1 / sqrt(h), 1)
    ls <- sum(log(s))
    nll_w <- function(pars) {
      pars[["b"]] <- bhat + s * pars[["b"]]
      nll(pars) - ls
    }
    tpl$b[] <- 0
    o <- RTMB::MakeADFun(nll_w, tpl, random = "b", map = map,
                         integrate = ns$quad_gk_integrate(), silent = TRUE,
                         inner.control = ictl)
  } else {
    o <- RTMB::MakeADFun(nll, tpl, random = "b", map = map,
                         integrate = ns$quad_gk_integrate(), silent = TRUE,
                         inner.control = ictl)
  }
  t0 <- proc.time()[[3]]
  v <- -o$fn(q)
  g <- max(abs(o$gr(q)))
  c(ll = v, grad = g, secs = proc.time()[[3]] - t0)
}

for (des in designs) for (s in seeds) {
  g <- gen(des, s)
  fam <- g$fam
  f <- tryCatch(suppressMessages(suppressWarnings(
    frm(fo[[fam]], family = frm_fam[[fam]], data = g$d, quadrature = TRUE))),
    error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("ROW\t%s\t%d\tERROR\t%s\n", des, s,
                substr(conditionMessage(f), 1, 120)))
    next
  }
  beta <- unname(f$estimates$beta); sd <- exp(unname(f$estimates$theta[1]))
  shape <- if (fam == "gam") exp(unname(f$estimates$betad[1])) else NA
  rf <- ref_ll(fam, g$d, beta, sd, shape)
  raw <- rule_at(f, FALSE)
  wh <- rule_at(f, TRUE)
  cat(sprintf(paste0("ROW\t%s\t%d\tsd\t%.4g\tfit\t%.3e\traw\t%.3e\t",
                     "white\t%.3e\trawgr\t%.2e\twhitegr\t%.2e\n"),
              des, s, sd, as.numeric(logLik(f)) - rf, raw[["ll"]] - rf,
              wh[["ll"]] - rf, raw[["grad"]], wh[["grad"]]))
}
cat("END\n")
