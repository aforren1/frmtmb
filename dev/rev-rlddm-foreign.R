# REVIEW round 2, attack 2: B2's fix refuses all five eam families. A
# FOREIGN family is still accepted. Is that the point of the seam, or
# the same hole one package further out?
#
#   Rscript dev/rev-rlddm-foreign.R <lib>
#
# Seed 707. Constructed rather than argued: two foreign families that
# both pass the guard, identical in every respect except one line.
#
#   ok      calls ndt_apply() where it reads the non-decision time
#   broken  reads dpars[["ndt"]] directly, which is the consumer who
#           attached a bound and did not read the rest of the page
#
# A shifted GAMMA, y = ndt + Gamma(shape, rate) with shape above one,
# because the whole question is what a density does with `ndt` and
# nothing here needs a diffusion to ask it. Shape above one is not
# cosmetic: a shifted EXPONENTIAL has its maximum at the sample
# minimum, so both arms pin `ndt` at the floor and the two cannot be
# told apart. With shape 3 the density is zero at the shift and the
# maximum is interior, which is the case a real response-time family
# is in.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("lib:", lib, " eam", format(packageVersion("frmtmb.eam")), "\n\n")

shifted_exp <- function(reads_floor) {
  fam <- custom_family(
    "shiftgam",
    dpars = c("rate", "shape", "ndt"),
    primary_dpars = "rate",
    links = list(rate = "log", shape = "log", ndt = "log"),
    accepts_aterms = c("weights", "ndt_group"),
    lpdf = function(y, dpars, aterms) {
      nd <- if (reads_floor) {
        ndt_apply(dpars, aterms, "shiftgam")
      } else {
        dpars[["ndt"]]
      }
      z <- y - nd
      sh <- dpars[["shape"]]
      rt <- dpars[["rate"]]
      (sh - 1) * log(z) - rt * z + sh * log(rt) - lgamma(sh)
    },
    init_dpars = list(rate = function(y, aterms) 3 / mean(y),
                      shape = function(y, aterms) 3,
                      ndt = function(y, aterms) 0.5 * min(y)),
    family_finalize = function(fam, y, aterms) {
      bd <- ndt_bound(y, aterms, NULL, "shiftgam")
      if (!is.null(ndt_bound_of(fam))) return(fam)
      ndt_bound_attach(fam, bd)
    })
  ndt_bound_attach(fam, ndt_bound_pending(NULL, "shiftgam"))
}

set.seed(707)
ng <- 4L
nper <- 300L
shift <- c(a = 0.20, b = 0.28, c = 0.36, d = 0.44)
rate <- 20
g <- rep(names(shift), each = nper)
y <- shift[g] + stats::rgamma(ng * nper, shape = 3, rate = rate)
d <- data.frame(y = as.numeric(y), g = factor(g))
own <- as.numeric(tapply(d$y, d$g, min))
names(own) <- levels(d$g)
cat("true shifts   :", paste(format(shift, digits = 4), collapse = " "),
    "\n")
cat("group floors  :", paste(format(own, digits = 6), collapse = " "),
    "\n")
cat("global min(y) :", format(min(d$y), digits = 6), "\n\n")

cat("-- does the guard let a foreign family through? --\n")
for (nm in c("ok", "broken")) {
  f <- tryCatch(shifted_exp(identical(nm, "ok")),
                error = function(e) e)
  cat(sprintf("  shiftexp(%-6s) attach at construction: %s\n", nm,
              if (inherits(f, "error"))
                paste("REFUSED:", substr(conditionMessage(f), 1, 60))
              else "ACCEPTED"))
}

one <- d[match(levels(d$g), as.character(d$g)), , drop = FALSE]
res <- list()
for (nm in c("ok", "broken")) {
  cat("\n--", nm, "--\n")
  out <- tryCatch({
    fit <- frm(bf(y | ndt_group(g) ~ 1, shape ~ 1, ndt ~ 0 + g),
               family = shifted_exp(identical(nm, "ok")), data = d)
    fr <- as.numeric(suppressWarnings(stats::predict(
      fit, newdata = one, dpar = "ndt", type = "response")))
    tm <- as.numeric(suppressWarnings(ndt_time(fit, newdata = one)))
    cat("  fit RETURNED, logLik",
        format(as.numeric(stats::logLik(fit)), digits = 9), "\n")
    cat("  predict(dpar='ndt', type='response'):",
        paste(format(fr, digits = 5), collapse = " "), "\n")
    cat("  ndt_time()                         :",
        paste(format(tm, digits = 5), collapse = " "), "\n")
    cat("  the truths                         :",
        paste(format(as.numeric(shift), digits = 5), collapse = " "),
        "\n")
    cat("  max abs error of ndt_time vs truth :",
        format(max(abs(tm - shift)), digits = 5), "\n")
    cat("  max abs error of predict  vs truth :",
        format(max(abs(fr - shift)), digits = 5), "\n")
    list(fr = fr, tm = tm)
  }, error = function(e) {
    cat("  REFUSED:", substr(conditionMessage(e), 1, 200), "\n")
    NULL
  })
  res[[nm]] <- out
}

cat("\n-- the same family with a SCALAR bound, no ndt_group --\n")
for (nm in c("ok", "broken")) {
  out <- tryCatch({
    fit <- frm(bf(y ~ 1, shape ~ 1, ndt ~ 1),
               family = shifted_exp(identical(nm, "ok")), data = d)
    format(as.numeric(suppressWarnings(ndt_time(fit)))[1L], digits = 6)
  }, error = function(e) paste("REFUSED:", substr(conditionMessage(e),
                                                  1, 60)))
  cat(sprintf("  %-6s ndt_time = %s   (global floor %s)\n", nm, out,
              format(min(d$y), digits = 6)))
}
