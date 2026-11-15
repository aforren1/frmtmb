# Short check of lane priors2's fix round, item 1: five ordinary fits with
# no bounds (gaussian, glmm, mixture, ordinal, nonlinear) on fresh data,
# bitwise against 0.82.0. One line per fit: code, objective, logLik and
# every outer parameter in hex.
# Rscript priors2-short-bitwise.R <lane|base> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/rellib-r22" else
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(23)
n <- 200
d <- data.frame(x = rnorm(n), z = runif(n), g = factor(rep(1:25, 8)))
d$y <- -0.5 + 0.8 * d$x - 0.4 * d$z + rnorm(25, sd = 0.7)[d$g] +
  rnorm(n, sd = 1.1)
k <- rbinom(n, 1, 0.35)
d$ym <- ifelse(k == 1, rnorm(n, -1.5 + 0.5 * d$x, 0.8),
               rnorm(n, 2, 1.2))
d$o <- factor(cut(d$y + rlogis(n), c(-Inf, -1, 0, 1.5, Inf),
                  labels = FALSE), ordered = TRUE)
d$xn <- runif(n, 0, 4)
d$ynl <- 5 * (1 - exp(-0.7 * d$xn)) + rnorm(n, sd = 0.3)
cases <- list(
  gaussian = function() frm(y ~ x + z, data = d),
  glmm = function() frm(y ~ x + z + (1 + x | g), data = d),
  mixture = function() frm(ym ~ x, data = d,
                           family = mixture(gaussian(), gaussian())),
  ordinal = function() frm(o ~ x + z, data = d, family = cumulative()),
  nonlinear = function() frm(bf(ynl ~ a * (1 - exp(-k * xn)), a ~ 1, k ~ 1,
                                nl = TRUE), data = d,
                             prior = set_prior("normal(4, 3)", nlpar = "a") +
                               set_prior("normal(1, 1)", nlpar = "k"))
)
hx <- function(v) paste(sprintf("%a", as.numeric(v)), collapse = ",")
out <- character(0)
for (nm in names(cases)) {
  r <- tryCatch(suppressWarnings(suppressMessages(cases[[nm]]())),
                error = function(e) e)
  line <- if (inherits(r, "error")) {
    paste(nm, "ERROR", gsub("[\t\n]", " ", conditionMessage(r)), sep = "\t")
  } else {
    paste(nm, r$opt$convergence, hx(r$opt$objective),
          hx(as.numeric(logLik(r))), hx(r$opt$par), sep = "\t")
  }
  cat(substr(line, 1, 120), "\n")
  out <- c(out, line)
}
writeLines(out, args[2])
