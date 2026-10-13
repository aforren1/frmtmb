# Reviewer, final check (r3): columns and nonlinear parameters named like
# specials, used inside other specials and nested calls; and real
# specials still refused. Each formula is compared with the same model
# where the column has a plain name (logLik), when both fit.
#   Rscript dev/moparity-rev3-names.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3); n <- 300
d <- data.frame(x = sample(0:3, n, TRUE), w = sample(0:2, n, TRUE),
                z = rnorm(n), g = factor(sample(1:20, n, TRUE)))
d$mo <- d$x; d$s <- d$z; d$gp <- d$z; d$ar <- d$z; d$I <- d$z
d$y <- 1 + 0.4 * d$x + sin(d$z) + rnorm(20, 0, 0.5)[d$g] + rnorm(n)
d$ynl <- 2 * exp(-0.5 * abs(d$z)) + rnorm(n, 0, 0.2)
ll <- function(fo, data = d, fam = gaussian()) {
  r <- tryCatch(suppressWarnings(suppressMessages(
    frm(fo, family = fam, data = data))), error = function(e) e)
  if (inherits(r, "error")) paste("ERROR", substr(conditionMessage(r), 1, 110))
  else format(as.numeric(logLik(r)), digits = 12)
}
pair <- function(lab, fo, fo_plain) {
  a <- ll(fo); b <- if (is.null(fo_plain)) "-" else ll(fo_plain)
  same <- !grepl("ERROR", a) && identical(a, b)
  cat(sprintf("%-34s special-named: %s | plain: %s | same %s\n", lab, a, b, same))
}
pair("mo(mo)", bf(y ~ mo(mo)), bf(y ~ mo(x)))
pair("mo(mo) * z", bf(y ~ mo(mo) * z), bf(y ~ mo(x) * z))
pair("mo(x) + mo:z (column mo)", bf(y ~ mo(x) + mo:z), bf(y ~ mo(x) + x:z))
pair("s(gp)", bf(y ~ s(gp)), bf(y ~ s(z)))
pair("gp(s)", bf(y ~ gp(s)), bf(y ~ gp(z)))
pair("s(s)", bf(y ~ s(s)), bf(y ~ s(z)))
pair("mo(x) + (1 + mo | g)", bf(y ~ mo(x) + (1 + mo | g)), bf(y ~ mo(x) + (1 + x | g)))
pair("ar * z (column ar)", bf(y ~ ar * z), bf(y ~ s * z))
pair("I(mo) (column)", bf(y ~ I(mo^2)), bf(y ~ I(x^2)))
pair("I * z (column I)", bf(y ~ I * z), bf(y ~ s * z))
cat("\n-- nonlinear parameters named like specials\n")
pair("nl a, s", bf(ynl ~ a * exp(-s * abs(z)), a ~ 1, s ~ 1, nl = TRUE),
     bf(ynl ~ a * exp(-k * abs(z)), a ~ 1, k ~ 1, nl = TRUE))
pair("nl a, gp", bf(ynl ~ a * exp(-gp * abs(z)), a ~ 1, gp ~ 1, nl = TRUE),
     bf(ynl ~ a * exp(-k * abs(z)), a ~ 1, k ~ 1, nl = TRUE))
pair("nl a, mo", bf(ynl ~ a * exp(-mo * abs(z)), a ~ 1, mo ~ 1, nl = TRUE),
     bf(ynl ~ a * exp(-k * abs(z)), a ~ 1, k ~ 1, nl = TRUE))
cat("\n-- real specials, which must still be refused\n")
for (fo in list(bf(y ~ I(mo(x))), bf(y ~ log(mo(x) + 1)), bf(y ~ mo(x):mo(w)),
                bf(y ~ z * cs(g)), bf(y ~ z:s(z)), bf(y ~ (mo(x) | g)),
                bf(y ~ z + z:ar(g)), bf(y ~ mo(x) + (1 | mo(x))),
                bf(y ~ exp(gp(z))), bf(y ~ mo(x) + mo(x):s(z)))) {
  cat(sprintf("%-34s %s\n", deparse1(fo$formula), ll(fo)))
}
