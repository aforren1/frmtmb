# Claim 2, follow-up: the thres(gr = ) shape, and the hole the rank test
# leaves where X carries a zero PLACEHOLDER column (mo(), me(), mi()).
#   Rscript dev/csfactor-rev-rank2.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 300L)
fitq <- function(...) {
  tryCatch(suppressWarnings(suppressMessages(frm(...))),
           error = function(e) structure(list(m = short(e)),
                                         class = "revfail"))
}
rep1 <- function(tag, r) {
  if (inherits(r, "revfail")) {
    cat(sprintf("%-28s REFUSED  %s\n", tag, r$m)); return(invisible(NULL))
  }
  v <- tryCatch(suppressWarnings(vcov(r)), error = function(e) NULL)
  nan_se <- if (is.null(v)) NA_integer_ else sum(is.nan(sqrt(diag(v))))
  cat(sprintf("%-28s FITTED npar=%-3d logLik=%.7f NaN-se=%s\n", tag,
              length(r$opt$par), as.numeric(logLik(r)), nan_se))
  invisible(r)
}

set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), sdx = 0.2)
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
d$g <- factor(sample(c("p", "q"), n, TRUE))
d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
eta <- 0.5 * d$x + 0.3 * d$z
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.8 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
d$yof <- factor(d$yo, ordered = TRUE)

cat("\n### thres(gr = g) beside cs()\n")
rep1("thres(gr) + x + cs(f)",
     fitq(bf(yo | thres(gr = g) ~ x + cs(f)), family = sratio(), data = d))
rep1("thres(gr) + cs(x)",
     fitq(bf(yo | thres(gr = g) ~ cs(x)), family = sratio(), data = d))
rep1("thres(gr) + x + cs(x)",
     fitq(bf(yo | thres(gr = g) ~ x + cs(x)), family = sratio(), data = d))
rep1("thres(gr) + g + cs(g)",
     fitq(bf(yo | thres(gr = g) ~ g + cs(g)), family = sratio(), data = d))
rep1("thres(gr) + cs(g)",
     fitq(bf(yo | thres(gr = g) ~ cs(g)), family = sratio(), data = d))

cat("\n### mo(): a zero placeholder column in X, so the rank test is blind\n")
a <- rep1("cs(m) alone", fitq(bf(yo ~ cs(m)), family = sratio(), data = d))
b <- rep1("mo(m) + cs(m)", fitq(bf(yo ~ mo(m) + cs(m)), family = sratio(),
                                data = d))
cc <- rep1("mo(m) alone", fitq(bf(yo ~ mo(m)), family = sratio(), data = d))
if (!inherits(a, "revfail") && !inherits(b, "revfail")) {
  la <- as.numeric(logLik(a)); lb <- as.numeric(logLik(b))
  cat(sprintf("  logLik(cs(m)) - logLik(mo(m)+cs(m)) = %.3e",
              la - lb), "\n")
  cat("  df: ", attr(logLik(a), "df"), " vs ", attr(logLik(b), "df"), "\n")
  cat("  mo(m)+cs(m) coefficient table:\n")
  print(signif(fixef(b), 4))
}
cat("\n### the same shape from two different starts: a flat direction\n")
b1 <- fitq(bf(yo ~ mo(m) + cs(m)), family = sratio(), data = d)
if (!inherits(b1, "revfail")) {
  st <- b1$estimates
  st$b <- st$b + 0.75
  b2 <- fitq(bf(yo ~ mo(m) + cs(m)), family = sratio(), data = d,
             start = st)
  if (!inherits(b2, "revfail")) {
    cat(sprintf("  logLik A = %.9f  logLik B = %.9f  diff = %.3e\n",
                as.numeric(logLik(b1)), as.numeric(logLik(b2)),
                as.numeric(logLik(b1)) - as.numeric(logLik(b2))))
    ca <- unlist(fixef(b1)[, 1L]); cb <- unlist(fixef(b2)[, 1L])
    cat("  max abs coefficient difference between the two = ",
        sprintf("%.4g", max(abs(ca - cb))), "\n")
    cat("  mom: ", sprintf("%.6f", fixef(b1)["mom", 1L]), " vs ",
        sprintf("%.6f", fixef(b2)["mom", 1L]), "\n")
  }
}

cat("\n### me(): also a placeholder column\n")
rep1("me(x,sdx) + cs(x)",
     fitq(bf(yo ~ me(x, sdx) + cs(x)), family = sratio(), data = d))
rep1("me(x,sdx) + cs(z)",
     fitq(bf(yo ~ me(x, sdx) + cs(z)), family = sratio(), data = d))

cat("\n### mi(): a missing covariate imputed, beside cs() on it\n")
d2 <- d
d2$xm <- d2$x
d2$xm[sample.int(n, 30L)] <- NA
r <- fitq(bf(yo ~ mi(xm) + cs(z)) + bf(xm ~ z),
          family = list(sratio(), gaussian()), data = d2)
rep1("mi(xm) + cs(z)", r)
cat("\nDONE ", TAG, "\n")
