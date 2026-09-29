# Two cs() terms that are the SAME column under different spellings: the
# parser cannot dedupe those, so the rank test has to.
#   Rscript dev/csfactor-rev-dup.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 170L)
one <- function(lab, f) {
  r <- tryCatch(suppressWarnings(suppressMessages(f())),
                error = function(e) structure(list(m = short(e)),
                                              class = "revfail"))
  if (inherits(r, "revfail")) {
    cat(sprintf("%-26s REFUSED  %s\n", lab, r$m)); return(invisible(NULL))
  }
  v <- tryCatch(suppressWarnings(vcov(r)), error = function(e) NULL)
  nse <- if (is.null(v)) NA_integer_ else sum(is.nan(sqrt(diag(v))))
  cat(sprintf("%-26s FITTED   df=%-3d logLik=%.9f  NaN-se=%s\n", lab,
              attr(logLik(r), "df"), as.numeric(logLik(r)), nse))
  invisible(r)
}
set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n))
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
d$fch <- as.character(d$f)
eta <- 0.5 * d$x
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.8 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))

one("cs(x) alone", function() frm(bf(yo ~ cs(x)), family = sratio(),
                                  data = d))
one("cs(x) + cs(I(x))", function() frm(bf(yo ~ cs(x) + cs(I(x))),
                                       family = sratio(), data = d))
one("cs(x) + cs(2*x)", function() frm(bf(yo ~ cs(x) + cs(I(2 * x))),
                                      family = sratio(), data = d))
one("cs(f) alone", function() frm(bf(yo ~ cs(f)), family = sratio(),
                                  data = d))
one("cs(f) + cs(fch)", function() frm(bf(yo ~ cs(f) + cs(fch)),
                                      family = sratio(), data = d))
cat("\nDONE ", TAG, "\n")
