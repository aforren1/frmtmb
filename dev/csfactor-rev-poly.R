# Is the poly(x, 2) + cs(x) refusal right? If the x direction really is
# the flat one, then dropping it and keeping only the quadratic must
# reach the same maximum. Measured against the reference build's fit of
# the refused model.
#   Rscript dev/csfactor-rev-poly.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 160L)
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
d <- data.frame(x = rnorm(n), z = rnorm(n))
eta <- 0.5 * d$x
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.8 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))

one("cs(x) alone", function() frm(bf(yo ~ cs(x)), family = sratio(),
                                  data = d))
one("poly(x,2) + cs(x)", function() frm(bf(yo ~ poly(x, 2) + cs(x)),
                                        family = sratio(), data = d))
one("I(x^2) + cs(x)", function() frm(bf(yo ~ I(x^2) + cs(x)),
                                     family = sratio(), data = d))
one("poly(x,2) alone", function() frm(bf(yo ~ poly(x, 2)),
                                      family = sratio(), data = d))
one("s(x) + cs(x)", function() frm(bf(yo ~ s(x) + cs(x)),
                                   family = sratio(), data = d))
one("s(x) alone", function() frm(bf(yo ~ s(x)), family = sratio(),
                                 data = d))
cat("\nDONE ", TAG, "\n")
