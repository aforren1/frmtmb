# Claim 2: hunt false alarms in check_cs_identified(). Each shape is
# fitted and the outcome recorded; a refusal prints its message so the
# verdict "right to refuse" can be judged from the text.
#   Rscript dev/csfactor-rev-rank.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 300L)

set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n))
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
d$g <- factor(rep(1:30, each = 10))
d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$gr <- factor(sample(c("p", "q"), n, TRUE))
eta <- 0.5 * d$x + 0.3 * d$z
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.8 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))
d$yo2 <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))

run <- function(tag, fexpr, ...) {
  r <- tryCatch(suppressWarnings(suppressMessages(
         frm(fexpr, family = sratio(), data = d, ...))),
        error = function(e) structure(list(m = short(e)), class = "revfail"))
  if (inherits(r, "revfail")) {
    cat(sprintf("%-38s REFUSED  %s\n", tag, r$m))
  } else {
    cat(sprintf("%-38s FITTED   npar=%d  logLik=%.5f  [%s]\n", tag,
                length(r$opt$par), as.numeric(logLik(r)),
                paste(rownames(fixef(r)), collapse = " ")))
  }
  invisible(r)
}

cat("\n### grammar and interaction shapes\n")
run("cs(x:z)",              bf(yo ~ x + z + cs(x:z)))
run("cs(x*z)",              bf(yo ~ cs(x * z)))
run("cs(f:z)",              bf(yo ~ z + cs(f:z)))
run("f:z main + cs(f)",     bf(yo ~ f:z + cs(f)))
run("f*z main + cs(f)",     bf(yo ~ f * z + cs(f)))
run("cs(x)+cs(x:z)",        bf(yo ~ cs(x) + cs(x:z)))
run("cs(x)+cs(z) share none", bf(yo ~ cs(x) + cs(z)))
run("cs(x)+cs(x)  twice",   bf(yo ~ cs(x) + cs(x)))
run("cs(f)+cs(f:z)",        bf(yo ~ cs(f) + cs(f:z)))
run("cs(f)+cs(f)  twice",   bf(yo ~ cs(f) + cs(f)))

cat("\n### cs() beside a basis or a monotonic term on the SAME variable\n")
run("poly(x,2) + cs(x)",    bf(yo ~ poly(x, 2) + cs(x)))
run("poly(x,2) + cs(z)",    bf(yo ~ poly(x, 2) + cs(z)))
run("s(x) + cs(x)",         bf(yo ~ s(x) + cs(x)))
run("s(x) + cs(z)",         bf(yo ~ s(x) + cs(z)))
run("mo(m) + cs(x)",        bf(yo ~ mo(m) + cs(x)))
run("mo(m) + cs(m)",        bf(yo ~ mo(m) + cs(m)))
run("x + cs(log(exp(x)))",  bf(yo ~ x + cs(log(exp(x)))))
run("x + cs(scale(x))",     bf(yo ~ x + cs(scale(x))))
run("x + cs(2*x+1)",        bf(yo ~ x + cs(I(2 * x + 1))))

cat("\n### cs(f) with f in a random slope, and a random intercept\n")
run("cs(f) + (f|g)",        bf(yo ~ cs(f) + (f | g)))
run("f + cs(z) + (f|g)",    bf(yo ~ f + cs(z) + (f | g)))
run("cs(x) + (1+x||g)",     bf(yo ~ cs(x) + (1 + x || g)))

cat("\n### thres(gr = ) with cs()\n")
r <- tryCatch(suppressWarnings(suppressMessages(
       frm(bf(yo ~ x + cs(f)) + thres(gr = gr), family = sratio(),
           data = d))), error = function(e) short(e))
cat(if (is.character(r)) paste("thres(gr) + cs(f) REFUSED ", r) else
    sprintf("thres(gr) + cs(f) FITTED npar=%d logLik=%.5f [%s]",
            length(r$opt$par), as.numeric(logLik(r)),
            paste(rownames(fixef(r)), collapse = " ")), "\n")

cat("\n### multivariate ordinal with cs()\n")
r <- tryCatch(suppressWarnings(suppressMessages(
       frm(bf(yo ~ x + cs(f)) + bf(yo2 ~ z + cs(f)),
           family = sratio(), data = d))),
     error = function(e) short(e))
cat(if (is.character(r)) paste("mv cs REFUSED ", r) else
    sprintf("mv cs FITTED npar=%d logLik=%.5f", length(r$opt$par),
            as.numeric(logLik(r))), "\n")
r <- tryCatch(suppressWarnings(suppressMessages(
       frm(bf(yo ~ x + cs(x)) + bf(yo2 ~ z + cs(f)),
           family = sratio(), data = d))),
     error = function(e) short(e))
cat(if (is.character(r)) paste("mv one bad lp REFUSED ", r) else
    "mv one bad lp FITTED (should have been refused)", "\n")

cat("\n### near collinearity: at what scale does the rank test fire?\n")
# main effect x, cs() on a column that is x plus noise of size eps.
# The design is legitimate at every eps > 0; the question is the qr()
# tolerance, which is 1e-7 by default on unit-norm columns.
for (eps in 10^-(1:11)) {
  d$xn <- d$x + eps * rnorm(n)
  cnum <- kappa(scale(cbind(1, d$x, d$xn), center = FALSE,
                      scale = sqrt(colSums(cbind(1, d$x, d$xn)^2))),
                exact = TRUE)
  r <- tryCatch(suppressWarnings(suppressMessages(
         frm(bf(yo ~ x + cs(xn)), family = sratio(), data = d))),
        error = function(e) structure(list(m = short(e)), class = "revfail"))
  cat(sprintf("  eps=%-8.0e kappa=%-10.3g %s\n", eps, cnum,
              if (inherits(r, "revfail")) "REFUSED" else
                sprintf("FITTED logLik=%.6f", as.numeric(logLik(r)))))
}
cat("\n### the refusals that SHOULD happen\n")
run("x + cs(x)",            bf(yo ~ x + cs(x)))
run("I(x) + cs(x)",         bf(yo ~ I(x) + cs(x)))
run("f + cs(f)",            bf(yo ~ f + cs(f)))
run("x1+x2 + cs(x1+x2)",    bf(yo ~ x + z + cs(I(x + z))))
cat("\nDONE ", TAG, "\n")
