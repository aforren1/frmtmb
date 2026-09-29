# Punch round 2 (nits): the spelling the mo() refusal should suggest, and
# the pre-existing mo(m) + m defect, measured before either is written
# into a message or the backlog. Seed 1907, the reviewer's.
#   Rscript dev/csfactor-p2.R <lib> > dev/csfactor-log/p2.txt
args <- commandArgs(TRUE)
LIB <- if (length(args)) args[[1L]] else
  "C:/Users/adf44/source/r/wt-csfactor-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(system.file(package = "frmtmb")), "\n")
set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
d$yo <- sample(1:3, n, TRUE)

rep1 <- function(lab, f) {
  cat("\n--", lab, "--\n")
  r <- tryCatch({
    fit <- frm(f, family = sratio(), data = d)
    fe <- fixef(fit)
    ev <- tryCatch(eigen(vcov(fit), only.values = TRUE)$values,
                   error = function(e) NA_real_)
    sprintf(paste("FITTED df %d logLik %.9f | non-finite SEs %d of %d |",
                  "min vcov eigenvalue %.4g"),
            attr(logLik(fit), "df") %||% nrow(fe),
            as.numeric(logLik(fit)),
            sum(!is.finite(fe[, "Est.Error"])), nrow(fe), min(ev))
  }, error = function(e) paste("REFUSED:", conditionMessage(e)))
  cat(r, "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("\n== (2) the spelling that keeps the interaction ==\n")
rep1("cs(m) alone", bf(yo ~ cs(m)))
rep1("mo(m) + cs(m): refused", bf(yo ~ mo(m) + cs(m)))
rep1("mo(m) * z + cs(m)", bf(yo ~ mo(m) * z + cs(m)))
rep1("z + mo(m):z + cs(m)", bf(yo ~ z + mo(m):z + cs(m)))
rep1("mo(m):z + cs(m)", bf(yo ~ mo(m):z + cs(m)))

cat("\n== (4) the pre-existing mo(m) + m defect, no cs() anywhere ==\n")
d$mn <- factor(as.integer(d$m))       # the same variable as a plain factor
rep1("mo(m) + m", bf(yo ~ mo(m) + m))
rep1("mo(m) + mn (m as an unordered factor)", bf(yo ~ mo(m) + mn))
rep1("mo(m) alone", bf(yo ~ mo(m)))
rep1("m alone", bf(yo ~ m))
cat("\ndone\n")
