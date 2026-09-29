# Claim 4: the six brms_lp_check() rows again on a DIFFERENT seed and a
# different n, plus an ordered factor, so the recorded identity is not a
# property of seed 405. Same helper the suite uses.
#   Rscript dev/csfactor-rev-lp.R <lib> <cachedir>
args <- commandArgs(TRUE)
LIB <- args[[1L]]
CACHE <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-csfactor"
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = CACHE,
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
H <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(WT, "tests", "testthat", "helper-brms.R"), envir = H)
brms_lp_check <- H[["brms_lp_check"]]
options(frmtmb.brms_lp_report = TRUE)
cat("frmtmb", as.character(packageVersion("frmtmb")),
    "brms", as.character(packageVersion("brms")),
    "StanHeaders", as.character(packageVersion("StanHeaders")), "\n")

# a DIFFERENT construction from dev/csfactor-lp.R: seed 2029, n = 350,
# four levels rather than three, and a different slope
set.seed(2029)
n <- 350
x <- rnorm(n)
fc <- factor(sample(c("p", "q", "r", "s"), n, TRUE))
eff <- c(p = -0.8, q = 0.2, r = 1.1, s = -0.4)[as.character(fc)]
p1 <- plogis(-0.5 + 0.7 * x + eff)
p2 <- (1 - p1) * plogis(0.4 - eff)
u <- runif(n)
d <- data.frame(x = x, fc = fc,
                y = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
d$fch <- as.character(fc)
d$fnum <- factor(as.integer(fc))
d$fo <- factor(as.character(fc), levels = c("p", "q", "r", "s"),
               ordered = TRUE)
d$z <- rnorm(n)
cat("seed 2029  n =", n, " table(y) =", paste(table(d$y), collapse = "/"),
    " table(fc) =", paste(table(d$fc), collapse = "/"), "\n")

row <- function(lab, bform, fform, bfam, ffam) {
  cat("\n--", lab, "--\n")
  r <- tryCatch({
    fit <- frm(fform, family = ffam, data = d)
    cat("bcs names:", paste(grep("^bcs", variables(fit), value = TRUE),
                            collapse = " "), "\n")
    res <- brms_lp_check(bform, bfam, d, fit)
    sprintf("PASS const %.3e grad %.3e logLik %.10g",
            res$measured_const, res$max_grad, res$ours)
  }, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(r, "\n")
}

row("sratio, cs(fc) factor", brms::bf(y ~ x + cs(fc)), bf(y ~ x + cs(fc)),
    brms::sratio(), sratio())
row("cratio, cs(fc) factor", brms::bf(y ~ x + cs(fc)), bf(y ~ x + cs(fc)),
    brms::cratio(), cratio())
row("acat, cs(fc) factor", brms::bf(y ~ x + cs(fc)), bf(y ~ x + cs(fc)),
    brms::acat(), acat())
row("sratio, cs(fch) character", brms::bf(y ~ x + cs(fch)),
    bf(y ~ x + cs(fch)), brms::sratio(), sratio())
row("sratio, cs(fnum) numeric-looking levels", brms::bf(y ~ x + cs(fnum)),
    bf(y ~ x + cs(fnum)), brms::sratio(), sratio())
row("sratio, cs(fc) + cs(z) two terms", brms::bf(y ~ cs(fc) + cs(z)),
    bf(y ~ cs(fc) + cs(z)), brms::sratio(), sratio())
row("sratio, cs(fo) ORDERED factor", brms::bf(y ~ x + cs(fo)),
    bf(y ~ x + cs(fo)), brms::sratio(), sratio())
cat("\ndone\n")
