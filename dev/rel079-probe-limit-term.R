# fitted(newdata) above 3000 rows beside a gp(k =) or a s() term, one fit each, against
# references that hold no frmtmb reader of the latent's joint: the mean's closed form
# exp(eta + s^2 / 2), and the two parts of Est.Error in quadrature, the no-latent fit's
# Est.Error at the row and the latent's draw. Inside and past the range of x.
#   Rscript dev/rel079-probe-limit-term.R [rows]
n <- as.integer(commandArgs(TRUE)[1])
sde <- as.numeric(commandArgs(TRUE)[2])
if (is.na(sde)) sde <- 0.15
if (is.na(n)) n <- 3500L
.libPaths(c("C:/Users/adf44/source/r/rellib-r17", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")),
    format(file.mtime(system.file("DESCRIPTION", package = "frmtmb"))), "\n")
cat("limit", getOption("frmtmb.latres_se_rows", 3000L), " rows", n, "\n")
q <- function(x) signif(x, 4)
for (seed in c(7551, 7552)) {
  set.seed(seed)
  nt <- 50L
  ng <- ceiling(n / nt)
  d <- expand.grid(t = seq_len(nt), s = factor(seq_len(ng)))[seq_len(n), ]
  d$x <- stats::runif(n, -2, 2)
  e <- stats::ave(stats::rnorm(n), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
  d$y <- stats::rpois(n, exp(1 + sin(d$x) + sde * e))
  nd <- data.frame(t = nt + 1L, s = factor(c(1, 2), levels = levels(d$s)), x = c(0, 1.5))
  for (term in c("gp(x, k = 10)", "s(x)")) {
    f <- suppressWarnings(suppressMessages(
      frm(stats::as.formula(paste0("y ~ ", term, " + ar(time = t, gr = s, cov = TRUE)")), data = d, family = poisson())))
    f0 <- suppressWarnings(suppressMessages(
      frm(stats::as.formula(paste0("y ~ ", term)), data = d, family = poisson())))
    said <- character()
    t0 <- proc.time()[[3]]
    a <- withCallingHandlers(suppressWarnings(fitted(f, newdata = nd)), message = function(m) {
      said <<- c(said, class(m)[1L])
      invokeRestart("muffleMessage")
    })
    tt <- proc.time()[[3]] - t0
    b <- suppressWarnings(suppressMessages(fitted(f0, newdata = nd)))
    sl <- ns$latres_row_sd(f, "y", NULL)
    eta <- as.numeric(suppressWarnings(frm_linpred(f, newdata = nd, type = "link", dpar = "mu")))
    m <- exp(eta + sl^2 / 2)
    two <- sqrt(b[, "Est.Error"]^2 + (a[, "Estimate"] * sqrt(exp(sl^2) - 1))^2)
    cat(sprintf("seed %d %-14s blocks %d sl %s conv %d in %.1f s, messages: %s\n", seed, term,
                length(f$frame$re_blocks), q(sl), f$opt$convergence, tt, paste(unique(said), collapse = ",")))
    cat("  estimate", q(a[, "Estimate"]), " closed form", q(m), " ratio", q(a[, "Estimate"] / m), "\n")
    cat("  Est.Error", q(a[, "Est.Error"]), " two parts", q(two), " ratio", q(a[, "Est.Error"] / two),
        " no-latent", q(b[, "Est.Error"]), "\n")
    cat("  95%", q(a[, "Q2.5"]), "to", q(a[, "Q97.5"]), "\n")
  }
}
cat("PROBE END\n")
