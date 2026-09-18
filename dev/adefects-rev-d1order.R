source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260917)
d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 6),
                g2 = rep(c(1, 2, 1, 2), each = 6),
                t = c(1:6, 1:6, 1:6, 7:12))
d$x <- rnorm(24)
d$g <- interaction(d$g1, d$g2)
d$y <- d$x + as.numeric(stats::filter(rnorm(24), 0.6, "recursive"))
d$gf1 <- factor(d$g1); d$gf2 <- factor(d$g2)
d$g3 <- rep(c(1, 2), 12)
d$a <- c(rep("x_1", 12), rep("x", 12))
d$b <- c(rep("2", 12), rep("1_2", 12))

v <- function(txt) {
  f <- stats::as.formula(paste("y ~ x +", txt))
  r <- tryCatch(sprintf("OK %.5f", as.numeric(stats::logLik(
        suppressWarnings(suppressMessages(frm(f, d)))))),
        error = function(e) paste0("REFUSE: ",
          substr(gsub("[\r\n]+", " ", conditionMessage(e)), 1, 105)))
  cat(sprintf("%-42s %s\n", txt, r))
}

cat("== the claims the lane records for factor(g) and interaction()\n")
v("ar(t, factor(g), cov = TRUE)")
v("ar(t, g, cov = TRUE)")
v("ar(t, interaction(g1, g2), cov = TRUE)")

cat("\n== refusal ORDER: grammar before the cov = TRUE refusal\n")
v("ar(x + t, g)")
v("ar(t, g1/g2)")
v("ar(t, factor(g))")
v("ar(x + t, g1/g2)")
v("ar(t, g)")

cat("\n== crossing collision: paste(sep = '_') is not injective\n")
v("ar(t, a:b, cov = TRUE)")
cat("  a:b as frmtmb pastes it: ",
    paste(unique(paste(d$a, d$b, sep = "_")), collapse = " | "), "\n")
cat("  R's factor interaction:  ",
    paste(levels(droplevels(interaction(d$a, d$b))), collapse = " | "), "\n")

cat("\n== a crossed component with NA\n")
d2 <- d; d2$g2[3] <- NA
r <- tryCatch(frm(y ~ x + ar(t, g1:g2, cov = TRUE), d2),
              error = function(e) conditionMessage(e))
cat("  NA in one component: ", substr(gsub("[\r\n]+", " ", r), 1, 130), "\n")
