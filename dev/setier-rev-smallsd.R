# Reviewer of lane setier: small but identified variance components and
# strong but identified correlations must not get the boundary message.
# Against lme4's isSingular() on the same data (ML for lmer).
#   Rscript dev/setier-rev-smallsd.R <lib> <out.tsv> [seeds]
# Designs (seed 1..S):
#   s1  y ~ x + (1 | g), 100 groups of 10, group sd 0.15, sigma 1
#   s2  y ~ x + (1 | g), 300 groups of 20, group sd 0.08, sigma 1
#   c95 y ~ x + (1 + x | g), 50 groups of 10, sds 1 and 0.5, cor 0.95
#   c99 the same with cor 0.99
#   p1  poisson y ~ x + (1 | g), 60 groups of 10, group sd 0.1
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
out <- args[2]
S <- if (length(args) > 2) as.integer(args[3]) else 40L
cat("lib", find.package("frmtmb"), "\n")
gen <- function(des, s) {
  set.seed(s)
  if (des %in% c("s1", "s2")) {
    G <- if (des == "s1") 100 else 300
    n <- if (des == "s1") 10 else 20
    sdg <- if (des == "s1") 0.15 else 0.08
    d <- data.frame(g = factor(rep(1:G, each = n)), x = rnorm(G * n))
    d$y <- 1 + 0.5 * d$x + rnorm(G, 0, sdg)[d$g] + rnorm(G * n)
    list(d = d, ff = y ~ x + (1 | g), fam = gaussian())
  } else if (des %in% c("c95", "c99")) {
    r <- if (des == "c95") 0.95 else 0.99
    G <- 50
    d <- data.frame(g = factor(rep(1:G, each = 10)), x = rnorm(G * 10))
    S2 <- matrix(c(1, r * 0.5, r * 0.5, 0.25), 2)
    u <- matrix(rnorm(2 * G), G) %*% chol(S2)
    d$y <- 1 + 0.5 * d$x + u[d$g, 1] + u[d$g, 2] * d$x + rnorm(G * 10)
    list(d = d, ff = y ~ x + (1 + x | g), fam = gaussian())
  } else {
    d <- data.frame(g = factor(rep(1:60, each = 10)), x = rnorm(600))
    d$y <- rpois(600, exp(0.5 + 0.3 * d$x + rnorm(60, 0, 0.1)[d$g]))
    list(d = d, ff = y ~ x + (1 | g), fam = poisson())
  }
}
rows <- list()
for (des in c("s1", "s2", "c95", "c99", "p1")) for (s in seq_len(S)) {
  g <- gen(des, s)
  lm4 <- suppressMessages(suppressWarnings(
    if (g$fam$family == "gaussian") lmer(g$ff, data = g$d, REML = FALSE)
    else glmer(g$ff, data = g$d, family = g$fam)))
  m <- character(); w <- character()
  f <- withCallingHandlers(frm(g$ff, family = g$fam, data = g$d),
    warning = function(x) {w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")},
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  vc <- as.data.frame(VarCorr(lm4))
  sdl <- vc$sdcor[vc$grp == "g" & is.na(vc$var2)]
  corl <- vc$sdcor[vc$grp == "g" & !is.na(vc$var2)]
  lost <- frmtmb:::sdr_of(f)$se_lost
  rows[[length(rows) + 1L]] <- data.frame(
    des = des, seed = s, lme4_singular = isSingular(lm4),
    lme4_minsd = min(sdl), lme4_cor = if (length(corl)) corl[1] else NA,
    code = f$opt$convergence, theta = paste(signif(f$estimates$theta, 4),
                                            collapse = ","),
    boundary = sum(grepl("^Boundary", m)),
    se_warn = sum(grepl("Standard errors are not available", w)),
    other_warn = sum(!grepl("Standard errors are not available", w)),
    lost = paste(names(lost), lost, sep = ":", collapse = ";"))
}
X <- do.call(rbind, rows)
write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
for (des in unique(X$des)) {
  x <- X[X$des == des, ]
  cat(sprintf(paste0("%-4s n %d | lme4 singular %d | boundary msg %d (on ",
                     "lme4-singular %d, on others %d) | SE warn %d | code!=0 ",
                     "%d | min lme4 sd among non-singular %.3g\n"),
              des, nrow(x), sum(x$lme4_singular), sum(x$boundary > 0),
              sum(x$boundary > 0 & x$lme4_singular),
              sum(x$boundary > 0 & !x$lme4_singular), sum(x$se_warn > 0),
              sum(x$code != 0),
              suppressWarnings(min(x$lme4_minsd[!x$lme4_singular]))))
}
