# Reviewer: the parser change's blast radius. Run once per library:
#   Rscript famsimplex-rev-parse.R lane|base
# Writes dev/famsimplex-rev-parse-<arm>.rds and prints a log.
arm <- commandArgs(TRUE)[1]
libs <- c(lane = "C:/Users/adf44/source/r/wt-famsimplex-lib",
          base = character(0))
.libPaths(c(if (arm == "lane") libs[["lane"]],
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " lib:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev"
run <- function(lab, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  out <- if (inherits(r, "frmtmb_fit")) {
    list(ll = as.numeric(logLik(r)), coef = fixef(r)[, "Estimate"],
         se = fixef(r)[, "Est.Error"], warn = w)
  } else list(err = r, warn = w)
  cat("\n---", lab, "\n")
  if (!is.null(out$err)) cat(out$err, "\n") else {
    cat(sprintf("logLik %.10f\n", out$ll))
    print(round(out$coef, 5))
  }
  if (length(w)) cat("WARN:", unique(w), sep = "\n  ")
  out
}
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:15, length.out = n)))
E <- cbind(0, 0.3 + 0.5 * d$z, -0.2 + 0.6 * d$x)
P <- exp(E) / rowSums(exp(E))
d$f <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample(3, 1, prob = p))])
d$N <- 10L
d$Ym <- t(vapply(seq_len(n), function(i) rmultinom(1, 10, P[i, ])[, 1],
                 numeric(3)))
colnames(d$Ym) <- c("a", "b", "c")
cl <- rbinom(n, 1, 0.4)
d$y <- ifelse(cl == 1, 2 + 1 * d$x, -1 + 0.5 * d$z) + rnorm(n, 0, 0.7)
d$y2 <- 1 + 0.5 * d$x + rnorm(n)
res <- list()
res$cat_main <- run("categorical f ~ x",
  frm(bf(f ~ x), family = categorical(), data = d))
res$cat_muc <- run("categorical f ~ x, muc ~ z",
  frm(bf(f ~ x, muc ~ z), family = categorical(), data = d))
res$cat_mub <- run("categorical f ~ x, mub ~ z (new)",
  frm(bf(f ~ x, mub ~ z), family = categorical(), data = d))
res$cat_both <- run("categorical f ~ w, mub ~ z, muc ~ x (main unused)",
  frm(bf(f ~ w, mub ~ z, muc ~ x), family = categorical(), data = d))
res$cat_both_ref <- run("categorical f ~ 1, mub ~ z, muc ~ x",
  frm(bf(f ~ 1, mub ~ z, muc ~ x), family = categorical(), data = d))
res$cat_re <- run("categorical f ~ x + (1 | ID | g)",
  frm(bf(f ~ x + (1 | ID | g)), family = categorical(), data = d))
res$mn_main <- run("multinomial Ym | trials(N) ~ x",
  frm(bf(Ym | trials(N) ~ x), family = multinomial(), data = d))
res$mn_muc <- run("multinomial muc ~ z",
  frm(bf(Ym | trials(N) ~ x, muc ~ z), family = multinomial(), data = d))
res$mn_mub <- run("multinomial mub ~ z (new)",
  frm(bf(Ym | trials(N) ~ x, mub ~ z), family = multinomial(), data = d))
mix <- mixture(gaussian(), gaussian())
res$mix_main <- run("mixture y ~ x",
  frm(bf(y ~ x), family = mix, data = d))
res$mix_mu2 <- run("mixture y ~ x, mu2 ~ z",
  frm(bf(y ~ x, mu2 ~ z), family = mix, data = d))
res$mix_mu1 <- run("mixture y ~ x, mu1 ~ z (new)",
  frm(bf(y ~ x, mu1 ~ z), family = mix, data = d))
res$mix_both <- run("mixture y ~ 1, mu1 ~ x, mu2 ~ z (new)",
  frm(bf(y ~ 1, mu1 ~ x, mu2 ~ z), family = mix, data = d))
res$mix_both_w <- run("mixture y ~ w, mu1 ~ x, mu2 ~ z (main unused)",
  frm(bf(y ~ w, mu1 ~ x, mu2 ~ z), family = mix, data = d))
res$mix_fix1 <- run("mixture y ~ x, mu1 = 0 (new)",
  frm(bf(y ~ x, mu1 = 0), family = mix, data = d))
res$mix_ordered <- run("mixture(order = TRUE) y ~ x, mu1 ~ z",
  frm(bf(y ~ x, mu1 ~ z), family = mixture(gaussian(), gaussian(),
                                           order = "mu"), data = d))
res$mv <- run("mvbf(y2 ~ x, f ~ x categorical)",
  frm(mvbf(bf(y2 ~ x), bf(f ~ x, family = categorical())), data = d))
res$mv_mub <- run("mvbf(y2 ~ x, f ~ x, mub ~ z categorical) (new)",
  frm(mvbf(bf(y2 ~ x), bf(f ~ x, mub ~ z, family = categorical())),
      data = d))
res$gauss_mu <- run("gaussian y ~ x, mu ~ z (still refused?)",
  frm(bf(y ~ x, mu ~ z), family = gaussian(), data = d))
res$cat_mu <- run("categorical f ~ x, mu ~ z",
  frm(bf(f ~ x, mu ~ z), family = categorical(), data = d))
# an NA only in the main formula's variable that every location overrides
d2 <- d
d2$w[1:10] <- NA
res$cat_na_unused <- run("categorical f ~ w, mub ~ z, muc ~ x; w NA rows 1:10",
  {
    f <- frm(bf(f ~ w, mub ~ z, muc ~ x), family = categorical(), data = d2)
    cat("nobs:", nobs(f), "\n")
    f
  })
saveRDS(res, file.path(wt, paste0("famsimplex-rev-parse-", arm, ".rds")))
