# Reviewer recheck round 1: the whole-term special refusal and the twin
# and slash handling, on designs of the reviewer's own, brms 2.23.0 /
# base / lane. Lane and base fit each accepted design and record logLik.
#   Rscript dev/priorform-rev2-special.R brms|ref|lane    seed 20260916
# Writes dev/priorform-rev2-special-<mode>.tsv
mode <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(switch(mode, ref = "C:/Users/adf44/source/r/rellib-r3",
                   lane = "C:/Users/adf44/source/r/priorform-lib", NULL),
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (mode == "brms") suppressMessages(library(brms)) else suppressMessages(library(frmtmb))
set.seed(20260916)
ng <- 12; n <- 240
lev <- paste0("i", 1:ng)
d <- data.frame(g = factor(rep(lev, 20), levels = lev), h = factor(rep(1:4, each = 60)),
                x = runif(n, 0, 3), z = runif(n, 0, 3), w = rnorm(n))
d$g1 <- factor(sample(lev, n, TRUE), levels = lev)
d$g2 <- factor(sample(lev, n, TRUE), levels = lev)
d$wt1 <- runif(n); d$wt2 <- runif(n)
d$mz <- sample(0:3, n, TRUE)
d$xm <- d$x; d$xm[sample(n, 20)] <- NA
d$sdx <- 0.1
d$off <- rnorm(n, 0, .1)
d$fg <- factor(rep(c("a", "b"), n / 2))
d$y <- sin(d$z) + 0.5 * d$x + rnorm(ng, 0, .5)[d$g] + rnorm(n, 0, .5)
d$ynl <- 2 * exp(-0.4 * d$x) + rnorm(ng, 0, .2)[d$g] + rnorm(n, 0, .2)
A <- diag(ng); for (i in seq(1, ng - 1, 2)) A[i, i + 1] <- A[i + 1, i] <- 0.5
dimnames(A) <- list(lev, lev)
d2 <- list(A = A)
designs <- c(
  s_by = "bf(y ~ s(z, by = x))",
  s_by_factor = "bf(y ~ fg + s(z, by = fg))",
  t2_xz = "bf(y ~ t2(x, z))",
  gp_by = "bf(y ~ gp(x, by = fg))",
  me_int = "bf(y ~ me(x, sdx) * w)",
  mo_int = "bf(y ~ mo(mz) * w)",
  mo_colon = "bf(y ~ w:mo(mz))",
  mi_int = "bf(y ~ mi(xm) * w) + bf(xm | mi() ~ w) + set_rescor(FALSE)",
  offset_plus_s = "bf(y ~ s(z) + offset(off))",
  s_in_nlpar = "bf(ynl ~ a * exp(-b * x), a ~ 1 + s(z), b ~ 1, nl = TRUE)",
  s_twice_star = "bf(y ~ x + s(z) * s(z))",
  s_twice_plus = "bf(y ~ x + s(z) + s(z))",
  s_paren = "bf(y ~ x + (s(z)))",
  s_paren_sum = "bf(y ~ (x + s(z)))",
  poly = "bf(y ~ poly(x, 2) + s(z))",
  s_in_sigma = "bf(y ~ x, sigma ~ s(z))",
  s_sigma_int = "bf(y ~ x, sigma ~ w * s(z))",
  s_fs = "bf(y ~ s(z, g, bs = 'fs', k = 4))",
  mmc_bar = "bf(y ~ (1 + mmc(wt1, wt2) | mm(g1, g2)))",
  s_power = "bf(y ~ (x + s(z))^2)",
  te_alone = "bf(y ~ te(x, z))",
  s_and_re = "bf(y ~ s(z) + (1 | g))",
  # twins and nesting
  tw_gh_hg = "bf(y ~ (1 | g:h) + (1 | h:g))",
  tw_gh_slash = "bf(y ~ (1 | g:h) + (1 | g/h))",
  tw_slash_slash = "bf(y ~ (1 | g/h) + (1 | g/h))",
  tw_slash_rev = "bf(y ~ (1 | g/h) + (1 | h/g))",
  tw_dblbar = "bf(y ~ (x || g) + (x || g))",
  tw_dblbar_bar = "bf(y ~ (x || g) + (x | g))",
  tw_single = "bf(y ~ x + (x | g))",
  tw_twin = "bf(y ~ x + (x | g) + (x | g))",
  tw_space = "bf(y ~ x + (1+x|g) + (1 + x | g))",
  tw_order = "bf(y ~ x + (1 + x | g) + (x + 1 | g))",
  tw_mm_args = "bf(y ~ (1 | mm(g1, g2)) + (1 | mm(g1, g2, weights = cbind(wt1, wt2))))",
  tw_gr_plain = "bf(y ~ (1 | gr(g)) + (1 | g))",
  tw_gr_args = "bf(y ~ (1 | gr(g, cov = A)) + (1 | gr(g, cov = A, dist = 'student')))",
  tw_ID = "bf(y ~ (1 | p | g) + (1 | p | g))",
  tw_ID_q = "bf(y ~ (1 | p | g) + (1 | q | g))",
  tw_nlpar = "bf(ynl ~ a * exp(-b * x), a ~ 1 + (1 | g) + (1 | g), b ~ 1, nl = TRUE)"
)
one <- function(nm) {
  txt <- designs[[nm]]
  if (mode == "ref" || mode == "lane") {
    txt <- sub(" + bf(xm | mi() ~ w) + set_rescor(FALSE)", "", txt, fixed = TRUE)
    if (nm == "mi_int") return(c(nm, "skip", "", ""))
  }
  r <- tryCatch({
    form <- eval(parse(text = txt))
    if (mode == "brms") {
      suppressWarnings(suppressMessages(stancode(form, data = d, data2 = d2)))
      c("accept", "", "")
    } else {
      fit <- suppressWarnings(suppressMessages(frm(form + gaussian(), data = d, data2 = d2)))
      c("accept", sprintf("%.10f", as.numeric(logLik(fit))),
        paste(colnames(fit$frame$linpreds[[1]]$X), collapse = ","))
    }
  }, error = function(e) {
    m <- gsub("[\t\r\n]+", " ", conditionMessage(e))
    v <- if (grepl("Duplicated group-level", m)) "refuse-dup" else
      if (grepl("is invalid", m)) "refuse-term" else "error-other"
    c(v, "", substr(m, 1, 100))
  })
  c(nm, r)
}
out <- do.call(rbind, lapply(names(designs), one))
colnames(out) <- c("id", "verdict", "logLik", "detail")
write.table(out, sprintf("C:/Users/adf44/source/r/frmtmb-wt-priorform/dev/priorform-rev2-special-%s.tsv", mode),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("done", mode, "\n")
