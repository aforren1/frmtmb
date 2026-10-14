# Lane optstart: the wide design (dev/optrap-p2-wide1.R, n 101000, 50
# columns), fitted from starts given through `start =`, to see which
# part of the start sets the first run's iterations and its code.
# Arms (mu at least squares in all):
#   rms: sigma at the rms of the residuals, the group sd at the unit
#        (the lane's default start, given explicitly)
#   within: sigma at the rms of the residuals less their group means,
#        the group sd at the unit
#   anova: sigma as within, the group sd at the sd of the residuals'
#        group means less the within part
# Per fit: the first run alone (restarts = 0) and then the default.
#   Rscript dev/optstart-wide-arms.R <lib or ""> <seed> <arms>
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
s <- as.integer(a[2])
arms <- strsplit(a[3], ",")[[1]]
cat("lib", dirname(find.package("frmtmb")), "seed", s, "\n")
set.seed(s)
n <- 101000
p <- 49
Z <- matrix(rnorm(n * p), n, p)
colnames(Z) <- paste0("z", seq_len(p))
d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                       "+ (1 | g)"))
l4 <- as.numeric(logLik(lmer(fo, data = d, REML = FALSE)))
X <- model.matrix(lme4::nobars(fo), d)
q <- qr(X)
r <- qr.resid(q, d$y)
gm <- ave(r, d$g)
rw <- r - gm
sw <- sqrt(mean(rw^2))
ng <- mean(table(d$g))
sb <- sqrt(max(var(tapply(r, d$g, mean)) - sw^2 / ng, sw^2 / ng))
st <- list(
  rms = list(beta = qr.coef(q, d$y), betad = log(sqrt(mean(r^2))), theta = 0),
  within = list(beta = qr.coef(q, d$y), betad = log(sw), theta = 0),
  anova = list(beta = qr.coef(q, d$y), betad = log(sw), theta = log(sb)))
cat(sprintf("sigma: rms %.4f within %.4f | group sd anova %.4f\n",
            sqrt(mean(r^2)), sw, sb))
for (nm in arms) {
  for (rs in c(0L, 1L)) {
    w <- character()
    t0 <- proc.time()[[3]]
    f <- withCallingHandlers(suppressMessages(frm(fo, data = d,
      start = st[[nm]], control = frmtmb_control(restarts = rs))),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    cat(sprintf(paste0("%-7s restarts %d: code %d | logLik - lmer %.3g | ",
                       "iter %d | warnings %d | secs %.0f\n"),
                nm, rs, f$opt$convergence, as.numeric(logLik(f)) - l4,
                f$opt$iterations, length(w), proc.time()[[3]] - t0))
  }
}
