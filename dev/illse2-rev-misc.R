# Review of lane illse2: odds and ends.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-misc.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
capm <- function(expr) {
  w <- character(); ms <- character(); wc <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "i2_err")),
    warning = function(x) { w <<- c(w, conditionMessage(x))
      wc <<- c(wc, class(x)[1]); invokeRestart("muffleWarning") },
    message = function(x) { ms <<- c(ms, paste(class(x)[1], conditionMessage(x)))
      invokeRestart("muffleMessage") })
  list(v = v, w = w, m = ms, wc = wc)
}
# 1. ML on sleepstudy with Days + shift, against lme4 (ML)
data("sleepstudy", package = "lme4")
set.seed(1)
ssu <- sleepstudy[sort(sample(nrow(sleepstudy), 120)), ]
for (sh in c(1e5, 1e6)) {
  d <- ssu
  d$t <- d$Days + sh
  r <- capm(frm(Reaction ~ t + (1 | Subject), data = d))
  l <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = d, REML = FALSE)
  cat(sprintf("ML shift %g: code %d | logLik - lme4(Days) %+.3e | warn %s\n",
              sh, r$v$opt$convergence,
              as.numeric(logLik(r$v)) - as.numeric(logLik(l)),
              if (length(r$w)) substr(paste(r$w, collapse = " || "), 1, 300) else "-"))
}
# 2. an active bound: what the held coefficient reports
set.seed(1)
deg <- 6
d <- data.frame(x = runif(200, 1, 2))
d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
xn <- paste0("x", 1:deg)
for (k in 1:deg) d[[xn[k]]] <- d$x^k
fo <- as.formula(paste("y ~", paste0(xn, collapse = " + ")))
l0 <- lm(fo, data = d)
ub <- coef(l0)[["x6"]] - sqrt(diag(vcov(l0)))[["x6"]]
r <- capm(frm(fo, data = d, prior = set_prior("", class = "b", coef = "x6",
                                                ub = ub)))
f <- r$v
print(suppressWarnings(fixef(f)))
cat("warnings:", r$wc, "\n")
for (x in r$w) cat("  W:", substr(x, 1, 300), "\n")
for (x in r$m) cat("  M:", substr(x, 1, 300), "\n")
cat("se_lost:", lost_of(f), "\n")
s <- capm(summary(f))
for (x in s$m) cat("  summary M:", substr(x, 1, 300), "\n")
for (x in s$w) cat("  summary W:", substr(x, 1, 300), "\n")
# 3. glmmTMB REML vcov, if installed: conditional or joint?
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  set.seed(2)
  G <- 30
  ng <- sample(2:15, G, replace = TRUE)
  g <- factor(rep(seq_len(G), ng))
  sim <- data.frame(g = g, x = rnorm(length(g)), z = rnorm(G)[g])
  sim$y <- 1 + 0.5 * sim$x + 0.8 * sim$z + rnorm(G, 0, 0.7)[g] +
    rnorm(nrow(sim), 0, 1)
  gt <- glmmTMB::glmmTMB(y ~ x + z + (1 | g), data = sim, REML = TRUE)
  cat("glmmTMB REML SE:", format(sqrt(diag(vcov(gt)$cond)), digits = 12), "\n")
} else cat("glmmTMB not installed\n")
