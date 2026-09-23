# Punch round 1: reproduce B1, B2, M2 and M3 before touching the code,
# then rerun the same script after the fix. LIB picks the build.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
op <- options(digits = 12)
esc <- function(f) {
  e <- f$opt[["stationary_escape"]]
  if (is.null(e)) "-" else paste0(e[["starts"]], "/",
                                  format(e[["gain"]], digits = 5))
}
line <- function(lab, f, ...) {
  a <- f$estimates$betad
  a <- a[grep("^alpha", names(a))]
  cat(sprintf("%-26s ll %16.7f  conv %d  %-24s max|a| %9.3g  esc %s\n",
              lab, as.numeric(logLik(f)), f$opt$convergence,
              substr(f$opt$message, 1, 24), max(abs(a)), esc(f)))
}

cat("\n===== B1: a dpar design with no intercept\n")
set.seed(7)
n <- 300
grp <- factor(rep(1:3, length.out = n))
x <- rnorm(n)
y <- 1 + 0.5 * x + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5
dd <- data.frame(y = y, x = x, grp = grp)
for (fo in list(bf(y ~ x, sigma ~ 1, alpha ~ grp),
                bf(y ~ x, sigma ~ 1, alpha ~ 0 + grp),
                bf(y ~ x, sigma ~ 1, alpha ~ x),
                bf(y ~ x, sigma ~ 1, alpha ~ 0 + x))) {
  lab <- deparse(fo$forms[[3]]$formula %||% fo)[1]
  f <- try(suppressWarnings(frm(fo, family = skew_normal(), data = dd)),
           silent = TRUE)
  if (inherits(f, "try-error")) {
    cat(sprintf("%-26s ERROR\n", lab))
  } else line(lab, f)
}

cat("\n===== B2: an ATOMIC post$stationary\n")
mk <- function(stat) {
  custom_family(
    "b2", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    lpdf = function(y, dpars, aterms) {
      -0.5 * ((y - dpars[["mu"]]) / dpars[["sigma"]])^2 -
        log(dpars[["sigma"]]) - 0.5 * log(2 * pi)
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) stats::sd(y)),
    post = c(list(mean_fn = function(dpars, aterms) dpars[["mu"]]),
             if (!is.null(stat)) list(stationary = stat)))
}
set.seed(9)
d2 <- data.frame(y = rnorm(120), x = rnorm(120))
shapes <- list(none = NULL, character = "alpha", numeric = c(0, 1),
               named_num = c(alpha = 0), logical = TRUE, na = NA)
for (nm in names(shapes)) {
  r <- try(frm(bf(y ~ x, sigma ~ 1), family = mk(shapes[[nm]]), data = d2),
           silent = TRUE)
  cat(sprintf("  %-10s %s\n", nm,
              if (inherits(r, "try-error"))
                paste("ERROR:", sub("\n.*", "",
                                    attr(r, "condition")$message))
              else paste("ok, ll", format(as.numeric(logLik(r)),
                                          digits = 10))))
}

cat("\n===== M2: offset() in the mu predictor\n")
for (s in 1:6) {
  dd2 <- make_data(s, FALSE)
  names(dd2)[names(dd2) == "xs"] <- "o"
  f <- suppressWarnings(frm(bf(y ~ offset(o), sigma ~ 1, alpha ~ 1),
                            family = skew_normal(), data = dd2))
  ref <- max(vapply(c(2, -2), function(a) as.numeric(logLik(
    suppressWarnings(frm(bf(y ~ offset(o), sigma ~ 1, alpha ~ 1),
                         family = skew_normal(), data = dd2,
                         start = list(betad = c(0, a)))))), 0))
  cat(sprintf("  offset seed %d: ll %14.6f  best %14.6f  short %10.5f  esc %s\n",
              s, as.numeric(logLik(f)), ref, ref - as.numeric(logLik(f)),
              esc(f)))
}

cat("\n===== M3: the spurious convergence warning\n")
set.seed(8)
n <- 200
xb <- rnorm(n, 0, 1e4)
yb <- 0.001 * xb + RTMBdist::rskewnorm2(n, 0, 1.5, 4)
d3 <- data.frame(y = yb, x = xb)
cnt <- 0L
f3 <- withCallingHandlers(
  frm(bf(y ~ x, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = d3),
  warning = function(w) {
    cnt <<- cnt + 1L
    cat("   WARNING:", sub("\n.*", "", conditionMessage(w)), "\n")
    invokeRestart("muffleWarning")
  })
cat(sprintf("   ll %16.8f  conv %d  '%s'  warnings %d  esc %s\n",
            as.numeric(logLik(f3)), f3$opt$convergence, f3$opt$message,
            cnt, esc(f3)))
cat("   sn::selm ll ", format(sn::selm(y ~ x, family = "SN",
                                       data = d3)@logL), "\n")
options(op)
