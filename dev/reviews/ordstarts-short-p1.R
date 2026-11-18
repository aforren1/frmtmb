# Reviewer of lane ordstarts, build round. Probe 1 (questions 1, 3, 6):
# - weights scale: the same data with constant case weights 1, 1/n and
#   10 on seeds whose swap the lane kept (cum_sratio, no random effects,
#   n 1000): is the swap still kept, are the jittered starts still run?
# - RNG: .Random.seed and RNGkind before and after a fit that runs the
#   swap and its holds check; no .Random.seed before a fit.
# - the hint: on default fits that warn "run them with
#   frmtmb_control(starts = 10)", does starts = 10 run any start?
#   Rscript ordstarts-build-review-p1.R > ordstarts-short-p1.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
grab <- function(expr) {
  w <- character()
  mc <- character()
  v <- withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      mc <<- c(mc, class(x)[1])
      invokeRestart("muffleMessage")
    })
  list(v = v, w = w, mc = mc)
}
rec <- function(f) {
  ms <- f$opt$mixture_starts
  if (is.null(ms)) return("record NULL")
  paste(names(ms), signif(ms, 4), sep = "=", collapse = " ")
}
cls_of <- function(f) {
  cl <- ns$mixture_ord_start_class(f, "y")
  if (length(cl$step)) "step" else if (length(cl$boundary)) "flat" else
    "proper"
}
fam <- function() mixture(cumulative(), sratio())

cat("\n== 1. Weights scale (cum_sratio, no random effects, n 1000) ==\n")
n <- 1000
for (seed in c(1, 2, 4, 6, 10, 13)) {
  d <- mk(seed, n)
  ref <- grab(frm(bf(y ~ x + z), family = fam(), data = d))$v
  for (wv in c(1, 1 / n, 10)) {
    d$w <- wv
    a0 <- grab(frm(bf(y | weights(w) ~ x + z), family = fam(), data = d,
                   control = frmtmb_control(starts = 0)))
    dd <- grab(frm(bf(y | weights(w) ~ x + z), family = fam(), data = d))
    s10 <- grab(frm(bf(y | weights(w) ~ x + z), family = fam(), data = d,
                    control = frmtmb_control(starts = 10)))
    se0 <- ns$mixture_ord_thres_se(a0$v, "y", force = TRUE, lost_inf = TRUE)
    cat(sprintf(paste0("seed %2d w %-8s | A0 class %-6s thres_se(lost_inf) %s",
                       " lost %s | D %s | D par vs unweighted D max|d| %.3g",
                       " | D degen warn %s | starts=10: %s\n"),
                seed, format(signif(wv, 3)), cls_of(a0$v),
                paste(signif(se0, 3), collapse = ","),
                paste(attr(se0, "lost"), collapse = ","), rec(dd$v),
                max(abs(dd$v$opt$par - ref$opt$par)),
                any(grepl("degenerate boundary", dd$w, fixed = TRUE)),
                rec(s10$v)))
  }
  cat(sprintf("seed %2d unweighted D: %s ll %.4f\n", seed, rec(ref),
              as.numeric(logLik(ref))))
}

cat("\n== 3b. Column scale: the holds jitter is 0.1 in raw units ==\n")
for (seed in c(1, 2, 4, 6, 10, 13)) {
  d <- mk(seed, n)
  for (sc in c(1, 1000, 1e-3)) {
    ds <- d
    ds$x <- ds$x * sc
    r <- grab(frm(bf(y ~ x + z), family = fam(), data = ds,
                  control = frmtmb_control(autoscale = FALSE)))
    cat(sprintf("seed %2d x * %-6g autoscale off | ll %.4f | %s\n", seed, sc,
                as.numeric(logLik(r$v)), rec(r$v)))
  }
}

cat("\n== 1b. The test's sound fit (cum_cum, seed 3, n 500, weights 1/500)\n")
omx_wdata <- function(seed, n, w) {
  set.seed(seed)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5), w = w)
}
for (wv in c(1, 1 / 500)) {
  dw <- omx_wdata(3, 500, wv)
  r <- grab(frm(bf(y | weights(w) ~ x), family = mixture(cumulative(),
                                                         cumulative()),
                data = dw, control = frmtmb_control(starts = 10)))
  cat(sprintf("w %-8s class at the fit: %s | starts=10 record: %s | warnings: %d\n",
              format(signif(wv, 3)), cls_of(r$v), rec(r$v), length(r$w)))
}

cat("\n== 3. RNG state (cum_sratio n 1000 seed 13, swap kept) ==\n")
d <- mk(13, 1000)
set.seed(42)
s0 <- .Random.seed
k0 <- RNGkind()
f <- grab(frm(bf(y ~ x + z), family = fam(), data = d))$v
cat("record:", rec(f), "\n")
cat("seed identical after fit:", identical(.Random.seed, s0),
    "; kind identical:", identical(RNGkind(), k0), "\n")
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
s1 <- .Random.seed
f2 <- grab(frm(bf(y ~ x + z), family = fam(), data = d))$v
cat("L'Ecuyer: seed identical:", identical(.Random.seed, s1), "; kind:",
    RNGkind()[1], "; same par as Mersenne run:",
    identical(f$opt$par, f2$opt$par), "\n")
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
if (exists(".Random.seed", envir = globalenv())) {
  rm(".Random.seed", envir = globalenv())
}
f3 <- grab(frm(bf(y ~ x + z), family = fam(), data = d))$v
cat("no seed before: exists after:",
    exists(".Random.seed", envir = globalenv()), "; same par:",
    identical(f$opt$par, f3$opt$par), "\n")

cat("\n== 6. The hint: does starts = 10 run any start where it is named? ==\n")
for (cfg in list(c("cc", 300), c("cs", 300))) {
  fm <- if (cfg[1] == "cc") mixture(cumulative(), cumulative()) else fam()
  nh <- 0
  tab <- character()
  for (seed in 1:40) {
    d <- mk(seed, as.integer(cfg[2]))
    r <- grab(frm(bf(y ~ x + z), family = fm, data = d))
    hint <- any(grepl("frmtmb_control(starts = 10)", r$w, fixed = TRUE))
    if (!hint) next
    nh <- nh + 1
    s <- grab(frm(bf(y ~ x + z), family = fm, data = d,
                  control = frmtmb_control(starts = 10)))
    ms <- s$v$opt$mixture_starts
    ran <- if (is.null(ms)) 0 else ms[["starts"]]
    still <- any(grepl("degenerate boundary", s$w, fixed = TRUE))
    moved <- as.numeric(logLik(s$v)) - as.numeric(logLik(r$v))
    tab <- c(tab, sprintf("%s n%s seed %2d class %-6s | starts=10 ran %d, moved %+.4f, still warns %s, text: %s",
                          cfg[1], cfg[2], seed, cls_of(r$v), as.integer(ran),
                          moved, still,
                          if (any(grepl("Compare several starts", s$w)))
                            "Compare several starts" else if
                          (any(grepl("None of", s$w))) "None of n" else "-"))
  }
  cat(tab, sep = "\n")
  cat(cfg[1], cfg[2], ": fits whose warning names starts = 10:", nh,
      "; of them starts = 10 ran no start:",
      sum(grepl("ran 0,", tab)), "\n")
}
cat("DONE p1\n")
