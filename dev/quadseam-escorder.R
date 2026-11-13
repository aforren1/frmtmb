# Item 828, the order of the escape and the quadrature tape. Two arms on
# the lane build, in one process, interleaved per seed:
#   lane:   escape on the Laplace objective before the rule is anchored,
#           then the escape on the rule, then re-anchoring (quad_settle)
#   settle: the Laplace-side escape removed from quad_fit(), so the escape
#           runs on the rule frozen at the stall and re-anchoring repairs
# Data: wt-skewinit's acceptance stream (make_data(seed, FALSE) of
# dev/skewinit-lib.R, inlined) plus a 20-level (1 | g), started at the
# stall (`bad` start of dev/skewinit-punch1c.R). The reference is the
# same arm's default-start fit; gap = rescued - default logLik.
#   Rscript dev/quadseam-escorder.R <lib> [seeds]
# log dev/quadseam-log/escorder-lib.txt
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
S <- if (length(a) > 1) as.integer(a[2]) else 10L
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
qf_lane <- get("quad_fit", ns)
qf_settle <- qf_lane
b <- body(qf_settle)
drop <- vapply(seq_along(b), function(i) {
  any(grepl("escape_stationary", deparse(b[[i]])))
}, NA)
stopifnot(sum(drop) == 1L)
body(qf_settle) <- b[!drop]
use <- function(f) {
  unlockBinding("quad_fit", ns)
  assign("quad_fit", f, envir = ns)
  lockBinding("quad_fit", ns)
}
skew <- function(v) mean((v - mean(v))^3) / stats::sd(v)^3
make <- function(seed) {
  set.seed(seed)
  n <- 200
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs)
  set.seed(99)
  dd$g <- factor(rep(1:20, length.out = n))
  dd
}
fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
fit <- function(dd, start = NULL) {
  t0 <- proc.time()[["elapsed"]]
  w <- character()
  f <- withCallingHandlers(
    frm(fo, family = skew_normal(), data = dd, quadrature = TRUE,
        start = start),
    warning = function(x) {
      w <<- c(w, class(x)[1])
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  list(ll = as.numeric(logLik(f)), secs = proc.time()[["elapsed"]] - t0,
       esc = f$opt[["stationary_escape"]], warn = w)
}
for (seed in seq_len(S)) {
  dd <- make(seed)
  m3 <- skew(dd$y)
  bad <- list(betad = c(log(stats::sd(dd$y)), 2 * sign(m3) + 0.5 * m3))
  for (arm in c("lane", "settle")) {
    use(if (arm == "lane") qf_lane else qf_settle)
    d0 <- fit(dd)
    d1 <- fit(dd, bad)
    cat(sprintf(paste0("SEED %2d %-6s default %.9f rescued %.9f gap %10.3e",
                       " esc_gain %s secs %.2f %.2f warn %s\n"),
                seed, arm, d0$ll, d1$ll, d1$ll - d0$ll,
                if (is.null(d1$esc)) "-" else format(d1$esc[["gain"]],
                                                    digits = 5),
                d0$secs, d1$secs,
                paste(c(d0$warn, "|", d1$warn), collapse = ",")))
  }
}
use(qf_lane)
cat("END\n")
