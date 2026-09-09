# Round 2, item 5. The lane rejects the obvious fix for N1 with a
# number: "offsetting `pp` removes the NaN and returns a finite WRONG
# answer ... the coefficients come out 1/3 and 1/3 instead of 1/2 and
# 1/2 and the whole trajectory is two-thirds of its true height.
# Measured in dev/lincmt/lincmt-n1.R."
#
# dev/lincmt/lincmt-n1.R does NOT apply the offset. It measures the
# shipped NaN and then prints ordinary values under a heading that says
# "the offset must change nothing". So the two-thirds is not produced
# by the script that is cited for it. This applies the offset.
#
# Script path: dev/rev-lincmt-n1.R. No seed; every point constructed.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

tt <- c(0, 1, 4, 12, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                 addl = 2L)
P0 <- list(ke = 0, k12 = 0, k21 = 0, k13 = 0, k31 = 0, ka = 1, V = 1)
# with no elimination and no transfer, everything absorbed stays put
want <- vapply(tt, function(t) {
  s <- 0
  for (d in c(0, 12, 24)) if (d < t) s <- s + 100 * (1 - exp(-(t - d)))
  s
}, 0)

ns <- asNamespace("frmtmb.ode")
orig <- get("lincmt_disp", envir = ns)
mk_offset <- function(off) {
  src <- deparse(orig)
  hit <- grep("pp <- a1 - a2 * a2/3", src, fixed = TRUE)
  stopifnot(length(hit) == 1L)
  src[[hit]] <- paste0("    pp <- a1 - a2 * a2/3 - ",
                       format(off, digits = 17))
  fn <- eval(parse(text = paste(src, collapse = "\n")))
  environment(fn) <- ns
  fn
}
run <- function(tag) {
  v <- frm_lincmt(parms = P0, times = tt, ncmt = 3, depot = TRUE,
                  events = ev)
  d <- frmtmb.ode:::lincmt_disp(3L, 0, 0, 0, 0, 0)
  cf <- vapply(d[["coef"]], as.numeric, 0)
  ok <- want > 0
  cat(sprintf("%-24s finite %-5s  coef sum %8.5f  height ratio %8.5f\n",
              tag, all(is.finite(v)), sum(cf),
              if (any(ok)) mean(v[ok] / want[ok]) else NA_real_))
  cat("     coefficients ",
      paste(format(cf, digits = 6), collapse = "  "), "\n")
  cat("     value        ",
      paste(format(v, digits = 6), collapse = "  "), "\n")
  cat("     want         ",
      paste(format(want, digits = 6), collapse = "  "), "\n")
}
run("shipped")
unlockBinding("lincmt_disp", ns)
for (off in c(1e-150, 1e-300, .Machine$double.xmin)) {
  assign("lincmt_disp", mk_offset(off), envir = ns)
  run(paste("pp offset by", format(off, digits = 2)))
}
assign("lincmt_disp", orig, envir = ns)
stopifnot(identical(get("lincmt_disp", envir = ns), orig))

cat("\nordinary rates, before and after the 1e-150 offset:\n")
P1 <- list(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.05, k31 = 0.01,
           ka = 1.1, V = 10)
a <- frm_lincmt(parms = P1, times = tt, ncmt = 3, depot = TRUE,
                events = ev)
assign("lincmt_disp", mk_offset(1e-150), envir = ns)
b <- frm_lincmt(parms = P1, times = tt, ncmt = 3, depot = TRUE,
                events = ev)
assign("lincmt_disp", orig, envir = ns)
cat("  max relative change:",
    format(max(abs(a - b)) / max(abs(a)), digits = 4), "\n")
