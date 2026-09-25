# Reviewer 2, item 3: decision (b). Per-boundary defective functions
# against RWiener and WienR; their sum against the marginal; the whole
# censored objective against a likelihood written by hand with RWiener;
# dec() flips; trunc() refusal; NA dec(); rdm's left rows.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.eam); library(RWiener); library(WienR)
})
ns <- asNamespace("frmtmb.eam")
options(width = 160)

cat("== 1. per boundary against RWiener::pwiener and WienR::pWDM\n")
set.seed(3)
n <- 400
g <- data.frame(v = runif(n, -4, 4), a = runif(n, 0.5, 3), w = runif(n, 0.1, 0.9),
                u = exp(runif(n, log(0.02), log(5))))
g$t <- g$u * g$a^2
tau <- 0.2
lu <- ns$ddm_rt_lcdf_b(g$t, g$v, g$a, g$w, 1)
ll <- ns$ddm_rt_lcdf_b(g$t, g$v, g$a, g$w, 0)
lF <- ns$ddm_rt_lcdf2(g$t, g$v, g$a, g$w)$lF
rw_u <- mapply(function(t, v, a, w) pwiener(t + tau, a, tau, w, v, resp = "upper"),
               g$t, g$v, g$a, g$w)
rw_l <- mapply(function(t, v, a, w) pwiener(t + tau, a, tau, w, v, resp = "lower"),
               g$t, g$v, g$a, g$w)
wr_u <- pWDM(g$t + tau, "upper", g$a, g$v, g$w, tau)$value
wr_l <- pWDM(g$t + tau, "lower", g$a, g$v, g$w, tau)$value
cat(sprintf("upper: max |exp(lane) - RWiener| %.2e, - WienR %.2e; max rel to WienR %.2e\n",
            max(abs(exp(lu) - rw_u)), max(abs(exp(lu) - wr_u)),
            max(abs(exp(lu) - wr_u) / wr_u)))
cat(sprintf("lower: max |exp(lane) - RWiener| %.2e, - WienR %.2e; max rel to WienR %.2e\n",
            max(abs(exp(ll) - rw_l)), max(abs(exp(ll) - wr_l)),
            max(abs(exp(ll) - wr_l) / wr_l)))
cat(sprintf("does the upper label mean the same boundary? corr(lane up, WienR up) %.6f, (lane up, WienR low) %.6f\n",
            cor(exp(lu), wr_u), cor(exp(lu), wr_l)))
cat(sprintf("sum of the two defective functions against the marginal: max rel %.2e\n",
            max(abs(exp(ll) + exp(lu) - exp(lF)) / exp(lF))))

cat("\n== 2. the censored objective against RWiener by hand (dry run, no fit)\n")
set.seed(11)
d <- ddm_simulate(300, mu = 0.7, bs = 1.3, ndt = 0.25, bias = 0.5)
d$code <- 0L; d$y2 <- d$rt
r <- d$rt > 1.2; d$code[r] <- 1L; d$rt[r] <- 1.2
l <- d$code == 0L & d$rt < 0.42; d$code[l] <- -1L; d$rt[l] <- 0.42
pk <- which(d$code == 0L)[seq(3, sum(d$code == 0L), by = 4)]
lo <- floor(d$rt[pk] * 10) / 10
d$code[pk] <- 2L; d$y2[pk] <- lo + 0.1; d$rt[pk] <- pmax(lo, 0.42)
print(table(d$code))
o <- frm(bf(rt | dec(upper) + cens(code, y2) ~ 1, bs ~ 1, ndt = 0.25,
            bias = 0.5), family = wiener(max_ndt = 0.4), data = d, dry_run = "objective")
cat("parameters:", names(o$obj$par), "\n")
hand <- function(p, dd = d) {
  v <- p[1]; a <- exp(p[2]); w <- 0.5; t0 <- 0.25
  resp <- ifelse(dd$upper == 1, "upper", "lower")
  ll <- numeric(nrow(dd))
  for (i in seq_len(nrow(dd))) {
    y <- dd$rt[i]
    ll[i] <- switch(as.character(dd$code[i]),
      "0" = dwiener(y, a, t0, w, v, resp = resp[i], give_log = TRUE),
      "1" = log(1 - pwiener(y, a, t0, w, v, resp = "both")),
      "-1" = log(pwiener(y, a, t0, w, v, resp = resp[i])),
      "2" = log(pwiener(dd$y2[i], a, t0, w, v, resp = resp[i]) -
                  pwiener(y, a, t0, w, v, resp = resp[i])))
  }
  -sum(ll)
}
for (p in list(c(0.7, log(1.3)), c(0.2, log(0.9)), c(-0.5, log(2)), c(0, log(1.3)))) {
  f1 <- o$obj$fn(p); f2 <- hand(p)
  cat(sprintf("p = (%5.2f, %6.3f)  lane %.10f  RWiener hand %.10f  rel %.2e\n",
              p[1], p[2], f1, f2, abs(f1 - f2) / abs(f2)))
}

cat("\n== 3. dec() flips\n")
flip <- function(dd, codes) { k <- dd$code %in% codes; dd$upper[k] <- 1L - dd$upper[k]; dd }
p <- c(0.7, log(1.3))
for (cc in list(1L, -1L, 2L)) {
  o2 <- frm(bf(rt | dec(upper) + cens(code, y2) ~ 1, bs ~ 1, ndt = 0.25, bias = 0.5),
            family = wiener(max_ndt = 0.4), data = flip(d, cc), dry_run = "objective")
  cat(sprintf("flip dec on code %2d rows: fn %.12f vs %.12f  identical %s\n", cc,
              o2$obj$fn(p), o$obj$fn(p), identical(o2$obj$fn(p), o$obj$fn(p))))
}

cat("\n== 4. left / interval with trunc() refused; right with trunc() allowed\n")
try_msg <- function(expr) tryCatch({ expr; "NO ERROR" }, error = function(e)
  paste(class(e)[1], "|", conditionMessage(e)))
dl <- d[d$code != 2L, ]; dl <- dl[dl$rt < 3, ]
cat("left + trunc(ub): ", try_msg(frm(bf(rt | dec(upper) + cens(code) + trunc(ub = 3) ~ 1,
    bs ~ 1, ndt = 0.25, bias = 0.5), family = wiener(max_ndt = 0.4), data = dl, dry_run = "objective")), "\n")
cat("left + trunc(lb): ", try_msg(frm(bf(rt | dec(upper) + cens(code) + trunc(lb = 0.3) ~ 1,
    bs ~ 1, ndt = 0.25, bias = 0.5), family = wiener(max_ndt = 0.4), data = dl, dry_run = "objective")), "\n")
cat("interval + trunc: ", try_msg(frm(bf(rt | dec(upper) + cens(code, y2) + trunc(ub = 3) ~ 1,
    bs ~ 1, ndt = 0.25, bias = 0.5), family = wiener(max_ndt = 0.4), data = d, dry_run = "objective")), "\n")
dr <- d; dr$code[dr$code != 1L] <- 0L; dr$y2 <- NULL
cat("right + trunc(ub): ", try_msg(frm(bf(rt | dec(upper) + cens(code) + trunc(ub = 3) ~ 1,
    bs ~ 1, ndt = 0.25, bias = 0.5), family = wiener(max_ndt = 0.4), data = dr, dry_run = "objective")), "\n")
# right + trunc against RWiener by hand: P(y < T <= 3) / P(T <= 3) and f / P(T <= 3)
o3 <- frm(bf(rt | dec(upper) + cens(code) + trunc(ub = 3) ~ 1, bs ~ 1, ndt = 0.25,
             bias = 0.5), family = wiener(max_ndt = 0.4), data = dr, dry_run = "objective")
hand_tr <- function(p) {
  v <- p[1]; a <- exp(p[2]); t0 <- 0.25; w <- 0.5
  Fub <- pwiener(3, a, t0, w, v, resp = "both")
  resp <- ifelse(dr$upper == 1, "upper", "lower")
  ll <- vapply(seq_len(nrow(dr)), function(i) {
    if (dr$code[i] == 0L) dwiener(dr$rt[i], a, t0, w, v, resp = resp[i], give_log = TRUE)
    else log(Fub - pwiener(dr$rt[i], a, t0, w, v, resp = "both"))
  }, 0)
  -sum(ll - log(Fub))
}
for (p in list(c(0.7, log(1.3)), c(-0.3, log(1.8)))) {
  cat(sprintf("right + trunc(ub = 3): lane %.10f  hand %.10f  rel %.2e\n",
              o3$obj$fn(p), hand_tr(p), abs(o3$obj$fn(p) - hand_tr(p)) / abs(hand_tr(p))))
}

cat("\n== 5. a row with NA dec()\n")
dn <- d; dn$upper[c(5, 17)] <- NA
msgs <- character(0)
on <- withCallingHandlers(
  frm(bf(rt | dec(upper) + cens(code, y2) ~ 1, bs ~ 1, ndt = 0.25, bias = 0.5),
      family = wiener(max_ndt = 0.4), data = dn, dry_run = "objective"),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
  warning = function(w) { msgs <<- c(msgs, paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") })
cat("messages:", paste(msgs, collapse = " || "), "\n")
cat("codes of the NA rows:", d$code[c(5, 17)], "\n")
cat(sprintf("fn with NA rows dropped %.10f, hand on the data without them %.10f\n",
            on$obj$fn(p), hand(p, d[-c(5, 17), ])))

cat("\n== 6. rdm(): left rows do not read vint()\n")
set.seed(5)
dd <- rdm_simulate(300, v = c(3, 2), A = 0.5, k = 0.5, ndt = 0.2)
print(head(dd, 2))
dd$code <- 0L
k <- dd$rt < 0.45; dd$code[k] <- -1L; dd$rt[k] <- 0.45
nm_choice <- intersect(c("choice", "response", "winner"), names(dd))[1]
fr <- function(x) frm(stats::as.formula(sprintf(
  "rt | vint(%s) + cens(code) ~ 1", nm_choice)), family = rdm(2), data = x,
  dry_run = "objective")
o5 <- fr(dd)
dd2 <- dd; dd2[[nm_choice]][k] <- 3L - dd2[[nm_choice]][k]
o6 <- fr(dd2)
p5 <- o5$obj$par
cat(sprintf("left rows %d; flipping their winner: fn %.12f vs %.12f identical %s\n",
            sum(k), o5$obj$fn(p5), o6$obj$fn(p5), identical(o5$obj$fn(p5), o6$obj$fn(p5))))
