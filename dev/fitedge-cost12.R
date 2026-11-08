# Items 1 and 2: the cost of each fix, with the fix on and off in one
# process, order alternating by row. Off is the lane's build with the fix's
# entry patched out: gp_body_start() returning NULL (item 1, 0.79.1's start
# for a body), nl_all_bodies() returning no body (item 2, no ridge pass).
# Rows: dev/fitedge-gpbody-sweep.R's bodies A and C (all 48), and
# dev/fitedge-nl2-sweep.R's disc group (30). Args: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
patch <- function(nm, f) {
  unlockBinding(nm, ns)
  assign(nm, f, envir = ns)
  lockBinding(nm, ns)
}
o_body <- ns$gp_body_start
o_all <- ns$nl_all_bodies
timed <- function(on, set, ...) {
  if (set == "item1") patch("gp_body_start", if (on) o_body else function(...) NULL)
  if (set == "item2") patch("nl_all_bodies", if (on) o_all else function(frame) list())
  t0 <- proc.time()[[3]]
  r <- suppressMessages(fe_fit(...))
  t <- proc.time()[[3]] - t0
  patch("gp_body_start", o_body)
  patch("nl_all_bodies", o_all)
  list(t = t, ll = if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit)))
}
pair <- function(i, set, ...) {
  if (i %% 2L) {
    a <- timed(TRUE, set, ...); b <- timed(FALSE, set, ...)
  } else {
    b <- timed(FALSE, set, ...); a <- timed(TRUE, set, ...)
  }
  c(on = a$t, off = b$t, ll_on = a$ll, ll_off = b$ll)
}
curves <- c("sin", "wig", "exp", "step")
mk1 <- function(seed, n, ci, body) {
  set.seed(18000 + 100 * seed + n + ci)
  x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
  f <- switch(curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
              exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
  y <- switch(body,
    A = f + 0.3 * z + 0.1 * rnorm(n),
    C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n))
  data.frame(x, z, w, y)
}
nl1 <- list(A = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
            C = bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE))
res <- list()
i <- 0L
for (body in c("A", "C")) for (ci in 1:4) for (n in c(60L, 120L)) for (seed in 1:3) {
  i <- i + 1L
  res[[length(res) + 1L]] <- c(set = 1, pair(i, "item1", nl1[[body]],
                                             data = mk1(seed, n, ci, body)))
}
gen <- function(seed, slope) {
  set.seed(seed); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z
  d$yo <- ordered(cut(lat * exp(0.3 * d$z) + rlogis(n),
                      c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  d
}
mub <- list(exp = "b1 * exp(b2 * x)", sat = "b1 * (1 - exp(-b2 * x))",
            two = "b1 * exp(b2 * x) + b3 * z")
for (seed in 1:5) for (slope in c(0.6, -0.6)) for (bn in names(mub)) {
  i <- i + 1L
  p <- c(setdiff(all.vars(str2lang(mub[[bn]])), c("x", "z")), "d1", "d2")
  f <- bf(yo ~ 1) + nlf(as.formula(paste("mu ~", mub[[bn]]))) +
    nlf(disc ~ exp(d1 * d2 * z))
  for (q in p) f <- f + lf(as.formula(paste(q, "~ 1")))
  res[[length(res) + 1L]] <- c(set = 2, pair(i, "item2", f,
                                             family = cumulative(),
                                             data = gen(seed * 100 + 1, slope)))
}
m <- as.data.frame(do.call(rbind, res))
for (s in 1:2) {
  d <- m[m$set == s, ]
  cat(sprintf("COST item %d rows %d: fix on %.1f s, off %.1f s, ratio %.2f; ll on above off %d, below %d\n",
              s, nrow(d), sum(d$on), sum(d$off), sum(d$on) / sum(d$off),
              sum(d$ll_on > d$ll_off + 1e-3, na.rm = TRUE),
              sum(d$ll_on < d$ll_off - 1e-3, na.rm = TRUE)))
}
cat("DONE\n")
