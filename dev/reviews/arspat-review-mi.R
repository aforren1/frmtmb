# Review of lane arspat: mi() beside the ar() family, the dropped-row
# margin against an independent dense margin over the full time-indexed
# covariance (row and column deleted). Output to stdout.
source("dev/reviews/arspat-review-common.R")

mk <- function(seed = 41, G = 8, K = 7, tdf = NULL, times = NULL) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  if (!is.null(times)) d$t <- times[d$t]
  n <- nrow(d)
  d$x <- rnorm(n)
  d$x2 <- rnorm(n)
  d$s <- runif(n, 0.1, 0.5)
  d$s0 <- d$s
  d$s0[c(2, 9, 30)] <- 0
  e <- as.vector(vapply(seq_len(G), function(i)
    as.vector(arima.sim(list(ar = 0.5, ma = 0.3), K, sd = 0.8)), numeric(K)))
  if (!is.null(tdf)) e <- e * rep(sqrt(tdf / rchisq(G, tdf)), each = K)
  d$y <- 1 + 0.6 * d$x + e + rnorm(n, 0, d$s)
  d
}
# missing patterns (K = 7, rows of group j are (j-1)*7 + 1:7)
pats <- list(
  interior = c(3, 12, 26),
  consecutive = c(3, 4, 5, 18, 19),
  start_end = c(1, 7, 8, 21, 22),
  nearly_all = c(1:6, 15:21),          # group 1 keeps 1 row, group 3 none
  mixed = c(1, 2, 7, 10, 11, 12, 35, 49)
)
terms <- c(ar = "ar(t, g, cov = TRUE)", ma = "ma(t, g, cov = TRUE)",
           arma = "arma(t, g, cov = TRUE)", cosy = "cosy(t, g)")

run1 <- function(tm, pat, se = FALSE, het = FALSE, fam = gaussian(),
                 times = NULL, s0 = FALSE) {
  d <- mk(tdf = if (fam$family == "student") 4, times = times)
  d$y[pats[[pat]]] <- NA
  sv <- if (s0) "s0" else "s"
  lhs <- if (se) paste0("y | mi() + se(", sv, ", sigma = TRUE)") else
    "y | mi()"
  fo <- stats::as.formula(paste(lhs, "~ x +", terms[[tm]]))
  b <- if (het) bf(fo, sigma ~ x2) else bf(fo)
  msgs <- character(0)
  fit <- withCallingHandlers(frm(b, family = fam, data = d),
    message = function(m) {
      msgs <<- c(msgs, class(m)[1]); invokeRestart("muffleMessage")
    })
  dp <- ns$eval_dpars(fit)[[1L]]
  n <- nrow(d)
  mu <- rep_len(dp$mu, n)
  sig <- rep_len(dp$sigma, n)
  nat <- ac_nat(fit)
  ar <- unname(nat[grep("^ar", names(nat))])
  ma <- unname(nat[grep("^ma", names(nat))])
  cz <- if (tm == "cosy") unname(nat[1]) else NULL
  tix <- match(d$t, sort(unique(d$t)))
  ref <- dense_margin(d$y, mu, sig, tix, d$g, ar, ma,
                      se = if (se) d[[sv]], nu = if (fam$family == "student")
                        dp$nu[1], cosy = cz)
  ll <- as.numeric(logLik(fit))
  cat(sprintf("%-5s %-12s se=%d het=%d %-8s times=%s s0=%d conv=%d ll=%.8f ref=%.8f rel=%.2e msgs=%s\n",
              tm, pat, se, het, fam$family, if (is.null(times)) "1:K" else
                "irreg", s0, fit$opt$convergence, ll, ref,
              abs(ll - ref) / abs(ref), paste(unique(msgs), collapse = ",")))
  invisible(fit)
}

for (tm in names(terms)) for (pat in names(pats)) {
  chk(paste(tm, pat), run1(tm, pat))
}
for (tm in c("ar", "arma")) for (pat in c("consecutive", "start_end", "mixed")) {
  chk(paste(tm, pat, "se"), run1(tm, pat, se = TRUE))
  chk(paste(tm, pat, "het"), run1(tm, pat, het = TRUE))
  chk(paste(tm, pat, "se het"), run1(tm, pat, se = TRUE, het = TRUE))
  chk(paste(tm, pat, "student"), run1(tm, pat, fam = student()))
  chk(paste(tm, pat, "student se"), run1(tm, pat, se = TRUE, fam = student()))
}
irr <- c(0.5, 1.25, 2, 4, 4.5, 7.1, 9)
for (tm in c("ar", "ma")) for (pat in c("interior", "start_end")) {
  chk(paste(tm, pat, "irregular"), run1(tm, pat, times = irr))
}
chk("ar mixed se with zero entries", run1("ar", "mixed", se = TRUE, s0 = TRUE))
chk("arma mixed se zero student", run1("arma", "mixed", se = TRUE, s0 = TRUE,
                                       fam = student()))
