# Review of lane quadseam2, attack 1: termination and cost of the
# repeated stationary escape. Three skew-normal components (mixture and
# hmm) started with every alpha at 0, on skewed data and on symmetric
# data (where alpha = 0 is legitimate and a pass can "gain" by rounding).
# Per fit: starts, gain, logLik, fit seconds, alphas.
#   Rscript quadseam2-review-passes.R <base|lib>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.latent)})
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
pl <- new.env(); pl$n <- 0L
`%||%` <- function(a, b) if (is.null(a)) b else a
# count the passes: each pass after the first calls stationary_escapes()
# again from inside escape_stationary()
orig <- get("stationary_escapes", ns)
unlockBinding("stationary_escapes", ns)
assign("stationary_escapes", function(obj, opt, frame) {
  o <- orig(obj, opt, frame)
  if (length(o)) pl$n <- pl$n + 1L
  o
}, envir = ns)
sn_draw <- function(n, mu, sigma, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  mu + sigma * (dl * abs(rnorm(n)) + sqrt(1 - dl^2) * rnorm(n))
}
run <- function(lbl, s, fo, fam, d) {
  u <- frm(fo, family = fam, data = d, dry_run = "objective")
  bd <- u$estimates$betad
  al <- grep("^alpha", names(bd))
  bd[al] <- 0
  pl$n <- 0L
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(suppressWarnings(suppressMessages(
    frm(fo, family = fam, data = d, start = list(betad = bd)))),
    error = function(e) e)
  t1 <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("ROW %s %d ERROR %s\n", lbl, s, conditionMessage(f)))
    return(invisible())
  }
  t2 <- proc.time()[["elapsed"]]
  fr <- suppressWarnings(suppressMessages(frm(fo, family = fam, data = d)))
  t3 <- proc.time()[["elapsed"]] - t2
  esc <- f$opt$stationary_escape
  cat(sprintf(paste0("ROW %s %d starts %s gain %s nonempty_calls %d ll %.9f ",
                     "free %.9f fit %.2f free %.2f alpha %s\n"), lbl, s,
              esc[["starts"]] %||% "-", format(esc[["gain"]] %||% NA, digits = 4),
              pl$n, as.numeric(logLik(f)), as.numeric(logLik(fr)), t1, t3,
              paste(format(f$estimates$betad[al], digits = 3), collapse = ",")))
}
for (s in 1:3) {
  set.seed(200 + s)
  n <- 600
  z <- sample(3, n, replace = TRUE)
  mu <- c(-5, 0, 5)[z]; al <- c(5, -5, 5)[z]
  y <- vapply(seq_len(n), function(i) sn_draw(1, mu[i], 1, al[i]), 0)
  run("MIX3skew", s, bf(y ~ 1), mixture(skew_normal(), skew_normal(),
                                         skew_normal()), data.frame(y = y))
  y0 <- c(-5, 0, 5)[z] + rnorm(n)
  run("MIX3sym", s, bf(y ~ 1), mixture(skew_normal(), skew_normal(),
                                        skew_normal()), data.frame(y = y0))
  y2 <- ifelse(z > 1, rnorm(n, 3), rnorm(n, -3))
  run("MIX2sym", s, bf(y ~ 1), mixture(skew_normal(), skew_normal()),
      data.frame(y = y2))
  # hmm(2) symmetric: both alphas truly 0
  G <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
  hd <- do.call(rbind, lapply(1:20, function(id) {
    st <- integer(40); st[1] <- sample.int(2, 1)
    for (t in 2:40) st[t] <- sample.int(2, 1, prob = G[st[t - 1], ])
    data.frame(id = id, t = 1:40, y = ifelse(st == 1, rnorm(40, -3),
                                             rnorm(40, 3)))
  }))
  run("HMM2sym", s, bf(y ~ 1), hmm(K = 2, skew_normal(), time = t, group = id),
      hd)
}
cat("END\n")
