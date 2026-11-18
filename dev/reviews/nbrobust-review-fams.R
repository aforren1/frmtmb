# Review of lane nbrobust, check 2: every family that now reads nb_lpmf(),
# on every mean link it allows, against a 300-bit Rmpfr reference written
# from the textbook masses (no frmtmb code), on the tape (with .eta_) and
# off it (the post-fit route). Gradient reference: central differences of
# the mpfr value in the linear predictors.
#   Rscript dev/reviews/nbrobust-review-fams.R <base|lib>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages({
  library(frmtmb)
  library(Rmpfr)
})
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
P <- 300
eps <- .Machine$double.eps
m <- function(x) mpfr(x, P)
inv <- list(
  log = function(e) exp(e),
  identity = function(e) e,
  sqrt = function(e) e^2,
  softplus = function(e) log1p(exp(e)),
  squareplus = function(e) (e + sqrt(e^2 + 4)) / 2)
fwd <- list(log = log, identity = identity, sqrt = sqrt,
            softplus = function(mu) log(expm1(mu)),
            squareplus = function(mu) (mu^2 - 1) / mu)
lnb <- function(y, mu, n) {
  lgamma(y + n) - lgamma(n) - lgamma(m(y + 1)) +
    n * log(n / (n + mu)) + y * log(mu / (n + mu))
}
plogis_m <- function(e) 1 / (1 + exp(-e))
# reference log-likelihood (sum over rows) at linear predictors `e`
ref_ll <- function(fam, lk, y, e, d = NULL) {
  mu <- inv[[lk]](e[["mu"]])
  switch(fam,
    negbinomial = {
      n <- exp(e[["shape"]])
      if (!is.null(d)) {
        mu <- mu * d
        n <- n * d
      }
      sum(lnb(y, mu, n))
    },
    negbinomial2 = sum(lnb(y, mu, 1 / exp(e[["sigma"]]))),
    nbinom1 = sum(lnb(y, mu, mu / exp(e[["phi"]]))),
    geometric = {
      n <- m(1)
      if (!is.null(d)) {
        mu <- mu * d
        n <- n * d
      }
      sum(lnb(y, mu, n))
    },
    zinb = {
      n <- exp(e[["shape"]])
      z <- plogis_m(e[["zi"]])
      l <- lnb(y, mu, n)
      p0 <- n * log(n / (n + mu))
      r <- log(1 - z) + l
      r[y == 0] <- log(z + (1 - z) * exp(p0))
      sum(r)
    },
    hurdle = {
      n <- exp(e[["shape"]])
      h <- plogis_m(e[["hu"]])
      l <- lnb(y, mu, n)
      p0 <- n * log(n / (n + mu))
      r <- log(1 - h) + l - log(1 - exp(p0))
      r[y == 0] <- log(h)
      sum(r)
    })
}
make_fam <- function(fam, lk) {
  switch(fam,
    negbinomial = negbinomial(link = lk),
    negbinomial2 = brmsfamily("negbinomial2", link = lk),
    nbinom1 = nbinom1(link = lk),
    geometric = geometric(link = lk),
    zinb = zero_inflated_negbinomial(link = lk),
    hurdle = hurdle_negbinomial(link = lk))
}
dp_of <- function(f, e, tape) {
  out <- list()
  for (nm in names(e)) {
    out[[nm]] <- f$links[[nm]]$linkinv(e[[nm]])
    if (tape) out[[paste0(".eta_", nm)]] <- e[[nm]]
  }
  out
}
set.seed(11)
y <- c(0, 0, 1, 2, 3, 5, 9, 24, 25, 40, 130)
fams <- c("negbinomial", "negbinomial2", "nbinom1", "geometric", "zinb",
          "hurdle")
other <- list(negbinomial = "shape", negbinomial2 = "sigma",
              nbinom1 = "phi", geometric = character(0),
              zinb = c("shape", "zi"), hurdle = c("shape", "hu"))
rows <- list()
for (fam in fams) for (lk in names(inv)) for (rep in 1:8) for (rt in c(FALSE, TRUE)) {
  if (rt && !fam %in% c("negbinomial", "geometric")) next
  f <- make_fam(fam, lk)
  mu0 <- exp(stats::runif(1, log(0.05), log(60)))
  e <- c(mu = fwd[[lk]](mu0))
  lsh <- if (rep <= 4) stats::runif(1, -2, 6) else stats::runif(1, 15, 40)
  for (o in other[[fam]]) {
    e[[o]] <- switch(o, shape = lsh, sigma = -lsh,
                     phi = log(mu0) - lsh,
                     zi = stats::qlogis(stats::runif(1, 0.05, 0.9)),
                     hu = stats::qlogis(stats::runif(1, 0.05, 0.9)))
  }
  d <- if (rt) c(0.5, 2, 1, 3, 1, 1, 0.7, 1, 2, 1, 4) else NULL
  at <- if (rt) list(rate = d) else list()
  fn <- function(p) {
    ee <- as.list(p)
    names(ee) <- names(e)
    sum(f$lpdf(y, dp_of(f, ee, TRUE), at))
  }
  tp <- RTMB::MakeTape(fn, unname(e))
  v <- tp(unname(e))
  g <- as.numeric(tp$jacobian(unname(e)))
  ee <- as.list(e)
  voff <- sum(f$lpdf(y, dp_of(f, ee, FALSE), at))
  em <- lapply(ee, m)
  vr <- ref_ll(fam, lk, y, em, d)
  h <- m(2)^-90
  gr <- vapply(names(e), function(nm) {
    up <- em
    dn <- em
    up[[nm]] <- up[[nm]] + h
    dn[[nm]] <- dn[[nm]] - h
    as.numeric((ref_ll(fam, lk, y, up, d) - ref_ll(fam, lk, y, dn, d)) /
                 (2 * h))
  }, 0)
  unit <- eps * sum(y + mu0 * (if (rt) d else 1) + 1)
  rows[[length(rows) + 1L]] <- data.frame(
    fam = fam, link = lk, rate = rt, rep = rep, mu = mu0, lsh = lsh,
    v = v, voff = voff, vref = as.numeric(vr),
    absv = abs(v - as.numeric(vr)), absvoff = abs(voff - as.numeric(vr)),
    gmu = g[1], gmu_ref = gr[1],
    umu = abs(g[1] - gr[1]) / max(unit, abs(gr[1]) * eps),
    gsh = if (length(g) > 1) g[2] else NA,
    gsh_ref = if (length(g) > 1) gr[2] else NA,
    ush = if (length(g) > 1) abs(g[2] - gr[2]) / unit else NA,
    gz = if (length(g) > 2) g[3] else NA,
    uz = if (length(g) > 2) abs(g[3] - gr[3]) / unit else NA)
}
tab <- do.call(rbind, rows)
write.table(tab, paste0("dev/reviews/nbrobust-review-fams-", arm, ".tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nper family and link: max abs value error (tape, off tape), max",
    "gradient error in units of eps * sum(y + mu + 1)\n")
for (fam in fams) for (lk in names(inv)) for (big in c(FALSE, TRUE)) {
  s <- tab[tab$fam == fam & tab$link == lk & (tab$lsh > 10) == big, ]
  cat(sprintf(paste("%-12s %-10s lsh %-5s n %2d  val %9.3g off %9.3g",
                    " umu %8.3g ush %9.3g uz %8.3g\n"),
              fam, lk, if (big) ">10" else "<=10", nrow(s), max(s$absv),
              max(s$absvoff), max(s$umu), max(s$ush, na.rm = TRUE),
              max(c(s$uz, 0), na.rm = TRUE)))
}
cat("\nbad rows (value err > 1e-9, or any gradient units > 1000):\n")
bad <- tab[tab$absv > 1e-9 | tab$absvoff > 1e-9 | tab$umu > 1000 |
             (!is.na(tab$ush) & tab$ush > 1000) |
             (!is.na(tab$uz) & tab$uz > 1000), ]
cat("count", nrow(bad), "of", nrow(tab), "\n")
print(head(bad[, c("fam", "link", "rate", "mu", "lsh", "absv", "absvoff",
                   "umu", "ush", "uz")], 40))

# identity link at a nonpositive mean: what the density gives, on and off
# the tape, value and gradient
cat("\nidentity link at mu <= 0\n")
for (fam in fams) for (mu in c(-0.5, -1e-8, 0)) {
  f <- make_fam(fam, "identity")
  e <- c(mu = mu)
  for (o in other[[fam]]) {
    e[[o]] <- switch(o, shape = 1, sigma = -1, phi = 0, zi = 0, hu = 0)
  }
  fn <- function(p) {
    ee <- as.list(p)
    names(ee) <- names(e)
    sum(f$lpdf(c(0, 2), dp_of(f, ee, TRUE), list()))
  }
  r <- tryCatch({
    tp <- RTMB::MakeTape(fn, unname(e))
    c(tp(unname(e)), as.numeric(tp$jacobian(unname(e)))[1])
  }, error = function(err) conditionMessage(err))
  off <- tryCatch(suppressWarnings(
    f$lpdf(c(0, 2), dp_of(f, as.list(e), FALSE), list())),
    error = function(err) conditionMessage(err))
  cat(sprintf("%-12s mu %6g  tape value, d/deta_mu: %s   off-tape rows: %s\n",
              fam, mu, paste(format(r, digits = 5), collapse = " "),
              paste(format(off, digits = 5), collapse = " ")))
}

# NA in the response
cat("\nNA response row\n")
for (fam in fams) {
  f <- make_fam(fam, "log")
  e <- c(mu = 0)
  for (o in other[[fam]]) {
    e[[o]] <- switch(o, shape = 1, sigma = -1, phi = 0, zi = 0, hu = 0)
  }
  r <- tryCatch(f$lpdf(c(1, NA, 3), dp_of(f, as.list(e), FALSE), list()),
                error = function(err) paste("ERROR", conditionMessage(err)))
  cat(sprintf("%-12s %s\n", fam, paste(format(r, digits = 5), collapse = " ")))
}
