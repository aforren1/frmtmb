# Short check, item 2: a missing value of a registered term's response
# that another density reads stays latent; one nothing reads leaves.
#  B1/B2: y | mi() ~ x + fcor(M) (gaussian / student), z | mi() ~ mi(y) + x
#    (gaussian), rescor off. y missing at 3, 11, 19, 27 (z observed: read)
#    and 6, 14, 33 (z missing too: unread); z also missing at 9, 22.
#    Reference: Laplace over the read values of [the y density of the rows
#    other than the unread ones, M cut to them] times z's density
#    (gaussian: exact, so also the dense joint normal margin).
#  B3: z | mi() ~ mi(y) + x + fcor(M): a registered reader keeps every row.
#  B4: ar(cov = FALSE) reading y's lag (not a registered term).
#  B5: ar() and fcor() on one response.
# Lane writes the points to an rds; the base arm (SHORT_LIB="") reads them
# and evaluates its own objective there (base's gaussian route is all
# latent, Laplace exact on a quadratic).
source("dev/reviews/arspat2-short-common.R")
ptsf <- file.path(outdir, "read-points.rds")
P <- if (arm == "base") readRDS(ptsf) else list()
V <- list()
n <- 40
M <- 0.6^abs(outer(1:n, 1:n, "-"))
R1 <- c(3, 11, 19, 27)
R2 <- c(6, 14, 33)
R3 <- c(9, 22)
mkd <- function(df = Inf) {
  set.seed(7)
  d <- data.frame(x = rnorm(n), t = seq_len(n))
  sc <- if (is.finite(df)) sqrt(df / rchisq(1, df)) else 1
  d$y <- 1 + 0.5 * d$x + as.vector(t(chol(M)) %*% rnorm(n)) * sc
  d$z <- 0.3 + 0.8 * d$y - 0.4 * d$x + rnorm(n, 0, 0.7)
  d$y[c(R1, R2)] <- NA
  d$z[c(R2, R3)] <- NA
  d
}
# Laplace reference over the read values; p by label
ref_fn <- function(d, p, lab, student) {
  S <- setdiff(seq_len(n), R2)
  pos <- match(R1, S)
  MS <- M[S, S]
  LS <- t(chol(MS))
  ldM <- 2 * sum(log(diag(LS)))
  Li <- forwardsolve(LS, diag(length(S)))
  zo <- which(!is.na(d$z))
  g <- function(nm) unname(p[lab == nm])
  f <- function(q) {
    "[<-" <- RTMB::ADoverload("[<-")
    "c" <- RTMB::ADoverload("c")
    y0 <- d$y
    y0[is.na(y0)] <- 0
    yS <- RTMB::advector(y0[S])
    yS[pos] <- q$ylat
    sy <- exp(g("y_sigma_(Intercept)"))
    r <- (yS - g("y_(Intercept)") - g("y_x") * d$x[S]) / sy
    u <- Li %*% r
    k <- length(S)
    Q <- sum(u * u)
    lly <- if (student) {
      nu <- 1 + exp(g("y_nu_(Intercept)"))
      lgamma((nu + k) / 2) - lgamma(nu / 2) - 0.5 * k * log(nu * pi) -
        0.5 * ldM - k * log(sy) - 0.5 * (nu + k) * log1p(Q / nu)
    } else {
      -0.5 * k * log(2 * pi) - 0.5 * ldM - k * log(sy) - 0.5 * Q
    }
    yall <- RTMB::advector(y0)
    yall[R1] <- q$ylat
    mz <- g("z_(Intercept)") + g("z_miy") * yall[zo] + g("z_x") * d$x[zo]
    llz <- sum(RTMB::dnorm(d$z[zo], mz, exp(g("z_sigma_(Intercept)")),
                           log = TRUE))
    -(lly + llz)
  }
  obj <- RTMB::MakeADFun(f, list(ylat = numeric(length(R1))),
                         random = "ylat", silent = TRUE)
  obj$fn(numeric(0))
}
# the dense joint normal margin of (y_o, z_o), gaussian only
dense_fn <- function(d, p, lab) {
  g <- function(nm) unname(p[lab == nm])
  sy <- exp(g("y_sigma_(Intercept)"))
  sz <- exp(g("z_sigma_(Intercept)"))
  c1 <- g("z_miy")
  my <- g("y_(Intercept)") + g("y_x") * d$x
  mz <- g("z_(Intercept)") + c1 * my + g("z_x") * d$x
  Syy <- sy^2 * M
  Szz <- c1^2 * Syy + diag(sz^2, n)
  Syz <- c1 * Syy
  S <- rbind(cbind(Syy, Syz), cbind(t(Syz), Szz))
  v <- c(d$y, d$z)
  m <- c(my, mz)
  o <- !is.na(v)
  So <- S[o, o]
  r <- v[o] - m[o]
  L <- chol(So)
  w <- backsolve(L, r, transpose = TRUE)
  0.5 * (sum(o) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(w * w))
}
run_mv <- function(key, student) {
  d <- mkd(if (student) 4 else Inf)
  fam_y <- if (student) student() else gaussian()
  fit <- quiet(frm(bf(y | mi() ~ x + fcor(M), family = fam_y) +
                     bf(z | mi() ~ mi(y) + x, family = gaussian()) +
                     set_rescor(FALSE), data = d, data2 = list(M = M)))
  lab <- ns$outer_par_labels(fit)
  mm <- fit$frame$mi_map
  cat(sprintf("%s [%s]: y latent rows %s, y unread %s, z unread %s, length(miss) %d, conv %d, exact flag %s\n",
              key, arm, paste(mm$y$rows, collapse = ","),
              paste(mm$y$unread %||% integer(0), collapse = ","),
              paste(mm$z$unread %||% integer(0), collapse = ","),
              length(fit$estimates$miss), fit$opt$convergence,
              ns$mi_laplace_exact(fit$frame)))
  if (arm == "lane") {
    P[[key]] <<- lapply(0:2, function(k) {
      p <- pts(fit$opt$par, k)
      if (student) p[grep("nu", lab)] <- log(c(3, 6, 21)[k + 1L] - 1)
      p
    })
  }
  val <- vapply(P[[key]], function(p) fit$obj$fn(p), 0)
  V[[key]] <<- val
  rf <- vapply(P[[key]], function(p) ref_fn(d, p, lab, student), 0)
  cat(sprintf("  vs Laplace reference (unread rows cut out): rel %s\n",
              paste(sprintf("%.1e", reld(val, rf)), collapse = " ")))
  if (!student) {
    dn <- vapply(P[[key]], function(p) dense_fn(d, p, lab), 0)
    cat(sprintf("  vs dense joint normal margin: rel %s\n",
                paste(sprintf("%.1e", reld(val, dn)), collapse = " ")))
  }
  cat(sprintf("  logLik %.10f\n", as.numeric(logLik(fit))))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
run_mv("B1 gaussian", FALSE)
run_mv("B2 student", TRUE)

# B3: a registered reader keeps every row of y
d <- mkd()
f3 <- quiet(frm(bf(y | mi() ~ x + fcor(M)) + bf(z | mi() ~ mi(y) + x + fcor(M)) +
                  set_rescor(FALSE), data = d, data2 = list(M = M),
                family = gaussian()))
mm <- f3$frame$mi_map
if (arm == "lane") P$B3 <- lapply(0:2, function(k) pts(f3$opt$par, k))
V$B3 <- vapply(P$B3, function(p) f3$obj$fn(p), 0)
cat(sprintf("B3 registered reader [%s]: y latent %s, y unread %s, z latent %s, z unread %s, conv %d, logLik %.10f\n",
            arm, paste(mm$y$rows, collapse = ","),
            paste(mm$y$unread %||% integer(0), collapse = ","),
            paste(mm$z$rows, collapse = ","),
            paste(mm$z$unread %||% integer(0), collapse = ","),
            f3$opt$convergence, as.numeric(logLik(f3))))
# B4: ar(cov = FALSE) reads y's lag
d4 <- mkd()
f4 <- quiet(frm(y | mi() ~ x + ar(t, p = 1, cov = FALSE), data = d4,
                family = gaussian()))
if (arm == "lane") P$B4 <- lapply(0:2, function(k) pts(f4$opt$par, k))
V$B4 <- vapply(P$B4, function(p) f4$obj$fn(p), 0)
cat(sprintf("B4 ar(cov = FALSE) [%s]: latent %s, unread %s, logLik %.10f\n", arm,
            paste(f4$frame$mi_map$y$rows, collapse = ","),
            paste(f4$frame$mi_map$y$unread %||% integer(0), collapse = ","),
            as.numeric(logLik(f4))))
# B5: ar() and fcor() on one response
f5 <- tryCatch(quiet(frm(y | mi() ~ x + ar(t, p = 1, cov = FALSE) + fcor(M),
                         data = d4, data2 = list(M = M))),
               error = function(e) conditionMessage(e))
cat("B5 ar + fcor on one response [", arm, "]: ",
    if (is.character(f5)) substr(f5, 1, 160) else "FITTED", "\n", sep = "")
if (arm == "lane") saveRDS(P, ptsf)
saveRDS(V, file.path(outdir, sprintf("read-values-%s.rds", arm)))
if (arm == "base") {
  L <- readRDS(file.path(outdir, "read-values-lane.rds"))
  for (k in names(V)) {
    cat(sprintf("%-12s lane - base objective: %s\n", k,
                paste(sprintf("%+.2e", L[[k]] - V[[k]]), collapse = " ")))
  }
}
