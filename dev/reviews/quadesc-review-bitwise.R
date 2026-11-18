# Review of lane quadesc: in ONE process (no cross-core ulps), fit each case
# with 0.83.0's escape_stationary(), then with the lane's
# escape_stationary() and escape_track() parsed from the worktree's R/fit.R
# and put into the same namespace, then with 0.83.0's again. Every fit
# whose escape runs and no restart raises must be bitwise the same.
#   Rscript dev/reviews/quadesc-review-bitwise.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.latent)})
cat("lib:", dirname(find.package("frmtmb")), "latent:",
    dirname(find.package("frmtmb.latent")), "\n")
ns <- asNamespace("frmtmb")
ex <- parse("R/fit.R", keep.source = FALSE)
le <- new.env(parent = ns)
grab <- function(nm) {
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
          identical(e[[2]], as.name(nm))) {
      return(eval(e[[3]], le))
    }
  }
  stop("not found: ", nm)
}
le$escape_track <- grab("escape_track")
# the lane's escape finds escape_track() in its enclosure, a child of the
# namespace, since the locked namespace takes no new binding
lane <- list(escape_stationary = grab("escape_stationary"))
base <- list(escape_stationary = ns$escape_stationary)
stopifnot(!identical(deparse(lane$escape_stationary),
                     deparse(base$escape_stationary)))
use <- function(set) {
  for (nm in names(set)) {
    if (exists(nm, ns, inherits = FALSE)) unlockBinding(nm, ns)
    assign(nm, set[[nm]], envir = ns)
  }
}
raised <- 0L
o_opt <- ns$optimize_obj
cases <- list(
  sn_stall = function() {
    set.seed(1)
    n <- 200
    xs <- -abs(rnorm(n)) * 3
    y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
    dd <- data.frame(y = y, xs = xs)
    m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
    frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(), data = dd,
        start = list(beta = c(mean(dd$y), 0),
                     betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3)))
  },
  sn_group = function() {
    set.seed(20)
    g <- factor(rep(1:15, each = 10))
    u <- rnorm(15, sd = 0.8)
    d <- data.frame(g = g, x = rnorm(150))
    d$y <- 1 + 0.5 * d$x + u[g] + rnorm(150)
    frm(bf(y ~ x + (1 | g)), family = skew_normal(), data = d)
  },
  sn_mixture = function() {
    set.seed(101)
    n <- 400
    sn_draw <- function(n, mu, sigma, alpha) {
      dl <- alpha / sqrt(1 + alpha^2)
      mu + sigma * (dl * abs(stats::rnorm(n)) +
                      sqrt(1 - dl^2) * stats::rnorm(n))
    }
    z <- stats::rbinom(n, 1, 0.5)
    dd <- data.frame(y = ifelse(z == 1, sn_draw(n, -3, 1, 5),
                                sn_draw(n, 3, 1, -5)))
    mf <- mixture(skew_normal(), skew_normal())
    u <- frm(bf(y ~ 1), family = mf, data = dd, dry_run = "objective")
    bd <- u$estimates$betad
    bd[grep("^alpha", names(bd))] <- 0
    frm(bf(y ~ 1), family = mf, data = dd, start = list(betad = bd))
  },
  sn_hmm = function() {
    set.seed(1)
    rsn <- function(n, mu, sigma, a) {
      d <- a / sqrt(1 + a^2)
      mu + sigma * (d * abs(stats::rnorm(n)) +
                      sqrt(1 - d^2) * stats::rnorm(n))
    }
    G <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2, byrow = TRUE)
    dd <- do.call(rbind, lapply(seq_len(20), function(id) {
      s <- integer(40)
      s[1L] <- sample.int(2L, 1L)
      for (t in 2:40) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
      y <- ifelse(s == 1L, rsn(40, -3, 1, 5), rsn(40, 3, 1, -5))
      data.frame(id = id, t = seq_len(40), y = y)
    }))
    fam <- hmm(K = 2, skew_normal(), time = t, group = id)
    u <- frm(bf(y ~ 1), family = fam, data = dd, dry_run = "objective")
    bd <- u$estimates$betad
    bd[["alpha1_(Intercept)"]] <- 0
    frm(bf(y ~ 1), family = fam, data = dd, start = list(betad = bd))
  },
  sn_quad_symsmall = function() {
    set.seed(1001)
    d <- data.frame(x = rnorm(60), g = factor(rep(1:15, each = 4)))
    d$y <- 0.5 * d$x + rnorm(15, 0, 0.8)[d$g] + rnorm(60)
    frm(bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ 1), family = skew_normal(),
        data = d, quadrature = TRUE, start = list(betad = c(0, 0)))
  },
  gaussian_group = function() {
    set.seed(3)
    d <- data.frame(x = rnorm(100), g = factor(rep(1:10, each = 10)))
    d$y <- 1 + d$x + rnorm(10)[d$g] + rnorm(100)
    frm(y ~ x + (1 | g), data = d)
  }
)
run <- function(f) {
  raised <<- 0L
  r <- withCallingHandlers(f(), warning = function(w) {
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  r
}
for (nm in names(cases)) {
  use(base)
  a <- run(cases[[nm]])
  use(lane)
  b <- run(cases[[nm]])
  use(base)
  c2 <- run(cases[[nm]])
  esc <- a$opt$stationary_escape
  cat(sprintf(paste0("%-18s escape %s | base==base %s | lane==base: par %s",
                     " objective %s logLik %s code %d/%d\n"),
              nm, if (is.null(esc)) "none" else
                paste(format(esc, digits = 4), collapse = "/"),
              identical(a$opt$par, c2$opt$par),
              identical(a$opt$par, b$opt$par),
              identical(a$opt$objective, b$opt$objective),
              identical(as.numeric(logLik(a)), as.numeric(logLik(b))),
              a$opt$convergence, b$opt$convergence))
}
cat("END\n")
