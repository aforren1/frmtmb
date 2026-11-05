# Does a second round buy accuracy? On a scalar block, with the exact optimum
# from integrate(), the estimates' error in SE at importance_rounds = 1 and at
# the default. importance_rounds = 1 is what the mutants one_round and
# grad_always do. Log dev/impsamp-out/scan-rounds.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
sq <- function(e) suppressWarnings(suppressMessages(e))
exact_nll <- function(p, d) {
  s <- exp(p[[3]])
  eta <- p[[1]] + p[[2]] * d$x
  -sum(vapply(split(seq_len(nrow(d)), d$g), function(r) {
    f <- function(u) {
      vapply(u, function(ui) {
        exp(sum(dbinom(d$y[r], 1, plogis(eta[r] + ui), log = TRUE))) *
          dnorm(ui, 0, s)
      }, 0)
    }
    log(integrate(f, -8 * s, 8 * s, subdivisions = 400L)$value)
  }, 0))
}
for (sdu in c(2, 3, 4)) for (m in c(3L, 5L)) for (s in 1:5) {
  set.seed(1000L * s + m)
  G <- 40L
  g <- factor(rep(seq_len(G), each = m))
  x <- round(rnorm(G * m), 3)
  u <- rnorm(G, 0, sdu)
  d <- data.frame(y = rbinom(G * m, 1, plogis(-0.5 + 0.5 * x + u[g])), x = x,
                  g = g)
  f <- bf(y ~ x + (1 | g))
  r <- tryCatch({
    def <- sq(frm(f, family = bernoulli(), data = d, importance = 2000L))
    one <- sq(frm(f, family = bernoulli(), data = d, importance = 2000L,
                  control = frmtmb_control(importance_rounds = 1L)))
    ex <- optim(def$opt$par, exact_nll, d = d, method = "BFGS",
                control = list(reltol = 1e-9))
    se <- sqrt(diag(frmtmb:::sdr_of(def)$cov.fixed))
    list(z_def = max(abs(def$opt$par - ex$par) / se),
         z_one = max(abs(one$opt$par - ex$par) / se),
         rounds = def$importance$rounds, moves = def$importance$moves,
         gap = max(abs(one$opt$par - def$opt$par) / se),
         seed = 1000L * s + m)
  }, error = function(e) NULL)
  if (is.null(r)) next
  cat(sprintf("sd %g m %d seed %d: z default %.4f one round %.4f gap %.4f rounds %d moves %s\n",
              sdu, m, r$seed, r$z_def, r$z_one, r$gap, r$rounds,
              paste(format(r$moves, digits = 3), collapse = ",")))
}
