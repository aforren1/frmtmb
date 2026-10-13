# Reviewer: how far a simo MAP fit's optimum moves with the convention.
# Three conventions on one likelihood:
#   lane:    frm(prior = dirichlet(a)), the density of the stereographic
#            coordinates with that chart's Jacobian;
#   natural: the mode of the posterior density of the weights
#            themselves, log L + sum (a - 1) log w, no Jacobian (Stan's
#            optimizing() default, jacobian = FALSE);
#   softmax: the density of softmax coordinates with their Jacobian,
#            log L + sum a log w (the chart frm_sample() samples in).
# The ML fit's tape is the likelihood for the hand optima.
# Usage: Rscript dev/moparity-rev-dirichlet-conv.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")

run_case <- function(label, dm, fo) {
  cat("\n#####", label, "\n")
  ml <- frm(bf(fo), family = gaussian(), data = dm)
  obj <- ml$obj
  p0 <- ml$opt$par
  zi <- which(names(p0) == "zeta1")
  bi <- which(colnames(ml$frame$linpreds[[1]]$X) == "moinc")
  w_of <- function(p) frmtmb::mo_simplex(p[zi])
  fmt <- function(w) paste(formatC(w, format = "f", digits = 4), collapse = " ")
  cat("ML weights:", fmt(w_of(p0)), "| bsp", format(ml$estimates$beta[bi],
                                                   digits = 5), "\n")
  hand_opt <- function(neglp, start) {
    best <- NULL
    set.seed(11)
    for (k in 1:8) {
      st <- start + c(rep(0, length(start) - length(zi)),
                      if (k > 1) rnorm(length(zi), 0, 0.4) else 0 * zi)
      o <- tryCatch(nlminb(st, neglp,
                           control = list(eval.max = 5000, iter.max = 5000,
                                          rel.tol = 1e-14)),
                    error = function(e) NULL)
      if (!is.null(o) && is.finite(o$objective) &&
            (is.null(best) || o$objective < best$objective)) best <- o
    }
    best
  }
  for (a in list(1, 2, c(2, 1.5, 3), 0.5)) {
    alab <- if (length(a) == 1) format(a) else
      paste0("c(", paste(a, collapse = ", "), ")")
    pr <- set_prior(paste0("dirichlet(", alab, ")"), class = "simo",
                    coef = "moinc1")
    mp <- tryCatch(frm(bf(fo), family = gaussian(), data = dm, prior = pr),
                   error = function(e) e)
    aa <- rep(a, length.out = length(zi) + 1L)
    nat <- function(p) {
      w <- w_of(p)
      v <- obj$fn(p) - sum((aa - 1) * log(w))
      if (!is.finite(v)) 1e10 else v
    }
    # in softmax coordinates: the other parameters as they are, the
    # simplex as softmax(c(0, s))
    to_p <- function(q) {
      p <- q
      s <- q[zi]
      w <- exp(c(0, s)) / sum(exp(c(0, s)))
      p[zi] <- frmtmb::mo_coords(w)
      p
    }
    soft <- function(q) {
      s <- q[zi]
      w <- exp(c(0, s)) / sum(exp(c(0, s)))
      v <- obj$fn(to_p(q)) - sum(aa * log(w))
      if (!is.finite(v)) 1e10 else v
    }
    start_nat <- if (inherits(mp, "error")) p0 else mp$opt$par
    on <- hand_opt(nat, start_nat)
    q0 <- start_nat
    w0 <- frmtmb::mo_simplex(start_nat[zi])
    w0 <- pmax(w0, 1e-3); w0 <- w0 / sum(w0)
    q0[zi] <- log(w0[-1] / w0[1])
    os <- hand_opt(soft, q0)
    lab_l <- if (inherits(mp, "error")) {
      paste("REFUSED:", substr(conditionMessage(mp), 1, 70))
    } else {
      paste0(fmt(w_of(mp$opt$par)), " | bsp ",
             format(mp$estimates$beta[bi], digits = 5))
    }
    w_nat <- w_of(on$par)
    w_sm <- {
      s <- os$par[zi]
      exp(c(0, s)) / sum(exp(c(0, s)))
    }
    cat(sprintf("dirichlet(%s)\n  lane (stereo + Jac): %s\n  natural (no Jac):    %s | bsp %s\n  softmax + Jac:       %s | bsp %s\n",
                alab, lab_l, fmt(w_nat), format(on$par[which(names(p0) == "beta")][bi], digits = 5),
                fmt(w_sm), format(os$par[which(names(p0) == "beta")][bi], digits = 5)))
    if (!inherits(mp, "error")) {
      wl <- w_of(mp$opt$par)
      cat(sprintf("  max |w lane - w natural| %.4f | max |w lane - w softmax| %.4f\n",
                  max(abs(wl - w_nat)), max(abs(wl - w_sm))))
      # the log likelihood each optimum gives up against ML
      cat(sprintf("  logLik(ML) - logLik at: lane %.4f natural %.4f softmax %.4f\n",
                  -ml$opt$objective + obj$fn(mp$opt$par),
                  -ml$opt$objective + obj$fn(on$par),
                  -ml$opt$objective + obj$fn(to_p(os$par))))
    }
  }
}

# the lane's own MAP case (dev/moparity-dirichlet-map.R): ML on a face
set.seed(3)
n <- 120
dm <- data.frame(inc = sample(0:3, n, TRUE), z = rnorm(n))
dm$y <- 1 + 0.6 * c(0, 1, 1.4, 1.5)[dm$inc + 1] + 0.3 * dm$z + rnorm(n)
run_case("seed 3, n 120 (ML weight 3 at a face)", dm, y ~ mo(inc) + z)
# a well identified interior simplex
set.seed(4)
n <- 400
dm2 <- data.frame(inc = sample(0:3, n, TRUE), z = rnorm(n))
dm2$y <- 1 + 1.2 * c(0, 0.3, 0.7, 1)[dm2$inc + 1] + 0.3 * dm2$z + rnorm(n, 0, 0.5)
run_case("seed 4, n 400 (interior)", dm2, y ~ mo(inc) + z)
# a weakly identified one
set.seed(5)
n <- 60
dm3 <- data.frame(inc = sample(0:3, n, TRUE), z = rnorm(n))
dm3$y <- 1 + 0.2 * c(0, 0.3, 0.7, 1)[dm3$inc + 1] + rnorm(n)
run_case("seed 5, n 60 (weak)", dm3, y ~ mo(inc) + z)
