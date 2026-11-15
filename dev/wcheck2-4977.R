# Lane wcheck2, backlog 4977: a treatment contrast of sigma ~ 1 + fac
# under a normal(0, 1) prior on the Intercept (which frmtmb puts on the
# centered intercept Ic = b0 + means_X'b). With Ic held, the log sigma
# of cell j can fall by t / n_j while that of cell k rises by t / n_k;
# the objective then gains E_j / n_j and loses W_k / n_k per unit of t
# (E: the cell's excess, W: its weight). For each case: the check's
# verdict, then the Laplace objective along that exact direction from
# the start (the sigma coefficients moved by t v, every other parameter
# re-optimized), with the predicted slope per unit of t.
#   Rscript dev/wcheck2-4977.R <base|lane|lib> [two|three]
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
part <- if (length(a) > 1L) a[2] else "two"
P <- function(...) set_prior(..., dpar = "sigma")
pr_ic <- P("normal(0, 1)", class = "Intercept")

# the objective along v in the sigma coefficients, from the start
profile_dir <- function(ob, v, g) {
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  isd <- which(names(p) == "betad")[seq_along(v)]
  base <- p[isd]
  free <- setdiff(seq_along(p), isd)
  vapply(g, function(t) {
    p[isd] <- base + t * v
    fr <- function(r) {
      q2 <- p
      q2[free] <- r
      val <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(val)) 1e10 else val
    }
    op <- tryCatch(nlminb(p[free], fr, control = list(iter.max = 500,
                                                       eval.max = 800)),
                   error = function(e) NULL)
    if (is.null(op)) return(NA_real_)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
}

# one case: weights per cell, the direction (which cell falls, which
# rises), the predicted slope
run_case <- function(lab, d, fo, wcell, down, up, prior = pr_ic) {
  d$w <- unname(wcell[as.character(d$fac)])
  lev <- levels(d$fac)
  n <- nrow(d)
  nj <- table(d$fac)[lev]
  E <- tapply(d$w - 1, d$fac, sum)[lev]
  W <- tapply(d$w, d$fac, sum)[lev]
  pred <- E[[down]] / nj[[down]] - W[[up]] / nj[[up]]
  cat(sprintf("--- %s: weights %s; down %s (E/n %.4g), up %s (W/n %.4g): G %.4g\n",
              lab, paste(names(wcell), wcell, sep = "=", collapse = " "),
              down, E[[down]] / nj[[down]], up, W[[up]] / nj[[up]], pred))
  cat("  frm():", outcome(frm(fo, family = gaussian(), data = d,
                               prior = prior)), "\n")
  # log sigma of each cell: b0 for the reference level, b0 + b_l else.
  # s = per-cell change: -1/n_down, +1/n_up, 0 elsewhere; then b0 =
  # s_ref and b_l = s_l - s_ref
  s <- setNames(rep(0, length(lev)), lev)
  s[[down]] <- -1 / nj[[down]]
  s[[up]] <- 1 / nj[[up]]
  v <- c(s[[1L]], s[-1L] - s[[1L]])
  m <- c(1, as.numeric(nj[-1L]) / n)
  cat(sprintf("  direction (b0, contrasts) %s; Ic change %.3g\n",
              paste(format(v, digits = 4), collapse = " "), sum(m * v)))
  ob <- q(frm(fo, family = gaussian(), data = d, prior = prior,
              dry_run = "objective"))
  g <- c(0, 5, 10, 15, 20, 25) * max(nj)
  vals <- profile_dir(ob, v, g)
  show_profile(sprintf("  profile along t (predicted slope %.4g)", pred),
               g, vals)
}

if (part == "two") {
  # the review's data (dev/reviews/wcheck-review-center2.R, seed 20261111)
  set.seed(20261111)
  d <- data.frame(obs = factor(seq_len(80)), x = rnorm(80),
                  fac = factor(rep(c("a", "b"), 40)))
  d$y <- 0.5 + 0.4 * d$x + rnorm(80, 0, ifelse(d$fac == "a", 0.6, 1.2))
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  run_case("R1", d, fo, c(a = 0.8, b = 1.05), "b", "a")
  run_case("R3", d, fo, c(a = 0.3, b = 1.5), "b", "a")
  run_case("edge +", d, fo, c(a = 0.3, b = 1.4), "b", "a")
  run_case("edge -", d, fo, c(a = 0.3, b = 1.2), "b", "a")
  # the reference level falls, the contrast's cell rises
  run_case("ref down +", d, fo, c(a = 1.4, b = 0.3), "a", "b")
  run_case("ref down -", d, fo, c(a = 1.2, b = 0.3), "a", "b")
  # R2: no prior; the cell rule refuses it
  d$w <- ifelse(d$fac == "a", 0.8, 1.05)
  cat("--- R2 (no prior):", outcome(frm(fo, family = gaussian(), data = d)),
      "\n")
  # a student_t prior on the contrast is free; a normal one bounds it
  for (pb in c("student_t(3, 0, 2.5)", "normal(0, 1)")) {
    d$w <- ifelse(d$fac == "a", 0.3, 1.5)
    cat(sprintf("--- R3 with %s on b too: %s\n", pb,
                outcome(frm(fo, family = gaussian(), data = d,
                            prior = c(pr_ic, P(pb, class = "b"))))))
  }
}

if (part == "three") {
  set.seed(20261113)
  n <- 90
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                  fac = factor(rep(c("a", "b", "c"), 30)))
  d$y <- 0.5 + 0.4 * d$x + rnorm(n)
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  # two contrasts, b0 held
  run_case("c down, b up, +", d, fo, c(a = 1, b = 0.3, c = 1.4), "c", "b")
  run_case("c down, b up, -", d, fo, c(a = 1, b = 0.3, c = 1.2), "c", "b")
  # the reference level falls
  run_case("a down, b up, +", d, fo, c(a = 1.4, b = 0.3, c = 1), "a", "b")
  run_case("a down, b up, -", d, fo, c(a = 1.2, b = 0.3, c = 1), "a", "b")
  # the reference level rises
  run_case("b down, a up, +", d, fo, c(a = 0.3, b = 1.4, c = 1), "b", "a")
  run_case("b down, a up, -", d, fo, c(a = 0.3, b = 1.2, c = 1), "b", "a")
}
cat("DONE\n")
