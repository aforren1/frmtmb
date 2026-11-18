# Short review of lane alapse, fix round 1, check 3: tied data.
# Part 1: 10 integer-y (rpois, x in 0:3) and 10 rounded-y (y to one
#   decimal, x in 0:4) data sets, n = 60, quantiles 0.25, 0.5, 0.75, 0.9,
#   unit and U(0.5, 2) weights: ald_vertex_try() (and ald_vertex_ties()
#   alone) against quantreg::rq.wfit(method = "br"), objective relative
#   difference; whether plain ald_vertex() failed (the ties path ran).
# Part 2: the proof under an adversarial perturbation: ald_vertex_ties()
#   with ald_tie_eps raised to 1e-3, 0.05, 0.3 and 1 (a copy of the
#   function evaluated in an environment that redefines the constant;
#   package code untouched), so the perturbed optimum's basis is often
#   not optimal for the data as they are. Every certified b must be at
#   rq's objective; every non-optimal basis must be rejected.
# Part 3: near-ties below the perturbation: y = integer + k * 1e-11.
# Part 4: frm() with a fixed quantile on the part-1 data: code, global,
#   objective against rq at the fit's quantile.
.libPaths(c("C:/Users/adf44/source/r/wt-alapse-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(quantreg))
options(frmtmb.notices = FALSE, width = 150)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
obj <- function(X, y, w, p, b) {
  r <- y - drop(X %*% b)
  sum(w * r * (p - (r < 0)))
}
rqobj <- function(X, y, w, p) {
  f <- suppressWarnings(rq.wfit(X, y, tau = p, weights = w, method = "br"))
  obj(X, y, w, p, f$coefficients)
}
mk <- function(kind, s, n = 60) {
  set.seed(s)
  if (kind == "int") {
    x <- sample(0:3, n, TRUE)
    y <- as.numeric(rpois(n, 2 + x))
  } else {
    x <- sample(0:4, n, TRUE)
    y <- round(1 + 0.4 * x + rnorm(n), 1)
  }
  list(X = cbind(1, x), y = y)
}

cat("== part 1\n")
p1 <- NULL
for (kind in c("int", "round")) for (s in 1:10) {
  dd <- mk(kind, 100 + s)
  for (p in c(0.25, 0.5, 0.75, 0.9)) for (wt in c("unit", "unif")) {
    n <- length(dd$y)
    set.seed(s + 1000)
    w <- if (wt == "unit") rep(1, n) else runif(n, 0.5, 2)
    b0 <- stats::lm.fit(dd$X, dd$y)$coefficients
    plain <- ns$ald_vertex(dd$X, dd$y, w, rep(p, n), b0)
    v <- ns$ald_vertex_try(dd$X, dd$y, w, rep(p, n), b0)
    vt <- ns$ald_vertex_ties(dd$X, dd$y, w, rep(p, n), b0)
    ref <- rqobj(dd$X, dd$y, w, p)
    rel <- function(v) if (is.null(v)) NA else
      (obj(dd$X, dd$y, w, p, v$b) - ref) / ref
    p1 <- rbind(p1, data.frame(kind, s, p, wt, plain_null = is.null(plain),
                               try_ok = !is.null(v), ties_ok = !is.null(vt),
                               rel_try = rel(v), rel_ties = rel(vt),
                               margin = if (is.null(v)) NA else v$margin))
  }
}
print(utils::head(p1[order(-abs(p1$rel_try)), ], 8), digits = 3)
cat(sprintf(paste0("PART1 n %d | plain ald_vertex NULL %d | try certified %d",
                   " | ties certified %d | max |rel| try %.3g | max |rel| ",
                   "ties %.3g | margin 0 count %d\n"),
            nrow(p1), sum(p1$plain_null), sum(p1$try_ok), sum(p1$ties_ok),
            max(abs(p1$rel_try), na.rm = TRUE),
            max(abs(p1$rel_ties), na.rm = TRUE),
            sum(p1$margin <= 1e-12, na.rm = TRUE)))

cat("== part 2\n")
for (eps in c(1e-3, 0.05, 0.3, 1)) {
  e <- new.env(parent = ns)
  e$ald_tie_eps <- eps
  vt_mod <- ns$ald_vertex_ties
  environment(vt_mod) <- e
  n_nonopt <- 0L; n_nonopt_cert <- 0L; n_opt <- 0L; n_opt_cert <- 0L
  n_noperturbed <- 0L; worst_cert <- 0
  for (kind in c("int", "round")) for (s in 1:10) {
    dd <- mk(kind, 100 + s)
    n <- length(dd$y)
    for (p in c(0.25, 0.5, 0.75, 0.9)) for (wt in c("unit", "unif")) {
      set.seed(s + 1000)
      w <- if (wt == "unit") rep(1, n) else runif(n, 0.5, 2)
      b0 <- stats::lm.fit(dd$X, dd$y)$coefficients
      # the basis the perturbed search finds, judged on the data as they are
      scale <- 1 + max(abs(dd$y))
      pat <- ((seq_len(n) * 0.6180339887498949) %% 1) - 0.5
      yp <- dd$y + eps * scale * pat
      vp <- ns$ald_vertex(dd$X, yp, w, rep(p, n), b0)
      ref <- rqobj(dd$X, dd$y, w, p)
      cert <- vt_mod(dd$X, dd$y, w, rep(p, n), b0)
      if (!is.null(cert)) {
        worst_cert <- max(worst_cert,
                          (obj(dd$X, dd$y, w, p, cert$b) - ref) / ref)
      }
      if (is.null(vp)) {
        n_noperturbed <- n_noperturbed + 1L
        next
      }
      bh <- solve(dd$X[vp$h, ], dd$y[vp$h])
      gap <- (obj(dd$X, dd$y, w, p, bh) - ref) / ref
      if (gap > 1e-12) {
        n_nonopt <- n_nonopt + 1L
        n_nonopt_cert <- n_nonopt_cert + !is.null(cert)
      } else {
        n_opt <- n_opt + 1L
        n_opt_cert <- n_opt_cert + !is.null(cert)
      }
    }
  }
  cat(sprintf(paste0("PART2 eps %-6g | perturbed basis not optimal on the ",
                     "data: %d, certified %d | optimal: %d, certified %d | ",
                     "perturbed search NULL %d | worst certified rel %.3g\n"),
              eps, n_nonopt, n_nonopt_cert, n_opt, n_opt_cert, n_noperturbed,
              worst_cert))
}

cat("== part 3\n")
n_cert <- 0L; worst <- 0; n_all <- 0L
for (s in 1:20) {
  set.seed(300 + s)
  x <- sample(0:3, 60, TRUE)
  y <- as.numeric(rpois(60, 2 + x)) + sample(-3:3, 60, TRUE) * 1e-11
  X <- cbind(1, x)
  for (p in c(0.25, 0.5, 0.75)) {
    n_all <- n_all + 1L
    v <- ns$ald_vertex_try(X, y, rep(1, 60), rep(p, 60), c(0, 0))
    ref <- rqobj(X, y, rep(1, 60), p)
    if (!is.null(v)) {
      n_cert <- n_cert + 1L
      worst <- max(worst, (obj(X, y, rep(1, 60), p, v$b) - ref) / ref)
    }
  }
}
cat(sprintf("PART3 near-ties 1e-11: n %d certified %d worst certified rel %.3g\n",
            n_all, n_cert, worst))

cat("== part 4\n")
p4 <- NULL
for (kind in c("int", "round")) for (s in 1:10) {
  dd <- mk(kind, 100 + s)
  d <- data.frame(x = dd$X[, 2], y = dd$y)
  for (p in c(0.25, 0.5, 0.9)) {
    wc <- character()
    f <- withCallingHandlers(
      frm(bf(y ~ x, quantile = p), family = asym_laplace(), data = d),
      warning = function(z) {
        wc <<- c(wc, class(z)[1]); invokeRestart("muffleWarning")
      }, message = function(z) invokeRestart("muffleMessage"))
    b <- ns$fixef_by_dpar(f)$mu
    ref <- rqobj(dd$X, dd$y, rep(1, nrow(d)), p)
    p4 <- rbind(p4, data.frame(kind, s, p, code = f$opt$convergence,
                               global = isTRUE(f$opt$ald_cert$global),
                               ties = isTRUE(f$opt$ald_cert$ok),
                               warn = paste(wc, collapse = ","),
                               rel = (obj(dd$X, dd$y, 1, p, b) - ref) / ref))
  }
}
print(table(p4$warn, p4$code))
cat(sprintf("PART4 n %d | code 0 %d | global %d | max rel %.3g | min rel %.3g\n",
            nrow(p4), sum(p4$code == 0), sum(p4$global), max(p4$rel),
            min(p4$rel)))
cat("DONE\n")
