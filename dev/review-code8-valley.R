# Review of lane code8, items 1 and 2: what conv_headroom() and the
# gradient verdict report at points short of the maximum on the floor of
# a curved valley (y = A exp(B x), x far from 0, natural A), where a
# straight step leaves the valley. Each point is the profile maximum over
# the other parameters at B = B_hat + m SE(B), so its true shortfall is
# known.
#   RV_ARM=before|after Rscript dev/review-code8-valley.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
res <- list()
for (x0 in c(10, 30)) for (seed in 1:3) {
  set.seed(seed)
  n <- 60
  d <- data.frame(x = runif(n, x0, x0 + 1))
  B0 <- 3 / x0
  d$y <- 2 * exp(B0 * d$x) + rnorm(n, 0, 1)
  r <- rfit(bf(y ~ A * exp(B * x), A ~ 1, B ~ 1, nl = TRUE), data = d,
            start = NULL)
  f <- r$f
  if (r$err) { cat("fit error", conditionMessage(f), "\n"); next }
  p <- f$opt$par
  nm <- ns$outer_par_names(f)
  iA <- match("A_(Intercept)", nm); iB <- match("B_(Intercept)", nm)
  obj <- f$obj
  f0 <- obj$fn(p)
  H <- ns$obj_he_cached(obj, p)
  V <- solve(H)
  seB <- sqrt(V[iB, iB])
  rho <- V[iA, iB] / sqrt(V[iA, iA] * V[iB, iB])
  cat(sprintf("x0 %g seed %d code %d f0 %.10f A %.4g B %.6g seB %.3g rho %.12f conds %s\n",
              x0, seed, f$opt$convergence, f0, p[iA], p[iB], seB, rho,
              cond_txt(r)))
  xbar <- mean(d$x)
  for (m in c(0.25, 0.5, 1, 2, 3, 5, 10, 30)) {
    for (sg in c(-1, 1)) {
      Bm <- p[iB] + sg * m * seB
      st <- p
      st[iB] <- Bm
      st[iA] <- p[iA] * exp(-(Bm - p[iB]) * xbar)
      oth <- setdiff(seq_along(p), iB)
      fo <- function(q) { z <- st; z[oth] <- q; obj$fn(z) }
      go <- function(q) { z <- st; z[oth] <- q; obj$gr(z)[oth] }
      o <- nlminb(st[oth], fo, go, control = list(rel.tol = 1e-15,
                                                   x.tol = 1e-15,
                                                   iter.max = 1000,
                                                   eval.max = 2000))
      P <- st; P[oth] <- o$par
      short <- obj$fn(P) - f0
      g <- drop(obj$gr(P))
      f2 <- f
      f2$cache <- new.env()
      f2$opt$par <- P
      f2$opt$convergence <- 1L
      f2$opt$message <- "false convergence (8)"
      Hp <- ns$obj_he_cached(obj, P)
      Dp <- sqrt(abs(diag(Hp)))
      gsmax <- max(abs(g / Dp))
      pred <- tryCatch(ns$newton_room(Hp, g, 1e-3), error = function(e) NA)
      room <- tryCatch(ns$conv_headroom(f2, g), error = function(e) NA)
      gv <- tryCatch(quiet(ns$grad_verdict(f2, g, f$control)),
                     error = function(e) list(headroom = NA, warn = NA))
      res[[length(res) + 1L]] <- data.frame(x0 = x0, seed = seed,
        m = sg * m, short = short, gsmax = gsmax, pred = pred, room = room,
        gv_head = gv$headroom %||% NA, gv_warn = gv$warn %||% NA)
    }
  }
}
X <- do.call(rbind, res)
options(width = 200)
print(X, digits = 4)
write.table(X, sprintf("dev/review-code8-out/valley-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("room <= 2e-5 while short > 2e-5:",
    sum(X$room <= 2e-5 & X$short > 2e-5, na.rm = TRUE), "\n")
cat("room <= 1e-3 while short > 1e-3:",
    sum(X$room <= 1e-3 & X$short > 1e-3, na.rm = TRUE), "\n")
cat("grad warn FALSE while short > 1e-3 (gradient tripped):",
    sum(!X$gv_warn & X$short > 1e-3 & !is.na(X$gv_head), na.rm = TRUE), "\n")
cat("DONE\n")
