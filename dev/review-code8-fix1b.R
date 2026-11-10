# Review of lane code8, item 3: the negative direction of fixture 1's
# Richardson Hessian (theta_3, theta_4 held) at the stop. Which
# parameters does it load on, and does the objective fall along it?
#   RV_ARM=after Rscript dev/review-code8-fix1b.R
source("dev/review-code8-common.R")
options(width = 180)
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
f <- quiet(frm(form, family = student(), data = d))
p <- f$opt$par
nm <- ns$outer_par_names(f)
obj <- f$obj
i3 <- match("theta_3", nm); i4 <- match("theta_4", nm)
keep <- setdiff(seq_along(p), c(i3, i4))
at <- function(P, lab) {
  gk <- function(q) { z <- P; z[keep] <- q; drop(obj$gr(z))[keep] }
  for (meth in c("Richardson", "simple")) {
    J <- numDeriv::jacobian(gk, P[keep], method = meth)
    H <- (J + t(J)) / 2
    D <- sqrt(abs(diag(H)))
    e <- eigen(H / outer(D, D), symmetric = TRUE)
    k <- length(e$values)
    cat(lab, meth, "smallest ev:", signif(e$values[(k - 2):k], 4),
        " asym:", signif(max(abs(J - t(J))), 3), "\n")
    for (kk in (k - 1):k) {
      v <- e$vectors[, kk]
      o <- order(-abs(v))[1:5]
      cat("   ev", signif(e$values[kk], 3), ":",
          paste(sprintf("%s %.3f", nm[keep][o], v[o]), collapse = ", "), "\n")
      dn <- numeric(length(P)); dn[keep] <- v / D
      f0 <- obj$fn(P)
      for (t in c(0.01, 0.1, 0.5, 1)) {
        cat(sprintf("      t %4.2f  df+ %.3e  df- %.3e\n", t,
                    obj$fn(P + t * dn) - f0, obj$fn(P - t * dn) - f0))
      }
    }
  }
}
at(p, "stop")
pb <- p; pb[i3] <- -ns$se_cor_clamp_theta()
at(pb, "clamp")
cat("DONE\n")
