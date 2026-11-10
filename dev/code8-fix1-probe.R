# Fixture 1: what the SE analysis and the headroom see at the stop, and
# the profile of the parameters at their edge.
#   CODE8_ARM=before|after Rscript dev/code8-fix1-probe.R
source("dev/code8-common.R")
options(width = 150)
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
f <- quiet(frm(form, family = student(), data = d))
p <- f$opt$par
nm <- fns$outer_par_names(f)
np <- length(p)
cat("objective", format(f$opt$objective, digits = 12), "\n")
cat("se_edge_cor:", nm[fns$se_edge_cor(f, np)], "\n")
H <- fns$fit_outer_hessian_at(f, p)
if (is.null(H)) H <- fns$fit_outer_hessian(f)$H
cat("Hessian finite:", all(is.finite(H)), " rows not finite:",
    nm[rowSums(!is.finite(H)) > 0], "\n")
cat("diag H:\n"); print(setNames(signif(diag(H), 4), nm))
cat("se_edge_sd:", nm[fns$se_edge_sd(f, H)], "\n")
an <- tryCatch(fns$se_fit_analysis(f, force = TRUE), error = function(e) e)
if (inherits(an, "error")) cat("analysis error", conditionMessage(an), "\n") else {
  cat("analysis tier", an$tier, "\n lost:\n"); print(an$lost)
}
g <- drop(f$obj$gr(p))
cat("conv_headroom", fns$conv_headroom(f, g), "\n")
# profile: theta_3 (the correlation's Cholesky parameter) and theta_4
prof <- function(j, vals) {
  out <- NULL
  for (v in vals) {
    q0 <- p; q0[j] <- v
    fn <- function(q) { z <- q0; z[-j] <- q; as.numeric(f$obj$fn(z)) }
    gr <- function(q) { z <- q0; z[-j] <- q; drop(f$obj$gr(z))[-j] }
    o <- quiet(nlminb(p[-j], fn, gr, control = list(iter.max = 500,
                                                    eval.max = 500)))
    out <- rbind(out, data.frame(par = nm[j], value = v,
                                 objective = o$objective, code = o$convergence,
                                 cor = if (j == 13) v / sqrt(1 + v^2) else NA))
  }
  out
}
pr <- rbind(prof(13, c(-0.5, -1, -2, -5, -10, -30, -100, p[13])),
            prof(14, c(-1, -3, -6, -10, -20, p[14])))
print(pr, digits = 12)
saveRDS(pr, "dev/code8-out/fix1-profile.rds")
cat("DONE\n")
