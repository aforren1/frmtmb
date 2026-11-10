# Fixture 1 of the ported brms suite: where does the fit stop, and why.
#   CODE8_ARM=before|after Rscript dev/code8-fix1.R
source("dev/code8-common.R")
options(width = 150)
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
r <- rfit(form, family = student(), data = d)
f <- r$f
cat("objective", format(f$opt$objective, digits = 12), "code",
    f$opt$convergence, f$opt$message, "\n")
cat("CONDITIONS\n"); cat(gsub(" \\|\\| ", "\n", cond_txt(r)), "\n")
p <- f$opt$par
nm <- fns$outer_par_names(f)
g <- drop(f$obj$gr(p))
print(data.frame(name = nm, par = p, grad = g), digits = 6)
si <- stop_info(r)
print(si)
sdr <- quiet(fns$sdr_of(f))
cat("se_lost:\n"); print(sdr$se_lost)
se <- sqrt(diag(as.matrix(sdr$cov.fixed)))
print(data.frame(name = nm, se = se), digits = 4)
fi <- quiet(fitted(f, dpar = "sigma"))
cat("fitted sigma head\n"); print(head(fi))
cat("all(fi > 0):", all(fi > 0), "\n")
cat("theta names:", names(f$estimates$theta), "\n")
print(f$estimates$theta)
cat("re blocks:\n")
for (bk in f$frame$re_blocks) cat(bk$term_label, " theta_idx",
                                   bk$theta_idx, " cov", bk$covstruct, "\n")
saveRDS(list(p = p, nm = nm), "dev/code8-out/fix1-stop.rds")
cat("DONE\n")
