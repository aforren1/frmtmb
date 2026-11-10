# Fixture 1, part 2: the profile of the edge parameters, reductions of the
# model in lme4 and glmmTMB beside frmtmb, and the alias of s(Age)'s
# unpenalized column with Age.
#   CODE8_ARM=before|after Rscript dev/code8-fix1-probe2.R
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
# the alias: s(Age)'s fixed column against the intercept and Age
X <- f$frame$linpreds[[1]]$X
cat("X columns:", colnames(X), "\n")
if (!is.null(X)) {
  X <- as.matrix(X)
  cat("column sds:", format(apply(X, 2, sd), digits = 4), "\n")
  q <- qr(X, tol = 1e-9)
  cat("rank of X:", q$rank, "of", ncol(X), "\n")
  fx <- grep("fx", colnames(X))
  if (length(fx)) {
    l <- lm(X[, fx] ~ d$Age)
    cat("s(Age).fx1 on Age: R^2 =", format(summary(l)$r.squared, digits = 17),
        "\n")
  }
}
prof <- function(j, vals) {
  out <- NULL
  for (v in vals) {
    q0 <- p; q0[j] <- v
    fn <- function(q) { z <- q0; z[-j] <- q; as.numeric(f$obj$fn(z)) }
    gr <- function(q) { z <- q0; z[-j] <- q; drop(f$obj$gr(z))[-j] }
    st <- p[-j]
    o <- tryCatch(quiet(nlminb(st, fn, gr, control = list(iter.max = 500,
                                                          eval.max = 500))),
                  error = function(e) quiet(nlminb(st, fn, control =
                    list(iter.max = 500, eval.max = 2000))))
    out <- rbind(out, data.frame(par = nm[j], value = v,
                                 objective = o$objective, code = o$convergence,
                                 cor = if (nm[j] == "theta_3")
                                   v / sqrt(1 + v^2) else NA))
  }
  out
}
j3 <- match("theta_3", nm)
pr <- prof(j3, c(-0.5, -1, -2, -5, -10, -30, -100, p[j3]))
print(pr, digits = 12)
cat("== reductions\n")
d$Expn <- as.numeric(d$Exp)
red <- list(
  r1 = count ~ Trt * Age + volume + offset(Age) + (1 + Trt | visit),
  r2 = count ~ Trt * Age + Expn + volume + offset(Age) + (1 + Trt | visit))
for (k in names(red)) {
  fo <- red[[k]]
  if (requireNamespace("lme4", quietly = TRUE)) {
    m <- quiet(lme4::lmer(fo, data = d, REML = FALSE))
    vc <- as.data.frame(lme4::VarCorr(m))
    cat(k, "lmer gaussian: logLik", format(as.numeric(logLik(m)), digits = 10),
        " singular", lme4::isSingular(m), " cor",
        format(vc$sdcor[vc$var1 == "(Intercept)" & !is.na(vc$var2)],
               digits = 6), "\n")
  }
  fg <- rfit(bf(fo), family = gaussian(), data = d)
  cat(k, "frm gaussian: logLik", format(ll_of(fg$f), digits = 10), " code",
      fg$f$opt$convergence, " theta", format(fg$f$estimates$theta, digits = 4),
      "\n   ", substr(cond_txt(fg), 1, 400), "\n")
  if (requireNamespace("glmmTMB", quietly = TRUE)) {
    g <- tryCatch(quiet(glmmTMB::glmmTMB(fo, dispformula = ~ Trt,
                                         family = glmmTMB::t_family(),
                                         data = d)), error = function(e) e)
    if (inherits(g, "error")) cat(k, "glmmTMB error", conditionMessage(g), "\n")
    else {
      vc <- glmmTMB::VarCorr(g)$cond$visit
      cat(k, "glmmTMB t: logLik", format(as.numeric(logLik(g)), digits = 10),
          " pdHess", g$sdr$pdHess, " cor", format(attr(vc, "correlation")[1, 2],
                                                  digits = 6),
          " conv", g$fit$convergence, g$fit$message, "\n")
      s <- summary(g)$coefficients$disp
      cat("   disp coefs:\n"); print(s)
    }
  }
  ft <- rfit(bf(fo, sigma ~ Trt), family = student(), data = d)
  cat(k, "frm student: logLik", format(ll_of(ft$f), digits = 10), " code",
      ft$f$opt$convergence, " theta", format(ft$f$estimates$theta, digits = 4),
      "\n   ", substr(cond_txt(ft), 1, 400), "\n")
  print(tryCatch(quiet(summary(ft$f)), error = function(e) conditionMessage(e)))
}
cat("DONE\n")
