# Probe C: what survives downstream of an nl-body ODE fit.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/lib",
            .libPaths()))
library(frmtmb)
library(RTMBode)
source("dev/ode/pk-common.R")

o <- readRDS("dev/ode/probeA2-fit.rds")
fit <- o$fit; d <- o$d

ck <- function(label, expr) {
  r <- try(force(expr), silent = TRUE)
  if (inherits(r, "try-error")) {
    cat("[FAIL] ", label, ": ", sub("\n$", "", as.character(r)), "\n", sep = "")
    invisible(NULL)
  } else {
    cat("[ok]   ", label, "\n", sep = "")
    invisible(r)
  }
}

cat("=== 1. truth recovery ===\n")
# The anchor is the SAMPLE truth (empirical RE draws), not the population
# constants: 12 subjects leave a visible Monte Carlo offset in the means.
emp <- data.frame(
  quantity = c("lka", "lke", "lV", "sd(lka)", "sd(lke)", "sigma"),
  pop_truth = c(PK_TRUTH$lka, PK_TRUTH$lke, PK_TRUTH$lV,
                PK_TRUTH$sd_lka, PK_TRUTH$sd_lke, PK_TRUTH$sigma),
  sample_truth = c(
    PK_TRUTH$lka + mean(tapply(d$b_ka, d$id, `[`, 1)),
    PK_TRUTH$lke + mean(tapply(d$b_ke, d$id, `[`, 1)),
    PK_TRUTH$lV,
    sd(tapply(d$b_ka, d$id, `[`, 1)), sd(tapply(d$b_ke, d$id, `[`, 1)),
    sd(d$conc - d$mu_true)),
  estimate = c(fixef(fit)$lka[[1]], fixef(fit)$lke[[1]], fixef(fit)$lV[[1]],
               sqrt(VarCorr(fit)[[1]][1, 1]), sqrt(VarCorr(fit)[[2]][1, 1]),
               exp(fixef(fit)$sigma[[1]])))
print(emp, digits = 4)
cat("\nranef vs simulated subject deviations (correlation):\n")
re <- ranef(fit)
cat("  lka:", cor(re[[1]][, 1], tapply(d$b_ka, d$id, `[`, 1)), "\n")
cat("  lke:", cor(re[[2]][, 1], tapply(d$b_ke, d$id, `[`, 1)), "\n")

cat("\n=== 2. sdreport / SEs ===\n")
ck("summary()", { s <- summary(fit); TRUE })
ck("vcov()", { v <- vcov(fit); cat("       dim ", paste(dim(v), collapse = "x"),
                                   "\n"); TRUE })
ck("confint(method='wald')", print(confint(fit, method = "wald")))
ck("confint(method='profile') on lka",
   print(confint(fit, parm = "lka_(Intercept)", method = "profile")))

cat("\n=== 3. predict / fitted ===\n")
p_in <- ck("predict() in-sample", predict(fit))
if (!is.null(p_in)) {
  cat("       max|predict - analytic truth| =",
      format(max(abs(p_in - d$mu_true)), digits = 3), "\n")
  cat("       cor(predict, conc) =", format(cor(p_in, d$conc), digits = 4), "\n")
}
ck("fitted()", { f <- fitted(fit); cat("       identical to predict:",
                                       isTRUE(all.equal(as.numeric(f),
                                                        as.numeric(p_in))),
                                       "\n"); TRUE })
ck("residuals()", { r <- residuals(fit); cat("       sd", sd(r), "\n"); TRUE })
ck("predict(dpar='lka')", { v <- predict(fit, dpar = "lka")
                            cat("       n unique:", length(unique(round(v, 8))),
                                "\n"); TRUE })
ck("predict(se.fit=TRUE)", predict(fit, se.fit = TRUE))

# newdata: dense grid for the observed subjects
nd <- do.call(rbind, lapply(levels(d$id), function(s)
  data.frame(id = factor(s, levels = levels(d$id)),
             time = seq(0.1, 12, length.out = 25), dose = PK_TRUTH$dose)))
p_nd <- ck("predict(newdata=) dense grid, existing subjects",
           predict(fit, newdata = nd))
if (!is.null(p_nd)) {
  cat("       range:", format(range(p_nd), digits = 4), "\n")
  # consistency: newdata restricted to the training rows must reproduce
  p_same <- predict(fit, newdata = d[, c("id", "time", "dose")])
  cat("       newdata==training reproduces in-sample:",
      isTRUE(all.equal(as.numeric(p_same), as.numeric(p_in),
                       tolerance = 1e-8)), "\n")
}
nd_new <- data.frame(id = factor("99"), time = c(0.5, 2, 6), dose = 100)
ck("predict(newdata=) NEW subject, allow_new_levels=TRUE",
   predict(fit, newdata = nd_new, allow_new_levels = TRUE))
ck("predict(newdata=) NEW subject, re.form=NA",
   print(predict(fit, newdata = nd_new, re.form = NA)))
# a newdata grid with only one row per subject: does the per-subject solve
# still work when the caller's grouping differs from training?
nd1 <- data.frame(id = factor(levels(d$id), levels = levels(d$id)),
                  time = 2, dose = 100)
ck("predict(newdata=) one row per subject", print(head(predict(fit, nd1), 3)))

cat("\n=== 4. simulate ===\n")
s <- ck("simulate(nsim=2)", simulate(fit, nsim = 2))
if (!is.null(s)) {
  cat("       dim", paste(dim(as.data.frame(s)), collapse = "x"),
      " col sds:", format(apply(as.data.frame(s), 2, sd), digits = 3), "\n")
}
ck("bootstrap-ish: simulate(nsim=1) refit", {
  ds <- d; ds$conc <- as.data.frame(simulate(fit, nsim = 1))[[1]]
  f2 <- frm(bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
               lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
              gaussian(), data = ds,
            start = list(beta = c(0, log(0.25), log(8))))
  cat("       refit lka/lke/lV:",
      format(c(fixef(f2)$lka[[1]], fixef(f2)$lke[[1]], fixef(f2)$lV[[1]]),
             digits = 4), "\n")
  TRUE
})

cat("\n=== 5. other post-processing ===\n")
ck("logLik/AIC", cat("       logLik", as.numeric(logLik(fit)), " AIC",
                     AIC(fit), "\n"))
ck("VarCorr print", { print(VarCorr(fit)); TRUE })
ck("coef()", { cc <- coef(fit); cat("       names:",
                                    paste(names(cc), collapse = ","), "\n")
               TRUE })
ck("anova/drop1-style: refit without BSV on lke", {
  f3 <- frm(bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
               lka ~ 1 + (1 | id), lke ~ 1, lV ~ 1, nl = TRUE) + gaussian(),
            data = d, start = list(beta = c(0, log(0.25), log(8))))
  cat("       dlogLik =", as.numeric(logLik(fit)) - as.numeric(logLik(f3)),
      "\n")
  TRUE
})
cat("PROBEC DONE\n")
