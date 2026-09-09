# ndt lane: does the per-group bound do what it says, on a small design?
# Seed 4242. Run: Rscript --vanilla dev/ndt-scripts/ndt-group-check.R
.libPaths(c("C:/Users/adf44/source/r/ndt-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

set.seed(4242)
ns <- 6L
nt <- 60L
u <- rnorm(ns, 0, 0.12)
s <- rep(seq_len(ns), each = nt)
cond <- rep(rep(0:1, each = nt / 2L), times = ns)
d <- ddm_simulate(ns * nt, mu = 0.4 + 0.9 * cond, bs = 1.4,
                  ndt = 0.25 * exp(u[s]), bias = 0.5)
d$s <- factor(s)
d$cond <- factor(cond, labels = c("a", "b"))

fit <- frm(bf(rt | dec(upper) + ndt_group(s) ~ cond, bs ~ 1,
              ndt ~ 1 + (1 | s), bias = 0.5),
           family = wiener(), data = d)
cat(sprintf("group fit  logLik %.6f  conv %d\n",
            as.numeric(logLik(fit)), fit$opt$convergence))

rsp <- frmtmb::single_response(fit)
fl <- fit$frame$aterm_values[[rsp$resp_name]]$ndt_floor
cat("ndt_floor length:", length(fl), " unique:", length(unique(fl)), "\n")
cat("per-group floors :", formatC(sort(unique(fl)), digits = 6,
                                  format = "g"), "\n")
cat("tapply(min)      :", formatC(sort(as.numeric(tapply(d$rt, d$s, min))),
                                  digits = 6, format = "g"), "\n")
cat("bound on family  :",
    formatC(rsp$family$ndt_bound$floors, digits = 6, format = "g"), "\n")

frac <- predict(fit, dpar = "ndt", type = "response")
cat("ndt fraction rng :", formatC(range(frac), digits = 6, format = "g"),
    "\n")
cat("ndt seconds rng  :", formatC(range(frac * fl), digits = 6,
                                  format = "g"), "\n")
cat("true ndt rng     :", formatC(range(0.25 * exp(u)), digits = 6,
                                  format = "g"), "\n")
cat("every row's ndt below its own floor:", all(frac * fl < fl), "\n")

# newdata: the bound must come from the FITTED data, not from newdata
one <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
p1 <- predict(fit, newdata = one, dpar = "ndt", type = "response")
cat("newdata fractions:", formatC(p1, digits = 6, format = "g"), "\n")
cat("fitted() on newdata (mean rt):",
    formatC(predict(fit, newdata = one, type = "response"), digits = 6,
            format = "g"), "\n")

# and a group the fit never saw is refused rather than given the
# global bound
bad <- one
levels(bad$s) <- c(levels(d$s), "99")[seq_len(nlevels(bad$s))]
bad$s <- factor(rep("99", nrow(bad)), levels = c(levels(d$s), "99"))
cat("unseen group ->",
    tryCatch({
      predict(fit, newdata = bad, type = "response")
      "NO ERROR (defect)"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 70))),
    "\n")

# max_ndt with ndt_group() is refused
cat("max_ndt + group ->",
    tryCatch({
      frm(bf(rt | dec(upper) + ndt_group(s) ~ cond, bs ~ 1, ndt ~ 1,
             bias = 0.5),
          family = wiener(max_ndt = 0.2), data = d)
      "NO ERROR (defect)"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 60))),
    "\n")

# a character grouping is refused
d2 <- d
d2$sc <- as.character(d2$s)
cat("character group ->",
    tryCatch({
      frm(bf(rt | dec(upper) + ndt_group(sc) ~ cond, bs ~ 1, ndt ~ 1,
             bias = 0.5), family = wiener(), data = d2)
      "NO ERROR (defect)"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 60))),
    "\n")

# gddm declares no ndt_group(), so frame assembly refuses it by name
cat("gddm + group ->",
    tryCatch({
      frm(bf(rt | vint(upper, cond) + ndt_group(s) ~ 1, bias = 0.5),
          family = gddm(), data = transform(d, cond = 1L))
      "NO ERROR (defect)"
    }, error = function(e) paste("refused:",
                                 substr(conditionMessage(e), 1, 80))),
    "\n")

# simulate() and the race families under a grouping
sm <- simulate(fit, nsim = 1, seed = 7)
cat("simulate() rows:", nrow(sm), " any below own floor:",
    sum(sm[[1L]] <= frac * fl), "\n")
