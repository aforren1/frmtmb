# Smoke and invariance checks for the brms link options added to
# frmtmb.eam's wiener(): link (mu), link_bs, link_ndt, link_bias.
# A bijective link on a SCALAR dpar reparameterizes the same model, so
# the maximized log-likelihood and the natural-scale estimate must agree
# with the default link's. Seed 20261013.
# Output: dev/famuni-out/wiener-links.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
seed <- 20261013L
set.seed(seed)
n <- 400
d <- ddm_simulate(n, mu = 0.8, bs = 1.3, ndt = 0.3, bias = 0.45)
d$x <- rnorm(n)
d$s <- factor(rep(1:20, each = n / 20))
cat("min rt", min(d$rt), "\n")
try_msg <- function(expr) {
  tryCatch({force(expr); "ok"}, error = function(e) conditionMessage(e))
}
cat("\n-- refusals and spellings\n")
for (lk in c("logit", "probit", "inverse", "sqrt")) {
  cat("link =", lk, ":", try_msg(wiener(link = lk)), "\n")
}
cat("link = log unquoted:", wiener(log)$links$mu$name, "\n")
cat("link = softplus unquoted:", wiener(softplus)$links$mu$name, "\n")
lk <- "squareplus"
cat("link = variable holding squareplus:", wiener(lk)$links$mu$name, "\n")
cat("link_bs = logit:", try_msg(wiener(link_bs = "logit")), "\n")
cat("link_bias = log:", try_msg(wiener(link_bias = "log")), "\n")
cat("link_ndt = logit:", try_msg(wiener(link_ndt = "logit")), "\n")
cat("link_ndt + max_ndt:", try_msg(wiener(link_ndt = "log", max_ndt = 0.2)),
    "\n")
cat("link_ndt + ndt_group:", try_msg(
  frm(bf(rt | dec(upper) + ndt_group(s) ~ 1), family = wiener(link_ndt = "log"),
      data = d)), "\n")

base <- frm(bf(rt | dec(upper) ~ x), family = wiener(), data = d)
ll0 <- as.numeric(logLik(base))
sp0 <- summary(base)$spec_pars
cat(sprintf("\ndefault links: logLik %.10f conv %d\n", ll0,
            base$opt$convergence))
print(sp0)
arms <- list(
  bs_softplus = wiener(link_bs = "softplus"),
  bs_squareplus = wiener(link_bs = "squareplus"),
  bs_identity = wiener(link_bs = "identity"),
  bias_identity = wiener(link_bias = "identity"),
  ndt_log = wiener(link_ndt = "log"),
  ndt_softplus = wiener(link_ndt = "softplus"),
  ndt_squareplus = wiener(link_ndt = "squareplus"),
  ndt_identity = wiener(link_ndt = "identity"),
  all_identity = wiener(link_bs = "identity", link_ndt = "identity",
                        link_bias = "identity"))
cat("\n-- scalar dpars under each link: logLik and estimates against default\n")
for (nm in names(arms)) {
  f <- tryCatch(frm(bf(rt | dec(upper) ~ x), family = arms[[nm]], data = d),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-15s ERROR %s\n", nm, conditionMessage(f))); next
  }
  sp <- summary(f)$spec_pars
  dz <- (sp[, "Estimate"] - sp0[rownames(sp), "Estimate"]) /
    sp0[rownames(sp), "Est.Error"]
  dse <- sp[, "Est.Error"] / sp0[rownames(sp), "Est.Error"]
  cat(sprintf(paste0("%-15s conv %d iters %3d dlogLik %.3e max|dEst|/SE ",
                     "%.3e SE ratio %.6f..%.6f links %s\n"),
              nm, f$opt$convergence, f$opt$iterations %||% NA_integer_,
              as.numeric(logLik(f)) - ll0, max(abs(dz)), min(dse), max(dse),
              paste(vapply(family(f)$links, function(l) l[["name"]], ""),
                    collapse = "/")))
}

cat("\n-- the support guard\n")
fam <- wiener(link_bias = "identity", link_bs = "identity",
              link_ndt = "identity")
ff <- frm(bf(rt | dec(upper) ~ 1), family = fam, data = d)
ffam <- family(ff)
args <- function(bs = 1.3, bias = 0.45, ndt = 0.3) {
  list(mu = 0.8, bs = bs, ndt = ndt, bias = bias)
}
av <- list(dec = d$upper[1:3])
cat("inside:", ffam$lpdf(d$rt[1:3], args(), av), "\n")
cat("default links at the same values:",
    family(base)$lpdf(d$rt[1:3], args(), av), "\n")
cat("bias 1.2:", ffam$lpdf(d$rt[1:3], args(bias = 1.2), av), "\n")
cat("bias 1.2 with the guard off (the default family's slot):",
    family(base)$lpdf(d$rt[1:3], args(bias = 1.2), av), "\n")
cat("bs -1:", ffam$lpdf(d$rt[1:3], args(bs = -1), av), "\n")
cat("ndt -0.05:", ffam$lpdf(d$rt[1:3], args(ndt = -0.05), av), "\n")

cat("\n-- link_ndt = log with formulas, as brms writes them\n")
f1 <- frm(bf(rt | dec(upper) ~ x, ndt ~ x), family = wiener(link_ndt = "log"),
          data = d)
print(summary(f1))
cat("ndt_time() head:", head(ndt_time(f1), 3), "\n")
f2 <- tryCatch(frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1 + (1 | s)),
                   family = wiener(link_ndt = "log"), data = d),
               error = function(e) e)
if (inherits(f2, "error")) cat("ndt ~ (1|s), log: ERROR", conditionMessage(f2),
                               "\n") else {
  cat("ndt ~ (1|s), log: conv", f2$opt$convergence, "logLik",
      as.numeric(logLik(f2)), "\n")
  print(frmtmb::VarCorr(f2))
  nd <- ndt_time(f2)
  cat("rows with ndt >= rt:", sum(nd >= d$rt), "\n")
}
f3 <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1, ndt ~ 1 + (1 | s)),
          family = wiener(), data = d)
cat("ndt ~ (1|s), ndt_group: conv", f3$opt$convergence, "logLik",
    as.numeric(logLik(f3)), "\n")

cat("\n-- link_ndt = log with right censoring at a 1.2 s deadline\n")
dc <- d
dc$cen <- ifelse(dc$rt > 1.2, "right", "none")
dc$rt <- pmin(dc$rt, 1.2)
cat("censored rows:", sum(dc$cen == "right"), "\n")
fc0 <- frm(bf(rt | dec(upper) + cens(cen) ~ 1), family = wiener(), data = dc)
fc1 <- frm(bf(rt | dec(upper) + cens(cen) ~ 1),
           family = wiener(link_ndt = "log"), data = dc)
cat(sprintf("default conv %d logLik %.10f; log conv %d dlogLik %.3e\n",
            fc0$opt$convergence, logLik(fc0), fc1$opt$convergence,
            as.numeric(logLik(fc1)) - as.numeric(logLik(fc0))))

cat("\n-- link_ndt = log with sv, contaminant, mixture\n")
f4 <- tryCatch(frm(bf(rt | dec(upper) ~ 1), family = wiener(link_ndt = "log",
                   variability = "sv"), data = d), error = function(e) e)
cat("sv:", if (inherits(f4, "error")) conditionMessage(f4) else
  sprintf("conv %d logLik %.6f", f4$opt$convergence, logLik(f4)), "\n")
f5 <- tryCatch(frm(bf(rt | dec(upper) ~ 1, lambda ~ 1),
                   family = wiener(link_ndt = "log", contaminant = TRUE,
                                   contaminant_range = c(0, 6)), data = d),
               error = function(e) e)
cat("contaminant:", if (inherits(f5, "error")) conditionMessage(f5) else
  sprintf("conv %d logLik %.6f", f5$opt$convergence, logLik(f5)), "\n")
f6 <- tryCatch(frm(bf(rt | dec(upper) ~ 1),
                   family = frmtmb::mixture(
                     wiener(link_ndt = "log", allow_unreachable = TRUE),
                     lognormal()), data = d), error = function(e) e)
cat("mixture:", if (inherits(f6, "error")) conditionMessage(f6) else
  sprintf("conv %d logLik %.6f", f6$opt$convergence, logLik(f6)), "\n")
