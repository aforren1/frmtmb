# rev-ndt: the backward-compatibility control, widened and sharpened.
#
# With no ndt_group() and no max_ndt the new parameterization is the old
# one written differently, so every number a user can reach must be
# BITWISE identical. This runs the same models on both builds and
# serializes the full doubles; rev-ndt-bc-compare.R does the comparison
# with identical() and a ulp count, because a printed digit is not a
# bit.
#
# REV_ARM = "new" | "old". Seed 4242 throughout, as the lane's
# ndt-smoke.R used.

arm <- Sys.getenv("REV_ARM", "new")
# REV_LIB picks the round: rev-ndt-lib was round one, rev-ndt-lib2
# is the build after the scalar bound went back into the link.
new_lib <- Sys.getenv("REV_LIB",
                      "C:/Users/adf44/source/r/rev-ndt-lib2")
ref_lib <- "C:/Users/adf44/source/r/reflib-r2"
usr_lib <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "new")) c(new_lib, ref_lib, usr_lib)
          else c(ref_lib, usr_lib))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm", arm, "eam", as.character(utils::packageVersion("frmtmb.eam")),
    "from", dirname(dirname(getNamespaceInfo("frmtmb.eam", "path"))),
    "\n")

res <- list()
grab <- function(nm, fit, dpars) {
  v <- list(
    logLik = as.numeric(stats::logLik(fit)),
    par = unname(fit$opt$par),
    par_names = names(fit$opt$par),
    conv = fit$opt$convergence,
    fixef = unlist(fixef(fit)),
    fitted = tryCatch(as.numeric(fitted(fit)),
                      error = function(e) conditionMessage(e)),
    resid = tryCatch(as.numeric(residuals(fit)),
                     error = function(e) conditionMessage(e)),
    vcov = tryCatch(as.numeric(vcov(fit)),
                    error = function(e) conditionMessage(e)))
  for (dp in dpars) {
    v[[paste0("resp_", dp)]] <- tryCatch(
      as.numeric(suppressWarnings(
        stats::predict(fit, dpar = dp, type = "response"))),
      error = function(e) paste("ERR:", conditionMessage(e)))
  }
  res[[nm]] <<- v
  cat(" ", nm, "logLik", sprintf("%.12f", v$logLik), "\n")
  invisible(NULL)
}

# ------------------------------------------------------------ wiener
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
res$min_rt <- min(d$rt)

fw <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = wiener(), data = d, se = TRUE)
grab("wiener", fw, c("mu", "bs", "ndt"))

fwm <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
           family = wiener(max_ndt = 0.15), data = d, se = TRUE)
grab("wiener_maxndt", fwm, c("ndt"))

for (vv in list("sv", "sz", "st", c("sv", "sz", "st"))) {
  nm <- paste0("wiener_", paste(vv, collapse = ""))
  fv <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1, bias = 0.5),
            family = wiener(variability = vv), data = d, se = FALSE)
  grab(nm, fv, c("ndt", if ("st" %in% vv) "st"))
}

# a random effect on ndt with NO ndt_group(): the old model exactly
fre <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt ~ 1 + (1 | s),
              bias = 0.5), family = wiener(), data = d, se = FALSE)
grab("wiener_re", fre, c("ndt"))
res$wiener_re_vc <- as.numeric(VarCorr(fre)[[1L]])

# A FIXED ndt in bf(). The constant is on the dpar's NATURAL scale, and
# the natural scale of `ndt` is what this change redefined, so a model
# that pinned it is the sharpest backward-compatibility case there is.
res$fixed_bare <- tryCatch({
  f <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt = 0.2, bias = 0.5),
           family = wiener(), data = d, se = FALSE)
  list(logLik = as.numeric(stats::logLik(f)),
       resp = as.numeric(suppressWarnings(
         stats::predict(f, dpar = "ndt", type = "response"))[1L]))
}, error = function(e) paste("ERR:", conditionMessage(e)))

res$fixed_maxndt <- tryCatch({
  f <- frm(bf(rt | dec(upper) ~ cond, bs ~ 1, ndt = 0.2, bias = 0.5),
           family = wiener(max_ndt = 0.25), data = d, se = FALSE)
  list(logLik = as.numeric(stats::logLik(f)),
       resp = as.numeric(suppressWarnings(
         stats::predict(f, dpar = "ndt", type = "response"))[1L]),
       fitted1 = as.numeric(fitted(f))[1L])
}, error = function(e) paste("ERR:", conditionMessage(e)))
str(res$fixed_bare)
str(res$fixed_maxndt)

# simulate(), through the wrapped sim slot
set.seed(11)
sm <- suppressWarnings(stats::simulate(fw, nsim = 2L))
res$sim_wiener <- as.numeric(unlist(lapply(sm, function(z)
  if (is.data.frame(z)) unlist(lapply(z, as.numeric)) else
    as.numeric(z))))

# --------------------------------------------------------------- rdm
set.seed(4242)
dr <- rdm_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.5, ndt = 0.2)
fr <- frm(bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
          family = rdm(2), data = dr, se = TRUE)
grab("rdm", fr, c("ndt", "k"))

# the censored path, which reaches the wrapped lccdf
dr2 <- dr
dr2$cn <- ifelse(dr2$rt > 0.8, "right", "none")
dr2$rtc <- pmin(dr2$rt, 0.8)
frc <- frm(bf(rtc | vint(choice) + cens(cn) ~ 1, v2 ~ 1, A ~ 1, k ~ 1,
              ndt ~ 1), family = rdm(2), data = dr2, se = FALSE)
grab("rdm_cens", frc, c("ndt"))

# --------------------------------------------------------------- lba
set.seed(4242)
dl <- lba_simulate(300, v = c(2.5, 1.5), A = 0.5, k = 0.4, ndt = 0.2)
fl <- frm(bf(rt | vint(choice) ~ 1, v2 ~ 1, A ~ 1, k ~ 1, ndt ~ 1),
          family = lba(2), data = dl, se = TRUE)
grab("lba", fl, c("ndt"))

# -------------------------------------------------------- wiener_gng
set.seed(4242)
dg <- wiener_gng_simulate(400, mu = 1.2, bs = 1.4, ndt = 0.25,
                          deadline = 1.5)
fg <- frm(bf(rt | dec(responded) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = wiener_gng(deadline = 1.5), data = dg, se = TRUE)
grab("wiener_gng", fg, c("ndt"))

# -------------------------------------------------------------- gddm
set.seed(5)
dq <- gddm_simulate(120, mu = 2, bs = 2.5, ndt = 0.25,
                    control = gddm_control(t_max = 2))
dq$cond <- 1L
fq <- frm(bf(rt | vint(upper, cond) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
          family = gddm(control = gddm_control(t_max = 2, dt = 0.05,
                                               ny = 51L)),
          data = dq, se = TRUE)
grab("gddm", fq, c("ndt"))

saveRDS(res, file.path("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev",
                       paste0("rev-ndt-bc-", arm, ".rds")))
cat("ARM", arm, "OK", length(res), "entries\n")
