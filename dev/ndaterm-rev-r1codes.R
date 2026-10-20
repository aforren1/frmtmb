# Reviewer of lane ndaterm, fix round 1: CE holding of vint()/thres()
# variables against brms (dev/ndaterm-rev-r1ce-brms.R), and the level
# tables through update(), a subset refit and saveRDS()/readRDS().
# Usage: Rscript dev/ndaterm-rev-r1codes.R <arm>   (mine | r9). Seed 7.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm", arm, find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste("ERROR:", conditionMessage(e)), class = "err")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [", label, "] WARN (", length(w), "):",
                     substr(unique(w)[1], 1, 120), "\n")
  if (inherits(r, "err")) cat("  [", label, "]", substr(r, 1, 300), "\n")
  invisible(r)
}
ok <- function(r) !inherits(r, "err") && !is.null(r)
set.seed(7)
n <- 200
d <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE),
                gnum = sample(1:3, n, TRUE))
d$y <- rbinom(n, d$N, plogis(0.3 + 0.5 * d$x))
d$yo <- sample(1:4, n, TRUE)
d$gf <- factor(d$gnum)
cbin <- custom_family(
  "rv_cbin", dpars = "mu", links = list(mu = "logit"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dbinom(y, aterms[["vint1"]], dpars[["mu"]], log = TRUE)
  },
  required_aterms = "vint1", type = "discrete",
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vint1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rbinom(n, aterms[["vint1"]], dpars[["mu"]])
  })
cat("\n== CE grid values (brms 2.23.0: vint N 10, thres gr first level)\n")
fv <- frm(bf(y | vint(N) ~ x), family = cbin, data = d)
r <- tr("CE vint", conditional_effects(fv, "x")[[1]])
if (ok(r)) cat("vint(N) CE epred grid N:", unique(r$N), "\n")
r <- tr("CE vint predict", conditional_effects(fv, "x", method = "predict",
                                               ndraws = 101)[[1]])
if (ok(r)) cat("vint(N) CE predict grid N:", unique(r$N),
               " estimate finite:", all(is.finite(r$estimate__)), "\n")
fvx <- frm(bf(y | vint(N) ~ x + N), family = cbin, data = d)
r <- tr("CE vint + predictor", conditional_effects(fvx, "x")[[1]])
if (ok(r)) cat("vint(N) with N a predictor, CE grid N:", unique(r$N), "\n")
ft <- frm(bf(yo | thres(gr = gf) ~ x), family = cumulative(), data = d)
r <- tr("CE thres gr factor", conditional_effects(ft, "x")[[1]])
if (ok(r)) cat("thres(gr = gf) CE grid gf:", as.character(unique(r$gf)),
               "\n")
ftn <- tr("thres gr numeric fit", frm(bf(yo | thres(gr = gnum) ~ x),
                                      family = cumulative(), data = d))
if (ok(ftn)) {
  r <- tr("CE thres gr numeric", conditional_effects(ftn, "x")[[1]])
  if (ok(r)) cat("thres(gr = gnum) CE grid gnum:", unique(r$gnum), "\n")
}

cat("\n== level tables: update(), subset refit, saveRDS()\n")
set.seed(7)
dw <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.3)
dw$x <- rnorm(300)
dw$side <- factor(ifelse(dw$upper == 1, "up", "down"),
                  levels = c("up", "down"))
ff <- frm(bf(rt | dec(side) ~ x, bias = 0.3), family = wiener(), data = dw)
nd <- data.frame(x = 0, side = c("down", "up"))
f0 <- fitted(ff, newdata = nd)[, "Estimate"]
cat("fit: fitted(down, up)", format(f0, digits = 5), "\n")
cat("aterm_codes on spec:", !is.null(ff$spec$responses[[1]]$aterm_codes),
    "\n")
fu <- tr("update formula", update(ff, . ~ . + I(x^2)))
if (ok(fu)) {
  r <- tr("update predict", fitted(fu, newdata = transform(nd, x = 0)))
  if (ok(r)) cat("update(): fitted(down, up)",
                 format(r[, "Estimate"], digits = 5), "\n")
}
fs <- tr("update to one-level subset", update(
  ff, data = dw[dw$side == "down", ]))
if (ok(fs)) {
  r <- tr("subset predict", fitted(fs, newdata = nd[1, ]))
  if (ok(r)) cat("update(data = only 'down' rows): fitted(down)",
                 format(r[, "Estimate"], digits = 5), "\n")
}
fs2 <- tr("frm on one-level subset", frm(bf(rt | dec(side) ~ x,
                                            bias = 0.3),
                                         family = wiener(),
                                         data = dw[dw$side == "down", ]))
if (ok(fs2)) cat("frm() on 'down' rows: coded as",
                 unique(fs2$frame$aterm_values[[1]]$dec), "\n")
tf <- tempfile(fileext = ".rds")
saveRDS(ff, tf)
fr <- readRDS(tf)
r <- tr("readRDS predict", fitted(fr, newdata = nd))
if (ok(r)) cat("readRDS(): fitted(down, up)",
               format(r[, "Estimate"], digits = 5), "\n")
set.seed(3)
r <- tr("readRDS predict draws", predict(fr, newdata = nd, ndraws = 200))
if (ok(r)) cat("readRDS(): predict(down, up)",
               format(r[, "Estimate"], digits = 5), "\n")
unlink(tf)
r <- tr("loo-style influence", frmtmb::influence_measures)
cat("done\n")
