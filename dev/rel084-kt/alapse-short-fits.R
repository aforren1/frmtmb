# rel084-kt copy of dev/reviews/alapse-short-fits.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Short review of lane alapse, fix round 1, checks 1, 4 and 5: which
# fits are proved global, which warn frmtmb_ald_local, check_se, the
# profile refusal, and bounds the proof must respect.
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(quantreg))
options(frmtmb.notices = FALSE, width = 150)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")

cap <- function(expr) {
  wc <- character(); w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(z) {
      wc <<- c(wc, class(z)[1]); w <<- c(w, conditionMessage(z))
      invokeRestart("muffleWarning")
    }, message = function(z) invokeRestart("muffleMessage"))
  list(v = v, wc = wc, w = w)
}
show <- function(lab, r) {
  if (inherits(r$v, "error")) {
    cat(sprintf("%-34s ERROR [%s] %s\n", lab, class(r$v)[1],
                substr(gsub("\\s+", " ", conditionMessage(r$v)), 1, 200)))
    if (length(r$wc)) cat(sprintf("%-34s   warnings before: %s\n", lab,
                                  paste(r$wc, collapse = ",")))
    return(invisible())
  }
  f <- r$v
  ce <- f$opt$ald_cert
  cat(sprintf("%-34s code %s | cert ok %-5s global %-5s | warn %s\n", lab,
              f$opt$convergence, isTRUE(ce$ok), isTRUE(ce$global),
              paste(r$wc, collapse = ",")))
  if (!is.null(ce$local)) cat(sprintf("%-34s   local: %s\n", lab, ce$local))
  if (!is.null(ce$why)) cat(sprintf("%-34s   why: %s\n", lab, ce$why))
}

set.seed(2000 + 109)
dB <- data.frame(x = rnorm(100))
dB$y <- 1 + 0.5 * dB$x + rnorm(100)
set.seed(7)
d <- data.frame(x = rnorm(200), z = rnorm(200), o = runif(200, -1, 1),
                w = runif(200, 0.2, 3))
d$y <- 1 + 0.5 * d$x + d$o + exp(0.4 * d$z) * rnorm(200)
d$y2 <- -1 + 0.3 * d$x + rnorm(200)
set.seed(9)
dz <- data.frame(x = rnorm(200))
dz$y <- ifelse(runif(200) < plogis(-1.5 + dz$x), 0,
               1 + 0.5 * dz$x + rnorm(200))

cat("== proved global expected\n")
show("fixq", cap(frm(bf(y ~ x, quantile = 0.3), family = asym_laplace(),
                     data = d)))
show("fixq weights", cap(frm(bf(y | weights(w) ~ x, quantile = 0.3),
                             family = asym_laplace(), data = d)))
rof <- cap(frm(bf(y ~ x + offset(o), quantile = 0.3),
               family = asym_laplace(), data = d))
show("fixq offset", rof)
if (!inherits(rof$v, "error")) {
  b <- ns$fixef_by_dpar(rof$v)$mu
  rq0 <- rq(I(y - o) ~ x, tau = 0.3, data = d)
  ob <- function(bb) { r <- d$y - d$o - drop(cbind(1, d$x) %*% bb)
    sum(r * (0.3 - (r < 0))) }
  cat(sprintf("fixq offset: rq objective rel diff %.3g\n",
              (ob(b) - ob(coef(rq0))) / ob(coef(rq0))))
}
show("estq offset weights", cap(frm(y | weights(w) ~ x + offset(o),
                                    family = asym_laplace(), data = d)))
show("zi fixq", cap(frm(bf(y ~ x, zi ~ x, quantile = 0.5),
                        family = zero_inflated_asym_laplace(), data = dz)))
show("zi estq", cap(frm(bf(y ~ x, zi ~ x),
                        family = zero_inflated_asym_laplace(), data = dz)))
show("zi identity link fixq",
     cap(frm(bf(y ~ x, quantile = 0.5),
             family = zero_inflated_asym_laplace(link_zi = "identity"),
             data = dz)))
show("sigma ~ 1 explicit, estq", cap(frm(bf(y ~ x, sigma ~ 1),
                                         family = asym_laplace(), data = dB)))

cat("== local expected (frmtmb_ald_local)\n")
show("sigma ~ z estq", cap(frm(bf(y ~ x, sigma ~ z), family = asym_laplace(),
                               data = d)))
show("sigma ~ z fixq", cap(frm(bf(y ~ x, sigma ~ z, quantile = 0.2),
                               family = asym_laplace(), data = d)))
show("quantile ~ z", cap(frm(bf(y ~ x, quantile ~ z),
                             family = asym_laplace(), data = d)))
show("multi-response (mvbind)", cap(frm(bf(mvbind(y, y2) ~ x),
                                        family = asym_laplace(), data = d)))
show("multi-response (bf + bf)",
     cap(frm(bf(y ~ x, quantile = 0.5) + bf(y2 ~ x, quantile = 0.5),
             family = asym_laplace(), data = d)))
show("prior density on sigma, fixq",
     cap(frm(bf(y ~ x, quantile = 0.3), family = asym_laplace(), data = d,
             prior = set_prior("student_t(3, 0, 2.5)", class = "sigma"))))
show("prior density on quantile, estq",
     cap(frm(y ~ x, family = asym_laplace(), data = dB,
             prior = set_prior("beta(2, 2)", class = "quantile"))))
show("prior density on b, fixq",
     cap(frm(bf(y ~ x, quantile = 0.3), family = asym_laplace(), data = d,
             prior = set_prior("normal(0, 10)", class = "b"))))
show("sigma softplus link, fixq",
     cap(frm(bf(y ~ x, quantile = 0.3),
             family = asym_laplace(link_sigma = "softplus"), data = d)))

cat("== check_se on a local fit (sigma ~ z)\n")
show("check_se ignore", cap(frm(bf(y ~ x, sigma ~ z), family = asym_laplace(),
                                data = d,
                                control = frmtmb_control(check_se = "ignore"))))
show("check_se stop", cap(frm(bf(y ~ x, sigma ~ z), family = asym_laplace(),
                              data = d,
                              control = frmtmb_control(check_se = "stop"))))
show("check_se ignore, global fixq",
     cap(frm(bf(y ~ x, quantile = 0.3), family = asym_laplace(), data = d,
             control = frmtmb_control(check_se = "ignore"))))

cat("== bounds on the quantile or sigma (B109: global p 0.7247, local 0.4197)\n")
bnd <- function(lab, expr) {
  r <- cap(expr)
  show(lab, r)
  if (!inherits(r$v, "error")) {
    f <- r$v
    dp <- ns$eval_dpars(f)[[1]]
    bx <- f$cache$bounds %||% ns$fit_outer_box(f)
    cat(sprintf("%-34s   p %.4f sigma %.4f logLik %.6f | par %s | lower %s | upper %s | search %s\n",
                lab, dp$quantile[1], dp$sigma[1], as.numeric(logLik(f)),
                paste(format(f$opt$par, digits = 4), collapse = " "),
                paste(format(rep_len(bx$lower %||% -Inf, length(f$opt$par)),
                             digits = 4), collapse = " "),
                paste(format(rep_len(bx$upper %||% Inf, length(f$opt$par)),
                             digits = 4), collapse = " "),
                paste(names(f$opt$ald_cert$search),
                      format(unlist(f$opt$ald_cert$search), digits = 4),
                      collapse = " ")))
  }
}
bnd("quantile ~ 1, uniform(-2, 0.3) Int",
    frm(bf(y ~ x, quantile ~ 1), family = asym_laplace(), data = dB,
        prior = set_prior("uniform(-2, 0.3)", class = "Intercept",
                          dpar = "quantile")))
bnd("quantile ~ 1, Int lb -2 ub 0.3",
    frm(bf(y ~ x, quantile ~ 1), family = asym_laplace(), data = dB,
        prior = set_prior("", class = "Intercept", dpar = "quantile",
                          lb = -2, ub = 0.3)))
bnd("quantile uniform(0.2, 0.6) class",
    frm(y ~ x, family = asym_laplace(), data = dB,
        prior = set_prior("uniform(0.2, 0.6)", class = "quantile",
                          lb = 0.2, ub = 0.6)))
bnd("sigma ~ 1, uniform(-5, -0.3) Int",
    frm(bf(y ~ x, sigma ~ 1), family = asym_laplace(), data = dB,
        prior = set_prior("uniform(-5, -0.3)", class = "Intercept",
                          dpar = "sigma")))

cat("== confint profile/uniroot\n")
fB <- cap(frm(y ~ x, family = asym_laplace(), data = dB))$v
cat("outer names:", ns$outer_par_names(fB), "\n")
for (m in c("profile", "uniroot")) for (pm in list("x", "(Intercept)",
                                                  "sigma_(Intercept)",
                                                  "quantile_(Intercept)")) {
  r <- cap(confint(fB, parm = pm, method = m))
  cat(sprintf("confint %-8s %-22s -> %s%s\n", m, paste(pm, collapse = ","),
              if (inherits(r$v, "error")) paste0("ERROR [", class(r$v)[1], "]")
              else paste(format(r$v, digits = 4), collapse = " "),
              if (length(r$wc)) paste0(" | warn ", paste(r$wc, collapse = ","))
              else ""))
}
fz <- cap(frm(bf(y ~ x, zi ~ x), family = zero_inflated_asym_laplace(),
              data = dz))$v
cat("zi outer names:", ns$outer_par_names(fz), "\n")
g <- glm(I(y == 0) ~ x, family = binomial(), data = dz)
cat("glm profile CI zi:\n"); print(suppressMessages(confint(g)))
for (m in c("profile", "uniroot")) for (pm in list("zi_x", "zi_(Intercept)",
                                                  c("zi_x", "x"))) {
  r <- cap(confint(fz, parm = pm, method = m))
  cat(sprintf("zi confint %-8s %-14s -> %s%s\n", m, paste(pm, collapse = ","),
              if (inherits(r$v, "error")) paste0("ERROR [", class(r$v)[1],
                                                 "] ", substr(conditionMessage(r$v), 1, 80))
              else paste(format(r$v, digits = 5), collapse = " "),
              if (length(r$wc)) paste0(" | warn ", paste(r$wc, collapse = ","))
              else ""))
}
cat("DONE\n")
