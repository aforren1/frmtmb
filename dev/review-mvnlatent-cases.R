# Reviewer's cases the lane did not test. Lane library first.
#   Rscript dev/review-mvnlatent-cases.R
.libPaths(c("C:/Users/adf44/source/r/wt-mvnlatent-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
cat("latent", dirname(find.package("frmtmb.latent")), "\n")
show <- function(label, expr) {
  r <- tryCatch({
    v <- force(expr)
    paste("OK", paste(utils::capture.output(utils::str(v, max.level = 1,
                                                        give.attr = FALSE)),
                      collapse = " ")[1])
  }, error = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(sprintf("%-55s %s\n", label, substr(r, 1, 220)))
}

## 1. loaded, NOT attached
set.seed(42)
n <- 200
cl <- rbinom(n, 1, 0.4)
x <- rnorm(n)
Y <- cbind(ifelse(cl == 1, 0, 3) + 0.5 * x + rnorm(n),
           ifelse(cl == 1, 0, 4) + rnorm(n))
dd <- data.frame(x = x)
dd$Y <- Y
cat("attached before:", "package:frmtmb.latent" %in% search(), "\n")
fam <- frmtmb.latent::mixture_mvn(K = 2, D = 2)
cat("loaded:", isNamespaceLoaded("frmtmb.latent"),
    " attached:", "package:frmtmb.latent" %in% search(), "\n")
show("frm(family = frmtmb.latent::mixture_mvn())",
     fit <- frm(bf(Y ~ x), family = frmtmb.latent::mixture_mvn(K = 2, D = 2),
                data = dd))
show("frm(bf() + frmtmb.latent::mixture_mvn())",
     fit2 <- frm(bf(Y ~ x) + frmtmb.latent::mixture_mvn(K = 2, D = 2),
                 data = dd))
cat("same logLik both spellings:",
    identical(as.numeric(logLik(fit)), as.numeric(logLik(fit2))), "\n")
show("logLik", logLik(fit))
show("predict()", predict(fit))
show("fitted()", fitted(fit))
show("predict(newdata=)", predict(fit, newdata = dd[1:5, ]))
show("mixture_probs()", mixture_probs(fit))
show("pp_mixture()", pp_mixture(fit))
show("simulate()", simulate(fit, nsim = 2, seed = 1))
show("frm_simulate()", frm_simulate(fit, nsim = 1, seed = 1))
show("conditional_effects(dpar = theta1)",
     conditional_effects(fit, dpar = "theta1"))
show("conditional_effects(method = predict) refuses",
     conditional_effects(fit, method = "predict"))
show("summary()", summary(fit))
show("confint()", confint(fit))
show("frm_compat(REML, mixture_mvn)", frm_compat("REML", "mixture_mvn")$status)
show("REML refusal", frm(bf(Y ~ 1), family = fam, data = dd, REML = TRUE))
show("update(fit) refit", fu <- update(fit))
if (exists("fu")) cat("update logLik identical:",
                      identical(as.numeric(logLik(fu)),
                                as.numeric(logLik(fit))), "\n")
show("update(fit, . ~ 1)", fu1 <- update(fit, formula. = . ~ 1))
show("update(fit, newdata = dd[1:150, ])", update(fit, newdata = dd[1:150, ]))
show("frm_allfit()", frm_allfit(fit))
show("anova(fit_x, fit_1)", if (exists("fu1")) anova(fu1, fit))

## 2. an unqualified call after library(), the user's ordinary spelling,
## saved and refitted in a session that never attaches frmtmb.latent
suppressMessages(library(frmtmb.latent))
fit3 <- frm(bf(Y ~ x) + mixture_mvn(K = 2, D = 2), data = dd)
f3 <- "C:/Users/adf44/source/r/frmtmb-wt-mvnlatent/dev/review-mvnlatent-out/fit3.rds"
saveRDS(list(fit = fit3, dd = dd), f3)
cat("saved fit3\n")

## 3. make_start() vs the moved test's mvn_start_vector()
mvn_start_vector <- function(fit, st) {
  tpl <- fit$frame$par_template
  miss <- setdiff(names(tpl)[lengths(tpl) > 0L], names(st))
  if (length(miss)) stop("the start leaves out ", toString(miss))
  v <- unlist(lapply(names(tpl), function(nm) {
    if (lengths(tpl[nm]) == 0L) numeric(0) else st[[nm]]
  }), use.names = FALSE)
  stopifnot(length(v) == length(fit$obj$par))
  v
}
set.seed(44)
n <- 200
g <- rbinom(n, 1, 0.5) + 1
Y <- cbind(ifelse(g == 1, 0, 3), ifelse(g == 1, 0, 3)) +
  matrix(rnorm(2 * n, sd = 0.9), n, 2)
d4 <- data.frame(row = seq_len(n))
d4$Y <- Y
ls <- log(0.9)
point <- list(
  EII = list(sigmaraw = ls),
  VII = list(sigmaraw1 = ls, sigmaraw2 = ls),
  EEI = list(sigmaraw = c(ls, ls)),
  VEI = list(sigmavol1 = ls, sigmavol2 = ls, sigmashape = 0),
  EVI = list(sigmavol = ls, sigmashape1 = 0, sigmashape2 = 0),
  VVI = list(sigmaraw1 = c(ls, ls), sigmaraw2 = c(ls, ls)),
  EEE = list(sigmaraw = c(ls, ls, 0)),
  VVV = list(sigmaraw1 = c(ls, ls, 0), sigmaraw2 = c(ls, ls, 0)))
for (m in names(point)) {
  st <- c(list(beta = c(0, 0, 3, 3), betad = 0.2), point[[m]])
  f <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = m), data = d4,
           start = st)
  a <- unlist(frmtmb:::make_start(f$frame, st), use.names = FALSE)
  b <- mvn_start_vector(f, st)
  cat(sprintf("start %s identical=%s  obj names order matches tpl=%s\n", m,
              identical(a, b),
              identical(unique(names(f$obj$par)),
                        names(f$frame$par_template)[
                          lengths(f$frame$par_template) > 0])))
}

## 4. us_sigma(): derivatives on the tape against finite differences,
## both arms of floored, d = 1 and d = 3
fd <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(length(f(x))))
}
for (cfg in list(list(d = 1, fl = FALSE), list(d = 1, fl = TRUE),
                 list(d = 3, fl = FALSE), list(d = 3, fl = TRUE))) {
  np <- cfg$d + cfg$d * (cfg$d - 1) / 2
  set.seed(cfg$d)
  th <- rnorm(np, 0, 0.5)
  f <- function(p) as.vector(us_sigma(p, cfg$d, floored = cfg$fl))
  tp <- RTMB::MakeTape(function(p) {
    S <- us_sigma(p, cfg$d, floored = cfg$fl)
    RTMB::advector(as.vector(S))
  }, numeric(np))
  J <- tp$jacobian(th)
  Jn <- fd(f, th)
  cat(sprintf("tape d=%d floored=%s value maxrel=%.2e jac maxrel=%.2e\n",
              cfg$d, cfg$fl,
              max(abs(tp(th) - f(th))) / max(abs(f(th))),
              max(abs(J - Jn)) / max(abs(Jn))))
}
# a wrong length on the tape: refused while taping, not a silent recycle
show("us_sigma wrong length on the tape",
     RTMB::MakeTape(function(p) sum(us_sigma(p, 2)), numeric(4)))
show("us_sigma(theta, 2L) integer d", us_sigma(c(0, 0, 0), 2L))
show("us_sigma(theta, d = 2 as double 2.0000000001)",
     us_sigma(c(0, 0, 0), 2 + 1e-12))
