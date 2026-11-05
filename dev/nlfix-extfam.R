# The audit of item 1 for the families the extensions register: does an
# nl = TRUE body reach the likelihood? Linear, nl and intercept-only maxima
# on one data set each. Also frmtmb's own custom-family route. Seed 7101.
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  for (p in c("survival", "latent", "eam", "learn", "spatial")) {
    library(paste0("frmtmb.", p), character.only = TRUE)
  }
})
ll <- function(expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    return(structure(NA_real_, err = conditionMessage(r)))
  }
  as.numeric(logLik(r))
}
show <- function(lab, lin, nl, int) {
  a <- ll(lin); b <- ll(nl); c0 <- ll(int)
  cat(sprintf("%-22s linear %10.4f nl %10.4f intercept %10.4f %s\n", lab, a, b,
              c0, if (!is.null(attr(b, "err"))) {
                paste("nl ERR:", substr(attr(b, "err"), 1, 170))
              } else if (!is.null(attr(a, "err"))) {
                paste("lin ERR:", substr(attr(a, "err"), 1, 170))
              } else ""))
}
set.seed(7101)
n <- 300
x <- rnorm(n)
d <- data.frame(x = x, id = factor(rep(1:30, each = 10)),
                tt = rep(1:10, times = 30))
d$time <- rexp(n, exp(0.5 * x))
d$ev <- rbinom(n, 1, 0.8)
d$cc <- 1 - d$ev
d$y <- ifelse(rep(rbinom(30, 1, 0.5), each = 10) == 1, -1, 1) + 0.5 * x +
  rnorm(n, 0, 0.5)
show("cox",
     frm(bf(time | cens(cc) ~ x), family = frmtmb.survival::cox(), data = d),
     frm(bf(time | cens(cc) ~ b1 * x, b1 ~ 1, nl = TRUE),
         family = frmtmb.survival::cox(), data = d),
     frm(bf(time | cens(cc) ~ 1), family = frmtmb.survival::cox(), data = d))
show("royston_parmar",
     frm(bf(time | cens(cc) ~ x), family = frmtmb.survival::royston_parmar(),
         data = d),
     frm(bf(time | cens(cc) ~ b1 * x, b1 ~ 1, nl = TRUE),
         family = frmtmb.survival::royston_parmar(), data = d),
     frm(bf(time | cens(cc) ~ 1), family = frmtmb.survival::royston_parmar(),
         data = d))
fam <- tryCatch(frmtmb.latent::hmm(gaussian(), 2), error = function(e) e)
if (!inherits(fam, "error")) {
  show("hmm", frm(bf(y ~ x), family = fam, data = d),
       frm(bf(y ~ b1 * x, b1 ~ 1, nl = TRUE), family = fam, data = d),
       frm(bf(y ~ 1), family = fam, data = d))
} else cat("hmm: constructor:", conditionMessage(fam), "\n")
d$rt <- 0.3 + rlnorm(n, -0.5, 0.4)
d$dc <- rbinom(n, 1, plogis(0.8 * x))
show("wiener",
     frm(bf(rt | dec(dc) ~ x), family = frmtmb.eam::wiener(), data = d),
     frm(bf(rt | dec(dc) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
         family = frmtmb.eam::wiener(), data = d),
     frm(bf(rt | dec(dc) ~ 1), family = frmtmb.eam::wiener(), data = d))
# a custom family that declares mu and reads it, and one that does not
fam_read <- frmtmb_family("nlfixread", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    RTMB::dnorm(y, dpars$mu, dpars$sigma, log = TRUE)
  })
fam_skip <- frmtmb_family("nlfixskip", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
    RTMB::dnorm(y, 0, dpars$sigma, log = TRUE)
  })
show("custom reads mu", frm(bf(y ~ x), family = fam_read, data = d),
     frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = fam_read,
         data = d),
     frm(bf(y ~ 1), family = fam_read, data = d))
show("custom skips mu", frm(bf(y ~ x), family = fam_skip, data = d),
     frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = fam_skip,
         data = d),
     frm(bf(y ~ 1), family = fam_skip, data = d))
b <- ll(frm(bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = fam_skip,
            data = d))
cat(attr(b, "err"), "\n")
cat("\nb ~ 0 is not a dropped design:\n")
print(ll(frm(bf(y ~ a + b * x, a ~ 1, b ~ 0, nl = TRUE), family = gaussian(),
             data = d)))
print(tryCatch(suppressMessages(as.character(brms::stancode(
  brms::bf(y ~ a + b * x, a ~ 1, b ~ 0, nl = TRUE), data = d))[1] != ""),
  error = function(e) conditionMessage(e)))
print(ll(frm(bf(y ~ 1), family = gaussian(), data = d)))
cat("\noffset-only nlpar is not a dropped design:\n")
d$o <- 0.5
print(ll(frm(bf(y ~ a + b * x, a ~ 1, b ~ 0 + offset(o), nl = TRUE),
             family = gaussian(), data = d)))
print(ll(frm(bf(y ~ 1 + offset(0.5 * x)), family = gaussian(), data = d)))
