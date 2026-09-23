# Blast radius of the two punch-round-1 start changes, which reach
# every family: make_start() now places an initializer's value into a
# dpar design that has NO intercept, and mu_residuals() now subtracts
# offset(). Both move the start, so both have to be shown not to move
# the optimum. Run on base and on the lane, then compare.
LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
TAG <- Sys.getenv("SKEWINIT_TAG", "fixed")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")

set.seed(31)
n <- 300
x <- rnorm(n)
g <- factor(rep(1:4, length.out = n))
o <- rnorm(n, 2, 0.5)
eta <- 0.4 + 0.8 * x
dd <- data.frame(
  x = x, g = g, o = o, logo = log(runif(n, 1, 5)),
  ygau = eta + rnorm(n, 0, 0.9),
  ypos = rgamma(n, shape = 3, scale = exp(eta) / 3),
  ycnt = rpois(n, exp(0.5 + 0.4 * x + log(runif(n, 1, 5)))),
  ybeta = rbeta(n, 2 * plogis(eta) * 5, 2 * (1 - plogis(eta)) * 5),
  ysn = eta + 1.2 * (0.9486833 * abs(rnorm(n)) + 0.3162278 * rnorm(n))
)
dd$ycnt_off <- dd$ycnt

jobs <- list(
  # dpar designs with NO intercept, across families
  gau_sig_noint   = list(bf(ygau ~ x, sigma ~ 0 + g), gaussian(), list()),
  gau_sig_noint_x = list(bf(ygau ~ x, sigma ~ 0 + x), gaussian(), list()),
  stu_sig_noint   = list(bf(ygau ~ x, sigma ~ 0 + g), student(), list()),
  gam_shape_noint = list(bf(ypos ~ x, shape ~ 0 + g), Gamma(link = "log"),
                         list()),
  beta_phi_noint  = list(bf(ybeta ~ x, phi ~ 0 + g), Beta(), list()),
  nb_shape_noint  = list(bf(ycnt ~ x, shape ~ 0 + g), negbinomial(), list()),
  sn_alpha_noint  = list(bf(ysn ~ x, sigma ~ 1, alpha ~ 0 + g),
                         skew_normal(), list()),
  sn_sig_noint    = list(bf(ysn ~ x, sigma ~ 0 + g, alpha ~ 1),
                         skew_normal(), list()),
  # a no-intercept design on mu itself
  gau_mu_noint    = list(bf(ygau ~ 0 + g), gaussian(), list()),
  pois_mu_noint   = list(bf(ycnt ~ 0 + g), poisson(), list()),
  # offset() models, across families
  gau_off         = list(bf(ygau ~ x + offset(o)), gaussian(), list()),
  gau_off_only    = list(bf(ygau ~ offset(o)), gaussian(), list()),
  pois_off        = list(bf(ycnt ~ x + offset(logo)), poisson(), list()),
  nb_off          = list(bf(ycnt ~ x + offset(logo)), negbinomial(), list()),
  gam_off         = list(bf(ypos ~ x + offset(o)), Gamma(link = "log"),
                         list()),
  stu_off         = list(bf(ygau ~ x + offset(o)), student(), list()),
  sn_off          = list(bf(ysn ~ x + offset(o), sigma ~ 1, alpha ~ 1),
                         skew_normal(), list()),
  sn_off_reml     = list(bf(ysn ~ x + offset(o), sigma ~ 1, alpha ~ 1),
                         skew_normal(), list(REML = TRUE)),
  # controls: no offset, intercept present, must be bitwise identical
  gau_plain       = list(bf(ygau ~ x), gaussian(), list()),
  sn_plain        = list(bf(ysn ~ x, sigma ~ 1, alpha ~ 1), skew_normal(),
                         list())
)

out <- list()
for (nm in names(jobs)) {
  j <- jobs[[nm]]
  r <- tryCatch({
    nw <- 0L
    f <- withCallingHandlers(
      do.call(frm, c(list(j[[1]], family = j[[2]], data = dd), j[[3]])),
      warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
    list(ll = as.numeric(logLik(f)), par = unname(as.numeric(f$opt$par)),
         se = unname(sqrt(diag(vcov(f)))), conv = f$opt$convergence,
         warns = nw, grad = max(abs(f$obj$gr(f$opt$par))))
  }, error = function(e) list(error = conditionMessage(e)))
  out[[nm]] <- r
  cat(sprintf("%-16s %s\n", nm,
              if (is.null(r$error)) sprintf("%s conv %d warn %d",
                                            format(r$ll, digits = 15),
                                            r$conv, r$warns)
              else paste("ERROR", r$error)))
}
p <- paste0("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/",
            "blast-", TAG, ".rds")
saveRDS(out, p)
cat("wrote", p, "\n")
