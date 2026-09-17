## Reviewer recheck, MAJOR 3 breadth: which columns does the lane store on
## the natural scale, per family, and does each name exist as a parameter
## in brms's own parse of the same model (brms::default_prior(), no
## compile)? Lane arm only. Data seed 63.
##   Rscript dev/brmsnames-rev2-families.R
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
set.seed(63)
n <- 300
d <- data.frame(x = rnorm(n))
d$yg <- 1 + 0.5 * d$x + rnorm(n)
d$ypos <- exp(0.3 + 0.2 * d$x + rnorm(n, 0, 0.4)) + 0.2
d$yp <- rpois(n, exp(0.5 + 0.3 * d$x))
d$yb <- plogis(rnorm(n, 0.2 * d$x, 0.8))
d$yb01 <- ifelse(runif(n) < 0.1, 1, d$yb)
d$ybin <- rbinom(n, 10, plogis(0.3 * d$x))
d$trials <- 10L
d$yang <- atan2(sin(d$x + rnorm(n, 0, 0.5)), cos(d$x + rnorm(n, 0, 0.5)))
d$ymix <- c(rnorm(n / 2, -2), rnorm(n / 2, 2)) + 0.3 * d$x
fams <- list(
  gaussian = list("yg", gaussian()),
  skew_normal = list("yg", skew_normal()),
  student = list("yg", student()),
  asym_laplace = list("yg", asym_laplace()),
  exgaussian = list("ypos", exgaussian()),
  shifted_lognormal = list("ypos", shifted_lognormal()),
  lognormal = list("ypos", lognormal()),
  Gamma = list("ypos", Gamma(link = "log")),
  weibull = list("ypos", weibull()),
  inverse.gaussian = list("ypos", inverse.gaussian(link = "log")),
  negbinomial = list("yp", negbinomial()),
  zero_inflated_negbinomial = list("yp", zero_inflated_negbinomial()),
  hurdle_lognormal = list("yp", hurdle_lognormal()),
  Beta = list("yb", Beta()),
  beta_binomial = list("ybin | trials(trials)", beta_binomial()),
  von_mises = list("yang", von_mises()),
  mixture = list("ymix", mixture(gaussian(), gaussian()))
)
for (nm in names(fams)) {
  F <- fams[[nm]]
  f <- stats::as.formula(paste(F[[1]], "~ x"))
  fit <- tryCatch(q(frm(bf(f), family = F[[2]], data = d)),
                  error = function(e) conditionMessage(e))
  if (is.character(fit)) {
    cat(sprintf("%-26s frm ERROR: %s\n", nm, substr(fit, 1, 100)))
    next
  }
  tab <- frmtmb:::brms_coef_table(fit)
  nat <- tab$brms[tab$natural]
  inv <- attr(tab, "linkinv")[tab$natural]
  fun <- attr(tab, "linkfun")[tab$natural]
  est <- unlist(fit$estimates[c("beta", "betad")])
  # natural value at the mode, through the stored inverse link
  vals <- vapply(seq_along(nat), function(k) {
    i <- which(tab$brms == nat[k])
    int <- c(fit$estimates$beta,
             fit$estimates$betad[setdiff(seq_along(fit$estimates$betad),
                                         fit$frame$betad_fixed_idx)])[i]
    rt <- fun[[k]](inv[[k]](int))
    c(inv[[k]](int), abs(rt - int))
  }, numeric(2))
  bfam <- tryCatch({
    bf2 <- if (nm == "mixture") brms::mixture(brms::gaussian(), brms::gaussian())
      else if (nm %in% c("Gamma", "inverse.gaussian")) F[[2]]
      else get(nm, asNamespace("brms"))()
    dp <- q(brms::default_prior(brms::bf(f), data = d, family = bf2))
    dp <- dp[!nzchar(dp$coef) & !nzchar(dp$group), ]
    unique(ifelse(nzchar(dp$dpar), paste0(dp$class, "|", dp$dpar), dp$class))
  }, error = function(e) paste("brms ERROR:", conditionMessage(e)))
  cat(sprintf("%-26s natural: %-40s value %-30s roundtrip %s\n", nm,
              paste(nat, collapse = ","),
              paste(signif(vals[1, ], 4), collapse = ","),
              paste(signif(vals[2, ], 2), collapse = ",")))
  cat(sprintf("%-26s   brms classes: %s\n", "", paste(bfam, collapse = " ")))
}
