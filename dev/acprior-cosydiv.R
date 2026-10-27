# Where do the cosy divergences come from? frm_sample()'s default on
# cosy is flat on (0, 1): a uniform density carried through the window
# logistic onto (-a, 1), with tmbstan's lower bound at the image of 0
# (lb + exp(u)). Arms on the cosy data of dev/acprior-agree.R, 4 chains
# of 2000, seeds 3 to 5, random inits, the same objective otherwise:
#   A  the shipped construction (window logistic, lb on the internal);
#   B  one logistic onto (0, 1): cosy = plogis(v), the internal
#      parameter th = qlogis((cosy + a) / (1 + a)), density on v the
#      logistic Jacobian, no bound (brms's own transform);
#   C  as A with a finite upper bound too, so tmbstan's logit replaces
#      its log transform.
# Log: dev/acprior-out/cosydiv.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-acprior"
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
src <- readLines(file.path(wt, "dev", "acprior-agree.R"))
i0 <- grep("^make_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
dc <- make_data()$cosy
fit <- frm(y ~ x + cosy(time = t, gr = g), data = dc)
ac <- fit$frame$autocor[[1]]
a <- 1 / (ac$d - 1)
rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
ent <- rp$ri$entries
cat("lower:", format(unlist(rp$ri$lower)), "\n")
nll <- build_objective(fit$frame)
other <- Filter(function(e) !identical(e$comp, "thetaac"), ent)
nlp_all <- neg_log_prior_fn(ent)
nlp_oth <- neg_log_prior_fn(other)
pars0 <- fit$estimates
pars0 <- pars0[names(fit$frame$par_template)]
div <- function(sf) sum(sapply(rstan::get_sampler_params(sf, inc_warmup = FALSE),
                               function(m) sum(m[, "divergent__"])))
cosy_of <- function(sf, arm) {
  m <- as.matrix(sf)
  v <- m[, grep("thetaac", colnames(m))]
  if (arm == "B") stats::plogis(v) else -a + (1 + a) * stats::plogis(v)
}
objA <- RTMB::MakeADFun(function(p) nll(p) + nlp_all(p), pars0,
                        silent = TRUE)
fnB <- function(p) {
  v <- p[["thetaac"]]
  cz <- 1 / (1 + exp(-v))
  q <- (cz + a) / (1 + a)
  th <- log(q) - log(1 - q)
  p[["thetaac"]] <- th
  # flat on cosy over (0, 1), read on v: log cosy + log(1 - cosy)
  nll(p) + nlp_oth(p) - (log(cz) + log(1 - cz))
}
pB <- pars0
pB$thetaac <- 0
objB <- RTMB::MakeADFun(fnB, pB, silent = TRUE)
nm <- names(objA$par)
lowA <- rep(-Inf, length(nm))
lowA[nm == "thetaac"] <- unlist(rp$ri$lower)
upC <- rep(Inf, length(nm))
upC[nm == "thetaac"] <- stats::qlogis((1 - 1e-6 + a) / (1 + a))
for (s in 3:5) {
  sA <- suppressWarnings(tmbstan::tmbstan(objA, lower = lowA,
    upper = rep(Inf, length(nm)), chains = 4, iter = 2000, seed = s,
    refresh = 0, init = "random"))
  sB <- suppressWarnings(tmbstan::tmbstan(objB, chains = 4, iter = 2000,
    seed = s, refresh = 0, init = "random"))
  sC <- suppressWarnings(tmbstan::tmbstan(objA, lower = lowA, upper = upC,
    chains = 4, iter = 2000, seed = s, refresh = 0, init = "random"))
  cat(sprintf("seed %d: A %d divergent, cosy mean %.4f | B %d, %.4f | C %d, %.4f\n",
              s, div(sA), mean(cosy_of(sA, "A")), div(sB),
              mean(cosy_of(sB, "B")), div(sC), mean(cosy_of(sC, "C"))))
}
