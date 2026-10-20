# Reviewer of lane sampleprior, fix round 1: set_prior()'s bound warning
# against brms 2.23.0's (validate_prior() on the same model), on cases
# the lane did not list; and whether frm() (MAP) adds any.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$wide
d$z <- rnorm(nrow(d))
d$t <- rep(1:10, 5)
kinds <- function(ws) {
  k <- c(if (any(grepl("lower bounded", ws))) "lb",
         if (any(grepl("upper bounded", ws))) "ub")
  if (length(k)) paste(k, collapse = "+") else "none"
}
catch <- function(expr) {
  ws <- character(0)
  r <- withCallingHandlers(tryCatch(expr, error = function(e)
    structure(conditionMessage(e), class = "err")),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  list(r = r, ws = ws)
}
cases <- list(
  list("beta(2,2) on b coef x", y ~ x + z, gaussian(),
       quote(set_prior("beta(2, 2)", class = "b", coef = "x")),
       quote(brms::set_prior("beta(2, 2)", class = "b", coef = "x"))),
  list("lognormal(0,1) on sd", y ~ x + (1 | g), gaussian(),
       quote(set_prior("lognormal(0, 1)", class = "sd")),
       quote(brms::set_prior("lognormal(0, 1)", class = "sd"))),
  list("uniform(-5,5) on Intercept", y ~ x, gaussian(),
       quote(set_prior("uniform(-5, 5)", class = "Intercept")),
       quote(brms::set_prior("uniform(-5, 5)", class = "Intercept"))),
  list("gamma(2,1) on b, lb = 0", y ~ x, gaussian(),
       quote(set_prior("gamma(2, 1)", class = "b", lb = 0)),
       quote(brms::set_prior("gamma(2, 1)", class = "b", lb = 0))),
  list("uniform(0,10) on b, lb = 0 only", y ~ x, gaussian(),
       quote(set_prior("uniform(0, 10)", class = "b", lb = 0)),
       quote(brms::set_prior("uniform(0, 10)", class = "b", lb = 0))),
  list("gamma(2,1) on Intercept dpar sigma", bf(y ~ x, sigma ~ x),
       gaussian(),
       quote(set_prior("gamma(2, 1)", class = "Intercept", dpar = "sigma")),
       quote(brms::set_prior("gamma(2, 1)", class = "Intercept",
                             dpar = "sigma"))),
  list("vector: gamma on b, normal on Intercept", y ~ x, gaussian(),
       quote(set_prior(c("gamma(2, 1)", "normal(0, 1)"),
                       class = c("b", "Intercept"))),
       quote(brms::set_prior(c("gamma(2, 1)", "normal(0, 1)"),
                             class = c("b", "Intercept")))),
  list("class lb row + coef gamma row", y ~ x + z, gaussian(),
       quote(set_prior("", class = "b", lb = 0) +
               set_prior("gamma(2, 1)", class = "b", coef = "x")),
       quote(brms::set_prior("", class = "b", lb = 0) +
               brms::set_prior("gamma(2, 1)", class = "b", coef = "x"))),
  list("uniform(-1,1) on ar", y ~ x + ar(time = t, gr = g), gaussian(),
       quote(set_prior("uniform(-1, 1)", class = "ar")),
       quote(brms::set_prior("uniform(-1, 1)", class = "ar"))),
  list("gamma(2,0.1) on nu (student)", y ~ x, frmtmb::brmsfamily("student"),
       quote(set_prior("gamma(2, 0.1)", class = "nu")),
       quote(brms::set_prior("gamma(2, 0.1)", class = "nu"))),
  list("uniform(0,10) on sigma", y ~ x, gaussian(),
       quote(set_prior("uniform(0, 10)", class = "sigma")),
       quote(brms::set_prior("uniform(0, 10)", class = "sigma"))),
  list("exponential(1) on b, ub = 5 only", y ~ x, gaussian(),
       quote(set_prior("exponential(1)", class = "b", ub = 5)),
       quote(brms::set_prior("exponential(1)", class = "b", ub = 5))),
  list("lognormal(0,1) on sd dpar sigma", bf(y ~ x, sigma ~ (1 | g)),
       gaussian(),
       quote(set_prior("lognormal(0, 1)", class = "sd", dpar = "sigma")),
       quote(brms::set_prior("lognormal(0, 1)", class = "sd",
                             dpar = "sigma")))
)
for (cs in cases) {
  fs <- catch(eval(cs[[4]]))
  fm <- if (inherits(fs$r, "err")) list(r = fs$r, ws = character(0)) else
    catch(frm(cs[[2]], data = d, family = cs[[3]], prior = fs$r))
  bfam <- if (inherits(cs[[3]], "family") &&
              identical(cs[[3]]$family, "gaussian")) gaussian() else
    brms::student()
  bform <- if (inherits(cs[[2]], "formula")) cs[[2]] else
    brms::bf(y ~ x, list(sigma ~ x, sigma ~ (1 | g))[[
      1 + grepl("1 | g", paste(deparse(cs[[2]]), collapse = ""), fixed = TRUE)]])
  if (grepl("dpar sigma", cs[[1]])) {
    bform <- if (grepl("sd", cs[[1]])) brms::bf(y ~ x, sigma ~ (1 | g)) else
      brms::bf(y ~ x, sigma ~ x)
  }
  bp <- catch(brms::validate_prior(eval(cs[[5]]), bform, data = d,
                                   family = bfam))
  cat(sprintf("%-40s set_prior: %-6s frm(): %-6s%s | brms: %-6s%s\n", cs[[1]],
              kinds(fs$ws), kinds(fm$ws),
              if (inherits(fm$r, "err")) paste0(" (frm ERROR ", substr(fm$r, 1, 50), ")") else "",
              kinds(bp$ws),
              if (inherits(bp$r, "err")) paste0(" (brms ERROR ", substr(bp$r, 1, 60), ")") else ""))
}
