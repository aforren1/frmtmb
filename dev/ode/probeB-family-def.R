# Custom family used by probe B: the whole PK likelihood, ODE solve
# included, lives in lpdf; times / subject ids / doses ride along as
# vreal()/vint() addition terms.
pk_family <- custom_family(
  "pk1cmt_oral",
  dpars = c("lka", "lke", "lV", "lsigma"),
  links = list(lka = "identity", lke = "identity", lV = "identity",
               lsigma = "identity"),
  primary_dpars = "lka",
  lpdf = function(y, dpars, aterms) {
    mu <- pk_ode(exp(dpars[["lka"]]), exp(dpars[["lke"]]),
                 exp(dpars[["lV"]]), aterms[["vreal1"]], aterms[["vint1"]],
                 aterms[["vreal2"]])
    RTMB::dnorm(y, mu, exp(dpars[["lsigma"]]), log = TRUE)
  },
  post = list(
    mean_fn = function(dpars, aterms)
      pk_ode(exp(dpars[["lka"]]), exp(dpars[["lke"]]), exp(dpars[["lV"]]),
             aterms[["vreal1"]], aterms[["vint1"]], aterms[["vreal2"]]),
    var_fn = function(dpars, aterms) exp(dpars[["lsigma"]])^2
  ),
  sim = function(dpars, aterms, n)
    rnorm(n, pk_ode(exp(dpars[["lka"]]), exp(dpars[["lke"]]),
                    exp(dpars[["lV"]]), aterms[["vreal1"]],
                    aterms[["vint1"]], aterms[["vreal2"]]),
          exp(dpars[["lsigma"]])),
  type = "continuous")
