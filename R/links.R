# AD-safe link registry. stats::make.link clamps at C level in ways the AD
# tape cannot see, so link functions are defined explicitly over plain
# arithmetic that RTMB overloads.
frmtmb_links <- list(
  identity = list(
    name    = "identity",
    linkfun = function(mu) mu,
    linkinv = function(eta) eta,
    mu_eta  = function(eta) rep(1, length(eta))
  ),
  log = list(
    name    = "log",
    linkfun = function(mu) log(mu),
    linkinv = function(eta) exp(eta),
    mu_eta  = function(eta) exp(eta)
  ),
  logit = list(
    name    = "logit",
    linkfun = function(mu) log(mu / (1 - mu)),
    linkinv = function(eta) 1 / (1 + exp(-eta)),
    mu_eta  = function(eta) {
      p <- 1 / (1 + exp(-eta))
      p * (1 - p)
    }
  ),
  cloglog = list(
    name    = "cloglog",
    linkfun = function(mu) log(-log(1 - mu)),
    linkinv = function(eta) 1 - exp(-exp(eta)),
    mu_eta  = function(eta) exp(eta - exp(eta))
  ),
  inverse = list(
    name    = "inverse",
    linkfun = function(mu) 1 / mu,
    linkinv = function(eta) 1 / eta,
    mu_eta  = function(eta) -1 / eta^2
  ),
  # log(x - 1): keeps student-t df above 1 (brms convention for nu)
  logm1 = list(
    name    = "logm1",
    linkfun = function(mu) log(mu - 1),
    linkinv = function(eta) 1 + exp(eta),
    mu_eta  = function(eta) exp(eta)
  ),
  # logit onto (1, 2): the tweedie power parameter's valid range
  power12 = list(
    name    = "power12",
    linkfun = function(mu) log((mu - 1) / (2 - mu)),
    linkinv = function(eta) 1 + 1 / (1 + exp(-eta)),
    mu_eta  = function(eta) {
      p <- 1 / (1 + exp(-eta))
      p * (1 - p)
    }
  )
)

#' Look up a link by name in the AD-safe link registry. An already
#' resolved link list passes through unchanged, so callers can accept
#' either a name or a custom link. An unknown name errors and lists the
#' available links.
#'
#' @noRd
get_link <- function(name) {
  if (is.list(name)) return(name)
  lk <- frmtmb_links[[name]]
  if (is.null(lk)) {
    stop("Unknown link: '", name, "'. Available links: ",
         paste(names(frmtmb_links), collapse = ", "), call. = FALSE)
  }
  lk
}
