# Reviewer of lane sampleprior, fix round 1: the bounds the natural-scale
# defaults hand tmbstan, per link, and none where the dpar has a formula.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
show <- function(lab, f, d, fam, expect) {
  r <- tryCatch({
    od <- cap_obj(f, d, fam)
    lo <- CAP$args$lower; up <- CAP$args$upper
    fin <- function(v) if (is.null(v)) "none" else {
      i <- which(is.finite(v))
      if (!length(i)) "none" else
        paste(sprintf("%s[%d]=%.6g", names(od$par)[i], i, v[i]), collapse = " ")
    }
    sprintf("lower %s | upper %s", fin(lo), fin(up))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 90)))
  cat(sprintf("%-38s %s   (expect: %s)\n", lab, r, expect))
}
my <- min(D$shift$y)
show("shifted_ln ndt, log link", y ~ x, D$shift,
     frmtmb::brmsfamily("shifted_lognormal"),
     sprintf("upper log(min_Y) = %.6g on ndt", log(my)))
show("shifted_ln ndt identity link", y ~ x, D$shift,
     frmtmb::brmsfamily("shifted_lognormal", link_ndt = "identity"),
     sprintf("lower 0, upper %.6g on ndt", my))
show("shifted_ln ndt ~ x", bf(y ~ x, ndt ~ x), D$shift,
     frmtmb::brmsfamily("shifted_lognormal"), "none")
show("student nu, logm1", y ~ x, D$real, frmtmb::brmsfamily("student"),
     "none")
show("student nu identity", y ~ x, D$real,
     frmtmb::brmsfamily("student", link_nu = "identity"), "lower 1 on nu")
show("student nu ~ x, identity link", bf(y ~ x, nu ~ x), D$real,
     frmtmb::brmsfamily("student", link_nu = "identity"), "none")
show("Gamma shape identity", y ~ x, D$wpos,
     frmtmb::brmsfamily("Gamma", link_shape = "identity"), "lower 0 on shape")
show("Gamma shape softplus", y ~ x, D$wpos,
     frmtmb::brmsfamily("Gamma", link_shape = "softplus"), "none")
show("zip zi identity", y ~ x, D$zcnt,
     frmtmb::brmsfamily("zero_inflated_poisson", link_zi = "identity"),
     "lower 0, upper 1 on zi")
show("zip zi ~ x, identity", bf(y ~ x, zi ~ x), D$zcnt,
     frmtmb::brmsfamily("zero_inflated_poisson", link_zi = "identity"),
     "none")
show("gaussian sigma identity", y ~ x, D$wide,
     frmtmb::brmsfamily("gaussian", link_sigma = "identity"),
     "lower 0 on sigma")
show("beta phi identity", y ~ x, D$unit,
     frmtmb::brmsfamily("beta", link_phi = "identity"), "lower 0 on phi")
show("mixture sigma2 identity", y ~ x, D$mix,
     frmtmb::mixture(gaussian(), frmtmb::brmsfamily("gaussian",
                                                    link_sigma = "identity")),
     "lower 0 on sigma2")
show("wiener (bounded ndt link)", rt | dec(upper) ~ x, D$wien,
     frmtmb.eam::wiener(), "none")
