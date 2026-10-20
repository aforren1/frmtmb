# Reviewer of lane sampleprior, fix round 1: is the shared-nu exemption
# narrow? A Student-t mvbf WITHOUT rescor (one nu per response), a
# gaussian rescor model (no nu), against brms 2.23.0 on the same priors.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$mv
try_frm <- function(lab, f, fam, pr) {
  r <- tryCatch({
    od <- cap_obj(f, d, fam, prior = pr)
    paste("captured,", length(od$par), "parameters")
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 110)))
  cat(sprintf("frmtmb %-40s %s\n", lab, r))
}
try_brm <- function(lab, f, fam, pr) {
  r <- tryCatch({
    v <- suppressWarnings(brms::validate_prior(pr, f, data = d, family = fam))
    rows <- v[v$class == "nu" & nzchar(v$prior), c("prior", "resp", "source")]
    paste("ok;", paste(apply(rows, 1, paste, collapse = "/"), collapse = "; "))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 110)))
  cat(sprintf("brms   %-40s %s\n", lab, r))
}
fs <- frmtmb::brmsfamily("student")
f_no <- bf(mvbind(y1, y2) ~ x) + set_rescor(FALSE)
b_no <- brms::bf(brms::mvbind(y1, y2) ~ x) + brms::set_rescor(FALSE)
f_g <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
b_g <- brms::bf(brms::mvbind(y1, y2) ~ x) + brms::set_rescor(TRUE)
try_frm("student, no rescor, defaults", f_no, fs, NULL)
try_frm("student, no rescor, nu without resp", f_no, fs,
        set_prior("gamma(2, 0.1)", class = "nu"))
try_brm("student, no rescor, nu without resp", b_no, brms::student(),
        brms::set_prior("gamma(2, 0.1)", class = "nu"))
try_frm("student, no rescor, nu resp = y1", f_no, fs,
        set_prior("gamma(2, 0.1)", class = "nu", resp = "y1"))
try_brm("student, no rescor, nu resp = y1", b_no, brms::student(),
        brms::set_prior("gamma(2, 0.1)", class = "nu", resp = "y1"))
try_frm("gaussian rescor, nu without resp", f_g, gaussian(),
        set_prior("gamma(2, 0.1)", class = "nu"))
try_brm("gaussian rescor, nu without resp", b_g, gaussian(),
        brms::set_prior("gamma(2, 0.1)", class = "nu"))
f_r <- bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE)
try_frm("student rescor, nu without resp", f_r, fs,
        set_prior("gamma(3, 0.1)", class = "nu"))
try_brm("student rescor, nu without resp", b_g, brms::student(),
        brms::set_prior("gamma(3, 0.1)", class = "nu"))
# a sigma prior without resp on the rescor model must still be refused
try_frm("student rescor, sigma without resp", f_r, fs,
        set_prior("exponential(1)", class = "sigma"))
try_brm("student rescor, sigma without resp", b_g, brms::student(),
        brms::set_prior("exponential(1)", class = "sigma"))
